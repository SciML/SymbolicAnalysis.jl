# Single read-only, memoized, bottom-up fold that replaces the Postwalk
# rebuild passes in `analyze`. Sign, curvature, and (when a manifold is given)
# geodesic curvature are computed once per node from the children's results and
# stored in an `IdDict` keyed by the node. No terms are reconstructed and no
# per-node metadata is written.
#
# The public `propagate_sign` / `propagate_curvature` / `propagate_gcurvature`
# walks keep writing metadata for callers that read annotated trees. While the
# fold is running, the getters consult the memo so rule logic that historically
# read child metadata (`getsign`, `find_curvature`, `increasing_if_positive`,
# …) sees the same values without a rebuild.
#
# The active fold context is stored in task-local storage (not a process-global
# `Ref`) so concurrent `analyze` calls on different tasks/threads cannot read
# each other's memos, and a `try`/`finally` clears the binding so nothing
# survives a throw.

"""
    NodeAnalysis

Per-node result of the read-only analysis fold. `sign` is ordinarily a
[`Sign`](@ref); atoms that store callable placeholders in the rule table
(e.g. `perspective`) may leave a `Function` here so that a later `::Sign`
assertion fails exactly as the metadata-based engine did.
"""
struct NodeAnalysis
    sign::Any
    curvature::Curvature
    gcurvature::Union{GCurvature, Nothing}
end

mutable struct FoldCtx
    memo::IdDict{Any, NodeAnalysis}
    M::Union{AbstractManifold, Nothing}
end

# Task-local key for the active fold. Must not be a process-global `Ref`: with
# several threads one Euclidean `analyze` was observed reading another thread's
# SPD memo and returning a false certificate.
const FOLD_TLS_KEY = :SymbolicAnalysis_fold_ctx

@inline function current_fold_ctx()
    return get(task_local_storage(), FOLD_TLS_KEY, nothing)
end

@inline function fold_lookup(ex)
    ctx = current_fold_ctx()
    ctx === nothing && return nothing
    return get(ctx.memo, ex, nothing)
end

function with_fold_ctx(f, ctx::FoldCtx)
    tls = task_local_storage()
    had_old = haskey(tls, FOLD_TLS_KEY)
    old = had_old ? tls[FOLD_TLS_KEY] : nothing
    tls[FOLD_TLS_KEY] = ctx
    try
        return f()
    finally
        if had_old
            tls[FOLD_TLS_KEY] = old
        else
            delete!(tls, FOLD_TLS_KEY)
        end
    end
end

# Recursively analyze `ex` bottom-up. Children that are symbols or calls are
# folded first so that getter lookups during this node's rule application hit
# the memo. Non-symbolic leaves (numbers, constant arrays) are not memoized;
# the existing getters already classify them.
function analyze_node!(ex, ctx::FoldCtx)
    hit = get(ctx.memo, ex, nothing)
    hit !== nothing && return hit

    # A bare `Vector{Num}` / `Matrix{Num}` (e.g. `conv` after Symbolics expands
    # it elementwise, or a broadcast that returned a Julia array) is not a
    # symbolic call. Postwalk never enters such containers, so elements stay
    # unannotated and `analyze` reads them through `getsign`/`getcurvature`'s
    # `AbstractArray` aggregators (which re-derive properties without metadata).
    # Do not pre-fold the elements into the memo — that would make composition
    # see child signs the Postwalk engine never wrote, and over-certify.
    if ex isa AbstractArray
        s = getsign(ex)
        c = getcurvature(ex)
        g = isnothing(ctx.M) ? nothing : getgcurvature(ex)
        props = NodeAnalysis(s, c, g)
        ctx.memo[ex] = props
        return props
    end

    if iscall(ex)
        for a in arguments(ex)
            if issym(a) || iscall(a)
                analyze_node!(a, ctx)
            end
            # Julia `AbstractArray` arguments are not entered by Postwalk either.
        end
    elseif !(issym(ex))
        # Numeric (or other non-symbolic) root/leaf: the Postwalk engine never
        # annotates these, and `analyze` reads them through the non-symbolic
        # `getsign`/`getcurvature` methods (`Positive`/`Affine` for a `Real`).
        # `node_sign` would wrongly return `AnySign` here.
        s = getsign(ex)
        c = getcurvature(ex)
        g = isnothing(ctx.M) ? nothing : getgcurvature(ex)
        props = NodeAnalysis(s, c, g)
        ctx.memo[ex] = props
        return props
    end

    # Match the historical pass order: signs (manifold-conditional when `M`
    # is set), then Euclidean curvature, then geodesic curvature. Callable
    # sign/curvature placeholders (e.g. `perspective`) are stored as-is so a
    # later `::Sign` assertion fails exactly as `getmetadata(...)::Sign` did.
    s = node_sign(ex, ctx.M)
    c = node_curvature(ex)
    g = isnothing(ctx.M) ? nothing : node_gcurvature(ex, ctx.M)
    props = NodeAnalysis(s, c, g)
    ctx.memo[ex] = props
    return props
end

"""
    analyze_fold(ex, M = nothing) -> NodeAnalysis

Read-only memoized fold of `ex`. `ex` must already be unwrapped and canonized
when called from [`analyze`](@ref).
"""
function analyze_fold(ex, M::Union{AbstractManifold, Nothing} = nothing)
    ctx = FoldCtx(IdDict{Any, NodeAnalysis}(), M)
    return with_fold_ctx(ctx) do
        analyze_node!(ex, ctx)
    end
end
