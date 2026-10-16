using Manifolds
using Symbolics: @register_symbolic
using LinearAlgebra

# @enum GSign GPositive GNegative GAnySign
"""
    GCurvature

Geodesic-curvature classification propagated by the gDCP analyzer.

# Values

- `GConvex`: geodesically convex.
- `GConcave`: geodesically concave.
- `GLinear`: geodesically affine.
- `GUnknownCurvature`: no supported geodesic-curvature rule applies.

These values are accepted by [`add_gdcprule`](@ref) and are returned in
`AnalysisResult.gcurvature` when manifold analysis is requested.
"""
@enum GCurvature GConvex GConcave GLinear GUnknownCurvature

"""
    GMonotonicity

Monotonicity classification for an argument of a geodesic DCP atom.

# Values

- `GIncreasing`: the atom is nondecreasing along the manifold argument.
- `GDecreasing`: the atom is nonincreasing along the manifold argument.
- `GAnyMono`: the atom's monotonicity is unrestricted for the rule.

These values are accepted by [`add_gdcprule`](@ref).
"""
@enum GMonotonicity GIncreasing GDecreasing GAnyMono

"""
    GDCPRule

Immutable descriptor for a registered geodesic DCP atom rule. Fields are
concretely typed so geodesic rule-table lookups do not box through `Any`.
"""
struct GDCPRule
    manifold::Any
    sign::Sign
    gcurvature::GCurvature
    gmonotonicity::Any
end

const gdcprules_dict = IdDict{Any, GDCPRule}()

"""
    add_gdcprule(f, manifold, sign, curvature, monotonicity)

Register a geodesic disciplined-convex-programming (gDCP) rule for the
symbolic atom `f` on a manifold type. This is the extension point for packages
that add atoms with known geodesic curvature.

Register the symbolic operation with `Symbolics.@register_symbolic` before
registering its rule. The function object passed here must be the operation
stored in the symbolic term.

# Arguments

- `f`: function used as the operation of the symbolic atom.
- `manifold`: manifold type on which the rule applies, such as
  `Manifolds.SymmetricPositiveDefinite` or `Manifolds.Lorentz`.
- `sign::Sign`: sign guaranteed by the atom.
- `curvature::GCurvature`: geodesic curvature of the atom.
- `monotonicity::GMonotonicity`: geodesic monotonicity of the atom. A tuple may
  be supplied for a multi-argument atom.

# Returns

The registered geodesic-rule descriptor.

# Examples

```jldoctest
julia> using SymbolicAnalysis, Symbolics, Manifolds

julia> test_gdcp_atom(p) = p[1];

julia> @register_symbolic test_gdcp_atom(p::Vector{Num})

julia> SymbolicAnalysis.add_gdcprule(
           test_gdcp_atom, Manifolds.Lorentz, SymbolicAnalysis.Positive,
           SymbolicAnalysis.GConvex, SymbolicAnalysis.GIncreasing
       );

julia> @variables p[1:3];

julia> result = analyze(test_gdcp_atom(p), Manifolds.Lorentz(2));

julia> result.gcurvature == SymbolicAnalysis.GConvex
true
```

This registration hook is public for packages extending SymbolicAnalysis. The
rule dictionary and `gdcprule` lookup function are implementation details and
should not be accessed directly.
"""
function add_gdcprule(f, manifold, sign, curvature, monotonicity)
    if !(monotonicity isa Tuple)
        monotonicity = (monotonicity,)
    end
    return gdcprules_dict[f] = makegrule(manifold, sign, curvature, monotonicity)
end
function makegrule(manifold, sign, curvature, monotonicity)
    return GDCPRule(manifold, sign, curvature, monotonicity)
end

"""
    constant_factors_ok(f, args)

Whether the constant factors of the gDCP atom call `f(args...)` meet the
preconditions under which the atom's registered rule holds. Atoms without such
preconditions accept any arguments; `hadamard_product` and `affine_map` add
methods.
"""
constant_factors_ok(f, args) = true

# Per-analysis memo of constant-factor checks. `valid` maps a call node to its
# `valid_constant_factors` result; `semidefinite` maps `(matrix, strict)` to the
# exact (semi)definiteness proof, keyed on the constant matrix's identity so
# trees rebuilt by the `propagate_*` walks still share it.
struct FactorCache
    valid::IdDict{Any, Bool}
    semidefinite::IdDict{Any, Bool}
end
FactorCache() = FactorCache(IdDict{Any, Bool}(), IdDict{Any, Bool}())

const FACTOR_CACHE_KEY = :SymbolicAnalysis_factor_cache

current_factor_cache() = get(task_local_storage(), FACTOR_CACHE_KEY, nothing)

function with_factor_cache(f)
    current_factor_cache() === nothing || return f()
    tls = task_local_storage()
    tls[FACTOR_CACHE_KEY] = FactorCache()
    try
        return f()
    finally
        delete!(tls, FACTOR_CACHE_KEY)
    end
end

"""
    valid_constant_factors(ex)

`false` when some call in `ex` (including `ex` itself) fails
`constant_factors_ok`. Such a node may not be a point of the manifold,
so no geodesic sign or curvature is derived for anything that contains it.
Within an analysis each node's result is computed once and reused by its
ancestors.
"""
function valid_constant_factors(ex, cache = current_factor_cache())
    iscall(ex) || return true
    if cache !== nothing
        hit = get(cache.valid, ex, nothing)
        hit === nothing || return hit
    end
    args = arguments(ex)
    ok = constant_factors_ok(operation(ex), args) &&
        all(a -> valid_constant_factors(a, cache), args)
    cache === nothing || (cache.valid[ex] = ok)
    return ok
end

hasgdcprule(f::Function) = haskey(gdcprules_dict, f)
hasgdcprule(f) = false
gdcprule(f, args...) = gdcprules_dict[f], args

setgcurvature(ex::Union{Symbolic, Num}, curv) = setmetadata(ex, GCurvature, curv)
setgcurvature(ex, curv) = ex
function getgcurvature(ex::Union{Symbolic, Num})
    v = unwrap(ex)
    hit = fold_lookup(v)
    if hit !== nothing
        g = hit.gcurvature
        return g === nothing ? GUnknownCurvature : g
    end
    if hasmetadata(ex, GCurvature)
        return getmetadata(ex, GCurvature)
    end
    return GUnknownCurvature
end
getgcurvature(ex) = GLinear
function hasgcurvature(ex::Union{Symbolic, Num})
    v = unwrap(ex)
    hit = fold_lookup(v)
    if hit !== nothing
        return hit.gcurvature !== nothing
    end
    return hasmetadata(ex, GCurvature)
end
hasgcurvature(ex) = ex isa Real

function mul_gcurvature(args)
    # Avoid allocations by not using findall
    non_constant_expr = nothing
    non_constant_count = 0
    constant_prod = one(Float64)
    for arg in args
        if issym(arg) || iscall(arg)
            non_constant_count += 1
            non_constant_expr = arg
            if non_constant_count > 1
                @warn "DGCP does not support multiple non-constant arguments in multiplication"
                return GUnknownCurvature
            end
        else
            constant_prod *= constval(arg)
        end
    end
    if non_constant_expr !== nothing
        curv = find_gcurvature(non_constant_expr)
        return if constant_prod < 0
            # flip
            curv == GConvex ? GConcave : curv == GConcave ? GConvex : curv
        else
            curv
        end
    end
    return GLinear
end

function add_gcurvature(args)
    # Avoid allocating intermediate arrays - check curvatures in one pass
    has_gconvex = false
    has_gconcave = false
    for arg in args
        curv = find_gcurvature(arg)
        if curv == GLinear
            continue
        elseif curv == GConvex
            has_gconvex = true
            if has_gconcave
                return GUnknownCurvature
            end
        elseif curv == GConcave
            has_gconcave = true
            if has_gconvex
                return GUnknownCurvature
            end
        else
            return GUnknownCurvature
        end
    end
    if has_gconvex
        return GConvex
    elseif has_gconcave
        return GConcave
    else
        return GLinear
    end
end

# A matrix argument's `Positive` sign and geodesic curvature on the SPD cone are
# in the Loewner order, so a Euclidean rule over a matrix argument keeps its
# monotonicity only when the atom is monotone in that order. `tr`, a full `sum`,
# `diag` and `quad_form`'s matrix slot are PSD-weighted linear functionals.
# On the PSD cone `|A[i,j]| ≤ √(A[i,i]·A[j,j])`, so `maximum` and `norm(A, Inf)`
# read only the diagonal; for an even integer `p` the gradient of `Σ A[i,j]^p` is
# a Hadamard power of `A`, PSD by the Schur product theorem. The entrywise 1-norm,
# `minimum`, `sum_largest`, column sums and nonlinear broadcasts are not monotone
# there (see test/dgp.jl).
const LOEWNER_LINEAR_READERS = (
    LinearAlgebra.tr, sum, LinearAlgebra.diag, adjoint, transpose, +, -, quad_form,
)

function loewner_monotonicity(f, args, monotonicity)
    mats = filter(a -> SymbolicUtils.symtype(a) <: AbstractMatrix, collect(args))
    isempty(mats) && return monotonicity
    signs_are(s) = all(a -> getsign(a) == s, mats)
    f in LOEWNER_LINEAR_READERS && return monotonicity
    if f isa SymbolicUtils.Mapreducer{typeof(identity)}
        f.dims isa Colon || return AnyMono
        f.reduce === Base.add_sum && return monotonicity
        f.reduce === max && return signs_are(Positive) ? monotonicity : AnyMono
        f.reduce === min && return signs_are(Negative) ? monotonicity : AnyMono
        return AnyMono
    end
    f === maximum && return signs_are(Positive) ? monotonicity : AnyMono
    f === minimum && return signs_are(Negative) ? monotonicity : AnyMono
    if f === LinearAlgebra.norm
        length(args) == 1 && return monotonicity
        p = constval(args[2])
        p isa Number && (p == Inf || (isinteger(p) && iseven(p) && p > 0)) &&
            return monotonicity
        return AnyMono
    end
    if f === broadcast
        g, xs = constval(args[1]), args[2:end]
        g in (+, -) && return monotonicity
        g === (*) && count(a -> !(constval(a) isa Number), xs) == 1 && return monotonicity
    end
    return AnyMono
end

function find_gcurvature(ex)
    # Checked before any stored curvature: a constant factor is a mutable matrix,
    # so metadata from an earlier analysis can certify a factor that is now invalid.
    iscall(ex) && !valid_constant_factors(ex) && return GUnknownCurvature
    if hasgcurvature(ex)
        return getgcurvature(ex)
    end
    if iscall(ex)
        f, args = operation(ex), arguments(ex)
        knowngcurv = false
        f_curvature = GUnknownCurvature
        f_monotonicity = (GAnyMono,)

        if hasgdcprule(f) && !any(iscall.(args))
            rule, args = gdcprule(f, args...)
            f_curvature = rule.gcurvature
            f_monotonicity = rule.gmonotonicity
            knowngcurv = true
        elseif f == LinearAlgebra.logdet
            if operation(args[1]) == conjugation ||
                    operation(args[1]) == LinearAlgebra.diag ||
                    operation(args[1]) === (+) ||
                    operation(args[1]) == affine_map ||
                    operation(args[1]) == hadamard_product
                return GConvex
            end
        elseif f == log &&
                iscall(args[1]) &&
                (operation(args[1]) == LinearAlgebra.tr || operation(args[1]) == quad_form)
            return GConvex
        elseif (f == schatten_norm || f == eigsummax) && operation(args[1]) == log
            return GConvex
        elseif f == sum_log_eigmax && hasdcprule(args[1])
            if dcprule(operation(args[1])) == Convex
                return GConvex
            else
                return GUnknownCurvature
            end
        elseif f == affine_map
            if args[1] == tr || args[1] == conjugation || args[1] == diag
                return GConvex
            else
                return GUnknownCurvature
            end
        elseif hasgdcprule(f) && any(iscall.(args))
            for i in eachindex(args)
                if iscall(args[i])
                    if operation(args[i]) == inv
                        rule, args = gdcprule(f, args...)
                        f_curvature = rule.gcurvature
                        f_monotonicity = if rule.gmonotonicity == GIncreasing
                            GDecreasing
                        elseif rule.gmonotonicity == GDecreasing
                            GIncreasing
                        else
                            GAnyMono
                        end
                        knowngcurv = true
                    elseif operation(args[i]) == broadcast
                        rule, args = gdcprule(f, args...)
                        f_curvature = rule.gcurvature
                        f_monotonicity = rule.gmonotonicity
                        knowngcurv = true
                    elseif operation(args[i]) == affine_map
                        rule, args = gdcprule(f, args...)
                        f_curvature = rule.gcurvature
                        f_monotonicity = rule.gmonotonicity
                        knowngcurv = true
                    elseif operation(args[i]) === (+) &&
                            count(!isconstarg, arguments(args[i])) <= 1
                        # A constant shift moves the atom by a constant, which leaves
                        # its geodesic curvature alone: `tr(X + C) = tr(X) + tr(C)`.
                        # Without this the rule is skipped and the Euclidean fallback
                        # below was what rescued `tr(X + γI)`.
                        rule, args = gdcprule(f, args...)
                        f_curvature = rule.gcurvature
                        f_monotonicity = rule.gmonotonicity
                        knowngcurv = true
                    end
                end
            end
        elseif f === (*)
            a1 = constval(args[1])
            if a1 isa Number && a1 > 0
                return find_gcurvature(args[2])
            elseif a1 isa Number && a1 < 0
                argscurv = find_gcurvature(args[2])
                if argscurv == GConvex
                    return GConcave
                elseif argscurv == GConcave
                    return GConvex
                else
                    return argscurv
                end
            else
                @warn "Disciplined Programming does not support multiple non-constant arguments in multiplication"
                return GUnknownCurvature
            end
        end

        # A Euclidean rule may supply the *composition* over an argument that already
        # carries a geodesic curvature — `distance(M, A, X)^2` gets its shape from
        # `^`'s Euclidean convex-increasing rule over a GConvex inner, soundly. It may
        # NOT classify an atom on the manifold: where every argument is GLinear the
        # Euclidean curvature is the only input and it implies nothing geodesically.
        # That path certified `eigmin(X)` as GConcave on the SPD cone from its
        # Euclidean `Concave`, but between A = [7.7517 1.132; 1.132 8.8903] and
        # B = [2.8936 0.3831; 0.3831 0.7551] the geodesic midpoint gives 2.3890
        # against a chord of 3.8712 — below it, so concavity is refuted.
        if !knowngcurv
            (hasdcprule(f) && any(a -> find_gcurvature(a) in (GConvex, GConcave), args)) ||
                return GUnknownCurvature
            rule, args = dcprule(f, args...)
            f_curvature = rule.curvature
            f_monotonicity = loewner_monotonicity(f, args, rule.monotonicity)
        end

        if f_curvature == Convex || f_curvature == Affine
            if all(enumerate(args)) do (i, arg)
                    arg_curv = find_gcurvature(arg)
                    m = get_arg_property(f_monotonicity, i, args)
                    # @show arg
                    if arg_curv == GConvex
                        m == Increasing
                    elseif arg_curv == GConcave
                        m == Decreasing
                    else
                        arg_curv == GLinear
                    end
                end
                return GConvex
            else
                return GUnknownCurvature
            end
        elseif f_curvature == Concave
            if all(enumerate(args)) do (i, arg)
                    arg_curv = find_gcurvature(arg)
                    m = get_arg_property(f_monotonicity, i, args)
                    if arg_curv == GConcave
                        m == Increasing
                    elseif arg_curv == GConvex
                        m == Decreasing
                    else
                        arg_curv == GLinear
                    end
                end
                return GConcave
            else
                return GUnknownCurvature
            end
        elseif f_curvature isa GCurvature
            return f_curvature
        else
            return GUnknownCurvature
        end
    elseif hasfield(typeof(ex), :val) && haskey(gdcprules_dict, operation(ex.val))
        f, args = operation(ex.val), arguments(ex.val)
        rule, args = gdcprule(f, args...)
        return rule.gcurvature
    else
        return GLinear
    end
    return GUnknownCurvature
end

# See `node_curvature`: geodesic curvature of a single node whose children the
# bottom-up walk has already annotated.
function node_gcurvature(ex)
    if iscall(ex)
        f = operation(ex)
        if f === (*)
            return mul_gcurvature(arguments(ex))
        elseif f === (+)
            return add_gcurvature(arguments(ex))
        end
    end
    return find_gcurvature(ex)
end

function propagate_gcurvature(ex, M::AbstractManifold)
    # Operate on the raw symbolic: on Symbolics v7 walking a `Num`/`Arr` wrapper
    # round-trips through wrap/unwrap and loses the gcurvature metadata that the
    # final `getgcurvature` reads. `analyze` already unwraps; do the same here so
    # the function is correct when called directly on a wrapped expression.
    ex = SymbolicUtils.unwrap(ex)
    # A geodesic rule's sign holds for a point on `M`, so the geodesic pass
    # re-derives signs under `M` rather than inheriting the Euclidean ones a bare
    # `propagate_sign` leaves behind. `find_gcurvature` falls back to the Euclidean
    # rule table for atoms with no geodesic rule and picks up `increasing_if_positive`
    # with it, so `distance(M, A, X)^2` needs `distance`'s manifold sign to compose.
    return with_factor_cache() do
        ex = propagate_sign(ex, M)
        Postwalk(x -> issym(x) || iscall(x) ? setgcurvature(x, node_gcurvature(x)) : x)(ex)
    end
end
