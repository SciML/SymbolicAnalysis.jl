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

hasgdcprule(f::Function) = haskey(gdcprules_dict, f)
hasgdcprule(f) = false
gdcprule(f, args...) = gdcprules_dict[f], args

# A rule applies only on the manifold type it was registered for. `distance` is
# registered for both SPD and Lorentz with the same Positive/GConvex/GAnyMono
# descriptor; Lorentz's `add_gdcprule` overwrites the SPD dict entry, so accept
# either supported type for that atom.
function gdcp_rule_applies(f, M)
    M === nothing && return false
    hasgdcprule(f) || return false
    rule = gdcprules_dict[f]
    M isa rule.manifold && return true
    return f === Manifolds.distance &&
        (M isa SymmetricPositiveDefinite || M isa Lorentz)
end

# A geodesic rule describes its atom at a point of the manifold, so it applies
# only to an argument that provably is one. `isometric_point` follows the
# variable through maps that send SPD geodesics to geodesics; `manifold_point`
# also admits PSD shifts and sums, which stay on the cone but bend geodesics,
# so they carry the rule's sign but not its geodesic curvature.
function fits_manifold(ex, M)
    M === nothing && return true
    T = SymbolicUtils.symtype(ex)
    if M isa SymmetricPositiveDefinite
        T <: AbstractMatrix || return false
    elseif M isa Lorentz
        T <: AbstractVector || return false
    else
        return false
    end
    sh = SymbolicUtils.shape(ex)
    sh isa AbstractVector || return true
    got = Tuple(map(length, sh))
    want = Manifolds.representation_size(M)
    got == want && return true
    return M isa Lorentz && got == (want[1] + 1,)
end

positive_const(x) = (v = constval(x); v isa Real && v > 0)

function symmetric_const(x)
    isconstarg(x) || return false
    v = constval(x)
    v isa Real && return true
    v = map(constval, v)
    return v isa AbstractMatrix && ishermitian(v)
end

# Exact test: every finite float is a rational, so converting to
# `Rational{BigInt}` and running a pivoted LDLᵀ decides PSD without rounding.
function psd_const(x)
    isconstarg(x) || return false
    v = map(constval, constval(x))
    v isa AbstractMatrix && size(v, 1) == size(v, 2) || return false
    all(y -> y isa Union{Integer, Rational, AbstractFloat} && isfinite(y), v) || return false
    A = Matrix{Rational{BigInt}}(v)
    A == transpose(A) || return false
    while !isempty(A)
        d = diag(A)
        any(<(0), d) && return false
        k = argmax(d)
        p = A[k, k]
        # A PSD matrix with a zero diagonal entry has a zero row there.
        iszero(p) && return iszero(A)
        r = [1:(k - 1); (k + 1):size(A, 1)]
        a = A[r, k]
        A = A[r, r] - a * transpose(a) / p
    end
    return true
end

function nonsingular_const(x)
    isconstarg(x) || return false
    v = map(constval, constval(x))
    return v isa AbstractMatrix && size(v, 1) == size(v, 2) &&
        LinearAlgebra.rank(float.(v)) == size(v, 1)
end

function point_kind(ex, M = nothing)::Int
    issym(ex) && return fits_manifold(ex, M) ? 2 : 0
    iscall(ex) || return 0
    M === nothing || M isa SymmetricPositiveDefinite || return 0
    ctx = current_fold_ctx()
    points = (ctx !== nothing && ctx.M === M) ? ctx.points : nothing
    if points !== nothing
        hit = get(points, ex, -1)
        hit != -1 && return hit
    end
    k = if !atom_dims_ok(ex, M)
        0
    else
        f, args = operation(ex), arguments(ex)
        if f === inv || f === adjoint || f === transpose
            point_kind(args[1], M)
        elseif f === conjugation
            nonsingular_const(args[2]) ? point_kind(args[1], M) : 0
        elseif f === (*)
            if count(!isconstarg, args) == 1 &&
                    all(a -> !isconstarg(a) || positive_const(a), args)
                idx = findfirst(!isconstarg, args)
                idx === nothing ? 0 : point_kind(args[idx], M)
            else
                0
            end
        elseif f === (+)
            if any(!isconstarg, args) &&
                    all(a -> isconstarg(a) ? psd_const(a) : point_kind(a, M) >= 1, args)
                1
            else
                0
            end
        else
            0
        end
    end
    if points !== nothing
        points[ex] = k
    end
    return k
end

manifold_point(ex, M = nothing; shifts::Bool = true) = point_kind(ex, M) >= (shifts ? 1 : 2)

isometric_point(ex, M = nothing) = manifold_point(ex, M; shifts = false)

spd_point(ex, M) = (M === nothing || M isa SymmetricPositiveDefinite) && isometric_point(ex, M)

# The SPD special cases treat every argument but the manifold point as a fixed
# parameter; one that depends on the point (`quad_form([sin(X[1,1]), 0], X)`)
# voids the rule.
consts_from(args, i) = all(isconstarg, args[i:end])

# Fixed data beside the point: a constant, or a free symbol of another kind than
# the points of `M` (the vector `y` of `log_quad_form(y, X)`). A symbol of the
# point's kind but the wrong shape (3×3 `Y` on SPD(2)) is neither.
function fixed_param(a, M)
    isconstarg(a) && return true
    issym(a) || return false
    T = SymbolicUtils.symtype(a)
    M isa SymmetricPositiveDefinite && return !(T <: AbstractMatrix)
    M isa Lorentz && return !(T <: AbstractVector)
    return M === nothing
end

function bare_arg_ok(a, M, role = nothing)
    if role === nothing
        return !issym(a) || fits_manifold(a, M) || fixed_param(a, M)
    elseif role === :spd_point || role === :lorentz_point
        return issym(a) && fits_manifold(a, M)
    elseif role === :spd_or_mat || role === :lorentz_or_vec
        return isconstarg(a) || (issym(a) && fits_manifold(a, M))
    elseif role === :mat || role === :cong_mat || role === :aff_mat ||
            role === :vec || role === :rect_mat || role === :rect_vec || role === :scalar
        return fixed_param(a, M)
    elseif role === :vec_or_vecs
        return fixed_param(a, M) || vector_collection_fits(a, nothing)
    elseif role === :cong_mat_or_mats
        return fixed_param(a, M) || matrix_collection_fits(a, nothing)
    elseif role === :fn
        return constval(a) isa Function
    elseif role === :manifold
        m = constval(a)
        return m isa AbstractManifold && (M === nothing || m == M)
    end
    return false
end

function sym_shape(a)
    if isconstarg(a)
        v = constval(a)
        return v isa AbstractArray ? size(v) : ()
    end
    issym(a) || iscall(a) || return nothing
    sh = SymbolicUtils.shape(a)
    sh isa AbstractVector && return Tuple(map(length, sh))
    if hasproperty(sh, :ndims) && sh.ndims >= 0
        return ntuple(_ -> nothing, sh.ndims)
    end
    T = SymbolicUtils.symtype(a)
    if T isa Type
        T <: Number && return ()
        if T <: AbstractArray && applicable(ndims, T)
            return ntuple(_ -> nothing, ndims(T))
        end
    end
    return nothing
end

function merge_shapes(s1, s2)
    (s1 === false || s2 === false) && return false
    s1 === nothing && return s2
    s2 === nothing && return s1
    length(s1) == length(s2) || return false
    for (a, b) in zip(s1, s2)
        (a !== nothing && b !== nothing && a != b) && return false
    end
    return map((a, b) -> a === nothing ? b : a, s1, s2)
end

function merge_bcast_shapes(s1, s2)
    (s1 === false || s2 === false) && return false
    s1 === nothing && return s2
    s2 === nothing && return s1
    s1 == () && return s2
    s2 == () && return s1
    nd = max(length(s1), length(s2))
    out = Vector{Union{Int, Nothing}}(undef, nd)
    for i in 1:nd
        a = i <= length(s1) ? s1[i] : 1
        b = i <= length(s2) ? s2[i] : 1
        if a !== nothing && b !== nothing && a != b && a != 1 && b != 1
            return false
        end
        out[i] = a === nothing ? (b == 1 ? nothing : b) :
            b === nothing ? (a == 1 ? nothing : a) : max(a, b)
    end
    return Tuple(out)
end

# Known rank and axis extents of `a` as a tuple of `Union{Int, Nothing}` whose
# length is the rank, or `nothing` when the rank itself is unknown.
function known_size(a, M = nothing)
    s = validated_shape(a, M)
    return s === false ? nothing : s
end

# Whether `a` has shape `want`, where a `nothing` entry in either `s` or `want`
# matches any length along that axis. The rank (`length(s)`) must match whenever
# it is known.
function shape_fits_val(s, want)
    s === false && return false
    s === nothing && return true
    return length(s) == length(want) &&
        all(((x, w),) -> x === nothing || w === nothing || x == w, zip(s, want))
end

shape_fits(a, want, M = nothing) = shape_fits_val(validated_shape(a, M), want)

# A constant vector of numeric vectors, each of length `n` (any length if `nothing`).
function vector_collection_fits(a, n)
    (issym(a) || iscall(a) || isconstarg(a)) && return false
    v = constval(a)
    v isa AbstractVector && !isempty(v) || return false
    return all(v) do y
        y isa AbstractVector && all(e -> constval(e) isa Number, y) &&
            (n === nothing || length(y) == n)
    end
end

function matrix_collection_fits(a, rdim, cdim = rdim)
    (issym(a) || iscall(a) || isconstarg(a)) && return false
    v = constval(a)
    v isa AbstractVector && !isempty(v) || return false
    m = cdim === nothing && first(v) isa AbstractMatrix ? size(first(v), 2) : cdim
    return all(v) do Y
        Y isa AbstractMatrix && all(e -> constval(e) isa Number, Y) &&
            (rdim === nothing || size(Y, 1) == rdim) &&
            (m === nothing || size(Y, 2) == m)
    end
end

# Per-atom argument roles: `:spd_point` / `:lorentz_point` mark the manifold
# variable slot, `:spd_or_mat` / `:lorentz_or_vec` accept a point or a constant
# array, and the remaining tags describe fixed parameter slots relative to the
# point dimension.
function atom_roles(f, args)
    if f === Manifolds.distance
        length(args) == 3 || return ()
        m = constval(args[1])
        m isa SymmetricPositiveDefinite && return (:manifold, :spd_or_mat, :spd_or_mat)
        m isa Lorentz && return (:manifold, :lorentz_or_vec, :lorentz_or_vec)
        return ()
    elseif f === affine_map
        isempty(args) && return ()
        op = constval(args[1])
        if op === conjugation && length(args) == 4
            return (:fn, :spd_point, :aff_mat, :cong_mat_or_mats)
        elseif op === hadamard_product && length(args) == 4
            return (:fn, :spd_point, :mat, :mat)
        end
        return ()
    end
    table = (
        LinearAlgebra.tr => ((:spd_point,),),
        sum => ((:spd_point,),),
        LinearAlgebra.diag => ((:spd_point,),),
        LinearAlgebra.eigmax => ((:spd_point,),),
        LinearAlgebra.logdet => ((:spd_point,),),
        inv => ((:spd_point,),),
        adjoint => ((:spd_point,),),
        eigsummax => ((:spd_point, :scalar),),
        schatten_norm => ((:spd_point,), (:spd_point, :scalar)),
        scalar_mat => ((:spd_point,), (:spd_point, :scalar)),
        sum_log_eigmax => ((:spd_point, :scalar), (:fn, :spd_point, :scalar)),
        sdivergence => ((:spd_or_mat, :spd_or_mat),),
        conjugation => ((:spd_point, :cong_mat),),
        hadamard_product => ((:spd_point, :mat),),
        quad_form => ((:vec, :spd_point),),
        log_quad_form => ((:vec_or_vecs, :spd_point),),
        lorentz_log_barrier => ((:lorentz_point,),),
        lorentz_homogeneous_quadratic => ((:mat, :lorentz_point),),
        lorentz_homogeneous_diagonal => ((:vec, :lorentz_point),),
        lorentz_nonhomogeneous_quadratic => ((:mat, :vec, :scalar, :lorentz_point),),
        lorentz_least_squares => ((:rect_mat, :rect_vec, :lorentz_point),),
        lorentz_transform => ((:mat, :lorentz_point),),
    )
    for (g, sigs) in table
        if f === g
            for r in sigs
                length(r) == length(args) && return r
            end
            return ()
        end
    end
    return nothing
end

shape_lead_dim(s) = (s === nothing || s === false || isempty(s)) ? nothing : s[1]
shape_col_dim(s) = (s === nothing || s === false || length(s) < 2) ? nothing : s[2]
lead_dim(a, M = nothing) = shape_lead_dim(validated_shape(a, M))
col_dim(a, M = nothing) = shape_col_dim(validated_shape(a, M))

function cong_out_dim(args, s_args)
    length(args) >= 4 || return nothing
    m = shape_lead_dim(s_args[3])
    m !== nothing && return m
    y = args[4]
    if issym(y) || iscall(y) || isconstarg(y)
        return shape_col_dim(s_args[4])
    end
    v = constval(y)
    if v isa AbstractVector && !isempty(v) && first(v) isa AbstractMatrix
        return size(first(v), 2)
    end
    return nothing
end

function validated_shape(a, M = nothing)
    if !iscall(a)
        isconstarg(a) && return sym_shape(a)
        issym(a) || return nothing
        s = sym_shape(a)
        if M isa SymmetricPositiveDefinite && s !== nothing && length(s) == 2
            fits_manifold(a, M) || return false
            return Manifolds.representation_size(M)
        elseif M isa Lorentz && s !== nothing && length(s) == 1
            fits_manifold(a, M) || return false
        end
        return s
    end
    ctx = current_fold_ctx()
    shapes = (ctx !== nothing && ctx.M === M) ? ctx.shapes : nothing
    if shapes !== nothing
        hit = get(shapes, a, :unvisited)
        hit !== :unvisited && return hit
    end
    f, args = operation(a), arguments(a)
    sh = if f === SymbolicUtils.array_literal
        sym_shape(a)
    else
        cs = validated_call_shape(f, args, M)
        cs === :fallback ? sym_shape(a) : cs
    end
    if shapes !== nothing
        shapes[a] = sh
    end
    return sh
end

function validated_call_shape(f, args, M = nothing)
    roles = atom_roles(f, args)
    if roles !== nothing
        isempty(roles) && return false
        eff_M = M === nothing && f === Manifolds.distance && !isempty(args) &&
            constval(args[1]) isa AbstractManifold ? constval(args[1]) : M
        has_spd = any(r -> r === :spd_point || r === :spd_or_mat, roles)
        has_lorentz = any(r -> r === :lorentz_point || r === :lorentz_or_vec, roles)
        if has_spd && eff_M !== nothing && !(eff_M isa SymmetricPositiveDefinite)
            return false
        elseif has_lorentz && eff_M !== nothing && !(eff_M isa Lorentz)
            return false
        end
        nargs = length(args)
        s1 = nargs >= 1 ? validated_shape(args[1], eff_M) : nothing
        s1 === false && return false
        s2 = nargs >= 2 ? validated_shape(args[2], eff_M) : nothing
        s2 === false && return false
        s3 = nargs >= 3 ? validated_shape(args[3], eff_M) : nothing
        s3 === false && return false
        s4 = nargs >= 4 ? validated_shape(args[4], eff_M) : nothing
        s4 === false && return false
        s_args = (s1, s2, s3, s4)
        d = point_dim(roles, args, s_args, eff_M)
        d === -1 && return false
        for i in eachindex(roles)
            slot_dims_ok(roles[i], args[i], s_args[i], d, args, s_args, eff_M) || return false
        end
        if f === conjugation
            return s2 !== nothing && length(s2) == 2 ? (s2[2], s2[2]) : (nothing, nothing)
        elseif f === affine_map
            op = constval(args[1])
            if op === conjugation
                m = cong_out_dim(args, s_args)
                return (m, m)
            elseif op === hadamard_product
                return merge_shapes(s2, s3)
            end
        elseif f === hadamard_product
            return merge_shapes(s1, s2)
        elseif f === scalar_mat
            length(args) == 1 && return s1
            k = constval(args[2])
            return k isa Integer ? (Int(k), Int(k)) : (nothing, nothing)
        elseif f === inv || f === adjoint
            return s1
        elseif f === LinearAlgebra.diag
            return s1 !== nothing && length(s1) == 2 ? (s1[1],) : (nothing,)
        elseif f === lorentz_transform
            d1 = s1 !== nothing && !isempty(s1) ? s1[1] : nothing
            return merge_shapes((d1,), s2)
        end
        return ()
    end
    if (f === (+) || f === (-)) && !isempty(args)
        acc = validated_shape(args[1], M)
        acc === false && return false
        for x in @view args[2:end]
            acc = merge_shapes(acc, validated_shape(x, M))
            acc === false && return false
        end
        return acc
    elseif f === (*) && !isempty(args)
        nonscalar = nothing
        many_nonscalar = false
        for x in args
            sx = validated_shape(x, M)
            sx === false && return false
            if sx != ()
                if nonscalar === nothing
                    nonscalar = sx
                else
                    many_nonscalar = true
                end
            end
        end
        many_nonscalar || return nonscalar === nothing ? () : nonscalar
        return :fallback
    elseif f === transpose && length(args) == 1
        s1 = validated_shape(args[1], M)
        s1 === false && return false
        return s1 !== nothing && length(s1) == 2 ? (s1[2], s1[1]) : s1
    elseif (f === log || f === exp || f === sqrt || f === abs || f === abs2) &&
            length(args) == 1
        return validated_shape(args[1], M)
    elseif f === broadcast && length(args) >= 2
        acc = validated_shape(args[2], M)
        acc === false && return false
        for x in @view args[3:end]
            acc = merge_bcast_shapes(acc, validated_shape(x, M))
            acc === false && return false
        end
        return acc
    end
    for x in args
        validated_shape(x, M) === false && return false
    end
    return :fallback
end

function point_dim(roles, args, s_args, M = nothing)
    has_spd = any(r -> r === :spd_point || r === :spd_or_mat, roles)
    has_lorentz = any(r -> r === :lorentz_point || r === :lorentz_or_vec, roles)
    d = nothing
    for i in eachindex(roles)
        r, a, s = roles[i], args[i], s_args[i]
        if r === :spd_point || r === :lorentz_point ||
                ((r === :spd_or_mat || r === :lorentz_or_vec) && !isconstarg(a))
            if r === :spd_point || r === :spd_or_mat
                s !== nothing && length(s) != 2 && return -1
                if s !== nothing
                    (s[1] !== nothing && s[2] !== nothing && s[1] != s[2]) && return -1
                    n = s[1] !== nothing ? s[1] : s[2]
                    if n !== nothing
                        d = n
                        break
                    end
                end
            else
                s !== nothing && length(s) != 1 && return -1
                if s !== nothing && s[1] !== nothing
                    d = s[1]
                    break
                end
            end
        end
    end
    if d === nothing && has_spd && M isa SymmetricPositiveDefinite
        d = Manifolds.representation_size(M)[1]
    end
    if d === nothing
        for i in eachindex(roles)
            r, a, s = roles[i], args[i], s_args[i]
            if r === :spd_or_mat || r === :lorentz_or_vec || r === :mat ||
                    r === :cong_mat || r === :vec
                n = shape_lead_dim(s)
            elseif r === :rect_mat
                n = shape_col_dim(s)
            elseif r === :vec_or_vecs
                n = (issym(a) || iscall(a) || isconstarg(a)) ? shape_lead_dim(s) :
                    (v = constval(a); v isa AbstractVector && !isempty(v) && first(v) isa AbstractVector ? length(first(v)) : nothing)
            elseif r === :cong_mat_or_mats
                n = (issym(a) || iscall(a) || isconstarg(a)) ? shape_lead_dim(s) :
                    (v = constval(a); v isa AbstractVector && !isempty(v) && first(v) isa AbstractMatrix ? size(first(v), 1) : nothing)
            else
                n = nothing
            end
            if n !== nothing
                d = n
                break
            end
        end
    end
    if has_lorentz && M isa Lorentz && d !== nothing
        want = Manifolds.representation_size(M)[1]
        (d == want || d == want + 1) || return -1
    end
    return d
end

function slot_dims_ok(role::Symbol, a, s, d, args, s_args, M = nothing)
    if role === :spd_point || role === :spd_or_mat || role === :mat
        return shape_fits_val(s, (d, d))
    elseif role === :cong_mat
        return shape_fits_val(s, (d, nothing))
    elseif role === :aff_mat
        m = cong_out_dim(args, s_args)
        return shape_fits_val(s, (m, m))
    elseif role === :cong_mat_or_mats
        m = cong_out_dim(args, s_args)
        (issym(a) || iscall(a) || isconstarg(a)) && return shape_fits_val(s, (d, m))
        return matrix_collection_fits(a, d, m)
    elseif role === :lorentz_point || role === :lorentz_or_vec || role === :vec
        return shape_fits_val(s, (d,))
    elseif role === :vec_or_vecs
        (issym(a) || iscall(a) || isconstarg(a)) && return shape_fits_val(s, (d,))
        return vector_collection_fits(a, d)
    elseif role === :rect_mat
        m = shape_lead_dim(s)
        return shape_fits_val(s, (m, d))
    elseif role === :rect_vec
        m = shape_lead_dim(s_args[1])
        return shape_fits_val(s, (m,))
    elseif role === :scalar
        return (issym(a) || iscall(a) || isconstarg(a)) && shape_fits_val(s, ())
    elseif role === :fn
        return constval(a) isa Function
    elseif role === :manifold
        m = constval(a)
        m isa AbstractManifold && (M === nothing || m == M) || return false
        d === nothing && return true
        want = Manifolds.representation_size(m)
        return m isa Lorentz ? (d == want[1] || d == want[1] + 1) : want == (d, d)
    end
    return false
end

# Whether each argument's known rank and dimensions fit its role in the atom,
# relative to the manifold `M` and the point argument.
atom_dims_ok(ex, M::Union{AbstractManifold, Nothing} = nothing) =
    !iscall(ex) || validated_shape(ex, M) !== false
atom_dims_ok(f, args, M = nothing) = validated_call_shape(f, args, M) !== false

function gdcp_arg_ok(f, a, M, role = nothing; allow_call::Bool = false)
    if iscall(a)
        allow_call || return false
        (
            role === nothing || role === :spd_point || role === :lorentz_point ||
                role === :spd_or_mat || role === :lorentz_or_vec
        ) || return false
        return geodesic_rule_arg(f, a, M)
    end
    return bare_arg_ok(a, M, role)
end

function gdcp_args_ok(f, args, M; allow_call::Bool = false)
    roles = atom_roles(f, args)
    roles === nothing && return all(a -> gdcp_arg_ok(f, a, M; allow_call), args)
    isempty(roles) && return false
    has_point = false
    for i in eachindex(roles)
        r, a = roles[i], args[i]
        gdcp_arg_ok(f, a, M, r; allow_call) || return false
        if (r === :spd_point || r === :lorentz_point || r === :spd_or_mat || r === :lorentz_or_vec) &&
                !isconstarg(a)
            has_point = true
        end
    end
    return has_point
end

gdcp_args_on_manifold(args, M) =
    all(a -> !(issym(a) || iscall(a)) || manifold_point(a, M), args)

function gdcp_args_on_manifold(f, args, M)
    atom_dims_ok(f, args, M) || return false
    roles = atom_roles(f, args)
    roles === nothing && return gdcp_args_on_manifold(args, M)
    isempty(roles) && return false
    has_point = false
    for i in eachindex(roles)
        r, a = roles[i], args[i]
        if r === :spd_point || r === :lorentz_point
            manifold_point(a, M) || return false
            has_point = true
        elseif r === :spd_or_mat || r === :lorentz_or_vec
            if !isconstarg(a)
                manifold_point(a, M) || return false
                has_point = true
            end
        elseif r === :manifold
            m = constval(a)
            (m isa AbstractManifold && (M === nothing || m == M)) || return false
        else
            (issym(a) || iscall(a)) && return false
        end
    end
    return has_point
end

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

function mul_gcurvature(args, M = nothing)
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
        curv = find_gcurvature(non_constant_expr, M)
        return if constant_prod < 0
            # flip
            curv == GConvex ? GConcave : curv == GConcave ? GConvex : curv
        else
            curv
        end
    end
    return GLinear
end

function add_gcurvature(args, M = nothing)
    # Avoid allocating intermediate arrays - check curvatures in one pass
    has_gconvex = false
    has_gconcave = false
    for arg in args
        curv = find_gcurvature(arg, M)
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

# Atoms that are a supremum of positive linear functionals of `X` (plus a
# constant), so `f(C + X)` and `f(B + Φ(X))` stay geodesically convex for any
# constant `C`, PSD `B` and positive linear `Φ`. `inv` and `logdet` are not:
# `tr(inv(X + I))` is `1/(1 + eᵗ)` along a scalar geodesic, which is not convex.
const SHIFT_INVARIANT_GATOMS = (LinearAlgebra.tr, sum, LinearAlgebra.diag, eigmax, eigsummax)

function positive_image(ex, M = nothing)
    spd_point(ex, M) && return true
    iscall(ex) || return false
    atom_dims_ok(ex, M) || return false
    g, gargs = operation(ex), arguments(ex)
    if g === conjugation || g === LinearAlgebra.diag || g === hadamard_product
        return spd_point(gargs[1], M) && consts_from(gargs, 2)
    elseif g === affine_map
        return spd_point(gargs[2], M) && consts_from(gargs, 3)
    end
    return false
end

# `logdet(B + Φ(X))` is geodesically convex for PSD `B` and positive linear `Φ`.
# Only these positive images qualify; isometric constructors (`inv`, `adjoint`,
# `transpose`, `*`) are absent so bare `logdet(X)` stays `GLinear` via its table rule.
function logdet_gconvex_arg(a, M = nothing)
    iscall(a) || return false
    atom_dims_ok(a, M) || return false
    g, gargs = operation(a), arguments(a)
    if g === (+)
        return count(!isconstarg, gargs) == 1 &&
            all(x -> isconstarg(x) ? psd_const(x) : positive_image(x, M), gargs)
    end
    return (
        g === conjugation || g === LinearAlgebra.diag || g === affine_map ||
            g === hadamard_product
    ) && positive_image(a, M)
end

function geodesic_rule_arg(f, a, M = nothing)
    iscall(a) && !atom_dims_ok(a, M) && return false
    isometric_point(a, M) && return true
    g = operation(a)
    any(h -> h === f, SHIFT_INVARIANT_GATOMS) || return false
    eig_shift = f === eigmax || f === eigsummax
    if g === affine_map
        return isometric_point(arguments(a)[2], M) && consts_from(arguments(a), 3)
    elseif g === (+)
        args = arguments(a)
        nonconst = filter(!isconstarg, args)
        return length(nonconst) == 1 && isometric_point(only(nonconst), M) &&
            (!eig_shift || all(x -> !isconstarg(x) || symmetric_const(x), args))
    elseif g === broadcast
        bargs = arguments(a)
        op = constval(bargs[1])
        return length(bargs) == 3 && (op === (+) || op === (-)) &&
            isometric_point(bargs[2], M) && isconstarg(bargs[3]) &&
            (!eig_shift || symmetric_const(bargs[3]))
    end
    return false
end

function find_gcurvature(ex, M = nothing)
    if hasgcurvature(ex)
        return getgcurvature(ex)
    end
    if iscall(ex)
        atom_dims_ok(ex, M) || return GUnknownCurvature
        f, args = operation(ex), arguments(ex)
        knowngcurv = false
        f_curvature = GUnknownCurvature
        f_monotonicity = (GAnyMono,)

        if gdcp_rule_applies(f, M) && !any(iscall, args) && gdcp_args_ok(f, args, M)
            rule, args = gdcprule(f, args...)
            f_curvature = rule.gcurvature
            f_monotonicity = rule.gmonotonicity
            knowngcurv = true
        elseif f == LinearAlgebra.logdet
            if logdet_gconvex_arg(args[1], M)
                return GConvex
            end
        elseif f == log &&
                iscall(args[1]) &&
                (
                (
                    operation(args[1]) == LinearAlgebra.tr &&
                        atom_dims_ok(args[1], M) &&
                        positive_image(arguments(args[1])[1], M)
                ) ||
                    (
                    operation(args[1]) == quad_form && fixed_param(arguments(args[1])[1], M) &&
                        atom_dims_ok(args[1], M) &&
                        positive_image(arguments(args[1])[2], M)
                )
            )
            return GConvex
        elseif (f == schatten_norm || f == eigsummax) && iscall(args[1]) &&
                operation(args[1]) == log && spd_point(arguments(args[1])[1], M) &&
                consts_from(args, 2)
            return GConvex
        elseif f == sum_log_eigmax && hasdcprule(args[1])
            if dcprule(operation(args[1])) == Convex && spd_point(args[2], M) && consts_from(args, 3)
                return GConvex
            else
                return GUnknownCurvature
            end
        elseif f == affine_map
            if (args[1] == conjugation || args[1] == hadamard_product) &&
                    spd_point(args[2], M) && consts_from(args, 3)
                return GConvex
            else
                return GUnknownCurvature
            end
        elseif gdcp_rule_applies(f, M) && any(iscall, args) &&
                gdcp_args_ok(f, args, M; allow_call = true)
            rule, args = gdcprule(f, args...)
            f_curvature = rule.gcurvature
            f_monotonicity = rule.gmonotonicity
            knowngcurv = true
        elseif f === (*)
            a1 = constval(args[1])
            if a1 isa Number && a1 > 0
                return find_gcurvature(args[2], M)
            elseif a1 isa Number && a1 < 0
                argscurv = find_gcurvature(args[2], M)
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
            (hasdcprule(f) && any(a -> find_gcurvature(a, M) in (GConvex, GConcave), args)) ||
                return GUnknownCurvature
            rule, args = dcprule(f, args...)
            f_curvature = rule.curvature
            f_monotonicity = rule.monotonicity
        end

        if f_curvature == Convex || f_curvature == Affine
            if all(enumerate(args)) do (i, arg)
                    arg_curv = find_gcurvature(arg, M)
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
                    arg_curv = find_gcurvature(arg, M)
                    m = f_monotonicity[i]
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
        gdcp_rule_applies(f, M) || return GLinear
        rule, args = gdcprule(f, args...)
        return rule.gcurvature
    else
        return GLinear
    end
    return GUnknownCurvature
end

# See `node_curvature`: geodesic curvature of a single node whose children the
# bottom-up walk has already annotated.
function node_gcurvature(ex, M = nothing)
    if iscall(ex)
        atom_dims_ok(ex, M) || return GUnknownCurvature
        f, args = operation(ex), arguments(ex)
        if f === (*)
            return mul_gcurvature(args, M)
        elseif f === (+)
            return add_gcurvature(args, M)
        end
    end
    return find_gcurvature(ex, M)
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
    ex = propagate_sign(ex, M)
    return Postwalk(x -> issym(x) || iscall(x) ? setgcurvature(x, node_gcurvature(x, M)) : x)(ex)
end
