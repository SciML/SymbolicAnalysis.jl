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
    v isa AbstractMatrix && size(v, 1) == size(v, 2) || return false
    A = float.(v)
    if eltype(A) <: Union{Float32, Float64, ComplexF32, ComplexF64}
        return LinearAlgebra.rank(A) == size(A, 1)
    end
    # Pivoted QR supports generic precision and a scale-relative rank threshold.
    R = LinearAlgebra.qr(A, LinearAlgebra.ColumnNorm()).R
    d = abs.(LinearAlgebra.diag(R))
    isempty(d) && return true
    return minimum(d) > size(A, 1) * eps(real(float(one(eltype(A))))) * maximum(d)
end

function manifold_point(ex, M = nothing; shifts::Bool = true)
    issym(ex) && return fits_manifold(ex, M)
    iscall(ex) || return false
    M === nothing || M isa SymmetricPositiveDefinite || return false
    f, args = operation(ex), arguments(ex)
    if f === inv || f === adjoint || f === transpose
        return manifold_point(args[1], M; shifts)
    elseif f === conjugation
        return manifold_point(args[1], M; shifts) && nonsingular_const(args[2])
    elseif f === (*)
        return count(!isconstarg, args) == 1 &&
            all(a -> isconstarg(a) ? positive_const(a) : manifold_point(a, M; shifts), args)
    elseif f === (+) && shifts
        return any(!isconstarg, args) &&
            all(a -> isconstarg(a) ? psd_const(a) : manifold_point(a, M; shifts), args)
    end
    return false
end

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

bare_arg_ok(a, M) = !issym(a) || fits_manifold(a, M) || fixed_param(a, M)

# Size of an argument when it is known, else `nothing`.
function known_size(a)
    if isconstarg(a)
        v = constval(a)
        return v isa AbstractArray ? size(v) : ()
    end
    issym(a) || iscall(a) || return nothing
    sh = SymbolicUtils.shape(a)
    sh isa AbstractVector || return nothing
    return Tuple(map(length, sh))
end

# Whether the known dimensions of an atom's parameters agree with its point
# argument; a parameter of unknown size is admitted.
function atom_dims_ok(f, args)
    sz = map(known_size, args)
    dim(i, k) = (s = sz[i]; s === nothing || length(s) < k ? nothing : s[k])
    agree(x, y) = x === nothing || y === nothing || x == y
    if f === quad_form || f === log_quad_form
        length(args) == 2 || return true
        return agree(dim(1, 1), dim(2, 1)) && agree(dim(1, 1), dim(2, 2))
    elseif f === lorentz_homogeneous_quadratic || f === lorentz_homogeneous_diagonal ||
            f === lorentz_nonhomogeneous_quadratic || f === lorentz_least_squares
        d = dim(length(args), 1)
        if f === lorentz_homogeneous_quadratic
            return agree(dim(1, 1), d) && agree(dim(1, 2), d)
        elseif f === lorentz_homogeneous_diagonal
            return agree(dim(1, 1), d)
        elseif f === lorentz_nonhomogeneous_quadratic
            return agree(dim(1, 1), d) && agree(dim(1, 2), d) && agree(dim(2, 1), d)
        else
            return agree(dim(1, 2), d) && agree(dim(1, 1), dim(2, 1))
        end
    end
    return true
end

gdcp_args_on_manifold(args, M) =
    all(a -> !(issym(a) || iscall(a)) || manifold_point(a, M), args)

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
            c = constval(arg)
            c isa Number || return GUnknownCurvature
            constant_prod *= c
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
    g = operation(a)
    if g === (+)
        return count(!isconstarg, arguments(a)) == 1 &&
            all(x -> isconstarg(x) ? psd_const(x) : positive_image(x, M), arguments(a))
    end
    return (
        g === conjugation || g === LinearAlgebra.diag || g === affine_map ||
            g === hadamard_product
    ) && positive_image(a, M)
end

function geodesic_rule_arg(f, a, M = nothing)
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
        f, args = operation(ex), arguments(ex)
        atom_dims_ok(f, args) || return GUnknownCurvature
        knowngcurv = false
        f_curvature = GUnknownCurvature
        f_monotonicity = (GAnyMono,)

        if gdcp_rule_applies(f, M) && !any(iscall, args) && all(a -> bare_arg_ok(a, M), args)
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
                (operation(args[1]) == LinearAlgebra.tr && positive_image(arguments(args[1])[1], M)) ||
                    (
                    operation(args[1]) == quad_form && fixed_param(arguments(args[1])[1], M) &&
                        atom_dims_ok(quad_form, arguments(args[1])) &&
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
            if (args[1] == tr || args[1] == conjugation || args[1] == diag) &&
                    spd_point(args[2], M) && consts_from(args, 3)
                return GConvex
            else
                return GUnknownCurvature
            end
        elseif gdcp_rule_applies(f, M) && any(iscall, args) &&
                all(a -> iscall(a) ? geodesic_rule_arg(f, a, M) : bare_arg_ok(a, M), args)
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
        f = operation(ex)
        if f === (*)
            return mul_gcurvature(arguments(ex), M)
        elseif f === (+)
            return add_gcurvature(arguments(ex), M)
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
