### DCP atom rules

add_dcprule(+, RealLine(), AnySign, Affine, Increasing)

"""
    dcprule(::typeof(-), x[, y])

Unary `-x` is decreasing, but binary `x - y` is *increasing* in `x` and
decreasing only in `y`. A single declared monotonicity is broadcast to every
argument by `get_arg_property`, so registering `Decreasing` flipped the first
slot too: `exp.(v) .- w`, whose Hessian is `diag(exp(vᵢ)) ⪰ 0`, certified as
`Concave`. Scalar `a - b` never reaches here — Symbolics rewrites it to
`a + (-1)*b` — but the broadcast form does.
"""
dcprule(::typeof(-), x) = makerule(RealLine(), AnySign, Affine, Decreasing), (x,)
function dcprule(::typeof(-), x, y)
    return makerule(RealLine(), AnySign, Affine, (Increasing, Decreasing)), (x, y)
end
hasdcprule(::typeof(-)) = true

add_dcprule(Base.Ref, RealLine(), AnySign, Affine, AnyMono)

"""
    dcprule(::typeof(dot), x, y)

`dot` is bilinear, so it is affine only when one side is constant: `dot(c, x)` is
affine in `x`, but `dot(x, x)` is the quadratic `‖x‖²`. Registering it as
unconditionally affine certified `dot(x, x)` as `Affine`, which is not a valid
certificate, so the curvature depends on the arguments like `^` and `norm`.
"""
function dcprule(::typeof(dot), x, y)
    return multilinear_rule(array_domain(RealLine()), (x, y))
end
hasdcprule(::typeof(dot)) = true

"""
    multilinear_rule(domain, args)

Shared by `dot`, `conv` and `kron`: each is linear in every argument separately
but not jointly, so each is affine only when all but one argument is constant.
`conv(x, x)` and `kron(X, X)` are quadratic in the same way `dot(x, x)` is.

The affine case declares `AnyMono` rather than `Increasing` because the constant
side holds the coefficients of the linear map and they may have either sign.
That only refuses a *curved* argument such as `kron(C, exp.(X))`; an affine one
composes without consulting monotonicity at all.
"""
function multilinear_rule(domain, args)
    domains = ntuple(_ -> domain, length(args))
    count(!isconstarg, args) > 1 &&
        return makerule(domains, AnySign, UnknownCurvature, AnyMono), args
    return makerule(domains, AnySign, Affine, AnyMono), args
end

"""
    dotsort(x, y)

Sorts `x` and `y` and returns the dot product of the sorted vectors.

# Arguments

    - `x::AbstractVector`: A vector.
    - `y::AbstractVector`: A vector.
"""
function dotsort(x::AbstractVector, y::AbstractVector)
    if length(x) != length(y)
        throw(DimensionMismatch("AbstractVectors must have same length"))
    end
    return dot(sort(x), sort(y))
end
Symbolics.@register_symbolic dotsort(x::AbstractVector, y::AbstractVector)

"""
    dcprule(::typeof(dotsort), x, y)

`dotsort` is a pointwise maximum of bilinear forms, so — like `dot` — it is a
convex atom only when one side is constant. With both sides symbolic it is
indefinite: for length-1 vectors it is literally `x[1]*y[1]`.
"""
function dcprule(::typeof(dotsort), x, y)
    args = (x, y)
    vec = array_domain(RealLine(), 1)
    isconstarg(x) || isconstarg(y) ||
        return makerule((vec, vec), AnySign, UnknownCurvature, AnyMono), args
    return makerule((vec, vec), AnySign, Convex, (AnyMono, increasing_if_positive ∘ minimum)),
        args
end
hasdcprule(::typeof(dotsort)) = true

add_dcprule(
    StatsBase.geomean,
    array_domain(HalfLine{Real, :open}(), 1),
    Positive,
    Concave,
    Increasing
)
add_dcprule(
    StatsBase.harmmean,
    array_domain(HalfLine{Real, :open}(), 1),
    Positive,
    Concave,
    Increasing
)

"""
    invprod(x::AbstractVector)

Returns the inverse of the product of the elements of `x`.

# Arguments

    - `x::AbstractVector`: A vector.
"""
function invprod(x::AbstractVector)
    if any(iszero(x))
        throw(DivideError())
    end
    return inv(prod(x))
end
Symbolics.@register_symbolic invprod(x::AbstractVector)

add_dcprule(invprod, array_domain(HalfLine{Real, :open}()), Positive, Convex, Decreasing)

# `eigmax`/`eigmin` build symbolic terms via Symbolics' own registration (a Base
# LinearAlgebra function belongs to Symbolics, not pirated here); we only attach
# the DCP curvature rule, exactly as for `tr`.
add_dcprule(eigmax, symmetric_domain(), AnySign, Convex, AnyMono)

add_dcprule(eigmin, symmetric_domain(), AnySign, Concave, AnyMono)

# The spectral norm is convex on every real matrix but not monotone in the
# Loewner order, so it composes only over affine arguments.
add_dcprule(opnorm, array_domain(RealLine(), 2), Positive, Convex, AnyMono)

"""
    eigsummax(m::Symmetric, k)

Returns the sum of the `k` largest eigenvalues of `m`.

# Arguments

    - `m::Symmetric`: A symmetric matrix.
    - `k::Int`: The number of largest eigenvalues to sum.
"""
function eigsummax(m::Symmetric, k::Int)
    if k < 1 || k > size(m, 1)
        throw(DomainError(k, "k must be between 1 and size(m, 1)"))
    end
    nrows = size(m, 1)
    return sum(eigvals(m, (nrows - k + 1):nrows))
end
Symbolics.@register_symbolic eigsummax(m::Matrix, k::Int)
add_dcprule(eigsummax, (array_domain(RealLine(), 2), RealLine()), AnySign, Convex, AnyMono)

"""
    eigsummin(m::Symmetric, k)

Returns the sum of the `k` smallest eigenvalues of `m`.

# Arguments

    - `m::Symmetric`: A symmetric matrix.
    - `k::Int`: The number of smallest eigenvalues to sum.
"""
function eigsummin(m::Symmetric, k::Int)
    if k < 1 || k > size(m, 1)
        throw(DomainError(k, "k must be between 1 and size(m, 1)"))
    end
    return sum(eigvals(m, 1:k))
end
Symbolics.@register_symbolic eigsummin(m::Matrix, k::Int)
add_dcprule(eigsummin, (array_domain(RealLine(), 2), RealLine()), AnySign, Concave, AnyMono)

add_dcprule(logdet, semidefinite_domain(), AnySign, Concave, AnyMono)

# `LogExpFunctions.logsumexp` on a symbolic vector must stay an unevaluated
# `logsumexp` term so the curvature pass can dispatch on it; Symbolics' own
# vector method expands it to `log(sum(exp, x))`, which erases the atom and
# leaves the expression UnknownCurvature. Each half of its
# `Union{AbstractVector{<:Num}, Arr}` signature is strictly more specific, so these
# extend rather than overwrite, and — unlike a `@register_symbolic ::Vector{Num}`
# form — emit no scalar `BasicSymbolic{SymReal}` method that could clobber the
# scalar `logsumexp`.
function LogExpFunctions.logsumexp(x::Symbolics.Arr)
    return Symbolics.wrap(
        SymbolicUtils.term(LogExpFunctions.logsumexp, SymbolicUtils.unwrap(x); type = Real)
    )
end
function LogExpFunctions.logsumexp(x::AbstractVector{<:Num})
    return Symbolics.wrap(
        SymbolicUtils.term(
            LogExpFunctions.logsumexp, map(SymbolicUtils.unwrap, x); type = Real
        )
    )
end

# `array_domain(RealLine())` (dimension-agnostic) not `..., 2)`: logsumexp is a
# vector reduction, so requiring `ndims == 2` never matched a 1-D vector.
add_dcprule(
    LogExpFunctions.logsumexp,
    array_domain(RealLine()),
    AnySign,
    Convex,
    Increasing
)

"""
    matrix_frac(x::AbstractVector, P::AbstractMatrix)

Returns the quadratic form `x' * P^{-1} * x`.

# Arguments

    - `x::AbstractVector`: A vector.
    - `P::AbstractMatrix`: A matrix.
"""
function matrix_frac(x::AbstractVector, P::AbstractMatrix)
    if length(x) != size(P, 1)
        throw(DimensionMismatch("x and P must have same length"))
    end
    return x' * inv(P) * x
end
Symbolics.@register_symbolic matrix_frac(x::AbstractVector, P::AbstractMatrix)
add_dcprule(
    matrix_frac,
    (array_domain(RealLine(), 1), definite_domain()),
    AnySign,
    Convex,
    AnyMono
)

add_dcprule(maximum, array_domain(RealLine()), AnySign, Convex, Increasing)

add_dcprule(minimum, array_domain(RealLine()), AnySign, Concave, Increasing)

# The matrix rules for `sqrt`, `log` and `inv` are stated in the Loewner order.
# Operator concavity licenses PSD-weighted functionals — `tr(f(X))`, or
# `sum(f(X)) = e'f(X)e` — but says nothing about a single entry `eᵢ'f(X)eⱼ`.
const LOEWNER_MATRIX_ATOMS = (sqrt, log, inv)

# Tests the *argument's* symtype, not the term's: `Postwalk` rebuilds each node
# through `maketerm`, which re-infers the symtype from the operation, so an
# annotated `sqrt(X)` reports `Real` rather than `Matrix{Real}` by the time the
# curvature pass reaches it. A leaf argument keeps its own symtype.
function is_loewner_matrix_term(x)
    iscall(x) || return false
    operation(x) in LOEWNER_MATRIX_ATOMS || return false
    args = arguments(x)
    return length(args) == 1 && SymbolicUtils.symtype(args[1]) <: AbstractArray{<:Any, 2}
end

"""
    dcprule(::typeof(minimum), x)

An entrywise `minimum` must not consume a Loewner-order matrix rule. The
smallest entry of a positive definite matrix can be an off-diagonal one, so
`minimum(sqrt(X))` and `minimum(log(X))` are indefinite over the SPD cone
(second differences of −0.274 and −0.753 at one base point, +0.462 and +1.305
at another) while both certified as `Concave`.

`maximum` needs no such guard: for a positive definite `M`,
`M[i,j] ≤ √(M[i,i]·M[j,j]) ≤ max(M[i,i], M[j,j])`, so the largest entry is
always on the diagonal and `maximum(inv(X))` is a max of convex functions.
"""
function dcprule(::typeof(minimum), x)
    dom = array_domain(RealLine())
    is_loewner_matrix_term(x) &&
        return makerule(dom, AnySign, UnknownCurvature, AnyMono), (x,)
    return makerule(dom, AnySign, Concave, Increasing), (x,)
end

# `norm(x, p)` is only a norm (hence convex) for p >= 1. For 0 < p < 1 the
# generalized "norm" (sum |x_i|^p)^(1/p) is concave on the nonnegative orthant
# and has no DCP curvature for sign-unknown arguments; for p <= 0 it is not
# DCP-representable. The rule depends on the value of `p`, so it is a
# specialized `dcprule` method like `^` rather than a static table entry.
function dcprule(::typeof(norm), x, p)
    pv = constval(p)
    args = (x, p)
    if !(pv isa Number)
        return makerule(array_domain(RealLine()), Positive, UnknownCurvature, AnyMono), args
    elseif pv >= 1
        return makerule(array_domain(RealLine()), Positive, Convex, increasing_if_positive), args
    elseif pv > 0
        curv = getsign(x) == Positive ? Concave : UnknownCurvature
        return makerule(array_domain(HalfLine()), Positive, curv, Increasing), args
    else
        return makerule(array_domain(RealLine()), Positive, UnknownCurvature, AnyMono), args
    end
end
function dcprule(::typeof(norm), x)
    return makerule(array_domain(RealLine()), Positive, Convex, increasing_if_positive), (x,)
end
hasdcprule(::typeof(norm)) = true


"""
    perspective(f::Function, x, s::Real)

Returns the perspective function `s * f(x / s)`.

# Arguments

    - `f::Function`: A function.
    - `x`: A Real.
    - `s::Real`: A positive Real.
"""
function perspective(f::Function, x, s::Real)
    if s < 0
        throw(DomainError(s, "s must be positive"))
    end
    if s == 0
        return zero(typeof(f(x)))
    end
    return s * f(x / s)
end
Symbolics.@register_symbolic perspective(f::Function, x, s::Real)
add_dcprule(
    perspective,
    (function_domain(), RealLine(), Positive),
    getsign,
    getcurvature,
    AnyMono
)

"""
    quad_form(x::AbstractVector, P::AbstractMatrix)

Returns the quadratic form `x' * P * x`.

# Arguments

    - `x::AbstractVector`: A vector.
    - `P::AbstractMatrix`: A matrix.
"""
function quad_form(x::AbstractVector, P::AbstractMatrix)
    if length(x) != size(P, 1)
        throw(DimensionMismatch("x and P must have same length"))
    end
    return x' * P * x
end
Symbolics.@register_symbolic quad_form(x::AbstractVector, P::AbstractMatrix)

"""
    dcprule(::typeof(quad_form), x, P)

`x'Px` is quadratic in `x` but **linear** in `P`, so it is convex only for a
constant positive definite `P`. Registering it as unconditionally `Convex`
certified both the both-symbolic form (indefinite: `quad_form([a], [b;;])` is
`a^2*b`, whose second difference along `(1, -1)` is `-2`) and a constant
indefinite `P` (`quad_form(x, [1 0; 0 -1])` is `x[1]^2 - x[2]^2`).

`P` is checked with `isposdef`, matching the `semidefinite_domain()` the rule
declares; a singular positive semidefinite `P` therefore gets no certificate.

Slot 1 is only `increasing_if_positive` when `P` is **entrywise** nonnegative.
`x'Px = Σ P[i,j]·x[i]·x[j]` is nondecreasing in each `x[i]` over the nonnegative
orthant only then; positive definiteness alone does not give it. With
`P = [1 -0.9; -0.9 1]` (which is positive definite) the composition certified
`quad_form(exp.(v), P)` as `Convex` while it is not — the midpoint exceeds the
chord by 1.92 between `v = [0.64, 2.0]` and `[1.64, 2.3]`.
"""
function dcprule(::typeof(quad_form), x, P)
    args = (x, P)
    dom = (array_domain(RealLine(), 1), semidefinite_domain())
    Pv = constval(P)
    if Pv isa AbstractMatrix && isconstarg(P) && issymmetric(Pv) && isposdef(Pv)
        mono1 = all(>=(0), Pv) ? increasing_if_positive : AnyMono
        return makerule(dom, Positive, Convex, (mono1, Increasing)), args
    end
    return makerule(dom, AnySign, UnknownCurvature, AnyMono), args
end
hasdcprule(::typeof(quad_form)) = true

function quad_over_lin(x::AbstractVector{<:Real}, y::Real)
    if getsign(y) == Negative
        throw(DomainError(y, "y must be positive"))
    end
    return sum(x .^ 2) / y
end

# On Symbolics v7 both scalar- and array-valued registered symbolics share the
# concrete type `BasicSymbolic{SymReal}`, so registering the vector form would
# generate the same `quad_over_lin(::BasicSymbolic{SymReal}, ::Real)` method as
# the scalar one below and collide. The scalar registration alone covers both.

"""
    quad_over_lin(x::Real, y::Real)

Returns the quadratic over linear form `x^2 / y`.

# Arguments

    - `x`: A Real or a vector.
    - `y::Real`: A positive Real.
"""
function quad_over_lin(x::Real, y::Real)
    if getsign(y) == Negative
        throw(DomainError(y, "y must be positive"))
    end
    return x^2 / y
end

Symbolics.@register_symbolic quad_over_lin(x::Real, y::Real)

add_dcprule(
    quad_over_lin,
    (array_domain(RealLine()), HalfLine{Real, :open}()),
    Positive,
    Convex,
    (increasing_if_positive, Decreasing)
)

add_dcprule(
    quad_over_lin,
    (RealLine(), HalfLine{Real, :open}()),
    Positive,
    Convex,
    (increasing_if_positive, Decreasing)
)

add_dcprule(sum, array_domain(RealLine(), 2), AnySign, Affine, Increasing)

"""
    sum_largest(x::AbstractMatrix, k)

Returns the sum of the `k` largest elements of `x`.

# Arguments

    - `x::AbstractMatrix`: A matrix.
    - `k::Int`: The number of largest elements to sum.
"""
function sum_largest(x::AbstractMatrix, k::Integer)
    return sum(sort(vec(x))[(end - k):end])
end
Symbolics.@register_symbolic sum_largest(x::AbstractMatrix, k::Integer)
add_dcprule(sum_largest, (array_domain(RealLine(), 2), ℤ), AnySign, Convex, Increasing)

"""
    sum_smallest(x::AbstractMatrix, k)

Returns the sum of the `k` smallest elements of `x`.

# Arguments

    - `x::AbstractMatrix`: A matrix.
    - `k::Int`: The number of smallest elements to sum.
"""
function sum_smallest(x::AbstractMatrix, k::Integer)
    return sum(sort(vec(x))[1:k])
end

Symbolics.@register_symbolic sum_smallest(x::AbstractArray, k::Integer)
add_dcprule(sum_smallest, (array_domain(RealLine(), 2), ℤ), AnySign, Concave, Increasing)

add_dcprule(tr, array_domain(RealLine(), 2), AnySign, Affine, Increasing)

"""
    trinv(x::AbstractMatrix)

Returns the trace of the inverse of `x`.

# Arguments

    - `x::AbstractMatrix`: A matrix.
"""
function trinv(x::AbstractMatrix)
    return tr(inv(x))
end
Symbolics.@register_symbolic trinv(x::AbstractMatrix)
add_dcprule(trinv, definite_domain(), Positive, Convex, AnyMono)

"""
    tv(x::AbstractVector{<:Real})

Returns the total variation of `x`, defined as `sum_i |x_{i+1} - x_i|`.

# Arguments

    - `x::AbstractVector`: A vector.
"""
function tv(x::AbstractVector{<:Real})
    return sum(abs.(x[2:end] - x[1:(end - 1)]))
end
Symbolics.@register_symbolic tv(x::AbstractVector) false
add_dcprule(tv, array_domain(RealLine(), 1), Positive, Convex, AnyMono)

"""
    tv(x::AbstractVector{<:AbstractMatrix})

Returns the total variation of `x`, defined as `sum_{i,j} |x_{k+1}[i,j] - x_k[i,j]|`.

# Arguments

    - `x::AbstractVector`: A vector of matrices.
"""
function tv(x::AbstractVector{<:AbstractMatrix})
    return sum(
        map(1:(size(x, 1) - 1)) do i
            map(1:(size(x, 2) - 1)) do j
                norm([x[k][i + 1, j] - x[k][i, j] for k in eachindex(x)])
            end
        end
    )
end
add_dcprule(tv, array_domain(array_domain(RealLine(), 2), 1), Positive, Convex, AnyMono)

add_dcprule(abs, ℂ, Positive, Convex, increasing_if_positive)

# `abs2(x) = |x|²` is a squared norm, so it shares `abs`'s shape: convex,
# nonnegative, and increasing only where its argument is.
add_dcprule(abs2, ℂ, Positive, Convex, increasing_if_positive)

add_dcprule(conj, ℂ, AnySign, Affine, AnyMono)

add_dcprule(exp, RealLine(), Positive, Convex, Increasing)
add_dcprule(exp2, RealLine(), Positive, Convex, Increasing)
add_dcprule(exp10, RealLine(), Positive, Convex, Increasing)
# `expm1` shares `exp`'s curvature but is negative below zero.
add_dcprule(expm1, RealLine(), AnySign, Convex, Increasing)
# `cosh` is even, so it increases only on the nonnegative half line.
add_dcprule(cosh, RealLine(), Positive, Convex, increasing_if_positive)
# `hypot(x, y)` is the Euclidean norm of `(x, y)`: convex and nonnegative, and
# increasing in each argument only where that argument is nonnegative.
add_dcprule(hypot, (RealLine(), RealLine()), Positive, Convex, increasing_if_positive)

"""
    dcprule(::typeof(exp), x)

`exp(X)` on a matrix is the matrix exponential, not the scalar law applied
pointwise: `expm` is neither operator convex nor operator monotone, so
`sum(exp(X))` and `tr(exp(X))` are indefinite even on the SPD cone. The
elementwise `exp.(X)` arrives through `broadcast` and keeps the scalar law via
`elementwise_dcprule`.
"""
function dcprule(::typeof(exp), x)
    SymbolicUtils.symtype(x) <: AbstractArray &&
        return makerule(array_domain(RealLine()), AnySign, UnknownCurvature, AnyMono), (x,)
    return @invoke dcprule(exp::Any, x)
end

Symbolics.@register_symbolic LogExpFunctions.xlogx(x::Real)
add_dcprule(xlogx, RealLine(), AnySign, Convex, AnyMono)

"""
    huber(x, M=1)

Returns the Huber loss function of `x` with threshold `M`.

# Arguments

    - `x::Real`: A Real.
    - `M::Real`: The threshold.
"""
function huber(x::Real, M::Real = 1)
    if M < 0
        throw(DomainError(M, "M must be positive"))
    end

    if abs(x) <= M
        return x^2
    else
        return 2 * M * abs(x) - M^2
    end
end
Symbolics.@register_symbolic huber(x::Real, M::Real)

"""
    dcprule(::typeof(huber), x, M)

Huber is convex in `x`, but for `abs(x) > M` it equals `2M*abs(x) - M^2`, which
is **concave** in `M` (`d²/dM² = -2`). The threshold must therefore be a
constant for the `Convex` certificate to hold.
"""
function dcprule(::typeof(huber), x, M)
    args = (x, M)
    dom = (RealLine(), HalfLine())
    isconstarg(M) ||
        return makerule(dom, AnySign, UnknownCurvature, AnyMono), args
    return makerule(dom, Positive, Convex, increasing_if_positive), args
end
hasdcprule(::typeof(huber)) = true

add_dcprule(imag, ℂ, AnySign, Affine, AnyMono)

add_dcprule(inv, HalfLine{Real, :open}(), Positive, Convex, Decreasing)
add_dcprule(log, HalfLine{Real, :open}(), AnySign, Concave, Increasing)
# `log2`/`log10` are `log` rescaled by a positive constant.
add_dcprule(log2, HalfLine{Real, :open}(), AnySign, Concave, Increasing)
add_dcprule(log10, HalfLine{Real, :open}(), AnySign, Concave, Increasing)

# Matrix-valued atoms (`log`, `inv`, `sqrt` of a symbolic matrix). On Symbolics
# v7 / SymbolicUtils v4 a symbolic matrix unwraps to `BasicSymbolic{SymReal}` —
# the same concrete type as a symbolic scalar — so a `@register_symbolic` macro
# would emit a `f(::BasicSymbolic{SymReal})` method that overwrites SymbolicUtils'
# own scalar `f` and aborts precompilation. Build the matrix term directly off the
# `Arr` wrapper instead via `SymbolicUtils.term`, which leaves the scalar methods
# untouched. Some scalar atoms (e.g. `sqrt`) register a `promote_shape` that
# rejects matrix shapes, which `term` would invoke unless the shape is supplied up
# front, so pass `shape` explicitly.
function matrix_atom(f, A::Symbolics.Arr)
    a = SymbolicUtils.unwrap(A)
    return Symbolics.wrap(
        SymbolicUtils.term(
            f, a; type = SymbolicUtils.symtype(a), shape = SymbolicUtils.shape(a)
        )
    )
end
function matrix_atom(f, A::AbstractMatrix{<:Num})
    a = SymbolicUtils.unwrap.(A)
    return Symbolics.wrap(
        SymbolicUtils.term(f, a; type = Matrix{Real}, shape = map(Base.OneTo, size(A)))
    )
end

# A symbolic matrix–matrix product returns a bare `BasicSymbolic`; re-wrap the
# 2-argument `Arr` product so it reaches the matrix-atom methods below (which
# dispatch on `Arr`), restoring the array shape. This 2-arg method is more
# specific than Symbolics' variadic `*`, so it does not overwrite it.
function Base.:*(x::Symbolics.Arr{<:Any, 2}, y::Symbolics.Arr{<:Any, 2})
    return Symbolics.wrap(SymbolicUtils.unwrap(x) * SymbolicUtils.unwrap(y))
end

# SymbolicUtils v4 gives `log` a matrix-permissive `promote_shape` but leaves
# `sqrt` scalar-only, so rewriting a matrix `sqrt` term (e.g. during the
# analysis walk's `maketerm`) throws "Invalid shapes for sqrt". Add a matrix
# rule for `sqrt`; dispatching on the concrete `ShapeVecT` (rather than the
# `ShapeT` union SymbolicUtils uses) makes this strictly more specific, so it
# extends rather than overwrites the existing method.
function SymbolicUtils.promote_shape(::typeof(sqrt), sh::SymbolicUtils.ShapeVecT)
    (length(sh) == 0 || length(sh) == 2) && return sh
    throw(ArgumentError("Invalid shapes for sqrt: $(sh)."))
end

Base.log(A::Symbolics.Arr) = matrix_atom(log, A)
Base.log(A::Matrix{Num}) = matrix_atom(log, A)
add_dcprule(log, array_domain(RealLine(), 2), Positive, Concave, Increasing)

add_dcprule(inv, semidefinite_domain(), AnySign, Convex, Decreasing)

Base.sqrt(A::Symbolics.Arr) = matrix_atom(sqrt, A)
add_dcprule(sqrt, semidefinite_domain(), Positive, Concave, Increasing)

add_dcprule(
    kldivergence,
    (array_domain(HalfLine{Real, :open}(), 1), array_domain(HalfLine{Real, :open}(), 1)),
    Positive,
    Convex,
    AnyMono
)

"""
    lognormcdf(x::Real)

Returns the log of the normal cumulative distribution function of `x`.

# Arguments

    - `x::Real`: A Real.
"""
function lognormcdf(x::Real)
    return logcdf(Normal(), x)
end
Symbolics.@register_symbolic lognormcdf(x::Real)
add_dcprule(lognormcdf, RealLine(), Negative, Concave, Increasing)

# `log1p` is negative only on `(-1, 0)`; it is positive on all of `(0, Inf)`, which
# is inside its own declared domain. A `Negative` sign here made
# `increasing_if_positive` hand `abs` a `Decreasing` slot, certifying the strictly
# concave `abs(log1p(a))` (for `a > 0`) as `Convex`.
add_dcprule(log1p, Interval{:open, :open}(-1, Inf), AnySign, Concave, Increasing)

# Symbolics rewrites these three into `log(1 + exp(x))`, `log(exp(x) + exp(y))` and
# `log(exp(x) - 1)`, each `concave of convex` and so uncertifiable as written.
Symbolics.@register_symbolic LogExpFunctions.log1pexp(x::Real)
add_dcprule(LogExpFunctions.log1pexp, RealLine(), Positive, Convex, Increasing)

Symbolics.@register_symbolic LogExpFunctions.logaddexp(x::Real, y::Real)
add_dcprule(
    LogExpFunctions.logaddexp, (RealLine(), RealLine()), AnySign, Convex, Increasing
)

# `logexpm1` inverts `log1pexp`, so it is real only for `x > 0`.
Symbolics.@register_symbolic LogExpFunctions.logexpm1(x::Real)
add_dcprule(LogExpFunctions.logexpm1, HalfLine{Real, :open}(), AnySign, Concave, Increasing)

add_dcprule(max, (RealLine(), RealLine()), AnySign, Convex, Increasing)
add_dcprule(min, (RealLine(), RealLine()), AnySign, Concave, Increasing)

# special cases which depend on arguments:

# The scalar power laws, applied to a scalar base or elementwise to an array one.
function power_rule(x, i)
    args = (x, i)
    if !(i isa Real)
        # The power laws below all need a numeric exponent. The remaining case
        # with a DCP curvature is a constant base: `c^g == exp(g*log(c))` is
        # convex for `c > 0`, monotone in the direction of `log(c)`.
        base = constval(x)
        if base isa Real && base > 0
            mono = base > 1 ? Increasing : base < 1 ? Decreasing : AnyMono
            return makerule(RealLine(), Positive, Convex, (AnyMono, mono)), args
        end
        return makerule(RealLine(), AnySign, UnknownCurvature, AnyMono), args
    end
    if isone(i)
        return makerule(RealLine(), AnySign, Affine, Increasing), args
    elseif isinteger(i) && iseven(i)
        return makerule(RealLine(), Positive, Convex, increasing_if_positive), args
    elseif isinteger(i) && isodd(i)
        return makerule(HalfLine(), Positive, Convex, Increasing), args
    elseif i >= 1
        return makerule(HalfLine(), Positive, Convex, Increasing), args
    elseif i > 0 && i < 1
        return makerule(HalfLine(), Positive, Concave, Increasing), args
    elseif i < 0
        return makerule(HalfLine{Float64, :closed}(), Positive, Convex, Increasing), args
    end
    return makerule(RealLine(), AnySign, UnknownCurvature, AnyMono), args
end

function dcprule(::typeof(^), x::Symbolic, i)
    # A literal exponent is wrapped as a constant `BasicSymbolic`, so
    # `isinteger`/`isone`/comparisons in `power_rule` would operate on a symbolic
    # and error; `constval` unwraps it to the underlying number (identity for an
    # already-numeric exponent).
    i = constval(i)
    # A matrix power is a different function from the scalar power laws:
    # `tr(X^2)` for an unconstrained `X` is indefinite (it contains cross terms
    # `X[i,j]*X[j,i]`), so "even integer power is convex" would certify a
    # non-convex expression. An elementwise power is not this case; it reaches
    # `power_rule` through the `broadcast` rule instead.
    if SymbolicUtils.symtype(x) <: AbstractArray
        return makerule(array_domain(RealLine()), AnySign, UnknownCurvature, AnyMono), (x, i)
    end
    return power_rule(x, i)
end
dcprule(::typeof(Base.literal_pow), f, x...) = dcprule(^, x...)

"""
    dcprule(::typeof(/), num, den)

`/` is bilinear in the same sense as `*`, so it is DCP only when one side is
constant. A nonzero constant denominator is an affine rescaling, flipping
curvature and sign when it is negative.

A constant numerator is the `inv` atom: `c/x` is convex-decreasing for `c > 0`
and concave-increasing for `c < 0` — but **only on `x > 0`**, and the
denominator must be *known* positive for the rule to fire. That is deliberately
stricter than the unchecked-precondition convention `log`, `sqrt` and `geomean`
use, and the asymmetry is the point: those are only wrong where the function
does not exist, whereas `1/x` exists for every nonzero `x` and is *concave*
below zero. Certifying it `Convex` for a sign-unknown argument would be a false
certificate over the function's own real domain — `f(-1) = -1`, `f(2) = 0.5`,
and the midpoint `f(0.5) = 2` is far above the chord `-0.25`.

Symbolics rewrites `inv(x)` to `/(1, x)` and `x^-n` to `/(1, x)^n`, so this rule
is also what gives those their curvature.
"""
function dcprule(::typeof(/), num, den)
    args = (num, den)
    nv = constval(num)
    dv = constval(den)
    unknown = makerule((RealLine(), RealLine()), AnySign, UnknownCurvature, AnyMono)
    if dv isa Number
        iszero(dv) && return unknown, args
        mono = dv > 0 ? Increasing : Decreasing
        return makerule((RealLine(), RealLine()), mul_sign(args), Affine, (mono, AnyMono)),
            args
    elseif nv isa Number && !iszero(nv) && known_positive(den)
        domain = (RealLine(), HalfLine{Real, :open}())
        return if nv > 0
            makerule(domain, Positive, Convex, (AnyMono, Decreasing)), args
        else
            makerule(domain, Negative, Concave, (AnyMono, Increasing)), args
        end
    end
    return unknown, args
end
hasdcprule(::typeof(/)) = true

# Positivity established rather than assumed: either a propagated sign says so, or
# the argument was declared on a domain inside the positive half line. `analyze`
# overwrites `Sign` metadata during propagation, so a `VarDomain` declaration is
# how a caller states the precondition durably.
function known_positive(x)
    getsign(x) == Positive && return true
    _has_vardomain(x) || return false
    return try
        issubset(getmetadata(x, VarDomain), HalfLine{Real, :open}())
    catch e
        e isa MethodError ? false : rethrow()
    end
end

hasdcprule(::typeof(^)) = true

add_dcprule(real, ℂ, AnySign, Affine, Increasing)

function rel_entr(x::Real, y::Real)
    if x < 0 || y < 0
        throw(DomainError((x, y), "x and y must be positive"))
    end
    if x == 0
        return 0
    end
    return x * log(x / y)
end
Symbolics.@register_symbolic rel_entr(x::Real, y::Real)
add_dcprule(
    rel_entr,
    (HalfLine{Real, :open}(), HalfLine{Real, :open}()),
    AnySign,
    Convex,
    (AnyMono, Decreasing)
)

add_dcprule(sqrt, HalfLine(), Positive, Concave, Increasing)

Symbolics.@register_symbolic LogExpFunctions.xexpx(x::Real)

"""
    dcprule(::typeof(xexpx), x)

`x*exp(x)` is convex only on `[-2, Inf)`: at `x = -4` the second difference is
`-0.0366`, at `x = -3` it is `-0.0498`, and the inflection sits exactly at `-2`.
So unlike `log` or `sqrt`, this atom **exists** outside its declared domain and is
concave there — the same shape that made `/` establish its precondition rather
than assume it, rather than the `log`/`sqrt` case where the declared domain is
only where the function is defined at all.

Registering it unguarded would have been a false certificate the moment the atom
became reachable: until this rule existed `xexpx(s)` traced to `*` and degraded,
which is the only reason the unenforced domain never bit.
"""
function dcprule(::typeof(xexpx), x)
    args = (x,)
    known_positive(x) ||
        return makerule(RealLine(), AnySign, UnknownCurvature, AnyMono), args
    return makerule(HalfLine(), Positive, Convex, Increasing), args
end
hasdcprule(::typeof(xexpx)) = true

dcprule(::typeof(conv), args...) = multilinear_rule(array_domain(RealLine(), 1), args)
hasdcprule(::typeof(conv)) = true

add_dcprule(cumsum, array_domain(RealLine()), AnySign, Affine, Increasing)

add_dcprule(diagm, array_domain(RealLine(), 1), AnySign, Affine, Increasing)

add_dcprule(diag, array_domain(RealLine(), 2), AnySign, Affine, Increasing)

add_dcprule(diff, array_domain(RealLine()), AnySign, Affine, Increasing)

add_dcprule(hcat, array_domain(array_domain(RealLine(), 1), 1), AnySign, Affine, Increasing)

dcprule(::typeof(kron), args...) = multilinear_rule(array_domain(RealLine(), 2), args)
hasdcprule(::typeof(kron)) = true

add_dcprule(reshape, array_domain(RealLine(), 2), AnySign, Affine, Increasing)

add_dcprule(tril, array_domain(RealLine(), 2), AnySign, Affine, Increasing)

add_dcprule(triu, array_domain(RealLine(), 2), AnySign, Affine, Increasing)

add_dcprule(vec, array_domain(RealLine(), 2), AnySign, Affine, Increasing)

add_dcprule(vcat, array_domain(array_domain(RealLine(), 1), 1), AnySign, Affine, Increasing)

"""
    elementwise_dcprule(f, args...)

The rule for `f` applied elementwise, as under a broadcast. Defaults to `f`'s own
rule. Atoms that mean something different on a matrix than pointwise — `^` is a
matrix power, `exp` is the matrix exponential — guard against the matrix meaning
in `dcprule` and restore the scalar law here.
"""
elementwise_dcprule(f, args...) = dcprule(f, args...)
elementwise_dcprule(::typeof(^), x, i) = power_rule(x, constval(i))
elementwise_dcprule(::typeof(exp), x) = @invoke dcprule(exp::Any, x)

"""
    dcprule(::typeof(broadcast), f, x...)

A broadcast carries its function as the first argument, so the rule is the one
`f` has when applied elementwise. Two cases are not a rule lookup: `*` has no
table entry (it is special-cased by the multiplication helpers, which pick out
the constant factor and flip the curvature on a negative one), and a broadcasted
function with no rule at all must degrade to `UnknownCurvature` rather than fail
the lookup.
"""
function dcprule(::typeof(broadcast), f, x...)
    # The broadcasted function is wrapped as a constant symbolic (e.g.
    # `broadcast(exp, z)` carries a symbolic `exp`); `constval` recovers the
    # underlying function (identity for a plain function).
    g = constval(f)
    if Symbol(g) == :*
        # Monotonicity in the non-constant factor is the direction of the
        # constant one, so that the composition check downstream reproduces the
        # flip `mul_curvature` has already applied.
        s = const_factor_sign(x)
        mono = s == Negative ? Decreasing : s == Positive ? Increasing : AnyMono
        return makerule(array_domain(RealLine()), mul_sign(x), mul_curvature(x), mono), x
    end
    hasdcprule(g) || return no_rule(array_domain(RealLine()), x)
    return elementwise_dcprule(g, x...)
end
hasdcprule(::typeof(broadcast)) = true

# add_dcprule(broadcast, (function_domain, array_domain(RealLine())), AnySign, Affine, (AnyMono, AnyMono))

add_dcprule(Base.adjoint, array_domain(RealLine(), 1), AnySign, Affine, Increasing)
# `Base.transpose`, not `LinearAlgebra.transpose`: on Julia 1.13 the owner is
# `Base` and the name is no longer public in `LinearAlgebra`, exactly as for
# `adjoint` just above. They are the same function on every supported version.
add_dcprule(Base.transpose, array_domain(RealLine(), 1), AnySign, Affine, Increasing)
add_dcprule(Base.getindex, array_domain(RealLine(), 1), AnySign, Affine, AnyMono)

# On Symbolics v7 / SymbolicUtils v4, reductions and maps over symbolic arrays
# trace to `SymbolicUtils.Mapreducer`/`SymbolicUtils.Mapper` operations rather
# than to `sum`/`map` themselves, so the static rule table never sees them.
# A plain sum — `mapreduce(identity, add_sum, x)`
# for any `dims`/`init` — delegates to the registered `sum` rule (an `init`
# only shifts by a constant, which preserves the affine composition), and
# `map(f, xs...)` delegates to `f`'s rule exactly like `broadcast`.
hasdcprule(::SymbolicUtils.Mapreducer{typeof(identity), typeof(Base.add_sum)}) = true
function dcprule(
        ::SymbolicUtils.Mapreducer{typeof(identity), typeof(Base.add_sum)}, args...
    )
    return dcprules_dict[sum], args
end

# `maximum`/`minimum` over a symbolic array reduce with `max`/`min`, so they
# trace to `Mapreducer{identity, max}` / `Mapreducer{identity, min}` and, like
# `sum`, never reach the `maximum`/`minimum` rules directly. Reducing over any
# `dims` preserves the convex-increasing / concave-increasing curvature, so
# delegate to the registered rule regardless of `dims`/`init`.
hasdcprule(::SymbolicUtils.Mapreducer{typeof(identity), typeof(max)}) = true
function dcprule(::SymbolicUtils.Mapreducer{typeof(identity), typeof(max)}, args...)
    return dcprules_dict[maximum], args
end

hasdcprule(::SymbolicUtils.Mapreducer{typeof(identity), typeof(min)}) = true
function dcprule(::SymbolicUtils.Mapreducer{typeof(identity), typeof(min)}, args...)
    # Through `dcprule` rather than the rule table, so the Loewner-order guard
    # applies — `minimum(sqrt(X))` reaches the reducer, never `minimum` itself.
    return dcprule(minimum, args...)
end

hasdcprule(op::SymbolicUtils.Mapper) = hasdcprule(op.f)
# `map` applies `op.f` elementwise, so it takes the elementwise rule for the same
# reason `broadcast` does.
dcprule(op::SymbolicUtils.Mapper, args...) = elementwise_dcprule(op.f, args...)

# A matrix or vector assembled from scalar expressions (`[u[1] u[2]; u[2] u[3]]`)
# traces to a `SymbolicUtils.array_literal` term whose first argument is the size
# tuple and whose remaining arguments are the elements in column-major order.
# Assembling an array is a linear map of its elements, so the construction is
# affine and monotonically increasing in each of them; without this rule the
# assembled matrix carries no curvature, which in turn blocks every atom taking a
# matrix argument (`logdet`, `eigmax`, ...) from composing.
add_dcprule(
    SymbolicUtils.array_literal, array_domain(RealLine()), AnySign, Affine, Increasing
)
