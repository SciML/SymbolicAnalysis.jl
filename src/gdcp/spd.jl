### DGCP Atoms

# `logdet` on SPD is geodesically linear, but not nonnegative: logdet(εI) → -∞
# as ε → 0 (e.g. logdet(0.01I₂) ≈ -9.21). Declaring `Positive` made
# `abs(logdet(X) + tr(X))` compose as GConvex while a geodesic midpoint refutes it.
add_gdcprule(
    LinearAlgebra.logdet,
    SymmetricPositiveDefinite,
    AnySign,
    GLinear,
    GIncreasing
)

"""
    conjugation(X, B)

Conjugation of a matrix `X` by a matrix `B` is defined as `B'X*B`.

# Arguments

    - `X::Matrix`: A symmetric positive definite matrix.
    - `B::Matrix`: A matrix.
"""
function conjugation(X, B)
    return B' * X * B
end

# Build an unevaluated matrix-valued term off a symbolic matrix argument. The
# `@register_array_symbolic`-generated `promote_shape` for these custom atoms is
# not `maketerm`-safe on Symbolics v7 (rewriting a nested expression such as
# `conjugation(inv(X), A)` recomputes a dimensionless `Array{T}` type and throws),
# so the symbolic methods are defined directly via `SymbolicUtils.term`, which
# round-trips cleanly through the rewriter.
function array_atom_term(f, X, args...; type = Matrix{Real})
    x = SymbolicUtils.unwrap(X)
    return Symbolics.wrap(
        SymbolicUtils.term(f, x, map(SymbolicUtils.unwrap, args)...; type = type)
    )
end

conjugation(X::Symbolics.Arr, B::AbstractMatrix) = array_atom_term(conjugation, X, B)

add_gdcprule(conjugation, SymmetricPositiveDefinite, Positive, GConvex, GIncreasing)

# Symbolics v7 builds a `tr` term natively for symbolic matrices, so no
# registration is needed here.
add_gdcprule(LinearAlgebra.tr, SymmetricPositiveDefinite, Positive, GConvex, GIncreasing)

add_gdcprule(sum, SymmetricPositiveDefinite, Positive, GConvex, GIncreasing)

add_gdcprule(adjoint, SymmetricPositiveDefinite, Positive, GLinear, GIncreasing)

"""
    scalar_mat(X, k=size(X, 1))

Scalar matrix of a symmetric positive definite matrix `X` is defined as `tr(X)*I(k)`.

# Arguments

    - `X::Matrix`: A symmetric positive definite matrix.
    - `k::Int`: The size of the identity matrix.
"""
function scalar_mat(X, k = size(X, 1))
    return tr(X) * I(k)
end

@register_symbolic scalar_mat(X::Matrix{Num}, k::Int)

add_gdcprule(scalar_mat, SymmetricPositiveDefinite, Positive, GConvex, GIncreasing)

add_gdcprule(LinearAlgebra.diag, SymmetricPositiveDefinite, Positive, GConvex, GIncreasing)

# """
#     pinching(X, Ps)

# Pinching of a symmetric positive definite matrix `X` by a set of symmetric positive definite matrices `Ps` is defined as `sum(Ps[i]*X*Ps[i])`.

# # Arguments
#     - `X::Matrix`: A symmetric positive definite matrix.
#     - `Ps::Vector`: A vector of symmetric positive definite matrices.
# """
# function pinching(X, Ps)
#     return sum(Ps[i]*X*Ps[i] for i in eachindex(Ps); dims = 1)
# end

# @register_symbolic pinching(X::Matrix{Num}, Ps::Vector{Union{Symbolics.Arr, Matrix{Num}}})

# add_gdcprule(pinching, SymmetricPositiveDefinite, Positive, GConvex, GIncreasing)

"""
    sdivergence(X, Y)

Symmetric divergence of two symmetric positive definite matrices `X` and `Y` is defined as `logdet((X+Y)/2) - 1/2*logdet(X*Y)`.

# Arguments

    - `X::Matrix`: A symmetric positive definite matrix.
    - `Y::Matrix`: A symmetric positive definite matrix.

The gDCP rule (`Positive`, `GConvex`) applies only when every constant argument is
a numeric positive definite matrix (see `pd_constant`); otherwise the analysis
returns `AnySign` and `GUnknownCurvature` for it and every expression containing it.
"""
function sdivergence(X, Y)
    return logdet((X + Y) / 2) - 1 / 2 * logdet(X * Y)
end

# `@register_symbolic sdivergence(X::Matrix{Num}, Y::Matrix)` is not used: on
# Symbolics v7 the macro expands a two-array signature into a combinatorial set of
# `Num`/`BasicSymbolic`/`Arr` wrapper methods that are mutually ambiguous (Aqua
# flags them). Build the unevaluated `sdivergence` term directly off the symbolic
# argument instead (the gDCP pass only needs `operation(ex) == sdivergence`),
# matching the `array_atom_term`/`conjugation` pattern above.
function sdivergence(X::Symbolics.Arr, Y::AbstractMatrix)
    return array_atom_term(sdivergence, X, Y; type = Real)
end
function sdivergence(X::AbstractMatrix, Y::Symbolics.Arr)
    return array_atom_term(sdivergence, X, Y; type = Real)
end
function sdivergence(X::Symbolics.Arr, Y::Symbolics.Arr)
    return array_atom_term(sdivergence, X, Y; type = Real)
end
add_gdcprule(sdivergence, SymmetricPositiveDefinite, Positive, GConvex, GIncreasing)

constant_factors_ok(::typeof(sdivergence), args) = all(a -> !isconstarg(a) || pd_constant(a), args)

# Symbolic geodesic distance must remain an unevaluated `distance` term so the
# gDCP pass can dispatch on `operation(ex) == Manifolds.distance`. The SPD and
# Lorentz registrations both collapse to the same all-`BasicSymbolic{SymReal}`
# method on Symbolics v7 and would collide, so the term is built directly off the
# symbolic point argument here.
function Manifolds.distance(
        M::Manifolds.SymmetricPositiveDefinite,
        X::AbstractMatrix,
        Y::Symbolics.Arr
    )
    return Symbolics.wrap(
        SymbolicUtils.term(Manifolds.distance, M, X, SymbolicUtils.unwrap(Y); type = Real)
    )
end
add_gdcprule(Manifolds.distance, SymmetricPositiveDefinite, Positive, GConvex, GAnyMono)

# @register_symbolic LinearAlgebra.exp(X::Union{Symbolics.Arr, Matrix{Num}})
# add_gdcprule(exp, SymmetricPositiveDefinite, Positive, GConvex, GIncreasing)

# add_gdcprule(sqrt, SymmetricPositiveDefinite, Positive, GConvex, GIncreasing)

add_gdcprule(
    SymbolicAnalysis.quad_form,
    SymmetricPositiveDefinite,
    Positive,
    GConvex,
    GIncreasing
)

add_gdcprule(
    LinearAlgebra.eigmax,
    SymmetricPositiveDefinite,
    Positive,
    GConvex,
    GIncreasing
)

"""
    log_quad_form(y, X)
    log_quad_form(ys, X)

Log of the quadratic form of a symmetric positive definite matrix `X` and a vector `y` is defined as `log(y'*X*y)` or for a vector of vectors `ys` as `log(sum(y'*X*y for y in ys))`.

# Arguments

    - `y::Vector`: A vector of `Number`s or a `Vector` of `Vector`s.
    - `X::Matrix`: A symmetric positive definite matrix.
"""
function log_quad_form(y::Vector{<:Number}, X::Matrix)
    return log(y' * X * y)
end

function log_quad_form(ys::Vector{<:Vector}, X::Matrix)
    return log(sum(y' * X * y for y in ys))
end

# See the `sdivergence` note: a `@register_symbolic` two-array registration is
# ambiguous on Symbolics v7. Build the `log_quad_form(y, X)` term directly off the
# symbolic argument, preserving the `(y, X)` argument order (the typical call site,
# e.g. `log_quad_form(x, inv(Σ))`, passes a concrete `y` and a symbolic `X`).
function _log_quad_form_term(y, X)
    return Symbolics.wrap(
        SymbolicUtils.term(
            log_quad_form, SymbolicUtils.unwrap(y), SymbolicUtils.unwrap(X); type = Real
        )
    )
end
log_quad_form(y::AbstractVector, X::Symbolics.Arr) = _log_quad_form_term(y, X)
log_quad_form(y::Symbolics.Arr, X::AbstractMatrix) = _log_quad_form_term(y, X)
log_quad_form(y::Symbolics.Arr, X::Symbolics.Arr) = _log_quad_form_term(y, X)
# `log(y'Xy)` is negative whenever `y'Xy < 1` (e.g. y = e₁, X = 0.01I ⇒ ≈ -4.61).
add_gdcprule(log_quad_form, SymmetricPositiveDefinite, AnySign, GConvex, GIncreasing)

add_gdcprule(inv, SymmetricPositiveDefinite, Positive, GConvex, GDecreasing)

# Matrix `log` (both `Arr` and `Matrix{Num}` forms) is defined in atoms.jl via
# `matrix_atom`; the `@register_array_symbolic` form collides with SymbolicUtils'
# scalar `log` on Symbolics v7 (see the note there).

add_gdcprule(eigsummax, SymmetricPositiveDefinite, Positive, GConvex, GIncreasing)

"""
    schatten_norm(X, p=2)

Schatten norm of a symmetric positive definite matrix `X`.

# Arguments

    - `X::Matrix`: A symmetric positive definite matrix.
    - `p::Int`: The p-norm.
"""
function schatten_norm(X::AbstractMatrix, p::Int = 2)
    return norm(eigvals(X), p)
end

@register_symbolic schatten_norm(X::Matrix{Num}, p::Int)
add_gdcprule(schatten_norm, SymmetricPositiveDefinite, Positive, GConvex, GIncreasing)

"""
    sum_log_eigmax(X, k)
    sum_log_eigmax(f, X, k)

Sum of the log of the maximum eigenvalues of a symmetric positive definite matrix `X`. If a function `f` is provided,
the sum is over `f` applied to the log of the eigenvalues.

# Arguments

    - `f::Function`: A function.
    - `X::Matrix`: A symmetric positive definite matrix.
    - `k::Int`: The number of eigenvalues to consider.
"""
function sum_log_eigmax(f::Function, X::AbstractMatrix, k::Int)
    nrows = size(X, 1)
    eigs = eigvals(X, (nrows - k + 1):nrows)
    return sum(f.(log.(eigs)))
end

@register_symbolic sum_log_eigmax(f::Function, X::Matrix{Num}, k::Int)

function sum_log_eigmax(X::AbstractMatrix, k::Int)
    nrows = size(X, 1)
    eigs = eigvals(X, (nrows - k + 1):nrows)
    return sum((log.(eigs)))
end

@register_symbolic sum_log_eigmax(X::Matrix{Num}, k::Int) false
# Without an outer `f`, this is `sum(log(λᵢ))` over the `k` largest eigenvalues,
# which is negative for λᵢ < 1 (e.g. X = 0.01I₂, k = 2 ⇒ ≈ -9.21).
add_gdcprule(sum_log_eigmax, SymmetricPositiveDefinite, AnySign, GConvex, GIncreasing)

"""
    psd_constant(B)
    pd_constant(B)

Whether `B` is a numeric, finite, real, square, symmetric matrix that is exactly
positive semidefinite (`psd_constant`) or positive definite (`pd_constant`).
Entries must be integers, rationals or floating-point numbers; each is converted
exactly to `Rational{BigInt}` and the test is a pivoted LDLᵀ elimination in exact
arithmetic, so a floating-point matrix that is PSD only up to rounding is refused.
Symbolic matrices return `false`.
"""
psd_constant(B) = exact_semidefinite(B, false)
pd_constant(B) = exact_semidefinite(B, true)

function exact_semidefinite(B, strict::Bool)
    real_constant_matrix(B) || return false
    Bv = map(constval, constval(B))
    all(x -> x isa Union{Integer, Rational, AbstractFloat}, Bv) || return false
    size(Bv, 1) == size(Bv, 2) || return false
    A = Matrix{Rational{BigInt}}(Bv)
    A == transpose(A) || return false
    while !isempty(A)
        d = diag(A)
        any(<(0), d) && return false
        k = argmax(d)
        p = A[k, k]
        # A PSD matrix with a zero diagonal entry has a zero row there.
        iszero(p) && return !strict && iszero(A)
        r = [1:(k - 1); (k + 1):size(A, 1)]
        a = A[r, k]
        A = A[r, r] - a * transpose(a) / p
    end
    return true
end

function real_constant_matrix(B)
    isconstarg(B) || return false
    Bv = constval(B)
    return Bv isa AbstractMatrix && all(x -> (y = constval(x); y isa Real && isfinite(y)), Bv)
end

"""
    affine_map(f, X, B, Y)
    affine_map(f, X, B, Ys)

Affine map, i.e., `B + f(X, Y)` or `B + sum(f(X, Y) for Y in Ys)` for a function `f` where `f` is a positive linear operator.

The gDCP rule (`Positive`, `GConvex`, `GIncreasing`), and the `logdet`/`log(tr(⋅))`
compositions over it, apply only when every constant factor makes the map positive:
`B` is a numeric positive semidefinite matrix (see `psd_constant`), `Y` is a
numeric real matrix for `conjugation`, and the Hadamard factor `Y` is numeric
positive semidefinite for `hadamard_product`. Then `X ↦ B + f(X, Y)` maps the SPD cone
into the PSD cone and is geodesically convex in the Löwner order, and
`logdet(B + f(X, Y))` is geodesically convex (Sra and Hosseini, "Conic geometric
optimization on the manifold of positive definite matrices", SIAM J. Optim. 2015).
Otherwise, including for symbolic factors, the analysis returns `AnySign` and
`GUnknownCurvature` for the map and for every expression containing it.

# Arguments

    - `f::Function`: One of the following functions: `conjugation`, `diag`, `tr` and `hadamard_product`.
    - `X::Matrix`: A symmetric positive definite matrix.
    - `B::Matrix`: A matrix.
    - `Y::Matrix`: A matrix.
    - `Ys::Vector{<:Matrix}`: A vector of matrices.
"""
function affine_map(f::typeof(conjugation), X::Matrix, B::Matrix, Y::Matrix)
    if !(LinearAlgebra.isposdef(B)) || !(eigvals(Symmetric(B), 1:1)[1] >= 0.0)
        throw(DomainError(B, "B must be positive semi-definite."))
    end
    return B + conjugation(X, Y)
end

function affine_map(f::typeof(conjugation), X::Matrix, B::Matrix, Ys::Vector{<:Matrix})
    if !(LinearAlgebra.isposdef(B)) || !(eigvals(Symmetric(B), 1:1)[1] >= 0.0)
        throw(DomainError(B, "B must be positive semi-definite."))
    end
    return B + sum(conjugation(X, Y) for Y in Ys)
end

function affine_map(f::Union{typeof(diag), typeof(tr)}, X::AbstractMatrix, B::AbstractMatrix)
    if !(LinearAlgebra.isposdef(B)) || !(eigvals(Symmetric(B), 1:1)[1] >= 0.0)
        throw(DomainError(B, "B must be positive semi-definite."))
    end
    return B + f(X)
end

add_gdcprule(affine_map, SymmetricPositiveDefinite, Positive, GConvex, GIncreasing)

function constant_factors_ok(::typeof(affine_map), args)
    g = constval(args[1])
    if g === conjugation
        return length(args) == 4 && psd_constant(args[3]) && real_constant_matrix(args[4])
    elseif g === hadamard_product
        return length(args) == 4 && psd_constant(args[3]) && psd_constant(args[4])
    elseif g === LinearAlgebra.tr || g === LinearAlgebra.diag
        return length(args) == 3 && psd_constant(args[3])
    end
    return false
end

"""
    hadamard_product(X, B)

Hadamard product or element-wise multiplication of a symmetric positive definite matrix `X` by a positive semi-definite matrix `B`.

# Arguments

    - `X::Matrix`: A symmetric positive definite matrix.
    - `B::Matrix`: A positive semi-definite matrix.

The gDCP rule (`Positive`, `GConvex`, `GIncreasing`) applies only when `B` is a
numeric positive semidefinite matrix (see `psd_constant`). Then `B ∘ X` is
positive semidefinite (Schur product theorem), `X ↦ B ∘ X` is a positive linear map
and hence geodesically convex in the Löwner order, and `logdet(B ∘ X)` is
geodesically convex (Sra and Hosseini 2015). For an indefinite, non-symmetric or
symbolic `B` none of this holds, and the analysis returns `AnySign` and
`GUnknownCurvature` for the product and for every expression containing it.
"""
function hadamard_product(X::AbstractMatrix, B::AbstractMatrix)
    if (!(LinearAlgebra.isposdef(B)) || !(eigvals(Symmetric(B), 1:1)[1] >= 0.0)) &&
            !(any(prod(r) == 0.0 for r in eachrow(B)))
        throw(DomainError(B, "B must be positive semi-definite and have no zero rows."))
    end
    return B .* X
end

function hadamard_product(X::Symbolics.Arr, B::AbstractMatrix)
    return array_atom_term(hadamard_product, X, B)
end

add_gdcprule(hadamard_product, SymmetricPositiveDefinite, Positive, GConvex, GIncreasing)

constant_factors_ok(::typeof(hadamard_product), args) = psd_constant(args[2])

function affine_map(f::typeof(hadamard_product), X::Matrix, Y::Matrix, B::Matrix)
    if !(LinearAlgebra.isposdef(B)) || !(eigvals(Symmetric(B), 1:1)[1] >= 0.0)
        throw(DomainError(B, "B must be positive semi-definite."))
    end
    return B + hadamard_product(X, Y)
end

# Symbolic `affine_map` over a symbolic matrix `X` must stay an unevaluated
# `affine_map` term (the gDCP rule and curvature pass dispatch on
# `operation(ex) == affine_map`). The `@register_array_symbolic` macro cannot
# express this on Symbolics v7: the function-typed first argument plus the
# `Matrix{Num}` argument collapse so the conjugation- and hadamard-flavoured
# registrations generate the same `BasicSymbolic{SymReal}` method and collide.
# Build the term directly off the symbolic `X` instead.
function affine_map_term(f, X, args...)
    x = SymbolicUtils.unwrap(X)
    uargs = map(SymbolicUtils.unwrap, args)
    return Symbolics.wrap(SymbolicUtils.term(affine_map, f, x, uargs...; type = Matrix{Real}))
end
function affine_map(f::typeof(conjugation), X::Symbolics.Arr, B::AbstractMatrix, Y::AbstractMatrix)
    return affine_map_term(f, X, B, Y)
end
function affine_map(f::Union{typeof(diag), typeof(tr)}, X::Symbolics.Arr, B::AbstractMatrix)
    return affine_map_term(f, X, B)
end
function affine_map(f::typeof(hadamard_product), X::Symbolics.Arr, Y::AbstractMatrix, B::AbstractMatrix)
    return affine_map_term(f, X, Y, B)
end
