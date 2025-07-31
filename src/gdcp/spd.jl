### DGCP Atoms

@register_symbolic LinearAlgebra.logdet(X::Matrix{Num})
add_gdcprule(
    LinearAlgebra.logdet,
    SymmetricPositiveDefinite,
    Positive,
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

@register_array_symbolic conjugation(X::Union{Symbolics.Arr, Matrix{Num}}, B::Matrix) begin
    size = (size(B, 2), size(B, 2))
end

add_gdcprule(conjugation, SymmetricPositiveDefinite, Positive, GConvex, GIncreasing)

@register_symbolic LinearAlgebra.tr(X::Union{Symbolics.Arr, Matrix{Num}})
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

@register_symbolic scalar_mat(X::Union{Symbolics.Arr, Matrix{Num}}, k::Int)

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
"""
function sdivergence(X, Y)
    return logdet((X + Y) / 2) - 1 / 2 * logdet(X * Y)
end

@register_symbolic sdivergence(X::Matrix{Num}, Y::Matrix)
add_gdcprule(sdivergence, SymmetricPositiveDefinite, Positive, GConvex, GIncreasing)

@register_symbolic Manifolds.distance(
    M::Manifolds.SymmetricPositiveDefinite,
    X::AbstractMatrix,
    Y::Union{Symbolics.Arr, Matrix{Num}}
)
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

@register_symbolic log_quad_form(y::Vector, X::Matrix{Num})
add_gdcprule(log_quad_form, SymmetricPositiveDefinite, Positive, GConvex, GIncreasing)

add_gdcprule(inv, SymmetricPositiveDefinite, Positive, GConvex, GDecreasing)

add_gdcprule(diag, SymmetricPositiveDefinite, Positive, GConvex, GIncreasing)

@register_array_symbolic Base.log(X::Matrix{Num}) begin
    size = (size(X, 1), size(X, 2))
end

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
add_gdcprule(sum_log_eigmax, SymmetricPositiveDefinite, Positive, GConvex, GIncreasing)

"""
    affine_map(f, X, B, Y)
    affine_map(f, X, B, Ys)

Affine map, i.e., `B + f(X, Y)` or `B + sum(f(X, Y) for Y in Ys)` for a function `f` where `f` is a positive linear operator.

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

@register_array_symbolic affine_map(
    conjf::typeof(conjugation),
    X::Matrix{Num},
    B::Matrix,
    Y::Union{Matrix, Vector{<:Matrix}}
) begin
    size = (size(B, 1), size(B, 2))
end

function affine_map(f::Union{typeof(diag), typeof(tr)}, X::AbstractMatrix, B::AbstractMatrix)
    if !(LinearAlgebra.isposdef(B)) || !(eigvals(Symmetric(B), 1:1)[1] >= 0.0)
        throw(DomainError(B, "B must be positive semi-definite."))
    end
    return B + f(X)
end

@register_array_symbolic affine_map(
    diagtrf::Union{typeof(diag), typeof(tr)},
    X::Matrix{Num},
    B::Matrix
) begin
    size = (size(B, 1), size(B, 2))
end false

add_gdcprule(affine_map, SymmetricPositiveDefinite, Positive, GConvex, GIncreasing)

"""
    hadamard_product(X, B)

Hadamard product or element-wise multiplication of a symmetric positive definite matrix `X` by a positive semi-definite matrix `B`.

# Arguments

    - `X::Matrix`: A symmetric positive definite matrix.
    - `B::Matrix`: A positive semi-definite matrix.
"""
function hadamard_product(X::AbstractMatrix, B::AbstractMatrix)
    if (!(LinearAlgebra.isposdef(B)) || !(eigvals(Symmetric(B), 1:1)[1] >= 0.0)) &&
       !(any(prod(r) == 0.0 for r in eachrow(B)))
        throw(DomainError(B, "B must be positive semi-definite and have no zero rows."))
    end
    return B .* X
end

@register_array_symbolic hadamard_product(X::Matrix{Num}, B::Matrix) begin
    size = (size(B, 1), size(B, 2))
end

add_gdcprule(hadamard_product, SymmetricPositiveDefinite, Positive, GConvex, GIncreasing)

function affine_map(f::typeof(hadamard_product), X::Matrix, Y::Matrix, B::Matrix)
    if !(LinearAlgebra.isposdef(B)) || !(eigvals(Symmetric(B), 1:1)[1] >= 0.0)
        throw(DomainError(B, "B must be positive semi-definite."))
    end
    return B + hadamard_product(X, Y)
end

@register_array_symbolic affine_map(
    hadamard_product::typeof(hadamard_product),
    X::Matrix{Num},
    Y::Matrix,
    B::Matrix
) begin
    size = (size(B, 1), size(B, 2))
end false
