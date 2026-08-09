# DGCP Atoms for Lorentz Model (Hyperbolic Space)
#
# This file implements geodesically convex atoms for the Lorentz model,
# a Cartan-Hadamard manifold of constant negative curvature.
# Based on the results from Ferreira, Németh, and Zhu (2022, 2023).

using Manifolds
using LinearAlgebra
using Symbolics: @register_symbolic

# See the SPD `distance` note in gdcp/spd.jl: build the term directly off the
# symbolic point so the SPD and Lorentz methods do not collide.
function Manifolds.distance(
        M::Manifolds.Lorentz,
        p::AbstractVector,
        q::Symbolics.Arr
    )
    return Symbolics.wrap(
        SymbolicUtils.term(Manifolds.distance, M, p, SymbolicUtils.unwrap(q); type = Real)
    )
end
add_gdcprule(Manifolds.distance, Manifolds.Lorentz, Positive, GConvex, GAnyMono)

"""
    lorentz_log_barrier(p) -> Real

Evaluate the log-barrier `-log(-1 - <a, p>_L)` for the Lorentz model, where
`a = (0, ..., 0, 1)` and `<., .>_L` is the Lorentzian inner product.

# Arguments

- `p::AbstractVector`: point on the Lorentz manifold. Its last coordinate must
  be greater than one for the barrier to be finite and real.

# Returns

- A real scalar containing the barrier value.

# Throws

- `DomainError`: if the last coordinate of `p` is outside the real logarithm's
  domain.

# Examples

```jldoctest
julia> using SymbolicAnalysis

julia> lorentz_log_barrier([0.0, 2.0]) == 0.0
true
```
"""
function lorentz_log_barrier(p::AbstractVector)
    # a = (0, ..., 0, 1), so the Lorentzian inner product
    # <a, p>_L = a1*p1 + ... + a_d*p_d - a_{d+1}*p_{d+1} reduces to -p[end].
    # The barrier is -log(-1 - <a, p>_L) = -log(-1 + p[end]).
    return -log(-1 + p[end])
end

@register_symbolic lorentz_log_barrier(p::Vector{Num})
add_gdcprule(lorentz_log_barrier, Manifolds.Lorentz, Positive, GConvex, GIncreasing)

"""
    lorentz_homogeneous_quadratic(A, p) -> Real

Evaluate the homogeneous quadratic `transpose(p) * A * p` on the Lorentz
model. The matrix must satisfy one of the implemented geodesic-convexity
conditions.

# Arguments

- `A::AbstractMatrix`: symmetric `(d + 1)` by `(d + 1)` coefficient matrix.
- `p::AbstractVector`: point on the Lorentz manifold.

# Returns

- A scalar containing the quadratic value.

# Throws

- `ArgumentError`: if `A` does not satisfy the geodesic-convexity conditions.

# Examples

```jldoctest
julia> using SymbolicAnalysis, LinearAlgebra

julia> A = Matrix{Float64}(I, 3, 3);

julia> lorentz_homogeneous_quadratic(A, [0.0, 0.0, 1.0])
1.0
```
"""
function lorentz_homogeneous_quadratic(A::AbstractMatrix, p::AbstractVector)
    d = size(A, 1) - 1

    # Extract the components from matrix A
    A_bar = A[1:d, 1:d]
    a_vec = A[1:d, d + 1]
    sigma = A[d + 1, d + 1]

    # Compute the minimum eigenvalue of A_bar
    lambda_min = minimum(eigvals(A_bar))

    # Check conditions from Theorem 21
    condition1 = isapprox(norm(a_vec), 0, atol = 1.0e-10) && sigma >= -lambda_min
    condition2 = sigma + lambda_min > 2 * sqrt(dot(a_vec, a_vec))

    if !(condition1 || condition2)
        throw(ArgumentError("Matrix A does not satisfy geodesic convexity conditions"))
    end

    return p' * A * p
end

@register_symbolic lorentz_homogeneous_quadratic(
    A::AbstractMatrix,
    p::Vector{Num}
)
add_gdcprule(lorentz_homogeneous_quadratic, Manifolds.Lorentz, Positive, GConvex, GAnyMono)

"""
    lorentz_homogeneous_diagonal(a, p) -> Real

Evaluate the diagonal quadratic `sum(a[i] * p[i]^2)` on the Lorentz model.
Geodesic convexity requires `minimum(a[1:end-1]) + a[end] >= 0`.

# Arguments

- `a::AbstractVector`: `(d + 1)` diagonal coefficients.
- `p::AbstractVector`: point on the Lorentz manifold.

# Returns

- A scalar containing the diagonal quadratic value.

# Throws

- `DimensionMismatch`: if `a` and `p` have different lengths.
- `ArgumentError`: if `a` does not satisfy the geodesic-convexity condition.

# Examples

```jldoctest
julia> using SymbolicAnalysis

julia> lorentz_homogeneous_diagonal([1.0, 2.0, 0.0], [0.0, 0.0, 1.0])
0.0
```
"""
function lorentz_homogeneous_diagonal(a::AbstractVector, p::AbstractVector)
    if length(a) != length(p)
        throw(DimensionMismatch("Vectors must have same length"))
    end

    if minimum(a[1:(end - 1)]) + a[end] < 0
        throw(
            ArgumentError(
                "For geodesic convexity, min(a[1:end-1]) + a[end] ≥ 0 is required",
            ),
        )
    end

    return sum(a .* p .^ 2)
end

@register_symbolic lorentz_homogeneous_diagonal(
    a::AbstractVector,
    p::Vector{Num}
)
add_gdcprule(lorentz_homogeneous_diagonal, Manifolds.Lorentz, Positive, GConvex, GAnyMono)

"""
    lorentz_nonhomogeneous_quadratic(A::AbstractMatrix, b::AbstractVector, c::Real, p::AbstractVector)

Computes the non-homogeneous quadratic function f(p) = p'Ap + b'p + c on the Lorentz model.
For geodesic convexity, p'Ap must be geodesically convex and b must be in the Lorentz cone L.

# Arguments

    - `A::AbstractMatrix`: A symmetric matrix in R^((d+1)×(d+1)).
    - `b::AbstractVector`: A vector in R^(d+1) which must be in the Lorentz cone.
    - `c::Real`: A constant term.
    - `p::AbstractVector`: A point on the Lorentz manifold.
"""
function lorentz_nonhomogeneous_quadratic(
        A::AbstractMatrix,
        b::AbstractVector,
        c::Real,
        p::AbstractVector
    )
    # Check if b is in the Lorentz cone
    b_head = b[1:(end - 1)]
    b_tail = b[end]

    if !(norm(b_head)^2 <= b_tail^2 && b_tail >= 0)
        throw(ArgumentError("Vector b must be in the Lorentz cone for geodesic convexity"))
    end

    # This call will check if A satisfies the geodesic convexity conditions
    homogeneous_part = lorentz_homogeneous_quadratic(A, p)
    affine_part = (Matrix(b') * p)
    return homogeneous_part + affine_part[1] + c
end

@register_symbolic lorentz_nonhomogeneous_quadratic(
    A::AbstractMatrix,
    b::AbstractVector,
    c::Real,
    p::Vector{Num}
)
add_gdcprule(lorentz_nonhomogeneous_quadratic, Manifolds.Lorentz, AnySign, GConvex, AnyMono)

"""
    lorentz_least_squares(X, y, p) -> Real

Evaluate the squared residual norm `sum(abs2, y - X * p)` on the Lorentz
model. The derived quadratic and linear terms must satisfy the implemented
geodesic-convexity conditions.

# Arguments

- `X::AbstractMatrix`: design matrix with `d + 1` columns.
- `y::AbstractVector`: response vector with one entry per row of `X`.
- `p::AbstractVector`: point on the Lorentz manifold.

# Returns

- A scalar containing the squared residual norm.

# Throws

- `ArgumentError`: if the derived homogeneous or linear term does not satisfy
  the geodesic-convexity conditions.
- `DimensionMismatch`: if the dimensions of `X`, `y`, and `p` are incompatible.

# Examples

```jldoctest
julia> using SymbolicAnalysis, LinearAlgebra

julia> X = Matrix{Float64}(I, 3, 3);

julia> lorentz_least_squares(X, [0.0, 0.0, -1.0], [0.0, 0.0, 1.0])
4.0
```
"""
function lorentz_least_squares(X::AbstractMatrix, y::AbstractVector, p::AbstractVector)
    A = X' * X      # Homogeneous quadratic coefficient
    b = -2 * X' * y # Linear coefficient
    c = y' * y      # Constant term

    # This call will check the geodesic convexity conditions for both
    # the homogeneous part (via lorentz_homogeneous_quadratic) and the linear term
    return lorentz_nonhomogeneous_quadratic(A, b, c, p)
end

# A `@register_symbolic` three-array registration explodes into a combinatorial,
# mutually-ambiguous set of wrapper methods on Symbolics v7 (Aqua flags them; see
# the `sdivergence` note in gdcp/spd.jl). Build the `lorentz_least_squares(X, y, p)`
# term directly off the symbolic point `p` instead (the gDCP pass only needs
# `operation(ex) == lorentz_least_squares`).
function lorentz_least_squares(X::AbstractMatrix, y::AbstractVector, p::Symbolics.Arr)
    return Symbolics.wrap(
        SymbolicUtils.term(
            lorentz_least_squares, X, y, SymbolicUtils.unwrap(p); type = Real
        )
    )
end
add_gdcprule(lorentz_least_squares, Manifolds.Lorentz, Positive, GConvex, AnyMono)

"""
    lorentz_transform(O, p) -> AbstractVector

Apply the orthochronous Lorentz transformation `O` to the point `p`.

# Arguments

- `O::AbstractMatrix`: element of the orthochronous Lorentz group.
- `p::AbstractVector`: point on the Lorentz manifold.

# Returns

- The transformed point `O * p`.

# Throws

- `ArgumentError`: if `O` does not preserve the Lorentz metric or the positive
  time direction.

# Examples

```jldoctest
julia> using SymbolicAnalysis, LinearAlgebra

julia> O = Matrix{Float64}(I, 3, 3);

julia> lorentz_transform(O, [0.0, 0.0, 1.0])
3-element Vector{Float64}:
 0.0
 0.0
 1.0
```
"""
function lorentz_transform(O::AbstractMatrix, p::AbstractVector)
    d = length(p) - 1
    J = Diagonal([ones(d)..., -1])

    # Check if O is in the Lorentz group
    if !isapprox(O' * J * O, J, rtol = 1.0e-10)
        throw(ArgumentError("Matrix is not in the Lorentz group"))
    end

    # Check if O preserves the positive time direction (orthochronous)
    if (O * [zeros(d)..., 1])[end] <= 0
        throw(ArgumentError("Matrix does not preserve the positive time direction"))
    end

    return O * p
end

@register_symbolic lorentz_transform(
    O::AbstractMatrix,
    p::Vector{Num}
)
# Not adding a rule since this preserves geodesic convexity but doesn't have a specific curvature

# Export functions
export lorentz_log_barrier, lorentz_homogeneous_quadratic
export lorentz_homogeneous_diagonal, lorentz_least_squares, lorentz_transform
