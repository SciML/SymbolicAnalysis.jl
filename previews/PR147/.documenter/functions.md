---
---

# Special functions {#Special-functions}

Since some atoms are not available in the base language or other packages we have implemented them here.
<details class='jldocstring custom-block' open>
<summary><a id='SymbolicAnalysis.dotsort-Tuple{AbstractVector, AbstractVector}' href='#SymbolicAnalysis.dotsort-Tuple{AbstractVector, AbstractVector}'><span class="jlbinding">SymbolicAnalysis.dotsort</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
dotsort(x, y)
```


Sorts `x` and `y` and returns the dot product of the sorted vectors.

**Arguments**

```julia
- `x::AbstractVector`: A vector.
- `y::AbstractVector`: A vector.
```



<Badge type="info" class="source-link" text="source"><a href="https://github.com/SciML/SymbolicAnalysis.jl" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='SymbolicAnalysis.eigsummax-Tuple{LinearAlgebra.Symmetric, Int64}' href='#SymbolicAnalysis.eigsummax-Tuple{LinearAlgebra.Symmetric, Int64}'><span class="jlbinding">SymbolicAnalysis.eigsummax</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
eigsummax(m::Symmetric, k)
```


Returns the sum of the `k` largest eigenvalues of `m`.

**Arguments**

```julia
- `m::Symmetric`: A symmetric matrix.
- `k::Int`: The number of largest eigenvalues to sum.
```



<Badge type="info" class="source-link" text="source"><a href="https://github.com/SciML/SymbolicAnalysis.jl" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='SymbolicAnalysis.eigsummin-Tuple{LinearAlgebra.Symmetric, Int64}' href='#SymbolicAnalysis.eigsummin-Tuple{LinearAlgebra.Symmetric, Int64}'><span class="jlbinding">SymbolicAnalysis.eigsummin</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
eigsummin(m::Symmetric, k)
```


Returns the sum of the `k` smallest eigenvalues of `m`.

**Arguments**

```julia
- `m::Symmetric`: A symmetric matrix.
- `k::Int`: The number of smallest eigenvalues to sum.
```



<Badge type="info" class="source-link" text="source"><a href="https://github.com/SciML/SymbolicAnalysis.jl" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='SymbolicAnalysis.huber' href='#SymbolicAnalysis.huber'><span class="jlbinding">SymbolicAnalysis.huber</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
huber(x, M=1)
```


Returns the Huber loss function of `x` with threshold `M`.

**Arguments**

```julia
- `x::Real`: A Real.
- `M::Real`: The threshold.
```



<Badge type="info" class="source-link" text="source"><a href="https://github.com/SciML/SymbolicAnalysis.jl" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='SymbolicAnalysis.invprod-Tuple{AbstractVector}' href='#SymbolicAnalysis.invprod-Tuple{AbstractVector}'><span class="jlbinding">SymbolicAnalysis.invprod</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
invprod(x::AbstractVector)
```


Returns the inverse of the product of the elements of `x`.

**Arguments**

```julia
- `x::AbstractVector`: A vector.
```



<Badge type="info" class="source-link" text="source"><a href="https://github.com/SciML/SymbolicAnalysis.jl" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='SymbolicAnalysis.lognormcdf-Tuple{Real}' href='#SymbolicAnalysis.lognormcdf-Tuple{Real}'><span class="jlbinding">SymbolicAnalysis.lognormcdf</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
lognormcdf(x::Real)
```


Returns the log of the normal cumulative distribution function of `x`.

**Arguments**

```julia
- `x::Real`: A Real.
```



<Badge type="info" class="source-link" text="source"><a href="https://github.com/SciML/SymbolicAnalysis.jl" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='SymbolicAnalysis.matrix_frac-Tuple{AbstractVector, AbstractMatrix}' href='#SymbolicAnalysis.matrix_frac-Tuple{AbstractVector, AbstractMatrix}'><span class="jlbinding">SymbolicAnalysis.matrix_frac</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
matrix_frac(x::AbstractVector, P::AbstractMatrix)
```


Returns the quadratic form `x' * P^{-1} * x`.

**Arguments**

```julia
- `x::AbstractVector`: A vector.
- `P::AbstractMatrix`: A matrix.
```



<Badge type="info" class="source-link" text="source"><a href="https://github.com/SciML/SymbolicAnalysis.jl" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='SymbolicAnalysis.perspective-Tuple{Function, Any, Real}' href='#SymbolicAnalysis.perspective-Tuple{Function, Any, Real}'><span class="jlbinding">SymbolicAnalysis.perspective</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
perspective(f::Function, x, s::Real)
```


Returns the perspective function `s * f(x / s)`.

**Arguments**

```julia
- `f::Function`: A function.
- `x`: A Real.
- `s::Real`: A positive Real.
```



<Badge type="info" class="source-link" text="source"><a href="https://github.com/SciML/SymbolicAnalysis.jl" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='SymbolicAnalysis.quad_form-Tuple{AbstractVector, AbstractMatrix}' href='#SymbolicAnalysis.quad_form-Tuple{AbstractVector, AbstractMatrix}'><span class="jlbinding">SymbolicAnalysis.quad_form</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
quad_form(x::AbstractVector, P::AbstractMatrix)
```


Returns the quadratic form `x' * P * x`.

**Arguments**

```julia
- `x::AbstractVector`: A vector.
- `P::AbstractMatrix`: A matrix.
```



<Badge type="info" class="source-link" text="source"><a href="https://github.com/SciML/SymbolicAnalysis.jl" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='SymbolicAnalysis.quad_over_lin-Tuple{Real, Real}' href='#SymbolicAnalysis.quad_over_lin-Tuple{Real, Real}'><span class="jlbinding">SymbolicAnalysis.quad_over_lin</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
quad_over_lin(x::Real, y::Real)
```


Returns the quadratic over linear form `x^2 / y`.

**Arguments**

```julia
- `x`: A Real or a vector.
- `y::Real`: A positive Real.
```



<Badge type="info" class="source-link" text="source"><a href="https://github.com/SciML/SymbolicAnalysis.jl" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='SymbolicAnalysis.sum_largest-Tuple{AbstractMatrix, Integer}' href='#SymbolicAnalysis.sum_largest-Tuple{AbstractMatrix, Integer}'><span class="jlbinding">SymbolicAnalysis.sum_largest</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
sum_largest(x::AbstractMatrix, k)
```


Returns the sum of the `k` largest elements of `x`.

**Arguments**

```julia
- `x::AbstractMatrix`: A matrix.
- `k::Int`: The number of largest elements to sum.
```



<Badge type="info" class="source-link" text="source"><a href="https://github.com/SciML/SymbolicAnalysis.jl" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='SymbolicAnalysis.sum_smallest-Tuple{AbstractMatrix, Integer}' href='#SymbolicAnalysis.sum_smallest-Tuple{AbstractMatrix, Integer}'><span class="jlbinding">SymbolicAnalysis.sum_smallest</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
sum_smallest(x::AbstractMatrix, k)
```


Returns the sum of the `k` smallest elements of `x`.

**Arguments**

```julia
- `x::AbstractMatrix`: A matrix.
- `k::Int`: The number of smallest elements to sum.
```



<Badge type="info" class="source-link" text="source"><a href="https://github.com/SciML/SymbolicAnalysis.jl" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='SymbolicAnalysis.trinv-Tuple{AbstractMatrix}' href='#SymbolicAnalysis.trinv-Tuple{AbstractMatrix}'><span class="jlbinding">SymbolicAnalysis.trinv</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
trinv(x::AbstractMatrix)
```


Returns the trace of the inverse of `x`.

**Arguments**

```julia
- `x::AbstractMatrix`: A matrix.
```



<Badge type="info" class="source-link" text="source"><a href="https://github.com/SciML/SymbolicAnalysis.jl" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='SymbolicAnalysis.tv-Tuple{AbstractVector{<:AbstractMatrix}}' href='#SymbolicAnalysis.tv-Tuple{AbstractVector{<:AbstractMatrix}}'><span class="jlbinding">SymbolicAnalysis.tv</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
tv(x::AbstractVector{<:AbstractMatrix})
```


Returns the total variation of `x`, defined as `sum_{i,j} |x_{k+1}[i,j] - x_k[i,j]|`.

**Arguments**

```julia
- `x::AbstractVector`: A vector of matrices.
```



<Badge type="info" class="source-link" text="source"><a href="https://github.com/SciML/SymbolicAnalysis.jl" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='SymbolicAnalysis.tv-Tuple{AbstractVector{<:Real}}' href='#SymbolicAnalysis.tv-Tuple{AbstractVector{<:Real}}'><span class="jlbinding">SymbolicAnalysis.tv</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
tv(x::AbstractVector{<:Real})
```


Returns the total variation of `x`, defined as `sum_i |x_{i+1} - x_i|`.

**Arguments**

```julia
- `x::AbstractVector`: A vector.
```



<Badge type="info" class="source-link" text="source"><a href="https://github.com/SciML/SymbolicAnalysis.jl" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='SymbolicAnalysis.affine_map-Tuple{typeofSymbolicAnalysis.conjugation, Matrix, Matrix, Matrix}' href='#SymbolicAnalysis.affine_map-Tuple{typeofSymbolicAnalysis.conjugation, Matrix, Matrix, Matrix}'><span class="jlbinding">SymbolicAnalysis.affine_map</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
affine_map(f, X, B, Y)
affine_map(f, X, B, Ys)
```


Affine map, i.e., `B + f(X, Y)` or `B + sum(f(X, Y) for Y in Ys)` for a function `f` where `f` is a positive linear operator.

**Arguments**

```julia
- `f::Function`: One of the following functions: `conjugation`, `diag`, `tr` and `hadamard_product`.
- `X::Matrix`: A symmetric positive definite matrix.
- `B::Matrix`: A matrix.
- `Y::Matrix`: A matrix.
- `Ys::Vector{<:Matrix}`: A vector of matrices.
```



<Badge type="info" class="source-link" text="source"><a href="https://github.com/SciML/SymbolicAnalysis.jl" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='SymbolicAnalysis.conjugation-Tuple{Any, Any}' href='#SymbolicAnalysis.conjugation-Tuple{Any, Any}'><span class="jlbinding">SymbolicAnalysis.conjugation</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
conjugation(X, B)
```


Conjugation of a matrix `X` by a matrix `B` is defined as `B'X*B`.

**Arguments**

```julia
- `X::Matrix`: A symmetric positive definite matrix.
- `B::Matrix`: A matrix.
```



<Badge type="info" class="source-link" text="source"><a href="https://github.com/SciML/SymbolicAnalysis.jl" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='SymbolicAnalysis.hadamard_product-Tuple{AbstractMatrix, AbstractMatrix}' href='#SymbolicAnalysis.hadamard_product-Tuple{AbstractMatrix, AbstractMatrix}'><span class="jlbinding">SymbolicAnalysis.hadamard_product</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
hadamard_product(X, B)
```


Hadamard product or element-wise multiplication of a symmetric positive definite matrix `X` by a positive semi-definite matrix `B`.

**Arguments**

```julia
- `X::Matrix`: A symmetric positive definite matrix.
- `B::Matrix`: A positive semi-definite matrix.
```



<Badge type="info" class="source-link" text="source"><a href="https://github.com/SciML/SymbolicAnalysis.jl" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='SymbolicAnalysis.log_quad_form-Tuple{Vector{<:Number}, Matrix}' href='#SymbolicAnalysis.log_quad_form-Tuple{Vector{<:Number}, Matrix}'><span class="jlbinding">SymbolicAnalysis.log_quad_form</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
log_quad_form(y, X)
log_quad_form(ys, X)
```


Log of the quadratic form of a symmetric positive definite matrix `X` and a vector `y` is defined as `log(y'*X*y)` or for a vector of vectors `ys` as `log(sum(y'*X*y for y in ys))`.

**Arguments**

```julia
- `y::Vector`: A vector of `Number`s or a `Vector` of `Vector`s.
- `X::Matrix`: A symmetric positive definite matrix.
```



<Badge type="info" class="source-link" text="source"><a href="https://github.com/SciML/SymbolicAnalysis.jl" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='SymbolicAnalysis.scalar_mat' href='#SymbolicAnalysis.scalar_mat'><span class="jlbinding">SymbolicAnalysis.scalar_mat</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
scalar_mat(X, k=size(X, 1))
```


Scalar matrix of a symmetric positive definite matrix `X` is defined as `tr(X)*I(k)`.

**Arguments**

```julia
- `X::Matrix`: A symmetric positive definite matrix.
- `k::Int`: The size of the identity matrix.
```



<Badge type="info" class="source-link" text="source"><a href="https://github.com/SciML/SymbolicAnalysis.jl" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='SymbolicAnalysis.schatten_norm' href='#SymbolicAnalysis.schatten_norm'><span class="jlbinding">SymbolicAnalysis.schatten_norm</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
schatten_norm(X, p=2)
```


Schatten norm of a symmetric positive definite matrix `X`.

**Arguments**

```julia
- `X::Matrix`: A symmetric positive definite matrix.
- `p::Int`: The p-norm.
```



<Badge type="info" class="source-link" text="source"><a href="https://github.com/SciML/SymbolicAnalysis.jl" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='SymbolicAnalysis.sdivergence-Tuple{Any, Any}' href='#SymbolicAnalysis.sdivergence-Tuple{Any, Any}'><span class="jlbinding">SymbolicAnalysis.sdivergence</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
sdivergence(X, Y)
```


Symmetric divergence of two symmetric positive definite matrices `X` and `Y` is defined as `logdet((X+Y)/2) - 1/2*logdet(X*Y)`.

**Arguments**

```julia
- `X::Matrix`: A symmetric positive definite matrix.
- `Y::Matrix`: A symmetric positive definite matrix.
```



<Badge type="info" class="source-link" text="source"><a href="https://github.com/SciML/SymbolicAnalysis.jl" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='SymbolicAnalysis.sum_log_eigmax-Tuple{Function, AbstractMatrix, Int64}' href='#SymbolicAnalysis.sum_log_eigmax-Tuple{Function, AbstractMatrix, Int64}'><span class="jlbinding">SymbolicAnalysis.sum_log_eigmax</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
sum_log_eigmax(X, k)
sum_log_eigmax(f, X, k)
```


Sum of the log of the maximum eigenvalues of a symmetric positive definite matrix `X`. If a function `f` is provided, the sum is over `f` applied to the log of the eigenvalues.

**Arguments**

```julia
- `f::Function`: A function.
- `X::Matrix`: A symmetric positive definite matrix.
- `k::Int`: The number of eigenvalues to consider.
```



<Badge type="info" class="source-link" text="source"><a href="https://github.com/SciML/SymbolicAnalysis.jl" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='SymbolicAnalysis.lorentz_homogeneous_diagonal-Tuple{AbstractVector, AbstractVector}' href='#SymbolicAnalysis.lorentz_homogeneous_diagonal-Tuple{AbstractVector, AbstractVector}'><span class="jlbinding">SymbolicAnalysis.lorentz_homogeneous_diagonal</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
lorentz_homogeneous_diagonal(a, p) -> Real
```


Evaluate the diagonal quadratic `sum(a[i] * p[i]^2)` on the Lorentz model. Geodesic convexity requires `minimum(a[1:end-1]) + a[end] >= 0`.

**Arguments**
- `a::AbstractVector`: `(d + 1)` diagonal coefficients.
  
- `p::AbstractVector`: point on the Lorentz manifold.
  

**Returns**
- A scalar containing the diagonal quadratic value.
  

**Throws**
- `DimensionMismatch`: if `a` and `p` have different lengths.
  
- `ArgumentError`: if `a` does not satisfy the geodesic-convexity condition.
  

**Examples**

```julia
julia> using SymbolicAnalysis

julia> lorentz_homogeneous_diagonal([1.0, 2.0, 0.0], [0.0, 0.0, 1.0])
0.0
```



<Badge type="info" class="source-link" text="source"><a href="https://github.com/SciML/SymbolicAnalysis.jl" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='SymbolicAnalysis.lorentz_homogeneous_quadratic-Tuple{AbstractMatrix, AbstractVector}' href='#SymbolicAnalysis.lorentz_homogeneous_quadratic-Tuple{AbstractMatrix, AbstractVector}'><span class="jlbinding">SymbolicAnalysis.lorentz_homogeneous_quadratic</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
lorentz_homogeneous_quadratic(A, p) -> Real
```


Evaluate the homogeneous quadratic `transpose(p) * A * p` on the Lorentz model. The matrix must satisfy one of the implemented geodesic-convexity conditions.

**Arguments**
- `A::AbstractMatrix`: symmetric `(d + 1)` by `(d + 1)` coefficient matrix.
  
- `p::AbstractVector`: point on the Lorentz manifold.
  

**Returns**
- A scalar containing the quadratic value.
  

**Throws**
- `ArgumentError`: if `A` does not satisfy the geodesic-convexity conditions.
  

**Examples**

```julia
julia> using SymbolicAnalysis, LinearAlgebra

julia> A = Matrix{Float64}(I, 3, 3);

julia> lorentz_homogeneous_quadratic(A, [0.0, 0.0, 1.0])
1.0
```



<Badge type="info" class="source-link" text="source"><a href="https://github.com/SciML/SymbolicAnalysis.jl" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='SymbolicAnalysis.lorentz_least_squares-Tuple{AbstractMatrix, AbstractVector, AbstractVector}' href='#SymbolicAnalysis.lorentz_least_squares-Tuple{AbstractMatrix, AbstractVector, AbstractVector}'><span class="jlbinding">SymbolicAnalysis.lorentz_least_squares</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
lorentz_least_squares(X, y, p) -> Real
```


Evaluate the squared residual norm `sum(abs2, y - X * p)` on the Lorentz model. The derived quadratic and linear terms must satisfy the implemented geodesic-convexity conditions.

**Arguments**
- `X::AbstractMatrix`: design matrix with `d + 1` columns.
  
- `y::AbstractVector`: response vector with one entry per row of `X`.
  
- `p::AbstractVector`: point on the Lorentz manifold.
  

**Returns**
- A scalar containing the squared residual norm.
  

**Throws**
- `ArgumentError`: if the derived homogeneous or linear term does not satisfy the geodesic-convexity conditions.
  
- `DimensionMismatch`: if the dimensions of `X`, `y`, and `p` are incompatible.
  

**Examples**

```julia
julia> using SymbolicAnalysis, LinearAlgebra

julia> X = Matrix{Float64}(I, 3, 3);

julia> lorentz_least_squares(X, [0.0, 0.0, -1.0], [0.0, 0.0, 1.0])
4.0
```



<Badge type="info" class="source-link" text="source"><a href="https://github.com/SciML/SymbolicAnalysis.jl" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='SymbolicAnalysis.lorentz_log_barrier-Tuple{AbstractVector}' href='#SymbolicAnalysis.lorentz_log_barrier-Tuple{AbstractVector}'><span class="jlbinding">SymbolicAnalysis.lorentz_log_barrier</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
lorentz_log_barrier(p) -> Real
```


Evaluate the log-barrier `-log(-1 - <a, p>_L)` for the Lorentz model, where `a = (0, ..., 0, 1)` and `<., .>_L` is the Lorentzian inner product.

**Arguments**
- `p::AbstractVector`: point on the Lorentz manifold. Its last coordinate must be greater than one for the barrier to be finite and real.
  

**Returns**
- A real scalar containing the barrier value.
  

**Throws**
- `DomainError`: if the last coordinate of `p` is outside the real logarithm's domain.
  

**Examples**

```julia
julia> using SymbolicAnalysis

julia> lorentz_log_barrier([0.0, 2.0]) == 0.0
true
```



<Badge type="info" class="source-link" text="source"><a href="https://github.com/SciML/SymbolicAnalysis.jl" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='SymbolicAnalysis.lorentz_nonhomogeneous_quadratic-Tuple{AbstractMatrix, AbstractVector, Real, AbstractVector}' href='#SymbolicAnalysis.lorentz_nonhomogeneous_quadratic-Tuple{AbstractMatrix, AbstractVector, Real, AbstractVector}'><span class="jlbinding">SymbolicAnalysis.lorentz_nonhomogeneous_quadratic</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
lorentz_nonhomogeneous_quadratic(A::AbstractMatrix, b::AbstractVector, c::Real, p::AbstractVector)
```


Computes the non-homogeneous quadratic function f(p) = p'Ap + b'p + c on the Lorentz model. For geodesic convexity, p'Ap must be geodesically convex and b must be in the Lorentz cone L.

**Arguments**

```julia
- `A::AbstractMatrix`: A symmetric matrix in R^((d+1)×(d+1)).
- `b::AbstractVector`: A vector in R^(d+1) which must be in the Lorentz cone.
- `c::Real`: A constant term.
- `p::AbstractVector`: A point on the Lorentz manifold.
```



<Badge type="info" class="source-link" text="source"><a href="https://github.com/SciML/SymbolicAnalysis.jl" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='SymbolicAnalysis.lorentz_transform-Tuple{AbstractMatrix, AbstractVector}' href='#SymbolicAnalysis.lorentz_transform-Tuple{AbstractMatrix, AbstractVector}'><span class="jlbinding">SymbolicAnalysis.lorentz_transform</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
lorentz_transform(O, p) -> AbstractVector
```


Apply the orthochronous Lorentz transformation `O` to the point `p`.

**Arguments**
- `O::AbstractMatrix`: element of the orthochronous Lorentz group.
  
- `p::AbstractVector`: point on the Lorentz manifold.
  

**Returns**
- The transformed point `O * p`.
  

**Throws**
- `ArgumentError`: if `O` does not preserve the Lorentz metric or the positive time direction.
  

**Examples**

```julia
julia> using SymbolicAnalysis, LinearAlgebra

julia> O = Matrix{Float64}(I, 3, 3);

julia> lorentz_transform(O, [0.0, 0.0, 1.0])
3-element Vector{Float64}:
 0.0
 0.0
 1.0
```



<Badge type="info" class="source-link" text="source"><a href="https://github.com/SciML/SymbolicAnalysis.jl" target="_blank" rel="noreferrer">source</a></Badge>

</details>

