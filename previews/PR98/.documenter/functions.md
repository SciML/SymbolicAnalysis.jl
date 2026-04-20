
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
lorentz_homogeneous_diagonal(a::AbstractVector, p::AbstractVector)
```


Computes the homogeneous diagonal quadratic function `∑(a_i * p_i^2)`. For geodesic convexity, min(a_1,...,a_d) + a_{d+1} ≥ 0.

**Arguments**

```julia
- `a::AbstractVector`: A (d+1)-vector where min(a_1,...,a_d) + a_{d+1} ≥ 0.
- `p::AbstractVector`: A point on the Lorentz manifold.
```



<Badge type="info" class="source-link" text="source"><a href="https://github.com/SciML/SymbolicAnalysis.jl" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='SymbolicAnalysis.lorentz_homogeneous_quadratic-Tuple{AbstractMatrix, AbstractVector}' href='#SymbolicAnalysis.lorentz_homogeneous_quadratic-Tuple{AbstractMatrix, AbstractVector}'><span class="jlbinding">SymbolicAnalysis.lorentz_homogeneous_quadratic</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
lorentz_homogeneous_quadratic(A::AbstractMatrix, p::AbstractVector)
```


Computes the homogeneous quadratic function f(p) = p&#39;Ap on the Lorentz model. For geodesic convexity, A must satisfy one of the conditions in Theorem 21.

**Arguments**

```julia
- `A::AbstractMatrix`: A symmetric matrix in R^((d+1)×(d+1)).
- `p::AbstractVector`: A point on the Lorentz manifold.
```



<Badge type="info" class="source-link" text="source"><a href="https://github.com/SciML/SymbolicAnalysis.jl" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='SymbolicAnalysis.lorentz_least_squares-Tuple{AbstractMatrix, AbstractVector, AbstractVector}' href='#SymbolicAnalysis.lorentz_least_squares-Tuple{AbstractMatrix, AbstractVector, AbstractVector}'><span class="jlbinding">SymbolicAnalysis.lorentz_least_squares</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
lorentz_least_squares(X::AbstractMatrix, y::AbstractVector, p::AbstractVector)
```


Computes the least squares function `‖y - Xp‖²_2 = y'y - 2y'Xp + p'X'Xp` for the Lorentz model.

**Arguments**

```julia
- `X::AbstractMatrix`: A matrix in R^(n×(d+1)).
- `y::AbstractVector`: A vector in R^n.
- `p::AbstractVector`: A point on the Lorentz manifold.
```



<Badge type="info" class="source-link" text="source"><a href="https://github.com/SciML/SymbolicAnalysis.jl" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='SymbolicAnalysis.lorentz_log_barrier-Tuple{AbstractVector}' href='#SymbolicAnalysis.lorentz_log_barrier-Tuple{AbstractVector}'><span class="jlbinding">SymbolicAnalysis.lorentz_log_barrier</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
lorentz_log_barrier(a, p)
```


Computes the log-barrier function for the Lorentz model: `-log(-1 - <a, p>_L)`.

**Arguments**

```julia
- `a`: The vector (0, ..., 0, 1) in R^(d+1).
- `p`: A point on the Lorentz manifold.
```



<Badge type="info" class="source-link" text="source"><a href="https://github.com/SciML/SymbolicAnalysis.jl" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='SymbolicAnalysis.lorentz_nonhomogeneous_quadratic-Tuple{AbstractMatrix, AbstractVector, Real, AbstractVector}' href='#SymbolicAnalysis.lorentz_nonhomogeneous_quadratic-Tuple{AbstractMatrix, AbstractVector, Real, AbstractVector}'><span class="jlbinding">SymbolicAnalysis.lorentz_nonhomogeneous_quadratic</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
lorentz_nonhomogeneous_quadratic(A::AbstractMatrix, b::AbstractVector, c::Real, p::AbstractVector)
```


Computes the non-homogeneous quadratic function f(p) = p&#39;Ap + b&#39;p + c on the Lorentz model. For geodesic convexity, p&#39;Ap must be geodesically convex and b must be in the Lorentz cone L.

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
lorentz_transform(O::AbstractMatrix, p::AbstractVector)
```


Applies a Lorentz transform to a point on the Lorentz manifold. The matrix O must be an element of the orthochronous Lorentz group O⁺(1,d).

**Arguments**

```julia
- `O::AbstractMatrix`: An element of the orthochronous Lorentz group.
- `p::AbstractVector`: A point on the Lorentz manifold.
```



<Badge type="info" class="source-link" text="source"><a href="https://github.com/SciML/SymbolicAnalysis.jl" target="_blank" rel="noreferrer">source</a></Badge>

</details>

