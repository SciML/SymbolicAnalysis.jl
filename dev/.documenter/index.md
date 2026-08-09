---
---

# SymbolicAnalysis.jl {#SymbolicAnalysis.jl}

Symbolics-based function property propagation for optimization

SymbolicAnalysis is a package for implementing the Disciplined Programming approach to optimization. Testing convexity structure in nonlinear programs relies on verifying the convexity of objectives and constraints. [Disciplined Convex Programming (DCP)](https://dcp.stanford.edu/), is a framework for automating this verification task for a wide range of convex functions that can be decomposed into basic convex functions (atoms) using convexity-preserving compositions and transformations (rules).

This package aims to utilize expression graph rewriting and metadata propagation provided by Symbolics.jl, for analysis of relevant properties - limited right now to Euclidean Convexity and Geodesic Convexity on the Symmetric Positive Definite manifold. This package provides an easy to expand implementation of "atoms", that are functions that have known properties. This allows users to add atoms to the library more easily than the previous implementations [CVXPY](https://www.cvxpy.org/index.html) and [Convex.jl](https://github.com/jump-dev/Convex.jl).

## Installation {#Installation}

To install this package, run the following in the Julia REPL:

```julia
using Pkg
Pkg.add("SymbolicAnalysis")
```


## Usage {#Usage}

The main interface to this package is the `analyze` function.
<details class='jldocstring custom-block' open>
<summary><a id='SymbolicAnalysis.analyze' href='#SymbolicAnalysis.analyze'><span class="jlbinding">SymbolicAnalysis.analyze</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
analyze(ex, M = nothing) -> AnalysisResult
```


Analyze the symbolic expression `ex` and return its Euclidean curvature and sign. When `M` is supplied, also determine the geodesic curvature on that manifold.

**Arguments**
- `ex`: Symbolics expression to analyze.
  
- `M::Union{AbstractManifold, Nothing} = nothing`: optional manifold for geodesic curvature analysis. `SymmetricPositiveDefinite` and `Lorentz` manifolds are supported.
  

**Returns**

An `AnalysisResult` with the fields:
- `curvature::SymbolicAnalysis.Curvature`: Euclidean curvature of `ex`.
  
- `sign::SymbolicAnalysis.Sign`: inferred sign of `ex`.
  
- `gcurvature::Union{SymbolicAnalysis.GCurvature, Nothing}`: geodesic curvature when `M` is supplied, or `nothing` otherwise.
  

**Throws**
- `AssertionError`: if `M` is not a supported manifold.
  

**Examples**

```julia
julia> using SymbolicAnalysis, Symbolics

julia> @variables x;

julia> result = analyze(exp(x));

julia> result.curvature == SymbolicAnalysis.Convex
true

julia> result.gcurvature === nothing
true
```



<Badge type="info" class="source-link" text="source"><a href="https://github.com/SciML/SymbolicAnalysis.jl" target="_blank" rel="noreferrer">source</a></Badge>

</details>


## Lorentz model {#Lorentz-model}
<details class='jldocstring custom-block' open>
<summary><a id='SymbolicAnalysis.lorentz_log_barrier' href='#SymbolicAnalysis.lorentz_log_barrier'><span class="jlbinding">SymbolicAnalysis.lorentz_log_barrier</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



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
<summary><a id='SymbolicAnalysis.lorentz_homogeneous_quadratic' href='#SymbolicAnalysis.lorentz_homogeneous_quadratic'><span class="jlbinding">SymbolicAnalysis.lorentz_homogeneous_quadratic</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



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
<summary><a id='SymbolicAnalysis.lorentz_homogeneous_diagonal' href='#SymbolicAnalysis.lorentz_homogeneous_diagonal'><span class="jlbinding">SymbolicAnalysis.lorentz_homogeneous_diagonal</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



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
<summary><a id='SymbolicAnalysis.lorentz_least_squares' href='#SymbolicAnalysis.lorentz_least_squares'><span class="jlbinding">SymbolicAnalysis.lorentz_least_squares</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



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
<summary><a id='SymbolicAnalysis.lorentz_transform' href='#SymbolicAnalysis.lorentz_transform'><span class="jlbinding">SymbolicAnalysis.lorentz_transform</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



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

