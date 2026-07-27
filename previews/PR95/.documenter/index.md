
# SymbolicAnalysis.jl {#SymbolicAnalysis.jl}

Symbolics-based function property propagation for optimization

SymbolicAnalysis is a package for implementing the Disciplined Programming approach to optimization. Testing convexity structure in nonlinear programs relies on verifying the convexity of objectives and constraints. [Disciplined Convex Programming (DCP)](https://dcp.stanford.edu/), is a framework for automating this verification task for a wide range of convex functions that can be decomposed into basic convex functions (atoms) using convexity-preserving compositions and transformations (rules).

This package aims to utilize expression graph rewriting and metadata propagation provided by Symbolics.jl, for analysis of relevant properties - limited right now to Euclidean Convexity and Geodesic Convexity on the Symmetric Positive Definite manifold. This package provides an easy to expand implementation of &quot;atoms&quot;, that are functions that have known properties. This allows users to add atoms to the library more easily than the previous implementations [CVXPY](https://www.cvxpy.org/index.html) and [Convex.jl](https://github.com/jump-dev/Convex.jl).

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
analyze(ex)
analyze(ex, M)
```


Analyze the expression `ex` and return the curvature and sign of the expression. If a manifold `M` from [Manifolds.jl](https://juliamanifolds.github.io/Manifolds.jl/stable/) is provided, also return the geodesic curvature of the expression. Currently supports the `SymmetricPositiveDefinite` and `Lorentz` manifolds.

The returned `AnalysisResult` contains the following fields:
- `curvature::SymbolicAnalysis.Curvature`: The curvature of the expression.
  
- `sign::SymbolicAnalysis.Sign`: The sign of the expression.
  
- `gcurvature::Union{SymbolicAnalysis.GCurvature,Nothing}`: The geodesic curvature of the expression if `M` is provided. Otherwise, `nothing`.
  


<Badge type="info" class="source-link" text="source"><a href="https://github.com/SciML/SymbolicAnalysis.jl" target="_blank" rel="noreferrer">source</a></Badge>

</details>

