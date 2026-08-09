# SymbolicAnalysis.jl

Symbolics-based function property propagation for optimization

SymbolicAnalysis is a package for implementing the Disciplined Programming approach to optimization. Testing convexity structure in nonlinear programs relies on verifying the convexity of objectives and constraints. [Disciplined Convex Programming (DCP)](https://dcp.stanford.edu/), is a framework for automating this verification task for a wide range of convex functions that can be decomposed into basic convex functions (atoms) using convexity-preserving compositions and transformations (rules).

This package aims to utilize expression graph rewriting and metadata propagation provided by Symbolics.jl, for analysis of relevant properties - limited right now to Euclidean Convexity and Geodesic Convexity on the Symmetric Positive Definite manifold. This package provides an easy to expand implementation of "atoms", that are functions that have known properties. This allows users to add atoms to the library more easily than the previous implementations [CVXPY](https://www.cvxpy.org/index.html) and [Convex.jl](https://github.com/jump-dev/Convex.jl).

## Installation

To install this package, run the following in the Julia REPL:

```julia
using Pkg
Pkg.add("SymbolicAnalysis")
```

## Usage

The main interface to this package is the `analyze` function.

```@autodocs
Modules=[SymbolicAnalysis]
Pages=["SymbolicAnalysis.jl"]
```

## Lorentz model

```@docs
SymbolicAnalysis.lorentz_log_barrier
SymbolicAnalysis.lorentz_homogeneous_quadratic
SymbolicAnalysis.lorentz_homogeneous_diagonal
SymbolicAnalysis.lorentz_least_squares
SymbolicAnalysis.lorentz_transform
```
