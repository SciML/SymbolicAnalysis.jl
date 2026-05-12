# SymbolicAnalysis.jl Documentation

SymbolicAnalysis.jl is a Julia package for automated verification of convexity properties of symbolic mathematical expressions. It implements both classical Disciplined Convex Programming (DCP) and Disciplined Geodesically Convex Programming (DGCP), extending convexity verification to optimization problems on Riemannian manifolds such as symmetric positive definite matrices and hyperbolic space.

## Quick Start

```julia
using SymbolicAnalysis, Manifolds, Symbolics, LinearAlgebra

@variables X[1:5, 1:5]
M = SymmetricPositiveDefinite(5)
result = analyze(tr(inv(X)) + logdet(X), M)
result.gcurvature  # GConvex (geodesically convex, but not Euclidean convex)
```

## Documentation Index

### Tutorials

- **[DGCP Analysis Workflow](tutorials/dgcp_tutorial.md)** -- A step-by-step guide to using the DGCP framework. Covers defining symbolic variables and manifolds, running `analyze`, interpreting results (GConvex, GConcave, GLinear, GUnknownCurvature), composition rules, canonicalization, and troubleshooting. Includes examples for both the SPD and Lorentz manifolds.

- **[Conic Form Generation and MOI Bridge](tutorials/conic_form_tutorial.md)** -- How to transform DCP-verified expressions into standard conic form and solve them with MathOptInterface (MOI) or JuMP solvers. Covers the epigraph reformulation pipeline, supported cone types, and integration with solvers like SCS, COSMO, and Clarabel.

### Reference

- **[Worked Examples](examples.md)** -- Six complete worked examples of optimization problems on the SPD manifold: Karcher mean, Tyler's M-estimator, Brascamp-Lieb bound, maximum likelihood estimation, matrix square root via S-divergence, and regularized distance minimization. Each example includes the mathematical formulation, Julia verification code, and interpretation of results.

- **[Atoms Reference Table](atoms_table.md)** -- Complete reference table of all DCP and DGCP atoms supported by SymbolicAnalysis.jl. Lists every atom with its domain, sign, curvature, monotonicity, cone type, and literature reference. Organized by category: SPD manifold atoms, Lorentz manifold atoms, standard DCP atoms (affine, convex, concave), and power atoms.

### Guides

- **[Porting Guide (Python/Matlab)](porting_guide.md)** -- Practical instructions for reimplementing DGCP in Python (using SymPy) or Matlab (using the Symbolic Math Toolbox). Describes the four-stage analysis pipeline, provides complete code for atom registries, expression tree traversal, composition rule application, and DCP/DGCP curvature propagation in both languages.
