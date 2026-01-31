# DGCP Atoms Reference Table

> **Verification Note**: This document was verified against the source code on 2026-01-30.
> All atoms, curvatures, and monotonicities have been confirmed to match the implementations in:
> - `src/gdcp/spd.jl` (SPD manifold atoms)
> - `src/gdcp/lorentz.jl` (Lorentz manifold atoms)
> - `src/gdcp/gdcp_rules.jl` (GDCP rule infrastructure)
> - `src/rules.jl` (DCP rule infrastructure)

This document provides a comprehensive table of all DGCP (Disciplined Geodesically Convex Programming) atoms supported by SymbolicAnalysis.jl. These atoms form the building blocks for constructing and verifying geodesically convex expressions.

## SPD Manifold Atoms (Symmetric Positive Definite Matrices)

These atoms are defined on the manifold of symmetric positive definite matrices with the affine-invariant Riemannian metric.

### Scalar-Valued Atoms

| Atom | Domain | Sign | G-Curvature | Monotonicity | Source | Reference |
|------|--------|------|-------------|--------------|--------|-----------|
| `logdet(X)` | SPD | AnySign | GLinear | GIncreasing | Literature | Vishnoi (2018); Bacak (2014) |
| `tr(X)` | SPD | Positive | GConvex | GIncreasing | Literature | Vishnoi (2018) |
| `sum(X)` | SPD | Positive | GConvex | GIncreasing | New | - |
| `sdivergence(X, Y)` | SPD | Positive | GConvex | GIncreasing | Literature | Sra (2015) |
| `distance(M, X, Y)` | SPD | Positive | GConvex | GAnyMono | Literature | Bacak (2014); Bhatia (2007) |
| `quad_form(h, X)` | SPD | Positive | GConvex | GIncreasing | Literature | - |
| `eigmax(X)` | SPD | Positive | GConvex | GIncreasing | Literature | - |
| `log_quad_form(y, X)` | SPD | Positive | GConvex | GIncreasing | Literature | Wiesel (2012) |
| `eigsummax(X, k)` | SPD | Positive | GConvex | GIncreasing | Literature | Sra (2015) |
| `schatten_norm(X, p)` | SPD | Positive | GConvex | GIncreasing | Literature | Sra (2015) |
| `sum_log_eigmax(X, k)` | SPD | Positive | GConvex | GIncreasing | Literature | Sra (2015) |
| `sum_log_eigmax(f, X, k)` | SPD | Positive | GConvex | GIncreasing | Literature | Sra (2015) |

### Matrix-Valued Atoms

| Atom | Domain | Sign | G-Curvature | Monotonicity | Source | Reference |
|------|--------|------|-------------|--------------|--------|-----------|
| `conjugation(X, B)` | SPD | Positive | GConvex | GIncreasing | Literature | Vishnoi (2018) |
| `adjoint(X)` | SPD | Positive | GLinear | GIncreasing | New | - |
| `inv(X)` | SPD | Positive | GConvex | GDecreasing | Literature | Bhatia (2007) |
| `diag(X)` | SPD | Positive | GConvex | GIncreasing | Literature | Vishnoi (2018) |
| `scalar_mat(X, k)` | SPD | Positive | GConvex | GIncreasing | New | - |
| `hadamard_product(X, B)` | SPD | Positive | GConvex | GIncreasing | Literature | Vishnoi (2018) |
| `affine_map(f, X, B, Y)` | SPD | Positive | GConvex | GIncreasing | New | Based on Sra (2015) |

## Lorentz Model Atoms (Hyperbolic Space)

These atoms are defined on the Lorentz model of hyperbolic space, a Cartan-Hadamard manifold of constant negative curvature.

| Atom | Domain | Sign | G-Curvature | Monotonicity | Source | Reference |
|------|--------|------|-------------|--------------|--------|-----------|
| `distance(M, p, q)` | Lorentz | Positive | GConvex | GAnyMono | Literature | Bacak (2014) |
| `lorentz_log_barrier(p)` | Lorentz | Positive | GConvex | GIncreasing | Literature | Ferreira et al. (2022) |
| `lorentz_homogeneous_quadratic(A, p)` | Lorentz | Positive | GConvex | GAnyMono | Literature | Ferreira et al. (2022) |
| `lorentz_homogeneous_diagonal(a, p)` | Lorentz | Positive | GConvex | GAnyMono | Literature | Ferreira et al. (2022) |
| `lorentz_nonhomogeneous_quadratic(A, b, c, p)` | Lorentz | AnySign | GConvex | AnyMono | Literature | Ferreira et al. (2023) |
| `lorentz_least_squares(X, y, p)` | Lorentz | Positive | GConvex | AnyMono | Literature | Ferreira et al. (2023) |

## Standard DCP Atoms

These atoms follow standard Disciplined Convex Programming rules and are defined on Euclidean domains. They can be composed with DGCP atoms through scalar composition rules.

### Affine Atoms

| Atom | Domain | Sign | Curvature | Monotonicity | Source | Reference |
|------|--------|------|-----------|--------------|--------|-----------|
| `+` | Real | AnySign | Affine | Increasing | Literature | Grant & Boyd (2006) |
| `-` | Real | AnySign | Affine | Decreasing | Literature | Grant & Boyd (2006) |
| `dot(x, y)` | Real arrays | AnySign | Affine | Increasing | Literature | Grant & Boyd (2006) |
| `sum(x)` | Real arrays | AnySign | Affine | Increasing | Literature | Grant & Boyd (2006) |
| `tr(X)` | Real matrices | AnySign | Affine | Increasing | Literature | Grant & Boyd (2006) |
| `diag(X)` | Real matrices | AnySign | Affine | Increasing | Literature | Grant & Boyd (2006) |
| `diagm(x)` | Real vectors | AnySign | Affine | Increasing | Literature | Grant & Boyd (2006) |
| `vec(X)` | Real matrices | AnySign | Affine | Increasing | Literature | Grant & Boyd (2006) |
| `reshape(X)` | Real matrices | AnySign | Affine | Increasing | Literature | Grant & Boyd (2006) |
| `hcat(...)` | Real vectors | AnySign | Affine | Increasing | Literature | Grant & Boyd (2006) |
| `vcat(...)` | Real vectors | AnySign | Affine | Increasing | Literature | Grant & Boyd (2006) |
| `kron(A, B)` | Real matrices | AnySign | Affine | Increasing | Literature | Grant & Boyd (2006) |
| `triu(X)` | Real matrices | AnySign | Affine | Increasing | Literature | Grant & Boyd (2006) |
| `cumsum(x)` | Real arrays | AnySign | Affine | Increasing | Literature | Grant & Boyd (2006) |
| `diff(x)` | Real arrays | AnySign | Affine | Increasing | Literature | Grant & Boyd (2006) |
| `conv(x, y)` | Real vectors | AnySign | Affine | AnyMono | Literature | Grant & Boyd (2006) |
| `real(z)` | Complex | AnySign | Affine | Increasing | Literature | Grant & Boyd (2006) |
| `imag(z)` | Complex | AnySign | Affine | AnyMono | Literature | Grant & Boyd (2006) |
| `conj(z)` | Complex | AnySign | Affine | AnyMono | Literature | Grant & Boyd (2006) |
| `adjoint(x)` | Real vectors | AnySign | Affine | Increasing | Literature | Grant & Boyd (2006) |

### Convex Atoms

| Atom | Domain | Sign | Curvature | Monotonicity | Source | Reference |
|------|--------|------|-----------|--------------|--------|-----------|
| `abs(x)` | Complex | Positive | Convex | increasing_if_positive | Literature | Grant & Boyd (2006) |
| `exp(x)` | Real | Positive | Convex | Increasing | Literature | Grant & Boyd (2006) |
| `huber(x, M)` | Real | Positive | Convex | increasing_if_positive | Literature | Grant & Boyd (2006) |
| `inv(x)` | Positive Real | Positive | Convex | Decreasing | Literature | Grant & Boyd (2006) |
| `inv(X)` | Semidefinite | AnySign | Convex | Decreasing | Literature | Grant & Boyd (2006) |
| `xlogx(x)` | Real | AnySign | Convex | AnyMono | Literature | Grant & Boyd (2006) |
| `logistic(x)` | Real | Positive | Convex | Increasing | Literature | Grant & Boyd (2006) |
| `max(x, y)` | Real | AnySign | Convex | Increasing | Literature | Grant & Boyd (2006) |
| `maximum(x)` | Real arrays | AnySign | Convex | Increasing | Literature | Grant & Boyd (2006) |
| `norm(x, p)` | Real arrays, p >= 1 | Positive | Convex | increasing_if_positive | Literature | Grant & Boyd (2006) |
| `dotsort(x, y)` | Real vectors | AnySign | Convex | varying | New | - |
| `eigmax(X)` | Symmetric | AnySign | Convex | AnyMono | Literature | Grant & Boyd (2006) |
| `eigsummax(X, k)` | Symmetric | AnySign | Convex | AnyMono | New | - |
| `logsumexp(X)` | Real arrays | AnySign | Convex | Increasing | Literature | Grant & Boyd (2006) |
| `matrix_frac(x, P)` | Real vector, PD | AnySign | Convex | AnyMono | Literature | Grant & Boyd (2006) |
| `quad_form(x, P)` | Real vector, PSD | Positive | Convex | (increasing_if_positive, Increasing) | Literature | Grant & Boyd (2006) |
| `quad_over_lin(x, y)` | Real, Positive | Positive | Convex | (increasing_if_positive, Decreasing) | Literature | Grant & Boyd (2006) |
| `sum_largest(X, k)` | Real matrices | AnySign | Convex | Increasing | Literature | Grant & Boyd (2006) |
| `trinv(X)` | Positive definite | Positive | Convex | AnyMono | Literature | Grant & Boyd (2006) |
| `tv(x)` | Real vectors | Positive | Convex | AnyMono | Literature | Grant & Boyd (2006) |
| `invprod(x)` | Positive Real | Positive | Convex | Decreasing | New | - |
| `rel_entr(x, y)` | Positive Real | AnySign | Convex | (AnyMono, Decreasing) | Literature | Grant & Boyd (2006) |
| `kldivergence(p, q)` | Positive vectors | Positive | Convex | AnyMono | Literature | Grant & Boyd (2006) |
| `xexpx(x)` | Positive | Positive | Convex | Increasing | Literature | Grant & Boyd (2006) |
| `perspective(f, x, s)` | varies | varies | varies | AnyMono | Literature | Grant & Boyd (2006) |

### Concave Atoms

| Atom | Domain | Sign | Curvature | Monotonicity | Source | Reference |
|------|--------|------|-----------|--------------|--------|-----------|
| `log(x)` | Positive Real | AnySign | Concave | Increasing | Literature | Grant & Boyd (2006) |
| `log(X)` | Real matrices | Positive | Concave | Increasing | Literature | Grant & Boyd (2006) |
| `log1p(x)` | x > -1 | Negative | Concave | Increasing | Literature | Grant & Boyd (2006) |
| `sqrt(x)` | Non-negative | Positive | Concave | Increasing | Literature | Grant & Boyd (2006) |
| `sqrt(X)` | Semidefinite | Positive | Concave | Increasing | Literature | Grant & Boyd (2006) |
| `logdet(X)` | Semidefinite | AnySign | Concave | AnyMono | Literature | Grant & Boyd (2006) |
| `lognormcdf(x)` | Real | Negative | Concave | Increasing | New | - |
| `min(x, y)` | Real | AnySign | Concave | Increasing | Literature | Grant & Boyd (2006) |
| `minimum(x)` | Real arrays | AnySign | Concave | Increasing | Literature | Grant & Boyd (2006) |
| `eigmin(X)` | Symmetric | AnySign | Concave | AnyMono | Literature | Grant & Boyd (2006) |
| `eigsummin(X, k)` | Symmetric | AnySign | Concave | AnyMono | New | - |
| `geomean(x)` | Positive vectors | Positive | Concave | Increasing | Literature | Grant & Boyd (2006) |
| `harmmean(x)` | Positive vectors | Positive | Concave | Increasing | Literature | Grant & Boyd (2006) |
| `sum_smallest(X, k)` | Real matrices | AnySign | Concave | Increasing | Literature | Grant & Boyd (2006) |

### Power Atoms

The power function `x^p` has curvature that depends on the exponent:

| Condition | Domain | Sign | Curvature | Monotonicity | Source |
|-----------|--------|------|-----------|--------------|--------|
| `p = 1` | Real | AnySign | Affine | Increasing | Literature |
| `p` even integer | Real | Positive | Convex | increasing_if_positive | Literature |
| `p` odd integer | Non-negative | Positive | Convex | Increasing | Literature |
| `p >= 1` | Non-negative | Positive | Convex | Increasing | Literature |
| `0 < p < 1` | Non-negative | Positive | Concave | Increasing | Literature |
| `p < 0` | Positive | Positive | Convex | Increasing | Literature |

## References

- Bacak, M. (2014). *Convex Analysis and Optimization in Hadamard Spaces*. De Gruyter.
- Bhatia, R. (2007). *Positive Definite Matrices*. Princeton University Press.
- Boyd, S. & Vandenberghe, L. (2004). *Convex Optimization*. Cambridge University Press.
- Ferreira, O.P., Nemeth, S.Z. & Zhu, J. (2022). Convexity of sets and quadratic functions on the hyperbolic space. *Journal of Optimization Theory and Applications*.
- Ferreira, O.P., Nemeth, S.Z. & Zhu, J. (2023). Convexity of non-homogeneous quadratic functions on the hyperbolic space. *Journal of Optimization Theory and Applications*.
- Grant, M. & Boyd, S. (2006). Disciplined Convex Programming. In *Global Optimization: From Theory to Implementation*, Springer.
- Sra, S. (2015). Conic Geometric Optimization on the Manifold of Positive Definite Matrices. *SIAM Journal on Optimization*.
- Vishnoi, N.K. (2018). Geodesic Convex Optimization: Differentiation on Manifolds, Geodesics, and Convexity. *arXiv preprint*.
- Wiesel, A. (2012). Geodesic convexity and covariance estimation. *IEEE Transactions on Signal Processing*.

## Notes

- **Domain abbreviations**: SPD = Symmetric Positive Definite matrices, Lorentz = Lorentz model of hyperbolic space, Real = real numbers, PD = Positive Definite, PSD = Positive Semi-Definite
- **Sign**: Indicates the sign of the function output (Positive, Negative, AnySign)
- **G-Curvature**: GConvex = geodesically convex, GConcave = geodesically concave, GLinear = both g-convex and g-concave
- **Monotonicity**: GIncreasing/GDecreasing = increasing/decreasing with respect to the Lowner order for matrix arguments, GAnyMono = monotonicity unknown or not applicable
- **Source**: "Literature" indicates the atom's g-convexity was established in prior work; "New" indicates atoms introduced or adapted in SymbolicAnalysis.jl

## Usage

To use these atoms in SymbolicAnalysis.jl:

```julia
using SymbolicAnalysis
using Symbolics
using LinearAlgebra
using Manifolds

# Define symbolic matrix
@variables X[1:3, 1:3]

# Create expression using atoms
expr = logdet(X) + tr(X)

# Analyze geodesic convexity
M = SymmetricPositiveDefinite(3)
result = analyze(expr, M)
println(result.gcurvature)  # GConvex
```
