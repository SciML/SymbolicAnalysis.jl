# Atoms for DCP and DGCP

This page is intended to be a reference for the atoms that are currently implemented in this package with their respective properties. As much as possible atoms are created with functions from base, standard libraries and popular packages, but we also inherit a few functions from the CVX family of packages such as `quad_form`, `quad_over_lin` etc. and also introduce some new functions in this package. Description of all such special functions implemented in this package is available in the [Special functions](@ref) section of the documentation.

## DCP Atoms

| Atom | Domain | Sign | Curvature | Monotonicity |
|:-----|:-------|:-----|:----------|:-------------|
| dot | (array_domain(ℝ), array_domain(ℝ)) | AnySign | Affine | Increasing |
| dotsort | (array_domain(ℝ, 1), array_domain(ℝ, 1)) | AnySign | Convex | (AnyMono, increasing_if_positive ∘ minimum) |
| StatsBase.geomean | array_domain(HalfLine{Real,:open}(), 1) | Positive | Concave | Increasing |
| StatsBase.harmmean | array_domain(HalfLine{Real,:open}(), 1) | Positive | Concave | Increasing |
| invprod | array_domain(HalfLine{Real,:open}()) | Positive | Convex | Decreasing |
| eigmax | symmetric_domain() | AnySign | Convex | AnyMono |
| eigmin | symmetric_domain() | AnySign | Concave | AnyMono |
| eigsummax | (array_domain(ℝ, 2), ℝ) | AnySign | Convex | AnyMono |
| eigsummin | (array_domain(ℝ, 2), ℝ) | AnySign | Concave | AnyMono |
| logdet | semidefinite_domain() | AnySign | Concave | AnyMono |
| LogExpFunctions.logsumexp | array_domain(ℝ, 2) | AnySign | Convex | Increasing |
| matrix_frac | (array_domain(ℝ, 1), definite_domain()) | AnySign | Convex | AnyMono |
| maximum | array_domain(ℝ) | AnySign | Convex | Increasing |
| minimum | array_domain(ℝ) | AnySign | Concave | Increasing |
| norm | (array_domain(ℝ), Interval{:closed, :open}(1, Inf)) | Positive | Convex | increasing_if_positive |
| norm | (array_domain(ℝ), Interval{:closed, :open}(0, 1)) | Positive | Convex | increasing_if_positive |
| perspective(f, x, s) | (function_domain(), ℝ, Positive) | Same as f | Same as f | AnyMono |
| quad_form | (array_domain(ℝ, 1), semidefinite_domain()) | Positive | Convex | (increasing_if_positive, Increasing) |
| quad_over_lin | (array_domain(ℝ), HalfLine{Real,:open}()) | Positive | Convex | (increasing_if_positive, Decreasing) |
| quad_over_lin | (ℝ, HalfLine{Real,:open}()) | Positive | Convex | (increasing_if_positive, Decreasing) |
| sum | array_domain(ℝ, 2) | AnySign | Affine | Increasing |
| sum_largest | (array_domain(ℝ, 2), ℤ) | AnySign | Convex | Increasing |
| sum_smallest | (array_domain(ℝ, 2), ℤ) | AnySign | Concave | Increasing |
| tr | array_domain(ℝ, 2) | AnySign | Affine | Increasing |
| trinv | definite_domain() | Positive | Convex | AnyMono |
| tv | array_domain(ℝ, 1) | Positive | Convex | AnyMono |
| tv | array_domain(array_domain(ℝ, 2), 1) | Positive | Convex | AnyMono |
| abs | ℂ | Positive | Convex | increasing_if_positive |
| conj | ℂ | AnySign | Affine | AnyMono |
| exp | ℝ | Positive | Convex | Increasing |
| xlogx | ℝ | AnySign | Convex | AnyMono |
| huber | (ℝ, HalfLine()) | Positive | Convex | increasing_if_positive |
| imag | ℂ | AnySign | Affine | AnyMono |
| inv | HalfLine{Real,:open}() | Positive | Convex | Decreasing |
| log | HalfLine{Real,:open}() | AnySign | Concave | Increasing |
| log | array_domain(ℝ, 2) | Positive | Concave | Increasing |
| inv | semidefinite_domain() | AnySign | Convex | Decreasing |
| sqrt | semidefinite_domain() | Positive | Concave | Increasing |
| kldivergence | (array_domain(HalfLine{Real,:open}, 1), array_domain(HalfLine{Real,:open}, 1)) | Positive | Convex | AnyMono |
| lognormcdf | ℝ | Negative | Concave | Increasing |
| log1p | Interval{:open,:open}(-1, Inf) | Negative | Concave | Increasing |
| logistic | ℝ | Positive | Convex | Increasing |
| max | (ℝ, ℝ) | AnySign | Convex | Increasing |
| min | (ℝ, ℝ) | AnySign | Concave | Increasing |
| ^(x, i) | See below | See below | See below | See below |
| real | ℂ | AnySign | Affine | Increasing |
| rel_entr | (HalfLine{Real,:open}(), HalfLine{Real,:open}()) | AnySign | Convex | (AnyMono, Decreasing) |
| sqrt | HalfLine() | Positive | Concave | Increasing |
| xexpx | HalfLine | Positive | Convex | Increasing |
| conv | (array_domain(ℝ, 1), array_domain(ℝ, 1)) | AnySign | Affine | AnyMono |
| cumsum | array_domain(ℝ) | AnySign | Affine | Increasing |
| diagm | array_domain(ℝ, 1) | AnySign | Affine | Increasing |
| diag | array_domain(ℝ, 2) | AnySign | Affine | Increasing |
| diff | array_domain(ℝ) | AnySign | Affine | Increasing |
| kron | (array_domain(ℝ, 2), array_domain(ℝ, 2)) | AnySign | Affine | Increasing |

### Special Cases for ^(x, i)

| Condition on i | Domain | Sign | Curvature | Monotonicity |
|:---------------|:-------|:-----|:----------|:-------------|
| i = 1 | ℝ | AnySign | Affine | Increasing |
| i is even integer | ℝ | Positive | Convex | increasing_if_positive |
| i is odd integer | HalfLine() | Positive | Convex | Increasing |
| i ≥ 1 | HalfLine() | Positive | Convex | Increasing |
| 0 < i < 1 | HalfLine() | Positive | Concave | Increasing |
| i < 0 | HalfLine{Float64,:closed}() | Positive | Convex | Increasing |

## DGCP Atoms (Symmetric Positive Definite)

| Atom | Sign | Geodesic Curvature | Monotonicity |
|:-----|:-----|:-------------------|:-------------|
| LinearAlgebra.logdet | Positive | GLinear | GIncreasing |
| conjugation | Positive | GConvex | GIncreasing |
| LinearAlgebra.tr | Positive | GConvex | GIncreasing |
| sum | Positive | GConvex | GIncreasing |
| adjoint | Positive | GLinear | GIncreasing |
| scalar_mat | Positive | GConvex | GIncreasing |
| LinearAlgebra.diag | Positive | GConvex | GIncreasing |
| sdivergence | Positive | GConvex | GIncreasing |
| Manifolds.distance | Positive | GConvex | GAnyMono |
| SymbolicAnalysis.quad_form | Positive | GConvex | GIncreasing |
| LinearAlgebra.eigmax | Positive | GConvex | GIncreasing |
| log_quad_form | Positive | GConvex | GIncreasing |
| inv | Positive | GConvex | GDecreasing |
| diag | Positive | GConvex | GIncreasing |
| eigsummax | Positive | GConvex | GIncreasing |
| schatten_norm | Positive | GConvex | GIncreasing |
| sum_log_eigmax | Positive | GConvex | GIncreasing |
| affine_map | Positive | GConvex | GIncreasing |
| hadamard_product | Positive | GConvex | GIncreasing |

## DGCP Atoms (Lorentz Model)

| Atom | Sign | Geodesic Curvature | Monotonicity |
|:-----|:-----|:-------------------|:-------------|
| lorentz_distance | Positive | GConvex | GAnyMono |
| lorentz_log_barrier | Positive | GConvex | GIncreasing |
| lorentz_homogeneous_quadratic | Positive | GConvex | GAnyMono |
| lorentz_homogeneous_diagonal | Positive | GConvex | GAnyMono |
| lorentz_least_squares | Positive | GConvex | GAnyMono |
| lorentz_transform | - | - | - |

Note: `lorentz_transform` does not have specific geodesic curvature properties by itself, but it preserves geodesic convexity when applied to geodesically convex functions.