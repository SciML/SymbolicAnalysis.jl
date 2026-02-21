# Worked Examples: Optimization on the SPD Manifold

This document presents worked examples of optimization problems on the
symmetric positive definite (SPD) manifold that can be verified using
SymbolicAnalysis.jl and its DGCP (Disciplined Geodesically Convex Programming)
framework.

Each example includes:
- A description of the problem and its applications
- The mathematical formulation
- Julia code for DGCP verification with SymbolicAnalysis.jl
- Interpretation of results

These problems are **geodesically convex** on the SPD manifold but
**not Euclidean convex**, meaning classical DCP tools (such as Convex.jl)
cannot verify them. DGCP extends convexity verification to this setting.

---

## 1. Karcher Mean (Frechet Mean on SPD)

### Problem Description

The Karcher mean (also called the Frechet mean) generalizes the arithmetic
mean to Riemannian manifolds. Given a set of SPD matrices
$A_1, \ldots, A_n \in \mathcal{S}_{++}^d$, the Karcher mean is the
minimizer of the sum of squared Riemannian distances:

$$\min_{X \in \mathcal{S}_{++}^d} \sum_{i=1}^{n} d^2(A_i, X)$$

where $d(A, X) = \|\log(A^{-1/2} X A^{-1/2})\|_F$ is the affine-invariant
Riemannian distance on SPD matrices. This problem arises in diffusion tensor
imaging, radar signal processing, and brain-computer interfaces.

The objective is geodesically convex on the SPD manifold (since the SPD
manifold with the affine-invariant metric is a Hadamard manifold, and
squared distance is g-convex on Hadamard manifolds), but it is not
Euclidean convex.

### Mathematical Formulation

$$f(X) = \sum_{i=1}^{n} \left\|\log\left(A_i^{-1/2} X A_i^{-1/2}\right)\right\|_F^2$$

### Julia Code

```julia
using SymbolicAnalysis
using Manifolds
using Symbolics
using LinearAlgebra
using Random

Random.seed!(42)

# Define symbolic matrix variable and manifold
@variables X[1:5, 1:5]
M = SymmetricPositiveDefinite(5)

# Generate sample SPD matrices
As = [let B = randn(5, 5); B * B' + I end for _ in 1:5]

# Construct the Karcher mean objective
objective = sum(Manifolds.distance(M, Ai, X)^2 for Ai in As)

# Verify with DGCP (manifold-aware analysis)
result = analyze(objective, M)
println("Geodesic curvature: ", result.gcurvature)  # GConvex

# Compare with Euclidean DCP analysis
result_eucl = analyze(objective)
println("Euclidean curvature: ", result_eucl.curvature)  # UnknownCurvature
```

### Interpretation

DGCP verifies the objective as `GConvex`, confirming it is geodesically
convex on SPD. Standard DCP returns `UnknownCurvature` because the
Riemannian distance function is not Euclidean convex. This verification
guarantees that any local minimizer found by a Riemannian optimization
algorithm is the global minimizer.

**Reference:** Karcher, H. (1977). Riemannian center of mass and mollifier
smoothing. *Communications on Pure and Applied Mathematics*.

---

## 2. Tyler's M-Estimator

### Problem Description

Tyler's M-estimator is a robust covariance estimator for heavy-tailed
distributions. Given data vectors $x_1, \ldots, x_k \in \mathbb{R}^d$,
the estimator minimizes:

$$\min_{X \in \mathcal{S}_{++}^d} \sum_{i=1}^{k} \log(x_i^\top X^{-1} x_i) + \frac{1}{d} \log\det(X)$$

This is the negative log-likelihood (up to constants) for a matrix-variate
elliptical distribution. The objective is geodesically convex on SPD but
not Euclidean convex, making it invisible to standard DCP analysis.

### Mathematical Formulation

$$f(X) = \sum_{i=1}^{k} \log(x_i^\top X^{-1} x_i) + \frac{1}{d} \log\det(X)$$

The first term uses the `log_quad_form` atom applied to the inverse, and
the second term uses `logdet`, which is geodesically linear on SPD.

### Julia Code

```julia
using SymbolicAnalysis
using Manifolds
using Symbolics
using LinearAlgebra
using Random

Random.seed!(42)

d = 5
@variables X[1:d, 1:d]
M = SymmetricPositiveDefinite(d)

# Generate random data vectors
xs = [randn(d) for _ in 1:5]

# Construct Tyler's M-estimator objective
objective = sum(SymbolicAnalysis.log_quad_form(xi, inv(X)) for xi in xs) +
            (1/d) * logdet(X)

# Verify with DGCP
result = analyze(objective |> Symbolics.unwrap, M)
println("Geodesic curvature: ", result.gcurvature)  # GConvex

# Euclidean analysis cannot verify this
println("Euclidean curvature: ", result.curvature)   # UnknownCurvature
```

### Interpretation

DGCP decomposes this expression as follows:
1. `log_quad_form(x, Y)` is a registered g-convex atom on SPD.
2. `inv(X)` reverses the monotonicity: since `log_quad_form` is g-increasing
   and `inv` is g-decreasing, their composition is g-convex.
3. `logdet(X)` is g-linear on SPD.
4. The sum of g-convex and g-linear terms is g-convex.

A human expert would need to manually verify each of these composition
steps. DGCP automates this process.

**Reference:** Tyler, D. E. (1987). A distribution-free M-estimator of
multivariate scatter. *Annals of Statistics*.

---

## 3. Brascamp-Lieb Bound

### Problem Description

The Brascamp-Lieb inequality is a fundamental result in analysis that
unifies several classical inequalities (Holder, Young, Loomis-Whitney).
Computing the Brascamp-Lieb constant involves optimizing over SPD matrices.
The dual formulation involves maximizing a geodesically concave function,
which is equivalent to minimizing a g-convex function:

$$\min_{X \in \mathcal{S}_{++}^d} \log\det(A^\top X A) - \log\det(X)$$

where $A$ is a given matrix.

### Mathematical Formulation

$$f(X) = \log\det(A^\top X A) - \log\det(X)$$

The first term is `logdet` composed with the `conjugation` atom
$A^\top X A$, and the second is a g-linear term.

### Julia Code

```julia
using SymbolicAnalysis
using Manifolds
using Symbolics
using LinearAlgebra
using Random

Random.seed!(42)

@variables X[1:5, 1:5]
M = SymmetricPositiveDefinite(5)

# A fixed matrix for the conjugation
A = randn(5, 5); A = A * A' + I

# Construct the Brascamp-Lieb bound objective
objective = logdet(SymbolicAnalysis.conjugation(X, A)) - logdet(X)

# Verify with DGCP
result = analyze(objective |> Symbolics.unwrap, M)
println("Geodesic curvature: ", result.gcurvature)  # GConvex
```

### Interpretation

DGCP verifies this by recognizing:
1. `conjugation(X, A) = A'XA` is a registered g-convex atom on SPD.
2. `logdet` composed with a g-convex function via special-case handling in
   `find_gcurvature`: `logdet(conjugation(...))` is detected as g-convex.
3. `-logdet(X)` is g-linear (negation of g-linear), so it is also g-linear.
4. The sum of g-convex and g-linear terms is g-convex.

**Reference:** Sra, S. and Hosseini, R. (2015). Conic geometric optimization
on the manifold of positive definite matrices. *SIAM Journal on Optimization*.

---

## 4. Maximum Likelihood Estimation on SPD

### Problem Description

Given $n$ observed covariance matrices $S_1, \ldots, S_n$ drawn from a
distribution on the SPD manifold, the maximum likelihood estimate of the
Frechet mean minimizes the sum of squared geodesic distances:

$$\min_{X \in \mathcal{S}_{++}^d} \sum_{i=1}^{n} d^2(X, S_i)$$

This is mathematically equivalent to the Karcher mean problem (Example 1),
but arises in a different context: statistical estimation. The problem
appears in covariance estimation for EEG data, financial time series,
and multivariate process control.

### Mathematical Formulation

$$\hat{\Sigma}_{\text{MLE}} = \arg\min_{X \in \mathcal{S}_{++}^d} \sum_{i=1}^{n} \left\|\log(S_i^{-1/2} X S_i^{-1/2})\right\|_F^2$$

### Julia Code

```julia
using SymbolicAnalysis
using Manifolds
using Symbolics
using LinearAlgebra
using Random

Random.seed!(42)

# Problem dimensions
n = 5          # matrix size
num_samples = 10  # number of observed covariance matrices

@variables X[1:n, 1:n]
M = SymmetricPositiveDefinite(n)

# Generate synthetic sample covariance matrices
function generate_samples(n, num_samples)
    samples = Matrix{Float64}[]
    A = randn(n, n); true_mean = A * A' + I
    for _ in 1:num_samples
        B = randn(n, n)
        push!(samples, B * B' + I)
    end
    return samples, true_mean
end

samples, true_mean = generate_samples(n, num_samples)

# Construct MLE objective
objective = sum(Manifolds.distance(M, S, X)^2 for S in samples)

# DGCP verification
dgcp_result = analyze(objective, M)
println("DGCP (geodesic): ", dgcp_result.gcurvature)  # GConvex

# DCP verification
dcp_result = analyze(objective)
println("DCP (Euclidean): ", dcp_result.curvature)     # UnknownCurvature
```

### Interpretation

The DGCP framework verifies this MLE objective as g-convex regardless of
the number of samples or the matrix dimension. This means:

1. The MLE problem has a **unique global minimizer** on the SPD manifold.
2. Any Riemannian optimization algorithm (gradient descent, conjugate
   gradient, trust regions) is guaranteed to converge to this global
   minimizer.
3. No Euclidean convexity-based tool can provide these guarantees, since
   the objective is non-convex in the Euclidean sense.

The verification scales well: DGCP analyzes the symbolic structure of the
expression tree, so the verification time depends on the number of distinct
terms, not on the numerical matrix size.

---

## 5. Matrix Square Root via S-Divergence

### Problem Description

The S-divergence (also called the symmetric Stein divergence or
Jensen-Bregman LogDet divergence) between two SPD matrices $X$ and $Y$ is:

$$S(X, Y) = \log\det\left(\frac{X + Y}{2}\right) - \frac{1}{2}\log\det(X Y)$$

The matrix geometric mean (or matrix square root) $\sqrt{A}$ can be
characterized as the minimizer of:

$$\min_{X \in \mathcal{S}_{++}^d} S(X, A) + S(X, I)$$

This problem is g-convex since each S-divergence term is g-convex in its
first argument and the sum of g-convex functions is g-convex.

### Mathematical Formulation

$$f(X) = S(X, A) + S(X, I) = \log\det\left(\frac{X + A}{2}\right) - \frac{1}{2}\log\det(XA) + \log\det\left(\frac{X + I}{2}\right) - \frac{1}{2}\log\det(X)$$

### Julia Code

```julia
using SymbolicAnalysis
using Manifolds
using Symbolics
using LinearAlgebra
using Random

Random.seed!(42)

@variables X[1:5, 1:5]
M = SymmetricPositiveDefinite(5)

# A fixed SPD matrix
A = randn(5, 5); A = A * A' + I

# Construct S-divergence objective
objective = SymbolicAnalysis.sdivergence(X, A) +
            SymbolicAnalysis.sdivergence(X, Matrix{Float64}(I(5)))

# Verify with DGCP
result = analyze(objective |> Symbolics.unwrap, M)
println("Geodesic curvature: ", result.gcurvature)  # GConvex
```

### Interpretation

DGCP recognizes `sdivergence` as a registered g-convex atom on the SPD
manifold. The sum of two g-convex terms is g-convex, so the overall
objective is verified automatically. The minimizer of this objective is
the matrix geometric mean $X^* = A^{1/2}$, providing a variational
characterization of the matrix square root.

**Reference:** Sra, S. (2016). Positive definite matrices and the
S-divergence. *Proceedings of the American Mathematical Society*.

---

## 6. Riemannian Distance Minimization with Regularization

### Problem Description

A common pattern in manifold optimization is to minimize a sum of squared
distances with a regularization term. For example, diagonal loading
regularization for robust covariance estimation:

$$\min_{X \in \mathcal{S}_{++}^d} \operatorname{tr}(X^{-1}) + \log\det(X) + \gamma \operatorname{tr}(X)$$

Each term has a known geodesic curvature on SPD:
- $\operatorname{tr}(X^{-1})$ is g-convex (trace of inverse)
- $\log\det(X)$ is g-linear
- $\operatorname{tr}(X)$ is g-convex

### Mathematical Formulation

$$f(X) = \operatorname{tr}(X^{-1}) + \log\det(X) + \gamma \operatorname{tr}(X)$$

### Julia Code

```julia
using SymbolicAnalysis
using Manifolds
using Symbolics
using LinearAlgebra

@variables X[1:5, 1:5]
M = SymmetricPositiveDefinite(5)

gamma = 0.5

# Construct the regularized objective
objective = tr(inv(X)) + logdet(X) + gamma * tr(X)

# Verify with DGCP
result = analyze(objective |> Symbolics.unwrap, M)
println("Geodesic curvature: ", result.gcurvature)  # GConvex
```

### Interpretation

DGCP verifies this by applying the composition rules:

| Term | Atom(s) | G-Curvature |
|---|---|---|
| `tr(inv(X))` | `tr` (g-convex, g-increasing) composed with `inv` (g-convex, g-decreasing) | GConvex |
| `logdet(X)` | `logdet` | GLinear |
| `gamma * tr(X)` | `tr` with positive scalar | GConvex |
| **Sum** | Sum of g-convex and g-linear | **GConvex** |

The scalar multiplication by `gamma > 0` preserves g-convexity. The sum
of g-convex and g-linear functions is g-convex.

**Reference:** Ledoit, O. and Wolf, M. (2004). A well-conditioned estimator
for large-dimensional covariance matrices. *Journal of Multivariate Analysis*.

---

## Summary

| Example | Objective | G-Curvature | Eucl. Curvature | Key Atoms |
|---|---|---|---|---|
| Karcher Mean | Sum of squared distances | GConvex | Unknown | `distance` |
| Tyler's M-Estimator | Log-quad forms + logdet | GConvex | Unknown | `log_quad_form`, `inv`, `logdet` |
| Brascamp-Lieb | logdet(conjugation) - logdet | GConvex | Unknown | `conjugation`, `logdet` |
| MLE on SPD | Sum of squared distances | GConvex | Unknown | `distance` |
| S-Divergence | Sum of S-divergences | GConvex | Unknown | `sdivergence` |
| Regularized Estimation | tr(inv) + logdet + tr | GConvex | Unknown | `tr`, `inv`, `logdet` |

All six problems are geodesically convex on the SPD manifold but cannot
be verified as Euclidean convex by classical DCP. DGCP provides automated
verification in milliseconds, replacing manual mathematical analysis that
can require significant expertise.
