# DGCP Analysis Workflow Tutorial

This tutorial covers the Disciplined Geodesically Convex Programming (DGCP) framework provided by SymbolicAnalysis.jl. You will learn how to define symbolic expressions on Riemannian manifolds, verify their geodesic convexity, and interpret the results.

## What is DGCP?

Disciplined Convex Programming (DCP) is a methodology for constructing and verifying convex optimization problems by composing a set of known convex atoms according to specific rules. DCP works in Euclidean space and can certify that an objective function is convex.

DGCP extends this idea to Riemannian manifolds. Many optimization problems arising in machine learning, statistics, and signal processing involve matrix-valued variables constrained to lie on manifolds such as the symmetric positive definite (SPD) matrices or hyperbolic space. These problems are often *geodesically convex* (g-convex) -- meaning they are convex along geodesics of the manifold -- even though they are non-convex in the Euclidean sense.

DGCP provides:

- A library of *geodesically convex atoms* (functions with known g-curvature on specific manifolds).
- *Composition rules* that propagate g-curvature through nested expressions.
- An `analyze` function that takes a symbolic expression and a manifold and returns the geodesic curvature classification.

When DGCP verifies an objective as g-convex, any Riemannian optimization solver is guaranteed to converge to the global optimum.

## Setup

```julia
using SymbolicAnalysis
using Manifolds
using Symbolics
using LinearAlgebra
```

## Defining Symbolic Variables

Use the `@variables` macro from Symbolics.jl to create symbolic matrix or vector variables:

```julia
# A 5x5 symbolic matrix (for SPD manifold problems)
@variables X[1:5, 1:5]

# A 3-element symbolic vector (for Lorentz manifold problems)
@variables p[1:3]
```

## Defining Manifolds

SymbolicAnalysis.jl currently supports two manifolds from Manifolds.jl:

```julia
# Symmetric Positive Definite matrices of size n x n
M_spd = SymmetricPositiveDefinite(5)

# Lorentz model of hyperbolic space (d-dimensional, (d+1)-dimensional ambient)
M_lor = Lorentz(2)  # 2D hyperbolic space in 3D ambient space
```

## Basic Analysis Workflow

The core function is `analyze(expression, manifold)`. It returns an `AnalysisResult` with three fields:

```julia
@variables X[1:5, 1:5]
M = SymmetricPositiveDefinite(5)

# Build a symbolic expression
expr = logdet(X)

# Analyze it on the SPD manifold
result = analyze(expr, M)

# Inspect the result
result.curvature   # Euclidean curvature: Convex, Concave, Affine, or UnknownCurvature
result.sign        # Sign: Positive, Negative, or AnySign
result.gcurvature  # Geodesic curvature: GConvex, GConcave, GLinear, or GUnknownCurvature
```

You can also call `analyze` without a manifold to get only the Euclidean DCP analysis:

```julia
result = analyze(expr)
result.curvature   # Euclidean curvature
result.sign        # Sign
result.gcurvature  # nothing (no manifold provided)
```

Internally, `analyze` performs these steps:
1. **Canonicalize** the expression (`canonize`) to rewrite it into DGCP-friendly forms.
2. **Propagate sign** information through the expression tree.
3. **Propagate Euclidean curvature** (DCP rules).
4. **Propagate geodesic curvature** (DGCP rules, only if a manifold is provided).

## Understanding Results

### Geodesic Curvature (`gcurvature`)

| Value | Meaning |
|---|---|
| `GConvex` | Verified as geodesically convex on the given manifold |
| `GConcave` | Verified as geodesically concave on the given manifold |
| `GLinear` | Verified as geodesically linear (both g-convex and g-concave) |
| `GUnknownCurvature` | Cannot be verified by DGCP composition rules |

### Euclidean Curvature (`curvature`)

| Value | Meaning |
|---|---|
| `Convex` | Verified as Euclidean convex by DCP |
| `Concave` | Verified as Euclidean concave by DCP |
| `Affine` | Verified as affine (both convex and concave) |
| `UnknownCurvature` | Cannot be verified by DCP |

A key insight is that many functions are `GConvex` on SPD but `UnknownCurvature` in the Euclidean sense. This is precisely the class of problems where DGCP adds value over classical DCP.

## SPD Manifold Examples

The SPD manifold has the richest set of DGCP atoms. Here are the main ones:

### logdet -- Geodesically Linear

```julia
@variables X[1:5, 1:5]
M = SymmetricPositiveDefinite(5)

expr = logdet(X)
result = analyze(expr, M)
# result.gcurvature == GLinear
```

`logdet` is the most fundamental atom on SPD. It is g-linear (both g-convex and g-concave), which means it can appear in both minimization and maximization objectives.

### tr(inv(X)) -- Geodesically Convex

```julia
expr = tr(inv(X))
result = analyze(expr, M)
# result.gcurvature == GConvex
```

The trace of the inverse is g-convex on SPD. The `inv` atom is g-convex with decreasing g-monotonicity, and `tr` is g-convex with increasing g-monotonicity, so their composition is g-convex.

### Riemannian Distance Squared

```julia
A = randn(5, 5); A = A * A' + I  # A fixed SPD matrix

expr = Manifolds.distance(M, A, X)^2
result = analyze(expr, M)
# result.gcurvature == GConvex
```

The squared Riemannian distance `d(A, X)^2` is g-convex in X. This is a fundamental result from Hadamard manifold theory. Note that this function is NOT Euclidean convex.

### Karcher (Frechet) Mean

The Karcher mean minimizes the sum of squared distances:

```julia
As = [let B = randn(5, 5); B * B' + I end for _ in 1:5]

expr = sum(Manifolds.distance(M, Ai, X)^2 for Ai in As)
result = analyze(expr, M)
# result.gcurvature == GConvex
```

This works because the sum of g-convex functions is g-convex.

### S-Divergence

```julia
A = randn(5, 5); A = A * A' + I

expr = SymbolicAnalysis.sdivergence(X, A)
result = analyze(expr, M)
# result.gcurvature == GConvex
```

The symmetric Stein divergence `S(X, Y) = logdet((X+Y)/2) - (1/2)*logdet(X*Y)` is g-convex in its first argument. It is used in matrix mean computations and covariance estimation.

### Conjugation

```julia
A = randn(5, 5); A = A * A' + I

expr = SymbolicAnalysis.conjugation(X, A)  # Computes A' * X * A
result = analyze(expr, M)
# result.gcurvature == GConvex
```

### Brascamp-Lieb Bound

```julia
expr = logdet(SymbolicAnalysis.conjugation(X, A)) - logdet(X)
result = analyze(expr, M)
# result.gcurvature == GConvex
```

This expression arises in the computation of Brascamp-Lieb constants. DGCP verifies it as g-convex because `logdet(conjugation(X, A))` is g-convex and `-logdet(X)` is g-convex (negation of a g-linear function).

### Tyler's M-Estimator

```julia
xs = [randn(5) for _ in 1:3]

expr = sum(SymbolicAnalysis.log_quad_form(x, inv(X)) for x in xs) +
       (1/5) * logdet(X)
result = analyze(expr, M)
# result.gcurvature == GConvex
```

Tyler's M-estimator objective is used for robust covariance estimation under heavy-tailed distributions. It is g-convex on SPD but not Euclidean convex.

### Spectral Functions

```julia
# Sum of k largest eigenvalues of log(X)
expr = SymbolicAnalysis.eigsummax(log(X), 2)
result = analyze(expr, M)
# result.gcurvature == GConvex

# Schatten norm of log(X)
expr = SymbolicAnalysis.schatten_norm(log(X), 3)
result = analyze(expr, M)
# result.gcurvature == GConvex
```

The `log` map pulls the SPD matrix back to the tangent space (symmetric matrices), and spectral functions like `eigsummax` and `schatten_norm` are convex on symmetric matrices.

### Additional Atoms

Other g-convex atoms on SPD include:

- `quad_form(x, X)` -- quadratic form `x' * X * x`
- `log_quad_form(x, X)` -- `log(x' * X * x)`
- `eigmax(X)` -- largest eigenvalue
- `scalar_mat(X)` -- `tr(X) * I`
- `diag(X)` -- diagonal extraction
- `hadamard_product(X, B)` -- element-wise product with a fixed PSD matrix B
- `affine_map(f, X, B, Y)` -- affine map `B + f(X, Y)` for positive linear operators
- `sum_log_eigmax(X, k)` -- sum of logs of k largest eigenvalues

## Lorentz Manifold Examples

The Lorentz model represents hyperbolic space. It is a Cartan-Hadamard manifold of constant negative curvature.

### Distance on Lorentz

```julia
M = Lorentz(2)  # 2D hyperbolic space
@variables p[1:3]

q = [0.0, 0.0, 1.0]  # A fixed point on the Lorentz model
expr = Manifolds.distance(M, q, p)
result = analyze(expr, M)
# result.gcurvature == GConvex
```

### Log-Barrier

```julia
expr = SymbolicAnalysis.lorentz_log_barrier(p)
result = analyze(expr, M)
# result.gcurvature == GConvex
```

The log-barrier function `-log(-1 - <a, p>_L)` is g-convex on the Lorentz model.

### Homogeneous Quadratic

```julia
# Matrix A must satisfy geodesic convexity conditions (Theorem 21)
A = [2.0 0.0 0.0; 0.0 2.0 0.0; 0.0 0.0 1.0]
expr = SymbolicAnalysis.lorentz_homogeneous_quadratic(A, p)
result = analyze(expr, M)
# result.gcurvature == GConvex
```

The matrix A must satisfy one of two conditions for geodesic convexity. The function checks these at construction time and throws an `ArgumentError` if they are not met.

### Diagonal Quadratic

```julia
a = [2.0, 2.0, 1.0]  # Must satisfy min(a[1:d]) + a[d+1] >= 0
expr = SymbolicAnalysis.lorentz_homogeneous_diagonal(a, p)
result = analyze(expr, M)
# result.gcurvature == GConvex
```

### Least Squares on Lorentz

```julia
X_data = [1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0]
y_data = [0.0, 0.0, -1.0]
expr = SymbolicAnalysis.lorentz_least_squares(X_data, y_data, p)
result = analyze(expr, M)
# result.gcurvature == GConvex
```

### Composing Lorentz Atoms

G-convex atoms can be combined on the Lorentz manifold:

```julia
q = [0.0, 0.0, 1.0]
expr = 2.0 * Manifolds.distance(M, q, p) + SymbolicAnalysis.lorentz_log_barrier(p)
result = analyze(expr, M)
# result.gcurvature == GConvex
```

This works because the sum of g-convex functions is g-convex, and a positive scalar multiple of a g-convex function is g-convex.

## Composition Rules

DGCP verifies expressions by propagating geodesic curvature through the expression tree. The composition rules mirror classical DCP but operate on g-curvature:

### Addition
- Sum of g-convex functions is g-convex.
- Sum of g-concave functions is g-concave.
- Sum of g-linear functions is g-linear.
- Mixing g-convex and g-concave in a sum produces `GUnknownCurvature`.

### Scalar Multiplication
- Positive constant times g-convex is g-convex.
- Negative constant times g-convex is g-concave (and vice versa).
- DGCP does not support multiplication of two non-constant symbolic expressions.

### Function Composition
For `f(g(x))` where `f` has known curvature and monotonicity:
- If `f` is convex and increasing, and `g` is g-convex, the composition is g-convex.
- If `f` is convex and decreasing, and `g` is g-concave, the composition is g-convex.
- If `f` is concave and increasing, and `g` is g-concave, the composition is g-concave.
- If `f` is concave and decreasing, and `g` is g-convex, the composition is g-concave.

### Inverse Composition
When `inv(X)` appears as an argument to a DGCP atom, the monotonicity is flipped (increasing becomes decreasing and vice versa), reflecting the order-reversing property of matrix inversion on SPD.

### DCP Fallback
If no DGCP-specific rule exists for a function but a DCP rule does, DGCP will use the DCP rule's curvature and monotonicity to propagate geodesic curvature through compositions. This means standard DCP-convex expressions are automatically handled by DGCP -- DGCP is a strict generalization of DCP.

## Canonicalization

Symbolic representation affects verifiability. Two mathematically equivalent expressions may have different DGCP outcomes depending on how they are written. SymbolicAnalysis provides canonicalization passes to rewrite expressions into DGCP-friendly forms.

### canonize(expr)

Applied automatically by `analyze`. Applies safe rewriting rules:

```julia
@variables X[1:5, 1:5]

# Double inverse simplification: inv(inv(X)) -> X
expr = inv(inv(X)) |> Symbolics.unwrap
canonical = SymbolicAnalysis.canonize(expr)
# Result: X

# log(det(X)) -> logdet(X)
# sum(diag(X)) -> tr(X)
# x'*A*x -> quad_form(x, A)
# B'*X*B -> conjugation(X, B)
```

### canonize_extended(expr)

More aggressive rewriting (not applied automatically):

```julia
# logdet(inv(X)) -> -logdet(X)
# log(a * b) -> log(a) + log(b)
expr = SymbolicAnalysis.canonize_extended(expr)
```

### is_canonical(expr)

Check whether an expression is already in canonical form:

```julia
expr1 = logdet(X) |> Symbolics.unwrap
SymbolicAnalysis.is_canonical(expr1)  # true

expr2 = inv(inv(X)) |> Symbolics.unwrap
SymbolicAnalysis.is_canonical(expr2)  # false
```

### equivalent_forms()

Returns documentation of known equivalent forms where one is DGCP-verifiable and the other is not:

```julia
forms = SymbolicAnalysis.equivalent_forms()
for f in forms
    println("Verifiable:     ", f.verifiable)
    println("Not verifiable: ", f.not_verifiable)
    println("Note:           ", f.note)
    println()
end
```

Key examples:

| Verifiable Form | Non-Verifiable Form | Note |
|---|---|---|
| `-logdet(X)` | `logdet(inv(X))` | Equivalent; use `canonize_extended` to transform |
| `2 * logdet(X)` | `logdet(X)^2` | NOT equivalent -- common mistake |
| `tr(inv(X))` | `sum(eigvals(inv(X)))` | Equivalent; use high-level atoms |

## When DGCP Returns GUnknownCurvature

`GUnknownCurvature` means the framework cannot verify the expression using its composition rules. This does NOT mean the function is not g-convex -- it means DGCP cannot prove it. Common causes:

### Product of Two Symbolic Expressions

```julia
@variables X[1:5, 1:5] Y[1:5, 1:5]

expr = sqrt(X * Y) |> Symbolics.unwrap
result = analyze(expr, M)
# result.gcurvature == GUnknownCurvature
```

DGCP does not support multiplication of two non-constant matrix variables.

### Sum of Matrix Variables

```julia
expr = (X + Y) |> Symbolics.unwrap
result = analyze(expr, M)
# result.gcurvature == GUnknownCurvature
```

Addition of two SPD matrix variables is not g-linear on SPD in general.

### Non-DGCP Compositions

```julia
expr = logdet(X)^2 |> Symbolics.unwrap
result = analyze(expr, M)
# result.gcurvature == GUnknownCurvature
```

Squaring a g-linear function does not preserve g-convexity. Note that `2 * logdet(X)` (which IS g-linear) is a different function from `logdet(X)^2`.

### Symbolic Non-Uniqueness

The same mathematical function can sometimes be written in forms that DGCP can or cannot verify. When you get `GUnknownCurvature`, try:

1. Use `canonize_extended(expr)` to apply additional rewriting rules.
2. Rewrite using high-level atoms (e.g., `distance` instead of manual eigenvalue formulas).
3. Consult `equivalent_forms()` for known problematic patterns.
4. Break the expression into simpler sub-expressions and verify them individually.

## DCP Fallback Behavior

When a function has no DGCP-specific rule but does have a standard DCP rule, DGCP uses the DCP classification to propagate geodesic curvature through compositions. This means:

```julia
@variables X[1:5, 1:5]
M = SymmetricPositiveDefinite(5)

# logdet is concave in DCP but g-linear on SPD
# tr(inv(X)) is convex in DCP and g-convex on SPD
# Their sum is not DCP-verifiable (convex + concave), but IS g-convex
expr = tr(inv(X)) + logdet(X) |> Symbolics.unwrap
result = analyze(expr, M)
# result.gcurvature == GConvex
```

This demonstrates that DGCP strictly generalizes DCP: problems that are not DCP-verifiable (because they mix convex and concave terms in Euclidean space) can still be DGCP-verified when all terms are g-convex on the manifold.

## Complete Example: Matrix Square Root via S-Divergence

Putting it all together, here is a complete example that defines a problem, verifies it with DGCP, and solves it with a Riemannian optimizer:

```julia
using SymbolicAnalysis
using Manifolds
using Symbolics
using LinearAlgebra
using Optimization
using OptimizationManopt

# Step 1: Define the problem
M = SymmetricPositiveDefinite(5)
A = randn(5, 5); A = A * A'  # Random SPD matrix

# Step 2: Verify with DGCP
@variables X[1:5, 1:5]
expr = SymbolicAnalysis.sdivergence(X, A) +
       SymbolicAnalysis.sdivergence(X, Matrix{Float64}(I(5)))
result = analyze(expr, M)
@assert result.gcurvature == SymbolicAnalysis.GConvex

# Step 3: Solve with a Riemannian optimizer (guaranteed global optimum)
f(X_val, p=nothing) = SymbolicAnalysis.sdivergence(X_val, A) +
                       SymbolicAnalysis.sdivergence(X_val, Matrix{Float64}(I(5)))

optf = OptimizationFunction(f, Optimization.AutoZygote())
prob = OptimizationProblem(optf, A / 2; manifold=M)
sol = solve(prob, GradientDescentOptimizer(), maxiters=1000)

# The minimizer is the matrix geometric mean sqrt(A)
@assert sqrt(A) ≈ sol.u rtol=1e-3
```

Because DGCP verified the objective as g-convex, we know the Riemannian gradient descent converges to the unique global minimizer.
