# Examples

Here are some examples demonstrating the use of the `analyze` function from the `SymbolicAnalysis` package.

## Basic Expression Analysis

```@example euclidean1
using SymbolicAnalysis, Symbolics, Symbolics.DomainSets

@variables x
ex1 = exp(x) - log(x)
result = analyze(ex1)
@show result.curvature
```

This example analyzes a simple expression `exp(x) - log(x)`, determining that it's convex and can have any sign.

## Analysis on Manifolds

We can perform DGCP analysis on the Symmetric Positive Definite (SPD) manifold by passing a manifold from [Manifolds.jl](https://juliamanifolds.github.io/Manifolds.jl/stable/) to the `analyze` function. We consider the Karcher mean problem which involves finding the geometric mean of SPD matrices:

```@example manifold1
using SymbolicAnalysis, Symbolics, Manifolds, LinearAlgebra

@variables X[1:5, 1:5]

M = SymmetricPositiveDefinite(5)

As = [rand(5, 5) for i in 1:5]
As = [As[i] * As[i]' for i in 1:5]  # Make them SPD

ex2 = sum(Manifolds.distance(M, As[i], X)^2 for i in 1:5)
result = analyze(ex2, M)
@show result.curvature
@show result.gcurvature
```

This analysis shows that the Karcher mean objective function is geodesically convex on the SPD manifold.

### Domain aware analysis

We can also assert the domain of the variable by assigning `VarDomain` metadata that takes a `Domain` from the [DomainSets.jl](https://juliaapproximation.github.io/DomainSets.jl/dev/) package.

```@example euclidean1
@variables x y

x = setmetadata(
    x,
    SymbolicAnalysis.VarDomain,
    OpenInterval(0,1),
)

y = setmetadata(
    y,
    SymbolicAnalysis.VarDomain,
    OpenInterval(0,1),
)

ex = SymbolicAnalysis.quad_over_lin(x - y, 1 - max(x, y))
result = analyze(ex)
@show result.curvature
```

This example analyzes a quadratic expression over a linear expression, showing that it's convex.

## Analysis on the Lorentz Manifold

We can also perform DGCP analysis on the Lorentz manifold, which is a model of hyperbolic space:

```@example lorentz1
using SymbolicAnalysis, Symbolics, Manifolds, LinearAlgebra

# Create a Lorentz manifold of dimension 2 (3D ambient space)
M = Lorentz(2)

# Define symbolic variables and fixed points
@variables p[1:3]
q = [0.0, 0.0, 1.0]  # A point on the Lorentz model

# Create a composite function from Lorentz atoms
ex = 2.0 * Manifolds.distance(M, q, p) + 
     SymbolicAnalysis.lorentz_log_barrier(p)

# Analyze the expression
result = analyze(ex, M)
@show result.gcurvature
```

This example shows that the sum of the Lorentz distance function and the log-barrier function is geodesically convex on the Lorentz manifold.