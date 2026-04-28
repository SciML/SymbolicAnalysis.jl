# Conic Form Generation and MOI Bridge Tutorial

This tutorial covers how to transform DCP-verified symbolic expressions into
standard conic form and solve them using MathOptInterface (MOI) or JuMP solvers.

## Overview

SymbolicAnalysis.jl can convert any DCP-compliant expression into a conic
formulation via **epigraph reformulation**. The key idea: every DCP atom
(exp, log, norm, etc.) has a corresponding MOI cone. When we walk the
expression tree bottom-up, each atom is replaced by an epigraph variable `t`
plus a cone constraint linking `t` to the atom's arguments. The result is a
linear objective over epigraph variables subject to cone constraints.

The pipeline is:

1. `to_conic_form(expr)` -- convert a symbolic expression to a `ConicFormulation`
2. `to_jump_model(cf)` -- convert to a JuMP model (for high-level modeling)
3. `to_moi_model(cf)` -- convert to a raw MOI model (for direct solver access)
4. Solve with any MOI-compatible solver (SCS, Mosek, ECOS, etc.)
5. `extract_solution(cf, model, var_map)` -- map solution back to original variables

## Basic Usage

### Setup

```julia
using SymbolicAnalysis
using Symbolics
using MathOptInterface
const MOI = MathOptInterface
```

### Converting a Simple Expression

```julia
@variables x

# exp(x) is convex, so this creates a minimization problem
cf = to_conic_form(exp(x) |> unwrap)
```

The returned `ConicFormulation` contains:

- `cf.objective_var` -- the top-level epigraph variable to optimize
- `cf.objective_sense` -- `:minimize` for convex, `:maximize` for concave
- `cf.constraints` -- vector of `ConeConstraint` objects
- `cf.variables` -- all variables (original + epigraph)
- `cf.original_variables` -- only the user's variables
- `cf.epigraph_map` -- maps epigraph variables to their source expressions

```julia
println(cf.objective_sense)       # :minimize
println(cf.original_variables)    # Set([:x])
println(length(cf.constraints))   # number of cone constraints
```

### Concave Expressions

For concave expressions, the system automatically sets the objective sense to
`:maximize`:

```julia
@variables x

cf = to_conic_form(log(x) |> unwrap)
println(cf.objective_sense)  # :maximize
```

## Inspecting Results with `print_conic_form`

The `print_conic_form` function provides a human-readable view of the
formulation:

```julia
@variables x

cf = to_conic_form(exp(x) |> unwrap)
print_conic_form(cf)
```

Output shows the objective, variables, and each constraint with its cone type
and affine expressions:

```
Conic Formulation:
  Objective: minimize _t1
  Original variables: x
  Epigraph variables: _t1
  Constraints (1):
    [1] exp: (x, 1, _t1) in ExponentialCone
        row 1: x
        row 2: 1.0
        row 3: _t1
```

You can also write to a file or buffer:

```julia
io = IOBuffer()
print_conic_form(cf; io=io)
output = String(take!(io))
```

## Examples by Atom

### Exponential Cone Atoms

**exp(x)**: `exp(x) <= t` is encoded as `(x, 1, t) in ExponentialCone`.

```julia
cf = to_conic_form(exp(x) |> unwrap)
exp_constraints = filter(c -> c.cone isa MOI.ExponentialCone, cf.constraints)
# exp_constraints[1] has 3 terms: (x, 1, t)
```

**log(x)**: `log(x) >= t` is encoded as `(t, 1, x) in ExponentialCone`.

```julia
cf = to_conic_form(log(x) |> unwrap)
# Objective sense is :maximize since log is concave
```

### Norm and Absolute Value

**abs(x)**: `|x| <= t` is encoded as `(t, x) in NormOneCone(2)`.

```julia
cf = to_conic_form(abs(x) |> unwrap)
norm_constraints = filter(c -> c.cone isa MOI.NormOneCone, cf.constraints)
```

**norm(x)**: `||x||_2 <= t` is encoded as `(t, x...) in SecondOrderCone`.

### RSOC Atoms

**sqrt(x)**: `sqrt(x) >= t` is encoded as `(x, 0.5, t) in RSOC(3)`, which
gives `2*x*0.5 >= t^2`, i.e., `t <= sqrt(x)`.

```julia
cf = to_conic_form(sqrt(x) |> unwrap)
rsoc = filter(c -> c.cone isa MOI.RotatedSecondOrderCone, cf.constraints)
# rsoc[1].terms[2].constant == 0.5  (the constant row)
```

**inv(x)**: `1/x <= t` is encoded as `(t, x, sqrt(2)) in RSOC(3)`, which
gives `2*t*x >= 2`, i.e., `t*x >= 1`.

### LP Atoms (max, min)

**max(x, y)**: Reformulated as LP constraints `t - x >= 0` and `t - y >= 0`.

```julia
@variables x y
cf = to_conic_form(max(x, y) |> unwrap)
nn = filter(c -> c.cone isa MOI.Nonnegatives, cf.constraints)
# Two Nonnegatives constraints
```

**min(x, y)**: Reformulated as `x - t >= 0` and `y - t >= 0`.

## Composite Expressions

The system handles composite DCP expressions by introducing epigraph variables
at each level:

```julia
@variables x

cf = to_conic_form((exp(x) + abs(x)) |> unwrap)
print_conic_form(cf)
```

This produces:
- An `ExponentialCone` constraint for `exp(x)`
- A `NormOneCone` constraint for `abs(x)`
- An equality constraint linking the sum to the objective variable

### Affine Flattening

Pure affine subexpressions are detected and flattened into a single equality
constraint, avoiding unnecessary epigraph variables:

```julia
@variables x y

cf = to_conic_form((2x + 3y + 5) |> unwrap)
# Only 1 epigraph variable and 1 equality constraint
println(length(setdiff(cf.variables, cf.original_variables)))  # 1
```

### Scaling and Constants

Multiplication by constants and addition of constants are handled directly:

```julia
cf = to_conic_form((2 * abs(x) - 1) |> unwrap)
# NormOneCone for abs(x), plus affine constraints for scaling
```

## Converting to a JuMP Model

`to_jump_model` converts a `ConicFormulation` into a JuMP `Model` that can be
solved with any compatible solver:

```julia
import JuMP

cf = to_conic_form(exp(x) |> unwrap)

# Without solver (for inspection)
model = to_jump_model(cf)
println(JuMP.num_variables(model))
println(JuMP.objective_sense(model))  # MIN_SENSE

# With solver
# using SCS
# model = to_jump_model(cf; solver=SCS.Optimizer)
# JuMP.optimize!(model)
```

The JuMP model contains:
- A `VariableRef` for each original and epigraph variable
- The objective set to Min or Max of the objective variable
- All cone constraints translated to JuMP constraint syntax

## Converting to an MOI Model

`to_moi_model` creates a raw `MOI.Utilities.Model{Float64}` for direct
solver access:

```julia
cf = to_conic_form(exp(x) |> unwrap)
moi_model, var_map = to_moi_model(cf)
```

The returned `var_map` is a `Dict{Symbol, MOI.VariableIndex}` mapping variable
names to their MOI indices. You can inspect the model:

```julia
# Check for ExponentialCone constraints
exp_ci = MOI.get(moi_model,
    MOI.ListOfConstraintIndices{
        MOI.VectorAffineFunction{Float64},
        MOI.ExponentialCone
    }())
println(length(exp_ci))  # >= 1
```

## Extracting Solutions

After solving an MOI model, use `extract_solution` to map results back to
original variable names:

```julia
# After solving:
# solution = extract_solution(cf, solved_model, var_map)
# solution[:x]  # optimal value of x
```

`extract_solution` returns a `Dict{Symbol, Float64}` containing only the
original (user) variables, not the epigraph auxiliaries.

## Supported Cone Types

The following table lists all atoms with their MOI cone mappings:

### Exponential Cone (`MOI.ExponentialCone`)

| Atom | Curvature | Reformulation |
|------|-----------|---------------|
| `exp(x)` | Convex | `(x, 1, t) in ExponentialCone` |
| `log(x)` | Concave | `(t, 1, x) in ExponentialCone` |
| `log1p(x)` | Concave | `(t, 1, 1+x) in ExponentialCone` |
| `logistic(x)` | Convex | Two `ExponentialCone` + one `Nonnegatives` |
| `xlogx(x)` | Convex | `(t, x, 1) in RelativeEntropyCone(3)` |
| `logsumexp(x)` | Convex | `ExponentialCone` (via decomposition) |
| `xexpx(x)` | Convex | `ExponentialCone` |

### Relative Entropy Cone (`MOI.RelativeEntropyCone`)

| Atom | Curvature | Reformulation |
|------|-----------|---------------|
| `xlogx(x)` | Convex | `(t, x, 1) in RelativeEntropyCone(3)` |
| `rel_entr(x, y)` | Convex | `(t, x, y) in RelativeEntropyCone(3)` |
| `kldivergence(p, q)` | Convex | `(t, p, q) in RelativeEntropyCone(3)` |

### Second Order Cone (`MOI.SecondOrderCone`)

| Atom | Curvature | Reformulation |
|------|-----------|---------------|
| `norm(x, 2)` | Convex | `(t, x...) in SecondOrderCone(n+1)` |
| `huber(x, M)` | Convex | Generic SOC (see Limitations) |

### Norm One Cone (`MOI.NormOneCone`)

| Atom | Curvature | Reformulation |
|------|-----------|---------------|
| `abs(x)` | Convex | `(t, x) in NormOneCone(2)` |
| `tv(x)` | Convex | `NormOneCone` |

### Rotated Second Order Cone (`MOI.RotatedSecondOrderCone`)

| Atom | Curvature | Reformulation |
|------|-----------|---------------|
| `sqrt(x)` | Concave | `(x, 0.5, t) in RSOC(3)` |
| `inv(x)` | Convex | `(t, x, sqrt(2)) in RSOC(3)` |
| `quad_over_lin(x, y)` | Convex | `(y/2, t, x) in RSOC(3)` |
| `x^2` | Convex | `(t, 0.5, x) in RSOC(3)` |
| `harmmean(x)` | Concave | `RSOC` |
| `invprod(x)` | Convex | `RSOC` |

### Power Cone (`MOI.PowerCone`)

| Atom | Curvature | Reformulation |
|------|-----------|---------------|
| `x^p` (p > 1) | Convex | `(t, 1, x) in PowerCone(1/p)` |
| `x^p` (0 < p < 1) | Concave | `(x, 1, t) in PowerCone(p)` |
| `x^p` (p < 0) | Convex | `(t, x, 1) in PowerCone(1/(1-p))` |

### Geometric Mean Cone (`MOI.GeometricMeanCone`)

| Atom | Curvature | Reformulation |
|------|-----------|---------------|
| `geomean(x)` | Concave | `(t, x...) in GeometricMeanCone(n+1)` |

### LP Reformulations (No Cone -- Linear Constraints)

| Atom | Curvature | Reformulation |
|------|-----------|---------------|
| `max(a, b)` | Convex | `t - a >= 0`, `t - b >= 0` |
| `min(a, b)` | Concave | `a - t >= 0`, `b - t >= 0` |
| `maximum(x)` | Convex | `t - x_i >= 0` for all i |
| `minimum(x)` | Concave | `x_i - t >= 0` for all i |

### PSD Cone (Registered but handled via generic fallback)

| Atom | Curvature | Cone Annotation |
|------|-----------|-----------------|
| `eigmax(X)` | Convex | `PositiveSemidefiniteConeTriangle` |
| `eigmin(X)` | Concave | `PositiveSemidefiniteConeTriangle` |
| `logdet(X)` | Concave | `LogDetConeTriangle` |
| `quad_form(x, P)` | Convex | `PositiveSemidefiniteConeTriangle` |
| `matrix_frac(x, P)` | Convex | `PositiveSemidefiniteConeTriangle` |
| `trinv(X)` | Convex | `PositiveSemidefiniteConeTriangle` |

## Introspection with `list_cone_annotations`

To see all registered atoms and their cone annotations:

```julia
annotations = list_cone_annotations()
for a in annotations
    println("$(a.atom): type=$(a.type), cone=$(a.cone)")
end
```

This returns a vector of named tuples with fields `atom`, `type` (`:DCP` or
`:GDCP`), `cone`, and either `curvature` or `gcurvature`.

## Thread Safety

`to_conic_form` is thread-safe. Each call creates its own local `ConicContext`
with no global mutable state. You can safely call it from multiple threads:

```julia
results = Vector{ConicFormulation}(undef, 4)
Threads.@threads for i in 1:4
    results[i] = to_conic_form(exp(x) |> unwrap)
end
```

## Known Limitations

1. **DCP compliance not checked.** `to_conic_form` does not verify that the
   input expression is DCP-compliant. If given a non-DCP expression, it may
   produce an incorrect formulation or error. Run `analyze(expr)` first to
   confirm DCP compliance.

2. **Huber loss.** The `huber(x, M)` atom falls through to a generic SOC
   constraint which is not a mathematically correct conic reformulation of the
   Huber loss. A proper decomposition into RSOC + LP constraints is not yet
   implemented.

3. **General division.** Division of two nonlinear expressions (`a/b` where
   both `a` and `b` are non-affine) cannot be correctly represented as a
   linear equality. The `constant / expr` case works correctly via RSOC.

4. **Vector KL divergence.** The `kldivergence` reformulation currently handles
   the scalar case. For element-wise vector KL divergence, the reformulation
   would need expansion to `RelativeEntropyCone(2n+1)`.

5. **Matrix-valued atoms.** Atoms like `eigmax`, `logdet`, `quad_form` have
   cone annotations registered (`PositiveSemidefiniteConeTriangle`,
   `LogDetConeTriangle`) but are handled through a generic fallback rather
   than specialized reformulations.

6. **Power edge cases.** The power atom `x^p` does not explicitly handle
   `p == 0` (constant) or `p == 1` (identity). These cases fall through to a
   generic handler, though they may be caught by affine detection earlier in
   the pipeline.
