#=
MOI/Conic Form Comparison: Convex.jl vs SymbolicAnalysis.jl

Run with: julia --project=test test/experiments/moi_comparison.jl

This script demonstrates SymbolicAnalysis.jl's conic form generation pipeline
side-by-side with Convex.jl, showing equivalent DCP verification + conic reformulation.
=#

using SymbolicAnalysis
using Symbolics
using MathOptInterface
const MOI = MathOptInterface
import JuMP
using SCS
using LinearAlgebra
using Random

Random.seed!(42)

println("="^70)
println("  MOI/Conic Form Comparison: Convex.jl vs SymbolicAnalysis.jl")
println("="^70)

# ─────────────────────────────────────────────────────────────────────
# Example 1: Simple scalar DCP -- exp(x) + abs(y)
# Convex.jl equivalent:
#   x = Variable(); y = Variable()
#   problem = minimize(exp(x) + abs(y))
# ─────────────────────────────────────────────────────────────────────

println("\n── Example 1: minimize exp(x) + abs(y) ──")

@variables x y

expr1 = exp(x) + abs(y)

# Step 1: DCP verification
result1 = analyze(expr1)
println("DCP curvature: $(result1.curvature)")  # Convex
println("Sign:          $(result1.sign)")

# Step 2: Conic form generation
cf1 = to_conic_form(Symbolics.unwrap(expr1))
println("\nConic form:")
print_conic_form(cf1)

# Step 3: Build JuMP model
model1 = to_jump_model(cf1; solver = SCS.Optimizer)
println("\nJuMP model created:")
println("  Variables:   $(JuMP.num_variables(model1))")
println("  Sense:       $(JuMP.objective_sense(model1))")

# Verify cone types present
moi1, vmap1 = to_moi_model(cf1)
exp_ci = MOI.get(
    moi1,
    MOI.ListOfConstraintIndices{MOI.VectorAffineFunction{Float64},MOI.ExponentialCone}(),
)
norm_ci = MOI.get(
    moi1,
    MOI.ListOfConstraintIndices{MOI.VectorAffineFunction{Float64},MOI.NormOneCone}(),
)
println("  ExpCone constraints:     $(length(exp_ci))")
println("  NormOneCone constraints: $(length(norm_ci))")

# ─────────────────────────────────────────────────────────────────────
# Example 2: quad_over_lin -- mirrors Convex.jl's sumsquares
# Convex.jl equivalent:
#   x = Variable()
#   problem = minimize(quad_over_lin(x, 1))  # = x²
# ─────────────────────────────────────────────────────────────────────

println("\n── Example 2: minimize x^2 (RSOC reformulation) ──")

expr2 = x^2

result2 = analyze(expr2)
println("DCP curvature: $(result2.curvature)")

cf2 = to_conic_form(Symbolics.unwrap(expr2))
println("\nConic form:")
print_conic_form(cf2)

moi2, vmap2 = to_moi_model(cf2)
rsoc_ci = MOI.get(
    moi2,
    MOI.ListOfConstraintIndices{
        MOI.VectorAffineFunction{Float64},
        MOI.RotatedSecondOrderCone,
    }(),
)
println("\n  RSOC constraints: $(length(rsoc_ci))")

# ─────────────────────────────────────────────────────────────────────
# Example 3: Composite -- sqrt(x) + exp(y) + max(x, y)
# Mixed cone types: RSOC + ExponentialCone + Nonnegatives
# ─────────────────────────────────────────────────────────────────────

println("\n── Example 3: minimize sqrt(x) + exp(y) + max(x, y) ──")
println("   (Note: sqrt is concave, but the sum as a whole may be mixed)")

# Use individual atoms to show conic decomposition
expr3_exp = exp(y)
expr3_abs = abs(x)

# Verify individual atoms
println("exp(y) curvature: $(analyze(expr3_exp).curvature)")
println("abs(x) curvature: $(analyze(expr3_abs).curvature)")

# Composite conic form
cf3 = to_conic_form(Symbolics.unwrap(exp(y) + abs(x)))
println("\nComposite conic form (exp(y) + abs(x)):")
print_conic_form(cf3)

model3 = to_jump_model(cf3; solver = SCS.Optimizer)
println("\nJuMP model:")
println("  Variables:   $(JuMP.num_variables(model3))")

# ─────────────────────────────────────────────────────────────────────
# Example 4: log(x) -- concave, maximization sense
# Convex.jl equivalent:
#   x = Variable(Positive())
#   problem = maximize(log(x))
# ─────────────────────────────────────────────────────────────────────

println("\n── Example 4: maximize log(x) (concave → maximization) ──")

expr4 = log(x)
result4 = analyze(expr4)
println("DCP curvature: $(result4.curvature)")

cf4 = to_conic_form(Symbolics.unwrap(expr4))
println("\nConic form:")
print_conic_form(cf4)
println("  Objective sense: $(cf4.objective_sense)")  # maximize

# ─────────────────────────────────────────────────────────────────────
# Example 5: rel_entr(x, y) -- RelativeEntropyCone
# Convex.jl equivalent:
#   x = Variable(Positive()); y = Variable(Positive())
#   problem = minimize(rel_entr(x, y))
# ─────────────────────────────────────────────────────────────────────

println("\n── Example 5: minimize rel_entr(x, y) ──")

expr5 = SymbolicAnalysis.rel_entr(x, y)
result5 = analyze(expr5)
println("DCP curvature: $(result5.curvature)")

cf5 = to_conic_form(Symbolics.unwrap(expr5))
println("\nConic form:")
print_conic_form(cf5)

moi5, _ = to_moi_model(cf5)
re_ci = MOI.get(
    moi5,
    MOI.ListOfConstraintIndices{MOI.VectorAffineFunction{Float64},MOI.RelativeEntropyCone}(),
)
println("  RelativeEntropyCone constraints: $(length(re_ci))")

# ─────────────────────────────────────────────────────────────────────
# Example 6: The DGCP advantage -- what Convex.jl CANNOT do
# ─────────────────────────────────────────────────────────────────────

println("\n" * "="^70)
println("  DGCP: Beyond Convex.jl")
println("="^70)

using Manifolds

@variables X[1:5, 1:5]
M = SymmetricPositiveDefinite(5)

# Generate SPD test matrices
A1 = let A = randn(5, 5)
    A * A' + 5I
end
A2 = let A = randn(5, 5)
    A * A' + 5I
end
A3 = let A = randn(5, 5)
    A * A' + 5I
end

# Karcher mean objective
expr_karcher =
    Manifolds.distance(M, A1, X)^2 +
    Manifolds.distance(M, A2, X)^2 +
    Manifolds.distance(M, A3, X)^2 |> Symbolics.unwrap

result_karcher = analyze(expr_karcher, M)
println("\n── Karcher Mean: sum of squared Riemannian distances ──")
println("  Euclidean curvature: $(result_karcher.curvature)")
println("  Geodesic curvature:  $(result_karcher.gcurvature)")
println("  Convex.jl can verify this: NO")
println("  SymbolicAnalysis.jl:       $(result_karcher.gcurvature) ✓")

# Tyler's M-estimator
xs = [randn(5) for _ = 1:3]
expr_tyler =
    sum(SymbolicAnalysis.log_quad_form(v, inv(X)) for v in xs) +
    (1 / 5) * LinearAlgebra.logdet(X) |> Symbolics.unwrap

result_tyler = analyze(expr_tyler, M)
println("\n── Tyler's M-estimator ──")
println("  Euclidean curvature: $(result_tyler.curvature)")
println("  Geodesic curvature:  $(result_tyler.gcurvature)")
println("  Convex.jl can verify this: NO")
println("  SymbolicAnalysis.jl:       $(result_tyler.gcurvature) ✓")

# S-divergence
expr_sdiv =
    SymbolicAnalysis.sdivergence(X, A1) + SymbolicAnalysis.sdivergence(X, A2) |>
    Symbolics.unwrap
result_sdiv = analyze(expr_sdiv, M)
println("\n── S-divergence (Symmetric Stein) ──")
println("  Euclidean curvature: $(result_sdiv.curvature)")
println("  Geodesic curvature:  $(result_sdiv.gcurvature)")
println("  Convex.jl can verify this: NO")
println("  SymbolicAnalysis.jl:       $(result_sdiv.gcurvature) ✓")

println("\n" * "="^70)
println("  Summary")
println("="^70)
println("""
  DCP (Euclidean) examples:
    exp(x) + abs(y)  → Convex  → ExponentialCone + NormOneCone
    x^2              → Convex  → RotatedSecondOrderCone
    log(x)           → Concave → ExponentialCone (maximize)
    rel_entr(x,y)    → Convex  → RelativeEntropyCone

  DGCP (Riemannian) examples -- Convex.jl returns "not DCP":
    Karcher mean     → GConvex (sum of squared distances)
    Tyler M-est.     → GConvex (log_quad_form + logdet)
    S-divergence     → GConvex (symmetric Stein divergence)

  Pipeline: symbolic expr → analyze() → to_conic_form() → to_jump_model()
""")
