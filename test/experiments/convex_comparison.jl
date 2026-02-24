#=
Direct comparison: Convex.jl vs SymbolicAnalysis.jl on the same problem

  minimize ||Ax - b||²  subject to  x >= 0

Run with: julia --project=test test/experiments/convex_comparison.jl
=#

using Random
Random.seed!(42)

m = 4; n = 5
A = randn(m, n); b = randn(m)

println("=" ^ 70)
println("  Problem: minimize ||Ax - b||²  s.t.  x >= 0")
println("  A is $m × $n, b is $m × 1")
println("=" ^ 70)

# ─────────────────────────────────────────────────────────────────────
# Convex.jl
# ─────────────────────────────────────────────────────────────────────

println("\n── Convex.jl ──")

using Convex, SCS

x_cvx = Variable(n)
problem = minimize(sumsquares(A * x_cvx - b), [x_cvx >= 0])
println("  problem is DCP: $(problem.head == :minimize)")
println("  number of variables: $n")
solve!(problem, SCS.Optimizer; silent = true)
println("  status:   $(problem.status)")
println("  optval:   $(problem.optval)")
println("  x*:       $(round.(vec(x_cvx.value), digits=6))")

# ─────────────────────────────────────────────────────────────────────
# SymbolicAnalysis.jl
# ─────────────────────────────────────────────────────────────────────

println("\n── SymbolicAnalysis.jl ──")

using SymbolicAnalysis
using Symbolics
using MathOptInterface
const MOI = MathOptInterface
import JuMP
using LinearAlgebra

# Use individual scalar symbolic variables (the conic form system
# operates on scalar expressions, not indexed arrays)
@variables x1 x2 x3 x4 x5
xvec = [x1, x2, x3, x4, x5]

# Build the same expression: ||Ax - b||²
residual = A * xvec - b
expr = sum(r^2 for r in residual)

# Step 1: DCP verification
result = analyze(expr)
println("  DCP curvature: $(result.curvature)")
println("  Sign:          $(result.sign)")

# Step 2: Conic form
cf = to_conic_form(Symbolics.unwrap(expr))
println("\n  Conic form summary:")
println("    Objective: $(cf.objective_sense) $(cf.objective_var)")
println("    Original variables: $(sort(collect(cf.original_variables)))")
println("    Epigraph variables: $(length(cf.variables) - length(cf.original_variables))")
println("    Constraints: $(length(cf.constraints))")

# Count cone types
cone_counts = Dict{String, Int}()
for c in cf.constraints
    cname = string(typeof(c.cone))
    cone_counts[cname] = get(cone_counts, cname, 0) + 1
end
for (cname, count) in sort(collect(cone_counts))
    println("      $cname: $count")
end

# Step 3: Build JuMP model and add constraint x >= 0
model = to_jump_model(cf; solver = SCS.Optimizer)

# Map original variable names to JuMP variables
all_vars = JuMP.all_variables(model)
jump_orig = Dict{Symbol, JuMP.VariableRef}()
for v in all_vars
    vname = Symbol(JuMP.name(v))
    if vname in cf.original_variables
        jump_orig[vname] = v
    end
end

# Add x >= 0 constraints
for vname in sort(collect(cf.original_variables))
    JuMP.@constraint(model, jump_orig[vname] >= 0)
end

JuMP.set_silent(model)
JuMP.optimize!(model)

println("\n  status:   $(JuMP.termination_status(model))")
println("  optval:   $(JuMP.objective_value(model))")
orig_names_sorted = sort(collect(cf.original_variables))
x_vals = [JuMP.value(jump_orig[vname]) for vname in orig_names_sorted]
println("  x*:       $(round.(x_vals, digits=6))")

# ─────────────────────────────────────────────────────────────────────
# Compare
# ─────────────────────────────────────────────────────────────────────

println("\n── Comparison ──")
cvx_val = problem.optval
sa_val = JuMP.objective_value(model)
println("  Convex.jl optval:            $(round(cvx_val, digits=8))")
println("  SymbolicAnalysis.jl optval:  $(round(sa_val, digits=8))")
println("  Difference:                  $(round(abs(cvx_val - sa_val), digits=10))")
