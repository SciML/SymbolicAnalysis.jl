#=
Generate complexity analysis plots for the MPC paper.
Produces:
  1. scaling_verification.pdf  -- verification time vs AST nodes (log-log) for 3 families
  2. phase_decomposition.pdf   -- stacked bar chart of phase fractions
  3. matrix_independence.pdf   -- verification time vs matrix dimension (flat line)

Run with: julia --project=test test/experiments/generate_complexity_plots.jl
=#

using SymbolicAnalysis
using Symbolics
using SymbolicUtils: iscall, arguments
using Manifolds
using LinearAlgebra
using Random
using Statistics
using Printf
using CairoMakie

Random.seed!(42)

# ============================================================================
# AST utilities
# ============================================================================

function count_ast_nodes(ex)
    ex = Symbolics.unwrap(ex)
    iscall(ex) || return 1
    return 1 + sum(count_ast_nodes(arg) for arg in arguments(ex); init = 0)
end

# ============================================================================
# Expression constructors
# ============================================================================

function make_karcher(m; n = 5)
    @variables X[1:n, 1:n]
    M = SymmetricPositiveDefinite(n)
    As = [
        let B = randn(n, n)
            B * B' + I
        end for _ = 1:m
    ]
    expr = sum(Manifolds.distance(M, Ai, X)^2 for Ai in As) |> Symbolics.unwrap
    return expr, M
end

function make_tyler(m; n = 5)
    @variables X[1:n, 1:n]
    M = SymmetricPositiveDefinite(n)
    xs = [randn(n) for _ = 1:m]
    expr =
        (
            sum(SymbolicAnalysis.log_quad_form(x, inv(X)) for x in xs) +
            (1 / n) * logdet(X)
        ) |> Symbolics.unwrap
    return expr, M
end

function make_scalar_dcp(m)
    @variables x[1:m]
    expr = sum(exp(x[i]) + log(x[i]) for i = 1:m) |> Symbolics.unwrap
    return expr
end

# ============================================================================
# Timing
# ============================================================================

const WARMUP = 5
const ITERS = 20

function time_min(f)
    for _ = 1:WARMUP
        f()
    end
    times = Vector{UInt64}(undef, ITERS)
    for i = 1:ITERS
        GC.gc(false)
        t0 = time_ns()
        f()
        t1 = time_ns()
        times[i] = t1 - t0
    end
    return minimum(times)
end

# ============================================================================
# Power-law fit
# ============================================================================

function fit_power_law(xs, ys)
    lx = log.(Float64.(xs))
    ly = log.(Float64.(ys))
    n = length(lx)
    mx, my = mean(lx), mean(ly)
    Sxx = sum((lx .- mx) .^ 2)
    Sxy = sum((lx .- mx) .* (ly .- my))
    Syy = sum((ly .- my) .^ 2)
    alpha = Sxy / Sxx
    log_c = my - alpha * mx
    SS_res = sum((ly .- (alpha .* lx .+ log_c)) .^ 2)
    R2 = 1.0 - SS_res / Syy
    return alpha, exp(log_c), R2
end

# ============================================================================
# PART 1: Verification time vs AST nodes
# ============================================================================

println("Running Part 1: Scaling verification...")

term_counts = [1, 2, 4, 8, 16, 32]

# Karcher (DGCP)
karcher_nodes = Int[]
karcher_times = Float64[]
for m in term_counts
    expr, M = make_karcher(m)
    nn = count_ast_nodes(expr)
    t_ns = time_min(() -> analyze(expr, M))
    push!(karcher_nodes, nn)
    push!(karcher_times, t_ns / 1e3)  # microseconds
    @printf("  Karcher m=%2d  nodes=%5d  time=%10.1f us\n", m, nn, t_ns / 1e3)
end

# Tyler (DGCP)
tyler_nodes = Int[]
tyler_times = Float64[]
for m in term_counts
    expr, M = make_tyler(m)
    nn = count_ast_nodes(expr)
    t_ns = time_min(() -> analyze(expr, M))
    push!(tyler_nodes, nn)
    push!(tyler_times, t_ns / 1e3)
    @printf("  Tyler  m=%2d  nodes=%5d  time=%10.1f us\n", m, nn, t_ns / 1e3)
end

# Scalar DCP
scalar_nodes = Int[]
scalar_times = Float64[]
for m in term_counts
    expr = make_scalar_dcp(m)
    nn = count_ast_nodes(expr)
    t_ns = time_min(() -> analyze(expr))
    push!(scalar_nodes, nn)
    push!(scalar_times, t_ns / 1e3)
    @printf("  Scalar m=%2d  nodes=%5d  time=%10.1f us\n", m, nn, t_ns / 1e3)
end

# Fit
alpha_k, c_k, R2_k = fit_power_law(karcher_nodes, karcher_times)
alpha_t, c_t, R2_t = fit_power_law(tyler_nodes, tyler_times)
alpha_s, c_s, R2_s = fit_power_law(scalar_nodes, scalar_times)

@printf("\nScaling exponents:\n")
@printf("  Karcher (DGCP): alpha=%.2f, R²=%.4f\n", alpha_k, R2_k)
@printf("  Tyler   (DGCP): alpha=%.2f, R²=%.4f\n", alpha_t, R2_t)
@printf("  Scalar  (DCP):  alpha=%.2f, R²=%.4f\n", alpha_s, R2_s)

# ---- Plot 1: Log-log scaling ----
fig1 = Figure(size = (500, 380), fontsize = 12)
ax1 = Axis(
    fig1[1, 1],
    xlabel = "AST node count (n)",
    ylabel = "Verification time (μs)",
    xscale = log10,
    yscale = log10,
    title = "Verification time vs. expression size",
)

scatter!(
    ax1,
    karcher_nodes,
    karcher_times,
    label = "Karcher mean (DGCP)",
    marker = :circle,
    markersize = 10,
    color = :steelblue,
)
scatter!(
    ax1,
    tyler_nodes,
    tyler_times,
    label = "Tyler M-est. (DGCP)",
    marker = :utriangle,
    markersize = 10,
    color = :firebrick,
)
scatter!(
    ax1,
    scalar_nodes,
    scalar_times,
    label = "Scalar DCP",
    marker = :diamond,
    markersize = 10,
    color = :forestgreen,
)

# Reference line: O(n)
ns_ref = range(
    minimum(vcat(karcher_nodes, tyler_nodes, scalar_nodes)),
    maximum(vcat(karcher_nodes, tyler_nodes, scalar_nodes)),
    length = 100,
)
# Use karcher fit as reference
lines!(
    ax1,
    collect(ns_ref),
    c_k .* collect(ns_ref) .^ alpha_k,
    linestyle = :dash,
    color = :gray60,
    label = @sprintf("O(n^{%.2f}) fit", alpha_k)
)

axislegend(ax1, position = :lt, framevisible = false, labelsize = 10)

save(
    "/Users/vaibhavdixit02/SymbolicAnalysis.jl/_MPC_v2__DGCP/figures/scaling_verification.pdf",
    fig1,
)
save(
    "/Users/vaibhavdixit02/SymbolicAnalysis.jl/_MPC_v2__DGCP/figures/scaling_verification.png",
    fig1,
    px_per_unit = 3,
)
println("\nSaved scaling_verification.pdf")

# ============================================================================
# PART 2: Phase decomposition
# ============================================================================

println("\nRunning Part 2: Phase decomposition...")

phase_term_counts = [2, 4, 8, 16, 32]
phase_data = []

for m in phase_term_counts
    expr, M = make_karcher(m)
    raw = Symbolics.unwrap(expr)
    nn = count_ast_nodes(raw)

    t_canon = time_min(() -> SymbolicAnalysis.canonize(raw))
    ex1 = SymbolicAnalysis.canonize(raw)

    t_sign = time_min(() -> SymbolicAnalysis.propagate_sign(ex1))
    ex2 = SymbolicAnalysis.propagate_sign(ex1)

    t_curv = time_min(() -> SymbolicAnalysis.propagate_curvature(ex2))
    ex3 = SymbolicAnalysis.propagate_curvature(ex2)

    t_gcurv = time_min(() -> SymbolicAnalysis.propagate_gcurvature(ex3, M))

    push!(
        phase_data,
        (
            m = m,
            nodes = nn,
            canon = t_canon / 1e3,
            sign = t_sign / 1e3,
            curv = t_curv / 1e3,
            gcurv = t_gcurv / 1e3,
        ),
    )

    total = (t_canon + t_sign + t_curv + t_gcurv) / 1e3
    @printf(
        "  m=%2d  nodes=%5d  canon=%6.1f  sign=%6.1f  curv=%6.1f  gcurv=%6.1f  total=%7.1f us\n",
        m,
        nn,
        t_canon / 1e3,
        t_sign / 1e3,
        t_curv / 1e3,
        t_gcurv / 1e3,
        total
    )
end

# Report DGCP/DCP ratio at largest
last = phase_data[end]
dcp_total = last.canon + last.sign + last.curv
dgcp_total = dcp_total + last.gcurv
@printf("\nAt m=%d (%d nodes):\n", last.m, last.nodes)
@printf("  DCP  (3 phases): %.1f us\n", dcp_total)
@printf("  DGCP (4 phases): %.1f us\n", dgcp_total)
@printf("  DGCP/DCP ratio:  %.2fx\n", dgcp_total / dcp_total)
@printf("  gcurvature fraction: %.1f%%\n", 100 * last.gcurv / dgcp_total)

# ---- Plot 2: Stacked bar chart ----
fig2 = Figure(size = (500, 380), fontsize = 12)
ax2 = Axis(
    fig2[1, 1],
    xlabel = "Number of composition terms (m)",
    ylabel = "Verification time (μs)",
    title = "Phase decomposition of DGCP verification",
    xticks = (1:length(phase_data), string.([d.m for d in phase_data])),
)

canon_vals = [d.canon for d in phase_data]
sign_vals = [d.sign for d in phase_data]
curv_vals = [d.curv for d in phase_data]
gcurv_vals = [d.gcurv for d in phase_data]

barplot!(
    ax2,
    repeat(1:length(phase_data), 4),
    vcat(canon_vals, sign_vals, curv_vals, gcurv_vals),
    stack = repeat(1:4, inner = length(phase_data)),
    color = repeat(
        [:steelblue, :forestgreen, :goldenrod, :firebrick],
        inner = length(phase_data),
    ),
)

# Manual legend
elem1 = PolyElement(color = :steelblue)
elem2 = PolyElement(color = :forestgreen)
elem3 = PolyElement(color = :goldenrod)
elem4 = PolyElement(color = :firebrick)
Legend(
    fig2[1, 2],
    [elem1, elem2, elem3, elem4],
    ["Canonicalize", "Sign prop.", "Curvature prop.", "G-curvature prop."],
    framevisible = false,
    labelsize = 10,
)

save(
    "/Users/vaibhavdixit02/SymbolicAnalysis.jl/_MPC_v2__DGCP/figures/phase_decomposition.pdf",
    fig2,
)
save(
    "/Users/vaibhavdixit02/SymbolicAnalysis.jl/_MPC_v2__DGCP/figures/phase_decomposition.png",
    fig2,
    px_per_unit = 3,
)
println("Saved phase_decomposition.pdf")

# ============================================================================
# PART 3: Matrix dimension independence
# ============================================================================

println("\nRunning Part 3: Matrix dimension independence...")

m_fixed = 4
dims = [3, 5, 8, 10, 15, 20, 30]

dim_nodes = Int[]
dim_times = Float64[]

for n in dims
    expr, M = make_karcher(m_fixed; n = n)
    nn = count_ast_nodes(expr)
    t_ns = time_min(() -> analyze(expr, M))
    push!(dim_nodes, nn)
    push!(dim_times, t_ns / 1e3)
    @printf("  n=%3d  nodes=%5d  time=%10.1f us\n", n, nn, t_ns / 1e3)
end

@printf(
    "\nNode count range: %d - %d (%.1fx variation)\n",
    minimum(dim_nodes),
    maximum(dim_nodes),
    maximum(dim_nodes) / minimum(dim_nodes)
)
@printf(
    "Time range: %.1f - %.1f us (%.1fx variation)\n",
    minimum(dim_times),
    maximum(dim_times),
    maximum(dim_times) / minimum(dim_times)
)

# ---- Plot 3: Matrix independence ----
fig3 = Figure(size = (500, 380), fontsize = 12)
ax3 = Axis(
    fig3[1, 1],
    xlabel = "Matrix dimension (p)",
    ylabel = "Verification time (μs)",
    title = "Verification time vs. matrix dimension (m = $m_fixed fixed)",
)

scatter!(ax3, dims, dim_times, marker = :circle, markersize = 12, color = :steelblue)
lines!(ax3, dims, dim_times, color = :steelblue, linewidth = 1.5)

# Add horizontal reference line at mean
mean_t = mean(dim_times)
hlines!(ax3, [mean_t], linestyle = :dash, color = :gray60, linewidth = 1)

save(
    "/Users/vaibhavdixit02/SymbolicAnalysis.jl/_MPC_v2__DGCP/figures/matrix_independence.pdf",
    fig3,
)
save(
    "/Users/vaibhavdixit02/SymbolicAnalysis.jl/_MPC_v2__DGCP/figures/matrix_independence.png",
    fig3,
    px_per_unit = 3,
)
println("Saved matrix_independence.pdf")

# ============================================================================
# Print summary for paper
# ============================================================================

println("\n" * "="^70)
println("SUMMARY FOR PAPER")
println("="^70)
println()
@printf("Scaling exponents (time ~ n^α):\n")
@printf("  Karcher mean (DGCP): α = %.2f, R² = %.4f\n", alpha_k, R2_k)
@printf("  Tyler M-est. (DGCP): α = %.2f, R² = %.4f\n", alpha_t, R2_t)
@printf("  Scalar DCP:          α = %.2f, R² = %.4f\n", alpha_s, R2_s)
println()
@printf("DGCP/DCP overhead ratio: %.2fx\n", dgcp_total / dcp_total)
@printf("G-curvature phase fraction: %.1f%%\n", 100 * last.gcurv / dgcp_total)
println()
@printf("Matrix dimension independence:\n")
@printf(
    "  Nodes: %d-%d across p=%d..%d (%.1fx)\n",
    minimum(dim_nodes),
    maximum(dim_nodes),
    minimum(dims),
    maximum(dims),
    maximum(dim_nodes) / minimum(dim_nodes)
)
@printf("  Time variation: %.1fx\n", maximum(dim_times) / minimum(dim_times))
println()
println("Figures saved to _MPC_v2__DGCP/figures/")
println("  scaling_verification.pdf")
println("  phase_decomposition.pdf")
println("  matrix_independence.pdf")
