"""
Empirical Scaling Analysis for SymbolicAnalysis.jl Verification Algorithms

This script provides rigorous empirical evidence for the O(n) time and space
complexity of the DCP/DGCP verification pipeline in SymbolicAnalysis.jl,
where n is the number of AST nodes in the input expression.

Methodology:
- Expressions with controlled AST node counts are constructed by varying the
  number of composition terms (e.g., Karcher mean with m distance terms).
- Matrix size is held constant (it does not affect AST size; matrices are
  numerical constants embedded in the expression tree).
- Each phase (canonize, propagate_sign, propagate_curvature, propagate_gcurvature)
  is timed separately to decompose overhead.
- Timing uses minimum-of-many-trials to remove GC and OS scheduling artifacts.
- Power-law curve fitting (time = c * n^alpha) on log-log data verifies the
  predicted linear scaling exponent alpha ~= 1.0.
- R^2 goodness-of-fit is reported.

Suitable for inclusion in a Mathematical Programming Computation (MPC) paper.
"""

using SymbolicAnalysis
using Symbolics
using SymbolicUtils: iscall, arguments, operation
using Manifolds
using LinearAlgebra
using Random
using Statistics
using Printf
import JuMP  # import to avoid @variables conflict with Symbolics

Random.seed!(42)

# ============================================================================
# Configuration
# ============================================================================

const WARMUP_ITERS = 3
const TIMING_ITERS = 15   # take minimum of this many trials
const MATRIX_DIM = 5      # fixed matrix dimension (does not affect AST size)

# ============================================================================
# AST Node Counting
# ============================================================================

"""
    count_ast_nodes(ex) -> Int

Count the total number of nodes (internal + leaf) in the expression tree.
"""
function count_ast_nodes(ex)
    ex = Symbolics.unwrap(ex)
    if !iscall(ex)
        return 1  # leaf: variable, number, or constant
    end
    return 1 + sum(count_ast_nodes(arg) for arg in arguments(ex); init = 0)
end

"""
    ast_depth(ex) -> Int

Maximum depth of the expression tree.
"""
function ast_depth(ex)
    ex = Symbolics.unwrap(ex)
    if !iscall(ex)
        return 1
    end
    args = arguments(ex)
    isempty(args) && return 1
    return 1 + maximum(ast_depth(arg) for arg in args)
end

# ============================================================================
# Controlled Expression Construction
# ============================================================================

"""
    make_karcher_expr(m; n=MATRIX_DIM) -> (expr, M)

Build a Karcher mean objective: sum_{i=1}^{m} d^2(A_i, X) on SPD(n).
The number of AST nodes scales linearly with m while matrix dimension n
is held constant. Returns the unwrapped expression and the manifold.
"""
function make_karcher_expr(m; n = MATRIX_DIM)
    @variables X[1:n, 1:n]
    M = SymmetricPositiveDefinite(n)
    As = [
        let B = randn(n, n)
                B * B' + I
        end for _ in 1:m
    ]
    expr = sum(Manifolds.distance(M, Ai, X)^2 for Ai in As) |> Symbolics.unwrap
    return expr, M
end

"""
    make_tyler_expr(m; n=MATRIX_DIM) -> (expr, M)

Build a Tyler M-estimator objective with m observation vectors.
"""
function make_tyler_expr(m; n = MATRIX_DIM)
    @variables X[1:n, 1:n]
    M = SymmetricPositiveDefinite(n)
    xs = [randn(n) for _ in 1:m]
    expr =
        (
        sum(SymbolicAnalysis.log_quad_form(x, inv(X)) for x in xs) +
            (1 / n) * logdet(X)
    ) |> Symbolics.unwrap
    return expr, M
end

"""
    make_scalar_dcp_expr(m) -> expr

Build a purely scalar DCP expression: sum of m terms exp(x_i) + log(x_i).
Each term adds a fixed number of AST nodes.
"""
function make_scalar_dcp_expr(m)
    @variables x[1:m]
    # Each term: exp(x_i) + log(x_i) contributes a fixed number of AST nodes
    expr = sum(exp(x[i]) + log(x[i]) for i in 1:m) |> Symbolics.unwrap
    return expr
end

# ============================================================================
# Timing Utilities
# ============================================================================

"""
    time_min(f; warmup=WARMUP_ITERS, iters=TIMING_ITERS) -> (min_ns, all_ns)

Time `f()` by running it `warmup` times (discarded), then `iters` times,
returning the minimum time in nanoseconds and the full vector of timings.
Uses `time_ns()` for sub-microsecond precision.
"""
function time_min(f; warmup = WARMUP_ITERS, iters = TIMING_ITERS)
    # Warmup
    for _ in 1:warmup
        f()
    end
    # Collect timings
    times = Vector{UInt64}(undef, iters)
    for i in 1:iters
        GC.gc(false)  # minor GC to reduce interference
        t0 = time_ns()
        f()
        t1 = time_ns()
        times[i] = t1 - t0
    end
    return minimum(times), times
end

"""
    time_with_alloc(f; warmup=WARMUP_ITERS, iters=TIMING_ITERS) -> (min_ns, min_alloc_bytes)

Time and measure allocations for `f()`.
"""
function time_with_alloc(f; warmup = WARMUP_ITERS, iters = TIMING_ITERS)
    for _ in 1:warmup
        f()
    end
    min_t = typemax(UInt64)
    min_alloc = typemax(Int)
    for _ in 1:iters
        GC.gc(false)
        alloc = @allocated begin
            t0 = time_ns()
            f()
            t1 = time_ns()
        end
        dt = t1 - t0
        if dt < min_t
            min_t = dt
            min_alloc = alloc
        end
    end
    return min_t, min_alloc
end

# ============================================================================
# Log-Log Linear Regression for Power-Law Fitting
# ============================================================================

"""
    fit_power_law(xs, ys) -> (alpha, log_c, R2)

Fit ys = c * xs^alpha by log-log OLS.  Returns the exponent alpha,
log(c), and the coefficient of determination R^2.
"""
function fit_power_law(xs, ys)
    lx = log.(Float64.(xs))
    ly = log.(Float64.(ys))
    n = length(lx)
    mx = sum(lx) / n
    my = sum(ly) / n
    Sxx = sum((lx .- mx) .^ 2)
    Sxy = sum((lx .- mx) .* (ly .- my))
    Syy = sum((ly .- my) .^ 2)
    alpha = Sxy / Sxx
    log_c = my - alpha * mx
    SS_res = sum((ly .- (alpha .* lx .+ log_c)) .^ 2)
    R2 = 1.0 - SS_res / Syy
    return alpha, log_c, R2
end

# ============================================================================
# Part 1: Verify O(n) Scaling of Full Verification
# ============================================================================

function run_part1_scaling()
    println("="^72)
    println("PART 1: Verify O(n) Scaling of Full Verification Pipeline")
    println("="^72)
    println()

    term_counts = [1, 2, 4, 8, 16, 32]

    # ---- Karcher mean (DGCP) ----
    println("1a. Karcher Mean Objective (DGCP), n=$MATRIX_DIM fixed")
    println("-"^60)

    karcher_nodes = Int[]
    karcher_times_ns = UInt64[]

    for m in term_counts
        expr, M = make_karcher_expr(m)
        nn = count_ast_nodes(expr)
        t_ns, _ = time_min(() -> analyze(expr, M))
        push!(karcher_nodes, nn)
        push!(karcher_times_ns, t_ns)
        @printf("  m=%2d  nodes=%5d  time=%10.1f us\n", m, nn, t_ns / 1.0e3)
    end

    alpha_k, _, R2_k = fit_power_law(karcher_nodes, karcher_times_ns)
    @printf("  Fit: time ~ n^%.3f   R^2 = %.4f\n", alpha_k, R2_k)
    println()

    # ---- Tyler M-estimator (DGCP) ----
    println("1b. Tyler M-Estimator Objective (DGCP), n=$MATRIX_DIM fixed")
    println("-"^60)

    tyler_nodes = Int[]
    tyler_times_ns = UInt64[]

    for m in term_counts
        expr, M = make_tyler_expr(m)
        nn = count_ast_nodes(expr)
        t_ns, _ = time_min(() -> analyze(expr, M))
        push!(tyler_nodes, nn)
        push!(tyler_times_ns, t_ns)
        @printf("  m=%2d  nodes=%5d  time=%10.1f us\n", m, nn, t_ns / 1.0e3)
    end

    alpha_t, _, R2_t = fit_power_law(tyler_nodes, tyler_times_ns)
    @printf("  Fit: time ~ n^%.3f   R^2 = %.4f\n", alpha_t, R2_t)
    println()

    # ---- Scalar DCP ----
    println("1c. Scalar DCP (sum of exp + log terms)")
    println("-"^60)

    scalar_nodes = Int[]
    scalar_times_ns = UInt64[]

    for m in term_counts
        expr = make_scalar_dcp_expr(m)
        nn = count_ast_nodes(expr)
        t_ns, _ = time_min(() -> analyze(expr))
        push!(scalar_nodes, nn)
        push!(scalar_times_ns, t_ns)
        @printf("  m=%2d  nodes=%5d  time=%10.1f us\n", m, nn, t_ns / 1.0e3)
    end

    alpha_s, _, R2_s = fit_power_law(scalar_nodes, scalar_times_ns)
    @printf("  Fit: time ~ n^%.3f   R^2 = %.4f\n", alpha_s, R2_s)
    println()

    println("Part 1 Summary:")
    println("-"^60)
    @printf("  Karcher (DGCP):   alpha = %.3f,  R^2 = %.4f\n", alpha_k, R2_k)
    @printf("  Tyler   (DGCP):   alpha = %.3f,  R^2 = %.4f\n", alpha_t, R2_t)
    @printf("  Scalar  (DCP):    alpha = %.3f,  R^2 = %.4f\n", alpha_s, R2_s)
    println("  Prediction: alpha ~= 1.0 (linear in AST node count)")
    println()

    return (
        karcher = (
            nodes = karcher_nodes,
            times_ns = karcher_times_ns,
            alpha = alpha_k,
            R2 = R2_k,
        ),
        tyler = (
            nodes = tyler_nodes,
            times_ns = tyler_times_ns,
            alpha = alpha_t,
            R2 = R2_t,
        ),
        scalar = (
            nodes = scalar_nodes,
            times_ns = scalar_times_ns,
            alpha = alpha_s,
            R2 = R2_s,
        ),
    )
end

# ============================================================================
# Part 2: Phase Decomposition
# ============================================================================

function run_part2_phase_decomposition()
    println("="^72)
    println("PART 2: Phase Decomposition of Verification Pipeline")
    println("="^72)
    println()
    println("Each phase is timed separately. DCP has 3 phases; DGCP adds a 4th.")
    println("The marginal cost of DGCP is one additional propagate_gcurvature pass.")
    println()

    term_counts = [1, 2, 4, 8, 16, 32]

    println("Karcher Mean, n=$MATRIX_DIM fixed")
    println("-"^72)
    @printf(
        "  %-4s  %6s  %10s  %10s  %10s  %10s  %10s\n",
        "m",
        "nodes",
        "canonize",
        "sign",
        "curvature",
        "gcurvature",
        "total"
    )
    @printf(
        "  %-4s  %6s  %10s  %10s  %10s  %10s  %10s\n",
        "",
        "",
        "(us)",
        "(us)",
        "(us)",
        "(us)",
        "(us)"
    )
    println("  " * "-"^68)

    phase_data = []

    for m in term_counts
        expr, M = make_karcher_expr(m)
        raw_expr = Symbolics.unwrap(expr)
        nn = count_ast_nodes(raw_expr)

        # Phase 1: canonize
        t_canon, _ = time_min() do
            SymbolicAnalysis.canonize(raw_expr)
        end
        ex1 = SymbolicAnalysis.canonize(raw_expr)

        # Phase 2: propagate_sign
        t_sign, _ = time_min() do
            SymbolicAnalysis.propagate_sign(ex1)
        end
        ex2 = SymbolicAnalysis.propagate_sign(ex1)

        # Phase 3: propagate_curvature
        t_curv, _ = time_min() do
            SymbolicAnalysis.propagate_curvature(ex2)
        end
        ex3 = SymbolicAnalysis.propagate_curvature(ex2)

        # Phase 4: propagate_gcurvature (DGCP only)
        t_gcurv, _ = time_min() do
            SymbolicAnalysis.propagate_gcurvature(ex3, M)
        end

        total = t_canon + t_sign + t_curv + t_gcurv

        @printf(
            "  m=%2d  %5d  %10.1f  %10.1f  %10.1f  %10.1f  %10.1f\n",
            m,
            nn,
            t_canon / 1.0e3,
            t_sign / 1.0e3,
            t_curv / 1.0e3,
            t_gcurv / 1.0e3,
            total / 1.0e3
        )

        push!(
            phase_data,
            (
                m = m,
                nodes = nn,
                canon_ns = t_canon,
                sign_ns = t_sign,
                curv_ns = t_curv,
                gcurv_ns = t_gcurv,
                total_ns = total,
            ),
        )
    end
    println()

    # Report phase fractions at largest size
    last = phase_data[end]
    dcp_total = last.canon_ns + last.sign_ns + last.curv_ns
    dgcp_total = last.total_ns
    @printf("Phase fractions at m=%d (%d nodes):\n", last.m, last.nodes)
    @printf("  canonize:            %5.1f%%\n", 100 * last.canon_ns / dgcp_total)
    @printf("  propagate_sign:      %5.1f%%\n", 100 * last.sign_ns / dgcp_total)
    @printf("  propagate_curvature: %5.1f%%\n", 100 * last.curv_ns / dgcp_total)
    @printf(
        "  propagate_gcurvature:%5.1f%%  <-- DGCP marginal cost\n",
        100 * last.gcurv_ns / dgcp_total
    )
    println()
    @printf("DCP-only time (3 phases):  %.1f us\n", dcp_total / 1.0e3)
    @printf("DGCP total   (4 phases):   %.1f us\n", dgcp_total / 1.0e3)
    @printf("DGCP / DCP ratio:          %.2fx\n", dgcp_total / dcp_total)
    println()

    # Fit each phase separately to check O(n)
    if length(phase_data) >= 3
        nodes_vec = [d.nodes for d in phase_data]
        println("Per-phase scaling exponents:")
        for (name, getter) in [
                ("canonize", d -> d.canon_ns),
                ("propagate_sign", d -> d.sign_ns),
                ("propagate_curvature", d -> d.curv_ns),
                ("propagate_gcurvature", d -> d.gcurv_ns),
            ]
            times_vec = [getter(d) for d in phase_data]
            if all(t -> t > 0, times_vec)
                alpha, _, R2 = fit_power_law(nodes_vec, times_vec)
                @printf("  %-24s alpha = %.3f,  R^2 = %.4f\n", name, alpha, R2)
            end
        end
    end
    println()

    return phase_data
end

# ============================================================================
# Part 3: Memory Scaling
# ============================================================================

function run_part3_memory()
    println("="^72)
    println("PART 3: Memory (Allocation) Scaling")
    println("="^72)
    println()

    term_counts = [1, 2, 4, 8, 16, 32]

    println("Karcher Mean (DGCP), n=$MATRIX_DIM fixed")
    println("-"^60)
    @printf("  %-4s  %6s  %12s  %12s\n", "m", "nodes", "time (us)", "alloc (KB)")
    println("  " * "-"^40)

    mem_nodes = Int[]
    mem_alloc = Int[]
    mem_time = UInt64[]

    for m in term_counts
        expr, M = make_karcher_expr(m)
        nn = count_ast_nodes(expr)
        t_ns, alloc = time_with_alloc(() -> analyze(expr, M))
        push!(mem_nodes, nn)
        push!(mem_alloc, alloc)
        push!(mem_time, t_ns)
        @printf("  m=%2d  %5d  %10.1f    %10.1f\n", m, nn, t_ns / 1.0e3, alloc / 1024)
    end
    println()

    if length(mem_nodes) >= 3
        alpha_m, _, R2_m = fit_power_law(mem_nodes, mem_alloc)
        @printf("Memory scaling: alloc ~ n^%.3f   R^2 = %.4f\n", alpha_m, R2_m)
        println("Prediction: alpha ~= 1.0 (linear in AST node count)")

        # Also report bytes per node
        bytes_per_node = [a / n for (a, n) in zip(mem_alloc, mem_nodes)]
        @printf(
            "Bytes per AST node: %.0f - %.0f (range)\n",
            minimum(bytes_per_node),
            maximum(bytes_per_node)
        )
    end
    println()

    return (nodes = mem_nodes, alloc_bytes = mem_alloc, times_ns = mem_time)
end

# ============================================================================
# Part 4: Conic Form Generation Scaling
# ============================================================================

function run_part4_conic()
    println("="^72)
    println("PART 4: Conic Form Generation Scaling")
    println("="^72)
    println()

    # Use scalar DCP expressions since to_conic_form operates on scalar DCP atoms
    term_counts = [1, 2, 4, 8, 16, 32]

    println("Scalar DCP expressions (sum of exp + log terms)")
    println("-"^72)
    @printf(
        "  %-4s  %6s  %12s  %10s  %12s\n",
        "m",
        "nodes",
        "conic (us)",
        "epi_vars",
        "constraints"
    )
    println("  " * "-"^56)

    conic_nodes = Int[]
    conic_times_ns = UInt64[]
    conic_epi = Int[]
    conic_cons = Int[]

    for m in term_counts
        expr = make_scalar_dcp_expr(m)
        nn = count_ast_nodes(expr)

        t_ns, _ = time_min() do
            SymbolicAnalysis.to_conic_form(expr)
        end

        cf = SymbolicAnalysis.to_conic_form(expr)
        n_epi = length(cf.variables) - length(cf.original_variables)
        n_con = length(cf.constraints)

        push!(conic_nodes, nn)
        push!(conic_times_ns, t_ns)
        push!(conic_epi, n_epi)
        push!(conic_cons, n_con)

        @printf("  m=%2d  %5d  %10.1f    %8d  %10d\n", m, nn, t_ns / 1.0e3, n_epi, n_con)
    end
    println()

    if length(conic_nodes) >= 3
        alpha_ct, _, R2_ct = fit_power_law(conic_nodes, conic_times_ns)
        @printf("Conic time scaling:   time ~ n^%.3f   R^2 = %.4f\n", alpha_ct, R2_ct)

        alpha_ce, _, R2_ce = fit_power_law(conic_nodes, conic_epi)
        @printf("Epigraph var scaling: vars ~ n^%.3f   R^2 = %.4f\n", alpha_ce, R2_ce)

        alpha_cc, _, R2_cc = fit_power_law(conic_nodes, conic_cons)
        @printf("Constraint scaling:   cons ~ n^%.3f   R^2 = %.4f\n", alpha_cc, R2_cc)
    end
    println()

    return (
        nodes = conic_nodes,
        times_ns = conic_times_ns,
        epi_vars = conic_epi,
        constraints = conic_cons,
    )
end

# ============================================================================
# Part 5: Comprehensive Data Table for Paper
# ============================================================================

function run_part5_summary_table(part1, part2, part3, part4)
    println("="^72)
    println("PART 5: Summary Data for Paper")
    println("="^72)
    println()

    println("Table 1: Verification Time vs AST Size (Karcher Mean, DGCP)")
    println("-"^60)
    @printf("  %6s  %10s  %10s  %10s\n", "nodes", "time(us)", "alloc(KB)", "us/node")
    println("  " * "-"^44)
    for i in eachindex(part1.karcher.nodes)
        nn = part1.karcher.nodes[i]
        t_us = part1.karcher.times_ns[i] / 1.0e3
        alloc_kb = i <= length(part3.alloc_bytes) ? part3.alloc_bytes[i] / 1024 : NaN
        @printf("  %5d  %10.1f  %10.1f  %10.3f\n", nn, t_us, alloc_kb, t_us / nn)
    end
    println()

    println("Table 2: Phase Decomposition at Largest Problem Size")
    println("-"^60)
    if !isempty(part2)
        last = part2[end]
        total = last.total_ns
        phases = [
            ("canonize", last.canon_ns),
            ("propagate_sign", last.sign_ns),
            ("propagate_curvature", last.curv_ns),
            ("propagate_gcurvature", last.gcurv_ns),
        ]
        @printf("  %-24s  %10s  %8s\n", "Phase", "Time(us)", "Fraction")
        println("  " * "-"^46)
        for (name, t) in phases
            @printf("  %-24s  %10.1f  %7.1f%%\n", name, t / 1.0e3, 100 * t / total)
        end
        dcp_only = last.canon_ns + last.sign_ns + last.curv_ns
        @printf(
            "  %-24s  %10.1f  %7.1f%%\n",
            "DCP total (3 phases)",
            dcp_only / 1.0e3,
            100 * dcp_only / total
        )
        @printf("  %-24s  %10.1f  %7.1f%%\n", "DGCP total (4 phases)", total / 1.0e3, 100.0)
        @printf("  DGCP/DCP ratio: %.2fx\n", total / dcp_only)
    end
    println()

    println("Table 3: Conic Form Generation Scaling")
    println("-"^60)
    @printf("  %6s  %10s  %8s  %11s\n", "nodes", "time(us)", "epi_vars", "constraints")
    println("  " * "-"^42)
    for i in eachindex(part4.nodes)
        @printf(
            "  %5d  %10.1f  %8d  %11d\n",
            part4.nodes[i],
            part4.times_ns[i] / 1.0e3,
            part4.epi_vars[i],
            part4.constraints[i]
        )
    end
    println()

    # Overall scaling exponents summary
    println("Table 4: Fitted Scaling Exponents (time ~ n^alpha)")
    println("-"^60)
    @printf("  %-30s  %8s  %8s\n", "Experiment", "alpha", "R^2")
    println("  " * "-"^50)
    @printf(
        "  %-30s  %8.3f  %8.4f\n",
        "Karcher (DGCP, full)",
        part1.karcher.alpha,
        part1.karcher.R2
    )
    @printf(
        "  %-30s  %8.3f  %8.4f\n",
        "Tyler (DGCP, full)",
        part1.tyler.alpha,
        part1.tyler.R2
    )
    @printf(
        "  %-30s  %8.3f  %8.4f\n",
        "Scalar (DCP, full)",
        part1.scalar.alpha,
        part1.scalar.R2
    )

    if length(part4.nodes) >= 3
        alpha_ct, _, R2_ct = fit_power_law(part4.nodes, part4.times_ns)
        @printf("  %-30s  %8.3f  %8.4f\n", "Conic form generation", alpha_ct, R2_ct)
    end
    if length(part3.nodes) >= 3
        alpha_m, _, R2_m = fit_power_law(part3.nodes, part3.alloc_bytes)
        @printf("  %-30s  %8.3f  %8.4f\n", "Memory allocation", alpha_m, R2_m)
    end
    println()

    # Log-log data points for plotting
    println("Log-Log Data (for external plotting):")
    println("-"^60)
    println("# Karcher DGCP: log(nodes), log(time_us)")
    for i in eachindex(part1.karcher.nodes)
        @printf(
            "  %.4f, %.4f\n",
            log(part1.karcher.nodes[i]),
            log(part1.karcher.times_ns[i] / 1.0e3)
        )
    end
    println("# Tyler DGCP: log(nodes), log(time_us)")
    for i in eachindex(part1.tyler.nodes)
        @printf(
            "  %.4f, %.4f\n",
            log(part1.tyler.nodes[i]),
            log(part1.tyler.times_ns[i] / 1.0e3)
        )
    end
    println("# Scalar DCP: log(nodes), log(time_us)")
    for i in eachindex(part1.scalar.nodes)
        @printf(
            "  %.4f, %.4f\n",
            log(part1.scalar.nodes[i]),
            log(part1.scalar.times_ns[i] / 1.0e3)
        )
    end
    return println()
end

# ============================================================================
# Part 6: Matrix Size Independence Check
# ============================================================================

function run_part6_matrix_independence()
    println("="^72)
    println("PART 6: Matrix Size Independence (Sanity Check)")
    println("="^72)
    println()
    println("Verification time should NOT depend on matrix dimension n,")
    println("because matrices are numerical constants in the AST.")
    println("We fix m=4 terms and vary n.")
    println()

    m_fixed = 4
    dims = [3, 5, 8, 10, 15]

    @printf("  %-4s  %6s  %6s  %10s\n", "n", "nodes", "depth", "time(us)")
    println("  " * "-"^34)

    independence_data = []

    for n in dims
        expr, M = make_karcher_expr(m_fixed; n = n)
        nn = count_ast_nodes(expr)
        dd = ast_depth(expr)
        t_ns, _ = time_min(() -> analyze(expr, M))

        @printf("  %3d  %5d  %5d  %10.1f\n", n, nn, dd, t_ns / 1.0e3)
        push!(independence_data, (n = n, nodes = nn, depth = dd, time_ns = t_ns))
    end
    println()

    nodes_vec = [d.nodes for d in independence_data]
    times_vec = [d.time_ns for d in independence_data]
    node_range = maximum(nodes_vec) - minimum(nodes_vec)
    time_range = maximum(times_vec) / minimum(times_vec)

    @printf(
        "Node count range: %d - %d (%.1fx)\n",
        minimum(nodes_vec),
        maximum(nodes_vec),
        maximum(nodes_vec) / minimum(nodes_vec)
    )
    @printf(
        "Time range: %.1f - %.1f us (%.1fx)\n",
        minimum(times_vec) / 1.0e3,
        maximum(times_vec) / 1.0e3,
        time_range
    )

    if maximum(nodes_vec) / minimum(nodes_vec) < 1.5
        println("Confirmed: AST node count is independent of matrix dimension.")
        println("Varying matrix size does NOT create larger verification problems.")
    end
    println()

    return independence_data
end

# ============================================================================
# Main Entry Point
# ============================================================================

function main()
    println()
    println("*"^72)
    println("  Empirical Scaling Analysis for SymbolicAnalysis.jl")
    println("  Verification Algorithm Complexity")
    println("*"^72)
    println()
    @printf(
        "Configuration: matrix_dim=%d, warmup=%d, timing_iters=%d\n",
        MATRIX_DIM,
        WARMUP_ITERS,
        TIMING_ITERS
    )
    println("Julia version: $(VERSION)")
    println("Timing method: minimum of $(TIMING_ITERS) trials (time_ns)")
    println()

    part1 = run_part1_scaling()
    part2 = run_part2_phase_decomposition()
    part3 = run_part3_memory()
    part4 = run_part4_conic()
    run_part6_matrix_independence()
    run_part5_summary_table(part1, part2, part3, part4)

    println("*"^72)
    println("  Analysis Complete")
    println("*"^72)

    return (part1 = part1, part2 = part2, part3 = part3, part4 = part4)
end

# Run if executed directly
if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
