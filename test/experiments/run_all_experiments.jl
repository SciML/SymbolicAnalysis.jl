"""
Run all experiments and save results (plots + CSV tables) to test/experiments/results/.
"""

using SymbolicAnalysis
using Manifolds
using Symbolics
using Symbolics: unwrap
using SymbolicUtils: iscall, arguments, operation
using LinearAlgebra
using Random
using Statistics
using Printf
using CairoMakie
using CSV
using DataFrames
using Test

Random.seed!(42)

const RESULTS_DIR = joinpath(@__DIR__, "results")
mkpath(RESULTS_DIR)

println("Results will be saved to: $RESULTS_DIR")
println()

#==============================================================================#
# Helpers from extended_benchmark.jl
#==============================================================================#

function count_ast_nodes(ex)
    ex = Symbolics.unwrap(ex)
    iscall(ex) || return 1
    return 1 + sum(count_ast_nodes(arg) for arg in arguments(ex); init=0)
end

function ast_depth(ex)
    ex = Symbolics.unwrap(ex)
    iscall(ex) || return 1
    args = arguments(ex)
    isempty(args) && return 1
    return 1 + maximum(ast_depth(arg) for arg in args)
end

function time_verification(f::Function, n_samples::Int=7)
    f()  # warmup
    times = [(@elapsed f()) for _ in 1:n_samples]
    return sort(times)[div(n_samples, 2) + 1]
end

#==============================================================================#
# 1. DCP vs DGCP Scope Comparison
#==============================================================================#

function run_and_save_scope()
    println("="^70)
    println("1. DCP vs DGCP Scope Comparison")
    println("="^70)

    @variables X[1:5, 1:5]
    M = SymmetricPositiveDefinite(5)
    A = let B = randn(5,5); B*B' + I end
    xs = [randn(5) for _ in 1:3]
    As = [let B = randn(5,5); B*B' + I end for _ in 1:3]

    cases = [
        ("logdet(X)",             logdet(X)),
        ("tr(inv(X))",            tr(inv(X))),
        ("distance²",            Manifolds.distance(M, A, X)^2),
        ("S-divergence",          SymbolicAnalysis.sdivergence(X, A)),
        ("logdet(A'X⁻¹A)",       logdet(SymbolicAnalysis.conjugation(inv(X), A))),
        ("Tyler M-Est",           sum(SymbolicAnalysis.log_quad_form(x, inv(X)) for x in xs) + (1/5)*logdet(X)),
        ("Karcher Mean",          sum(Manifolds.distance(M, Ai, X)^2 for Ai in As)),
    ]

    rows = []
    for (name, expr) in cases
        expr_u = unwrap(expr)
        r = analyze(expr_u, M)
        is_gcvx = r.gcurvature in (SymbolicAnalysis.GConvex, SymbolicAnalysis.GLinear)
        is_ecvx = r.curvature in (SymbolicAnalysis.Convex, SymbolicAnalysis.Affine)
        push!(rows, (Expression=name, DGCP=string(r.gcurvature),
                     EuclConvex=is_ecvx, GConvex=is_gcvx))
        println("  $name → DGCP=$(r.gcurvature), Eucl=$(r.curvature)")
    end

    df = DataFrame(rows)
    CSV.write(joinpath(RESULTS_DIR, "scope_comparison.csv"), df)
    println("  → Saved scope_comparison.csv")
    println()
    return df
end

#==============================================================================#
# 2. Timing Comparison (DCP vs DGCP overhead)
#==============================================================================#

function run_and_save_timing()
    println("="^70)
    println("2. DCP vs DGCP Timing Comparison")
    println("="^70)

    @variables X[1:5, 1:5]
    M = SymmetricPositiveDefinite(5)
    A = let B = randn(5,5); B*B' + I end

    cases = [
        ("logdet(X)",        logdet(X) |> unwrap,                     true),
        ("tr(X)",            tr(X) |> unwrap,                         true),
        ("tr(inv(X))",       tr(inv(X)) |> unwrap,                    true),
        ("-logdet(X)",       -logdet(X) |> unwrap,                    true),
        ("distance²",       Manifolds.distance(M, A, X)^2 |> unwrap, false),
        ("S-divergence",     SymbolicAnalysis.sdivergence(X, A) |> unwrap, false),
    ]

    rows = []
    for (name, expr, both) in cases
        dcp_t  = time_verification(7) do; analyze(expr); end
        dgcp_t = time_verification(7) do; analyze(expr, M); end
        overhead = dgcp_t / dcp_t
        push!(rows, (Function=name, DCP_us=dcp_t*1e6, DGCP_us=dgcp_t*1e6,
                     Overhead=overhead, BothVerify=both))
        println(@sprintf("  %-20s DCP=%8.1f us  DGCP=%8.1f us  overhead=%.2fx",
                name, dcp_t*1e6, dgcp_t*1e6, overhead))
    end

    df = DataFrame(rows)
    CSV.write(joinpath(RESULTS_DIR, "timing_comparison.csv"), df)

    # Plot: overhead bar chart
    n_funcs = nrow(df)
    xs = 1:n_funcs
    fig = Figure(size = (700, 400))
    ax = Axis(fig[1, 1],
        ylabel = "DGCP / DCP Overhead",
        xlabel = "Function",
        title = "DGCP vs DCP Verification Overhead",
        xticks = (collect(xs), df.Function),
        xticklabelrotation = pi / 7,
    )
    barplot!(ax, collect(xs), df.Overhead; color = :steelblue)
    hlines!(ax, [1.0]; linestyle = :dash, color = :red)
    ylims!(ax, 0, max(2.0, maximum(df.Overhead) * 1.2))
    save(joinpath(RESULTS_DIR, "timing_overhead.png"), fig)
    println("  → Saved timing_comparison.csv, timing_overhead.png")
    println()
    return df
end

#==============================================================================#
# 3. Scaling Analysis
#==============================================================================#

function run_and_save_scaling()
    println("="^70)
    println("3. Scaling Analysis")
    println("="^70)

    rows = []

    # Part A: vary matrix size
    println("  Part A: Varying matrix size (3 terms)")
    for n in [3, 5, 8, 10]
        @variables Xn[1:n, 1:n]
        M = SymmetricPositiveDefinite(n)
        As = [let B = randn(n,n); B*B' + I end for _ in 1:3]
        expr = sum(Manifolds.distance(M, Ai, Xn)^2 for Ai in As) |> unwrap
        dcp_t  = time_verification(5) do; analyze(expr); end
        dgcp_t = time_verification(5) do; analyze(expr, M); end
        push!(rows, (Problem="Karcher", MatrixSize=n, Terms=3,
                     DCP_us=dcp_t*1e6, DGCP_us=dgcp_t*1e6, Overhead=dgcp_t/dcp_t))
        println(@sprintf("    n=%2d: DCP=%8.1f us  DGCP=%8.1f us", n, dcp_t*1e6, dgcp_t*1e6))
    end

    # Part B: vary number of terms
    println("  Part B: Varying terms (n=5)")
    for nt in [1, 3, 5, 10]
        n = 5
        @variables Xn[1:n, 1:n]
        M = SymmetricPositiveDefinite(n)
        As = [let B = randn(n,n); B*B' + I end for _ in 1:nt]
        expr = sum(Manifolds.distance(M, Ai, Xn)^2 for Ai in As) |> unwrap
        dcp_t  = time_verification(5) do; analyze(expr); end
        dgcp_t = time_verification(5) do; analyze(expr, M); end
        push!(rows, (Problem="Karcher", MatrixSize=n, Terms=nt,
                     DCP_us=dcp_t*1e6, DGCP_us=dgcp_t*1e6, Overhead=dgcp_t/dcp_t))
        println(@sprintf("    terms=%2d: DCP=%8.1f us  DGCP=%8.1f us", nt, dcp_t*1e6, dgcp_t*1e6))
    end

    # Part C: Tyler's M-estimator varying vectors
    println("  Part C: Tyler's M-estimator (n=5)")
    for nv in [1, 3, 5, 8]
        n = 5
        @variables Xn[1:n, 1:n]
        M = SymmetricPositiveDefinite(n)
        xs = [randn(n) for _ in 1:nv]
        expr = (sum(SymbolicAnalysis.log_quad_form(x, inv(Xn)) for x in xs) +
                (1/n)*logdet(Xn)) |> unwrap
        dcp_t  = time_verification(5) do; analyze(expr); end
        dgcp_t = time_verification(5) do; analyze(expr, M); end
        push!(rows, (Problem="Tyler", MatrixSize=n, Terms=nv,
                     DCP_us=dcp_t*1e6, DGCP_us=dgcp_t*1e6, Overhead=dgcp_t/dcp_t))
        println(@sprintf("    vectors=%2d: DCP=%8.1f us  DGCP=%8.1f us", nv, dcp_t*1e6, dgcp_t*1e6))
    end

    df = DataFrame(rows)
    CSV.write(joinpath(RESULTS_DIR, "scaling_analysis.csv"), df)

    # Plot: scaling with terms (Karcher n=5)
    karcher_terms = filter(r -> r.Problem == "Karcher" && r.MatrixSize == 5, df)
    fig1 = Figure(size = (600, 400))
    ax1 = Axis(fig1[1, 1],
        xlabel = "Number of terms", ylabel = "Time (us)",
        title = "Verification Time vs Problem Size (Karcher, n=5)")
    scatterlines!(ax1, karcher_terms.Terms, karcher_terms.DCP_us;
        label = "DCP", marker = :circle, linewidth = 2)
    scatterlines!(ax1, karcher_terms.Terms, karcher_terms.DGCP_us;
        label = "DGCP", marker = :rect, linewidth = 2)
    axislegend(ax1; position = :lt)
    save(joinpath(RESULTS_DIR, "scaling_terms.png"), fig1)

    # Plot: scaling with matrix size (Karcher 3 terms)
    karcher_size = filter(r -> r.Problem == "Karcher" && r.Terms == 3, df)
    fig2 = Figure(size = (600, 400))
    ax2 = Axis(fig2[1, 1],
        xlabel = "Matrix size n", ylabel = "Time (us)",
        title = "Verification Time vs Matrix Size (Karcher, 3 terms)")
    scatterlines!(ax2, karcher_size.MatrixSize, karcher_size.DCP_us;
        label = "DCP", marker = :circle, linewidth = 2)
    scatterlines!(ax2, karcher_size.MatrixSize, karcher_size.DGCP_us;
        label = "DGCP", marker = :rect, linewidth = 2)
    axislegend(ax2; position = :lt)
    save(joinpath(RESULTS_DIR, "scaling_matrix_size.png"), fig2)

    println("  → Saved scaling_analysis.csv, scaling_terms.png, scaling_matrix_size.png")
    println()
    return df
end

#==============================================================================#
# 4. Extended Benchmark (AST complexity)
#==============================================================================#

function run_and_save_benchmark()
    println("="^70)
    println("4. Extended Benchmark (AST Complexity)")
    println("="^70)

    configs = [
        ("Tyler",        collect(5:5:30)),
        ("Karcher",      collect(25:25:150)),
        ("LogDet",       collect(50:50:400)),
        ("BrascampLieb", collect(5:5:30)),
    ]

    rows = []
    for (ptype, sizes) in configs
        for sz in sizes
            @variables Xb[1:sz, 1:sz]
            M = SymmetricPositiveDefinite(sz)

            expr = if ptype == "Tyler"
                xs = [randn(sz) for _ in 1:min(10, sz)]
                sum(SymbolicAnalysis.log_quad_form(x, inv(Xb)) for x in xs) + (1/sz)*logdet(Xb)
            elseif ptype == "Karcher"
                As = [let B = randn(sz,sz); B*B' + I end for _ in 1:5]
                sum(Manifolds.distance(M, Ai, Xb)^2 for Ai in As)
            elseif ptype == "LogDet"
                logdet(Xb)
            elseif ptype == "BrascampLieb"
                A = let B = randn(sz,sz); B*B' + I end
                logdet(SymbolicAnalysis.conjugation(Xb, A)) - logdet(Xb)
            end

            expr_u = unwrap(expr)

            # Warmup
            analyze(expr_u, M)
            analyze(expr_u, M)

            nodes = count_ast_nodes(expr_u)
            depth = ast_depth(expr_u)

            t = median([@elapsed(analyze(expr_u, M)) for _ in 1:5]) * 1000
            alloc = @allocated(analyze(expr_u, M))

            push!(rows, (Problem=ptype, Size=sz, Time_ms=t, Nodes=nodes,
                         Depth=depth, Memory_KB=alloc/1024))
            println(@sprintf("  %-15s %3dx%-3d  %.3f ms  %3d nodes  depth %d",
                    ptype, sz, sz, t, nodes, depth))
        end
    end

    df = DataFrame(rows)
    CSV.write(joinpath(RESULTS_DIR, "extended_benchmark.csv"), df)

    # Plot: time vs AST nodes by problem type
    fig = Figure(size = (700, 450))
    ax = Axis(fig[1, 1], title = "Verification Time vs AST Nodes",
        xlabel = "AST Nodes", ylabel = "Time (ms)")
    for (i, ptype) in enumerate(unique(df.Problem))
        sub = filter(r -> r.Problem == ptype, df)
        scatter!(ax, sub.Nodes, sub.Time_ms; label = ptype, markersize = 8)
    end
    axislegend(ax; position = :lt)
    save(joinpath(RESULTS_DIR, "benchmark_nodes_vs_time.png"), fig)

    # Plot: time vs matrix size by problem type
    fig2 = Figure(size = (700, 450))
    ax2 = Axis(fig2[1, 1], title = "Verification Time vs Matrix Size",
        xlabel = "Matrix Size n", ylabel = "Time (ms)")
    for (i, ptype) in enumerate(unique(df.Problem))
        sub = filter(r -> r.Problem == ptype, df)
        scatterlines!(ax2, sub.Size, sub.Time_ms;
            label = ptype, marker = :circle, linewidth = 2)
    end
    axislegend(ax2; position = :lt)
    save(joinpath(RESULTS_DIR, "benchmark_size_vs_time.png"), fig2)

    println("  → Saved extended_benchmark.csv, benchmark_nodes_vs_time.png, benchmark_size_vs_time.png")
    println()
    return df
end

#==============================================================================#
# 5. Expert Examples
#==============================================================================#

function run_and_save_expert()
    println("="^70)
    println("5. Expert Examples")
    println("="^70)

    @variables X[1:5, 1:5]
    M = SymmetricPositiveDefinite(5)

    A = let B = randn(5,5); B*B' + I end
    xs = [randn(5) for _ in 1:3]
    As = [let B = randn(5,5); B*B' + I end for _ in 1:3]

    cases = [
        ("Tyler M-Est",      "Hard",
         sum(SymbolicAnalysis.log_quad_form(x, inv(X)) for x in xs) + (1/5)*logdet(X)),
        ("Brascamp-Lieb",    "Hard",
         logdet(SymbolicAnalysis.conjugation(X, A)) - logdet(X)),
        ("S-Divergence Sum", "Medium",
         SymbolicAnalysis.sdivergence(X, A) + SymbolicAnalysis.sdivergence(X, Matrix(I(5)*1.0))),
        ("Karcher Mean",     "Hard",
         sum(Manifolds.distance(M, Ai, X)^2 for Ai in As)),
        ("Diagonal Loading", "Medium",
         tr(inv(X)) + logdet(X) + 0.1*tr(X)),
        ("Spectral Fn",      "Hard",
         SymbolicAnalysis.eigsummax(log(X), 2)),
    ]

    rows = []
    for (name, difficulty, expr) in cases
        expr_u = unwrap(expr)
        # Warmup
        analyze(expr_u, M)
        t_ms = (@elapsed analyze(expr_u, M)) * 1000
        r = analyze(expr_u, M)
        push!(rows, (Case=name, Difficulty=difficulty,
                     Result=string(r.gcurvature), Time_ms=t_ms))
        println(@sprintf("  %-20s [%s] → %s  (%.2f ms)", name, difficulty,
                r.gcurvature, t_ms))
    end

    df = DataFrame(rows)
    CSV.write(joinpath(RESULTS_DIR, "expert_examples.csv"), df)

    # Plot: expert verification times
    n_cases = nrow(df)
    ys = 1:n_cases
    colors = [d == "Hard" ? :firebrick : :orange for d in df.Difficulty]
    fig = Figure(size = (700, 400))
    ax = Axis(fig[1, 1],
        xlabel = "Time (ms)",
        title = "Expert-Level DGCP Verification Time",
        yticks = (collect(ys), df.Case),
        xticklabelrotation = pi / 9,
    )
    barplot!(ax, collect(ys), df.Time_ms; direction = :x, color = colors)
    save(joinpath(RESULTS_DIR, "expert_times.png"), fig)

    println("  → Saved expert_examples.csv, expert_times.png")
    println()
    return df
end

#==============================================================================#
# 6. MLE Experiment
#==============================================================================#

function run_and_save_mle()
    println("="^70)
    println("6. MLE Experiment")
    println("="^70)

    rows = []

    # Frechet Mean
    for (n, m) in [(3,3), (3,5), (3,10), (5,3), (5,5), (5,10)]
        @variables Xm[1:n, 1:n]
        M = SymmetricPositiveDefinite(n)
        samples = [let B = randn(n,n); B*B' + I end for _ in 1:m]
        expr = sum(Manifolds.distance(M, S, Xm)^2 for S in samples) |> unwrap

        dgcp_t = @elapsed dgcp_r = analyze(expr, M)
        dcp_t  = @elapsed dcp_r  = analyze(expr)

        is_gcvx = dgcp_r.gcurvature in (SymbolicAnalysis.GConvex, SymbolicAnalysis.GLinear)
        is_ecvx = dcp_r.curvature in (SymbolicAnalysis.Convex, SymbolicAnalysis.Affine)

        push!(rows, (Problem="Frechet", n=n, Samples=m,
                     GConvex=is_gcvx, EuclConvex=is_ecvx,
                     DGCP_s=dgcp_t, DCP_s=dcp_t))
        println(@sprintf("  Frechet n=%d m=%2d: DGCP=%.4fs DCP=%.4fs gcvx=%s",
                n, m, dgcp_t, dcp_t, is_gcvx))
    end

    # Tyler
    for (n, k) in [(3,3), (3,5), (5,3), (5,5)]
        @variables Xm[1:n, 1:n]
        M = SymmetricPositiveDefinite(n)
        xs = [randn(n) for _ in 1:k]
        expr = (sum(SymbolicAnalysis.log_quad_form(x, inv(Xm)) for x in xs) +
                (1/n)*logdet(Xm)) |> unwrap

        dgcp_t = @elapsed dgcp_r = analyze(expr, M)
        dcp_t  = @elapsed dcp_r  = analyze(expr)

        is_gcvx = dgcp_r.gcurvature in (SymbolicAnalysis.GConvex, SymbolicAnalysis.GLinear)
        is_ecvx = dcp_r.curvature in (SymbolicAnalysis.Convex, SymbolicAnalysis.Affine)

        push!(rows, (Problem="Tyler", n=n, Samples=k,
                     GConvex=is_gcvx, EuclConvex=is_ecvx,
                     DGCP_s=dgcp_t, DCP_s=dcp_t))
        println(@sprintf("  Tyler   n=%d k=%2d: DGCP=%.4fs DCP=%.4fs gcvx=%s",
                n, k, dgcp_t, dcp_t, is_gcvx))
    end

    df = DataFrame(rows)
    CSV.write(joinpath(RESULTS_DIR, "mle_experiment.csv"), df)

    # Plot: MLE verification times
    labels = df.Problem .* " n=" .* string.(df.n) .* " k=" .* string.(df.Samples)
    n_mle = nrow(df)
    xs = 1:n_mle
    fig = Figure(size = (800, 450))
    ax = Axis(fig[1, 1],
        ylabel = "Time (ms)",
        title = "MLE Verification Time (DGCP)",
        xticks = (collect(xs), labels),
        xticklabelrotation = pi / 5,
    )
    barplot!(ax, collect(xs), df.DGCP_s .* 1000; color = :steelblue)
    save(joinpath(RESULTS_DIR, "mle_times.png"), fig)

    println("  → Saved mle_experiment.csv, mle_times.png")
    println()
    return df
end

#==============================================================================#
# Run All
#==============================================================================#

println("="^70)
println("RUNNING ALL EXPERIMENTS")
println("="^70)
println()

scope_df    = run_and_save_scope()
timing_df   = run_and_save_timing()
scaling_df  = run_and_save_scaling()
bench_df    = run_and_save_benchmark()
expert_df   = run_and_save_expert()
mle_df      = run_and_save_mle()

println("="^70)
println("ALL EXPERIMENTS COMPLETE")
println("="^70)
println()
println("Results saved to: $RESULTS_DIR")
println("  CSV files:")
for f in filter(f -> endswith(f, ".csv"), readdir(RESULTS_DIR))
    println("    $f")
end
println("  Plots:")
for f in filter(f -> endswith(f, ".png"), readdir(RESULTS_DIR))
    println("    $f")
end
