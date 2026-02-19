"""
Experiment 3: Extended Verification Benchmarks with AST Complexity Metrics

This experiment extends the timing benchmarks to include symbolic complexity
metrics (AST node count, depth) to better understand verification performance.

Addresses:
- Reviewer 399: "symbolic complexity and verification time experiments"
- Reviewer 400: "Section 4.4 focuses exclusively on verification time for
  small to moderate-scale problem instances"
"""

using SymbolicAnalysis, Manifolds, LinearAlgebra
using Symbolics
using SymbolicUtils: iscall, arguments, operation
using Random
using Statistics
using Printf
using Test

Random.seed!(42)

#==============================================================================#
# AST Complexity Metrics
#==============================================================================#

"""
    count_ast_nodes(ex)

Count the total number of nodes in an expression tree.
Returns the number of operations + leaves in the symbolic expression.
"""
function count_ast_nodes(ex)
    ex = Symbolics.unwrap(ex)
    if !iscall(ex)
        return 1
    end
    return 1 + sum(count_ast_nodes(arg) for arg in arguments(ex); init=0)
end

"""
    ast_depth(ex)

Compute the maximum depth of an expression tree.
"""
function ast_depth(ex)
    ex = Symbolics.unwrap(ex)
    if !iscall(ex)
        return 1
    end
    args = arguments(ex)
    if isempty(args)
        return 1
    end
    return 1 + maximum(ast_depth(arg) for arg in args)
end

"""
    count_unique_operations(ex)

Count the number of unique operations in an expression.
"""
function count_unique_operations(ex)
    ops = Set{Any}()
    _collect_ops!(ops, ex)
    return length(ops)
end

function _collect_ops!(ops, ex)
    ex = Symbolics.unwrap(ex)
    if iscall(ex)
        push!(ops, operation(ex))
        for arg in arguments(ex)
            _collect_ops!(ops, arg)
        end
    end
end

#==============================================================================#
# Expression Generation
#==============================================================================#

function generate_test_data(size::Int, problem_type::String)
    if problem_type == "Tyler"
        A = randn(size, size)
        Sigma = A * A' + I
        xs = [randn(size) for _ in 1:min(10, size)]
        return (Sigma=Sigma, xs=xs)
    elseif problem_type == "Karcher"
        matrices = []
        for _ in 1:5
            A = randn(size, size)
            push!(matrices, A * A' + I)
        end
        return (matrices=matrices,)
    elseif problem_type == "LogDet"
        A = randn(size, size)
        return (matrix=A * A' + I,)
    elseif problem_type == "BrascampLieb"
        A = randn(size, size)
        A = A * A' + I
        return (A=A,)
    end
end

function create_expression(data, size::Int, problem_type::String)
    @variables X[1:size, 1:size]

    if problem_type == "Tyler"
        return sum(SymbolicAnalysis.log_quad_form(x, inv(X)) for x in data.xs) +
               (1/size) * logdet(X)
    elseif problem_type == "Karcher"
        M = SymmetricPositiveDefinite(size)
        return sum(Manifolds.distance(M, As, X)^2 for As in data.matrices)
    elseif problem_type == "LogDet"
        return logdet(X)
    elseif problem_type == "BrascampLieb"
        return logdet(SymbolicAnalysis.conjugation(X, data.A)) - logdet(X)
    end
end

#==============================================================================#
# Extended Benchmark with Complexity Metrics
#==============================================================================#

struct BenchmarkResult
    problem_type::String
    expression_name::String
    size::Int
    median_time_ms::Float64
    std_time_ms::Float64
    ast_nodes::Int
    ast_depth::Int
    unique_ops::Int
    memory_kb::Float64
end

function benchmark_with_complexity(problem_type::String, size::Int; n_samples=5)
    M = SymmetricPositiveDefinite(size)

    # Warmup
    for _ in 1:3
        test_data = generate_test_data(size, problem_type)
        expr = create_expression(test_data, size, problem_type)
        SymbolicAnalysis.analyze(expr, M)
    end

    # Benchmark with metrics
    times = Float64[]
    node_counts = Int[]
    depths = Int[]
    op_counts = Int[]
    allocations = Int[]

    for _ in 1:n_samples
        test_data = generate_test_data(size, problem_type)
        expr = create_expression(test_data, size, problem_type)

        # Measure complexity
        push!(node_counts, count_ast_nodes(expr))
        push!(depths, ast_depth(expr))
        push!(op_counts, count_unique_operations(expr))

        # Measure time and allocations
        alloc = @allocated begin
            time_ms = @elapsed(SymbolicAnalysis.analyze(expr, M)) * 1000
        end
        push!(times, time_ms)
        push!(allocations, alloc)
    end

    return BenchmarkResult(
        problem_type,
        problem_type,
        size,
        median(times),
        length(times) > 1 ? std(times) : 0.0,
        Int(median(node_counts)),
        Int(median(depths)),
        Int(median(op_counts)),
        median(allocations) / 1024,
    )
end

function run_extended_benchmark()
    println("="^70)
    println("EXPERIMENT 3: Extended DGCP Verification Benchmarks")
    println("="^70)
    println()
    println("Measuring verification time + symbolic complexity metrics")
    println()

    configs = [
        ("Tyler", "Tyler's M-Estimator", collect(5:5:30)),
        ("Karcher", "Karcher Mean", collect(25:25:150)),
        ("LogDet", "Log-Determinant", collect(50:50:400)),
        ("BrascampLieb", "Brascamp-Lieb", collect(5:5:30)),
    ]

    all_results = BenchmarkResult[]

    for (problem_type, expr_name, sizes) in configs
        println("\nBenchmarking: $expr_name")
        println("-"^50)

        for size in sizes
            print("  Size $(size)x$(size)... ")
            flush(stdout)

            try
                result = benchmark_with_complexity(problem_type, size, n_samples=5)
                push!(all_results, result)

                println(@sprintf("%.3f ms, %d nodes, depth %d, %d ops",
                    result.median_time_ms, result.ast_nodes,
                    result.ast_depth, result.unique_ops))

            catch e
                println("FAILED: $e")
            end
        end
    end

    return all_results
end

#==============================================================================#
# Complexity Analysis (text-based, no plotting dependencies)
#==============================================================================#

function run_complexity_analysis(results::Vector{BenchmarkResult})
    println()
    println("=" ^ 70)
    println("COMPLEXITY ANALYSIS")
    println("=" ^ 70)

    # Full results table
    println()
    println("Full Results Table:")
    println("-" ^ 90)
    println(rpad("Problem", 18), " | ",
            rpad("Size", 5), " | ",
            rpad("Time(ms)", 10), " | ",
            rpad("Nodes", 7), " | ",
            rpad("Depth", 6), " | ",
            rpad("Ops", 5), " | ",
            "Mem(KB)")
    println("-" ^ 90)

    for r in results
        println(
            rpad(r.problem_type, 18), " | ",
            rpad(string(r.size), 5), " | ",
            rpad(@sprintf("%.3f", r.median_time_ms), 10), " | ",
            rpad(string(r.ast_nodes), 7), " | ",
            rpad(string(r.ast_depth), 6), " | ",
            rpad(string(r.unique_ops), 5), " | ",
            @sprintf("%.1f", r.memory_kb)
        )
    end
    println("-" ^ 90)

    # Per-problem-type analysis
    problem_types = unique(r.problem_type for r in results)
    for ptype in problem_types
        pdata = filter(r -> r.problem_type == ptype, results)
        if length(pdata) < 2
            continue
        end

        println()
        println("$ptype:")
        println("  Size range: $(minimum(r.size for r in pdata)) - $(maximum(r.size for r in pdata))")
        println("  Node count range: $(minimum(r.ast_nodes for r in pdata)) - $(maximum(r.ast_nodes for r in pdata))")
        println("  Depth range: $(minimum(r.ast_depth for r in pdata)) - $(maximum(r.ast_depth for r in pdata))")
        println("  Time range: $(@sprintf("%.3f", minimum(r.median_time_ms for r in pdata))) - $(@sprintf("%.3f", maximum(r.median_time_ms for r in pdata))) ms")

        # Estimate scaling exponent via log-log linear regression
        if length(pdata) >= 3
            x = log.(Float64[r.ast_nodes for r in pdata])
            y = log.(Float64[r.median_time_ms for r in pdata])
            n = length(x)
            denom = n * sum(x .^ 2) - sum(x)^2
            if abs(denom) > 1e-10
                slope = (n * sum(x .* y) - sum(x) * sum(y)) / denom
                println("  Approximate scaling (time vs nodes): O(nodes^$(@sprintf("%.2f", slope)))")
            end

            # Depth-based scaling
            xd = log.(Float64[r.ast_depth for r in pdata])
            yd = y
            nd = length(xd)
            denomd = nd * sum(xd .^ 2) - sum(xd)^2
            if abs(denomd) > 1e-10
                sloped = (nd * sum(xd .* yd) - sum(xd) * sum(yd)) / denomd
                println("  Approximate scaling (time vs depth): O(depth^$(@sprintf("%.2f", sloped)))")
            end
        end
    end

    # Depth vs time table (grouped by depth)
    println()
    println("AST Depth vs Verification Time (all problems):")
    println("-" ^ 50)
    println(rpad("Depth", 8), " | ", rpad("Avg Time (ms)", 15), " | ", "Count")
    println("-" ^ 50)
    depths_seen = sort(unique(r.ast_depth for r in results))
    for d in depths_seen
        ddata = filter(r -> r.ast_depth == d, results)
        avg_time = mean(r.median_time_ms for r in ddata)
        println(rpad(string(d), 8), " | ",
                rpad(@sprintf("%.3f", avg_time), 15), " | ",
                string(length(ddata)))
    end
    println("-" ^ 50)
end

#==============================================================================#
# Main
#==============================================================================#

function main()
    println("Extended DGCP Verification Benchmark")
    println("Measuring symbolic complexity + verification time...")
    println()

    results = run_extended_benchmark()
    run_complexity_analysis(results)

    println()
    println("=" ^ 70)
    println("EXTENDED BENCHMARK COMPLETE")
    println("=" ^ 70)

    return results
end

#==============================================================================#
# Tests
#==============================================================================#

@testset "Extended Benchmark" begin
    @testset "AST Complexity Metrics" begin
        @variables X[1:3, 1:3]
        M = SymmetricPositiveDefinite(3)

        expr = logdet(X)
        @test count_ast_nodes(expr) >= 1
        @test ast_depth(expr) >= 1
        @test count_unique_operations(expr) >= 1

        A = randn(3, 3)
        A = A * A' + I
        expr2 = Manifolds.distance(M, A, X)^2
        @test count_ast_nodes(expr2) > count_ast_nodes(expr)
    end

    @testset "Benchmark Small Problem" begin
        result = benchmark_with_complexity("LogDet", 5, n_samples=3)
        @test result.median_time_ms > 0
        @test result.ast_nodes >= 1
        @test result.ast_depth >= 1
        @test result.memory_kb > 0
    end
end

# Run if executed directly
if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
