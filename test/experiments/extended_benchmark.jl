"""
Experiment 3: Extended Verification Benchmarks with AST Complexity Metrics

This experiment extends the timing benchmarks to include symbolic complexity
metrics (AST node count, depth) to better understand verification performance.

Addresses:
- Reviewer 399: "symbolic complexity and verification time experiments"
- Reviewer 400: "Section 4.4 focuses exclusively on verification time for 
  small to moderate-scale problem instances"
"""

using Plots, DataFrames, CSV, Statistics
using SymbolicAnalysis, Manifolds, LinearAlgebra
using Symbolics
using SymbolicUtils: iscall, arguments, operation
using Random

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
# Expression Generation (from original benchmark)
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

function benchmark_with_complexity(problem_type::String, size::Int; n_samples=5)
    """Benchmark with AST complexity metrics"""
    
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
    allocations = Int[]
    
    for _ in 1:n_samples
        test_data = generate_test_data(size, problem_type)
        expr = create_expression(test_data, size, problem_type)
        
        # Measure complexity
        push!(node_counts, count_ast_nodes(expr))
        push!(depths, ast_depth(expr))
        
        # Measure time and allocations
        alloc = @allocated begin
            time_ms = @elapsed(SymbolicAnalysis.analyze(expr, M)) * 1000
        end
        push!(times, time_ms)
        push!(allocations, alloc)
    end
    
    return (
        median_time_ms = median(times),
        median_nodes = median(node_counts),
        median_depth = median(depths),
        median_alloc_kb = median(allocations) / 1024,
        std_time_ms = std(times)
    )
end

function run_extended_benchmark()
    """Run extended benchmark with complexity metrics"""
    
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
    
    all_results = DataFrame(
        problem_type = String[],
        expression_name = String[],
        size = Int[],
        median_time_ms = Float64[],
        std_time_ms = Float64[],
        ast_nodes = Int[],
        ast_depth = Int[],
        memory_kb = Float64[]
    )
    
    for (problem_type, expr_name, sizes) in configs
        println("\nBenchmarking: $expr_name")
        println("-"^50)
        
        for size in sizes
            print("  Size $(size)×$(size)... ")
            flush(stdout)
            
            try
                result = benchmark_with_complexity(problem_type, size, n_samples=5)
                
                push!(all_results, (
                    problem_type = problem_type,
                    expression_name = expr_name,
                    size = size,
                    median_time_ms = result.median_time_ms,
                    std_time_ms = result.std_time_ms,
                    ast_nodes = Int(result.median_nodes),
                    ast_depth = Int(result.median_depth),
                    memory_kb = result.median_alloc_kb
                ))
                
                println("$(round(result.median_time_ms, digits=3)) ms, " *
                        "$(Int(result.median_nodes)) nodes, " *
                        "depth $(Int(result.median_depth))")
                
            catch e
                println("FAILED: $e")
            end
        end
    end
    
    return all_results
end

#==============================================================================#
# Plotting Functions
#==============================================================================#

function create_complexity_plots(results)
    """Create plots showing time vs complexity"""
    
    # Save results
    CSV.write("extended_benchmark_results.csv", results)
    println("\n✓ Results saved to: extended_benchmark_results.csv")
    
    # Filter successful results
    successful = filter(row -> !isnan(row.median_time_ms), results)
    
    if nrow(successful) == 0
        println("❌ No results to plot")
        return
    end
    
    # Plot 1: Time vs AST Node Count
    p1 = plot(
        title = "Verification Time vs Expression Complexity",
        xlabel = "AST Node Count",
        ylabel = "Time (ms)",
        legend = :topleft,
        grid = true,
        size = (700, 500),
        dpi = 300
    )
    
    colors = Dict(
        "Tyler's M-Estimator" => :blue,
        "Karcher Mean" => :red,
        "Log-Determinant" => :green,
        "Brascamp-Lieb" => :purple
    )
    markers = Dict(
        "Tyler's M-Estimator" => :circle,
        "Karcher Mean" => :square,
        "Log-Determinant" => :diamond,
        "Brascamp-Lieb" => :star5
    )
    
    for expr_name in unique(successful.expression_name)
        data = filter(row -> row.expression_name == expr_name, successful)
        if nrow(data) > 0
            scatter!(p1, data.ast_nodes, data.median_time_ms,
                    label = expr_name,
                    color = get(colors, expr_name, :gray),
                    marker = get(markers, expr_name, :circle),
                    markersize = 6)
        end
    end
    
    savefig(p1, "complexity_vs_time.png")
    println("✓ Plot saved: complexity_vs_time.png")
    
    # Plot 2: Time vs Matrix Size (by problem type)
    p2 = plot(
        title = "Verification Time vs Matrix Size",
        xlabel = "Matrix Size (n)",
        ylabel = "Time (ms, log scale)",
        legend = :topleft,
        grid = true,
        yscale = :log10,
        size = (700, 500),
        dpi = 300
    )
    
    for expr_name in unique(successful.expression_name)
        data = filter(row -> row.expression_name == expr_name, successful)
        if nrow(data) > 0
            plot!(p2, data.size, data.median_time_ms,
                  label = expr_name,
                  color = get(colors, expr_name, :gray),
                  marker = get(markers, expr_name, :circle),
                  linewidth = 2,
                  markersize = 5)
        end
    end
    
    savefig(p2, "size_vs_time.png")
    println("✓ Plot saved: size_vs_time.png")
    
    # Summary statistics
    println("\n" * "="^70)
    println("COMPLEXITY ANALYSIS SUMMARY")
    println("="^70)
    
    for expr_name in unique(successful.expression_name)
        data = filter(row -> row.expression_name == expr_name, successful)
        if nrow(data) > 0
            println("\n$expr_name:")
            println("  • Size range: $(minimum(data.size)) - $(maximum(data.size))")
            println("  • Node count range: $(minimum(data.ast_nodes)) - $(maximum(data.ast_nodes))")
            println("  • Time range: $(round(minimum(data.median_time_ms), digits=3)) - " *
                    "$(round(maximum(data.median_time_ms), digits=3)) ms")
            
            # Estimate scaling
            if nrow(data) >= 3
                # Simple linear regression on log-log
                x = log.(data.ast_nodes)
                y = log.(data.median_time_ms)
                n = length(x)
                slope = (n * sum(x .* y) - sum(x) * sum(y)) / (n * sum(x.^2) - sum(x)^2)
                println("  • Approximate scaling: O(n^$(round(slope, digits=2)))")
            end
        end
    end
end

#==============================================================================#
# Main
#==============================================================================#

function main()
    println("Extended DGCP Verification Benchmark")
    println("Measuring symbolic complexity + verification time...")
    println()
    
    results = run_extended_benchmark()
    create_complexity_plots(results)
    
    println("\n" * "="^70)
    println("EXTENDED BENCHMARK COMPLETE!")
    println("="^70)
    
    return results
end

# Run if executed directly
if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
