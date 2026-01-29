using Plots, DataFrames, CSV, Statistics
using SymbolicAnalysis, Manifolds, LinearAlgebra
using Random

"""
Simple, reliable DGCP benchmarking that extracts timing from your working approach
"""

Random.seed!(42)

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
    end
end

function warmup_and_benchmark(problem_type::String, size::Int; n_samples=10)
    """Warmup and benchmark with multiple samples"""
    
    M = SymmetricPositiveDefinite(size)
    
    # Warmup (5 runs)
    for _ in 1:5
        test_data = generate_test_data(size, problem_type)
        expr = create_expression(test_data, size, problem_type)
        SymbolicAnalysis.analyze(expr, M)
    end
    
    # Benchmark (multiple samples)
    times = Float64[]
    for _ in 1:n_samples
        test_data = generate_test_data(size, problem_type)
        expr = create_expression(test_data, size, problem_type)
        
        # Simple, reliable timing
        time_ms = @elapsed(SymbolicAnalysis.analyze(expr, M)) * 1000
        push!(times, time_ms)
    end
    
    return median(times)
end

function run_benchmark()
    """Run the benchmark and extract results"""
    
    println("="^60)
    println("DGCP VERIFICATION TIMING BENCHMARK") 
    println("="^60)
    
    # Problem configurations matching your ranges
    configs = [
        ("Tyler", "Tyler's M-Estimator", collect(5:5:40)),
        ("Karcher", "Karcher Mean", collect(25:25:200)), 
        ("LogDet", "Log-Determinant", collect(100:100:800))
    ]
    
    all_results = DataFrame(
        problem_type=String[],
        expression_name=String[],
        size=Int[], 
        median_time_ms=Float64[],
        success=Bool[]
    )
    
    for (problem_type, expr_name, sizes) in configs
        println("\n" * "="^50)
        println("BENCHMARKING: $expr_name")
        println("="^50)
        
        for size in sizes
            print("  Size $(size)×$(size)... ")
            flush(stdout)
            
            try
                median_time = warmup_and_benchmark(problem_type, size, n_samples=10)
                
                push!(all_results, (
                    problem_type=problem_type,
                    expression_name=expr_name,
                    size=size, 
                    median_time_ms=median_time,
                    success=true
                ))
                
                println("$(round(median_time, digits=3)) ms")
                
            catch e
                println("FAILED: $e")
                push!(all_results, (
                    problem_type=problem_type,
                    expression_name=expr_name,
                    size=size,
                    median_time_ms=NaN,
                    success=false
                ))
            end
        end
    end
    
    return all_results
end

function create_plots(results)
    """Create the performance plots"""
    
    # Save results
    CSV.write("dgcp_clean_benchmark_results.csv", results)
    println("\n✓ Results saved to: dgcp_clean_benchmark_results.csv")
    
    # Filter successful results
    successful = filter(row -> row.success, results)
    
    if nrow(successful) == 0
        println("❌ No successful results to plot")
        return
    end
    
    # Create individual plots
    expr_types = [
        ("Tyler's M-Estimator", :blue, :circle, "tyler_estimator_performance.png"),
        ("Karcher Mean", :red, :square, "karcher_mean_performance.png"),
        ("Log-Determinant", :green, :diamond, "logdet_performance.png")
    ]
    
    plots_created = []
    
    for (expr_name, color, marker, filename) in expr_types
        data = filter(row -> row.expression_name == expr_name, successful)
        
        if nrow(data) > 0
            # Determine if we need log scale
            use_log = expr_name == "Karcher Mean"
            
            p = plot(
                title="$expr_name Verification",
                xlabel="Matrix Size (n×n)",
                ylabel="Time (ms)",
                grid=true,
                legend=false,
                size=(600, 400),
                dpi=300,
                linewidth=4,
                markersize=8,
                guidefontsize=12,
                titlefontsize=14
            )
            
            if use_log
                plot!(p, yscale=:log10)
            end
            
            plot!(p, data.size, data.median_time_ms,
                  marker=marker,
                  color=color,
                  linewidth=4,
                  markersize=8)
            
            savefig(p, filename)
            push!(plots_created, p)
            println("✓ $expr_name plot saved: $filename")
        end
    end
    
    # Create combined plot if we have all three
    if length(plots_created) == 3
        combined = plot(plots_created..., 
                       layout=(1,3), 
                       size=(1200, 400),
                       plot_title="DGCP Performance Analysis")
        savefig(combined, "dgcp_three_panel.png")
        println("✓ Combined plot saved: dgcp_three_panel.png")
    end
    
    # Print summary
    println("\n" * "="^50)
    println("BENCHMARK SUMMARY")
    println("="^50)
    
    for expr_name in ["Tyler's M-Estimator", "Karcher Mean", "Log-Determinant"]
        data = filter(row -> row.expression_name == expr_name, successful)
        if nrow(data) > 0
            min_time = minimum(data.median_time_ms)
            max_time = maximum(data.median_time_ms) 
            mean_time = mean(data.median_time_ms)
            
            println("\n$expr_name:")
            println("  • $(nrow(data)) measurements")
            println("  • Range: $(round(min_time, digits=3))ms - $(round(max_time, digits=3))ms")
            println("  • Mean: $(round(mean_time, digits=3))ms")
        end
    end
end

# Main execution
function main()
    println("Simple DGCP Verification Benchmark")
    println("Measuring symbolic analysis time with reliable statistical sampling...")
    
    results = run_benchmark()
    create_plots(results)
    
    println("\n" * "="^50)
    println("BENCHMARK COMPLETE!")
    println("="^50)
end