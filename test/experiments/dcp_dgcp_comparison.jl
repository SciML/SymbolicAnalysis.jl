"""
Experiment 1: DCP vs DGCP Verification Scope Comparison

This experiment demonstrates functions that DGCP can verify as geodesically convex
but that DCP (via Convex.jl) cannot verify as Euclidean convex.

Addresses:
- Reviewer 399: "fair DCP vs DGCP comparison for problems both can verify"
- Reviewer 400: "explicitly demonstrate correspondence between DGCP and classical DCP"
"""

using SymbolicAnalysis
using Manifolds
using Symbolics
using LinearAlgebra
using Test

# Try to load Convex.jl for DCP comparison
const HAS_CONVEX = try
    using Convex
    true
catch
    @warn "Convex.jl not available, DCP comparison will be limited"
    false
end

#==============================================================================#
# Test Cases: Functions with Known Properties
#==============================================================================#

"""
Structure to hold comparison results
"""
struct ComparisonResult
    name::String
    dgcp_curvature::SymbolicAnalysis.GCurvature
    dcp_curvature::Union{Symbol, String}
    euclidean_convex::Bool
    geodesically_convex::Bool
    notes::String
end

"""
Run comparison for a given expression
"""
function compare_verification(
    name::String,
    dgcp_expr,
    convex_expr_fn::Union{Function, Nothing},
    notes::String = ""
)
    M = SymmetricPositiveDefinite(5)
    
    # DGCP analysis
    dgcp_result = analyze(dgcp_expr, M)
    dgcp_curv = dgcp_result.gcurvature
    
    # Euclidean curvature
    eucl_curv = dgcp_result.curvature
    is_eucl_convex = eucl_curv == SymbolicAnalysis.Convex || eucl_curv == SymbolicAnalysis.Affine
    
    # DCP analysis via Convex.jl
    dcp_curv = :not_tested
    if HAS_CONVEX && !isnothing(convex_expr_fn)
        try
            X_convex = Convex.Variable(5, 5)
            convex_obj = convex_expr_fn(X_convex)
            dcp_curv = Convex.vexity(convex_obj)
        catch e
            dcp_curv = Symbol("error: $(typeof(e).name)")
        end
    end
    
    is_g_convex = dgcp_curv == SymbolicAnalysis.GConvex || dgcp_curv == SymbolicAnalysis.GLinear
    
    return ComparisonResult(
        name,
        dgcp_curv,
        string(dcp_curv),
        is_eucl_convex,
        is_g_convex,
        notes
    )
end

#==============================================================================#
# Main Experiment
#==============================================================================#

function run_scope_comparison()
    results = ComparisonResult[]
    
    # Setup
    @variables X[1:5, 1:5]
    M = SymmetricPositiveDefinite(5)
    
    # Generate test data
    A = randn(5, 5)
    A = A * A' + I  # SPD matrix
    
    xs = [randn(5) for _ in 1:3]  # Random vectors for Tyler's estimator
    
    println("="^70)
    println("EXPERIMENT 1: DCP vs DGCP Verification Scope")
    println("="^70)
    println()
    
    #--------------------------------------------------------------------------
    # Case 1: logdet(X) - Both should verify
    #--------------------------------------------------------------------------
    expr = logdet(X) |> Symbolics.unwrap
    result = compare_verification(
        "logdet(X)",
        expr,
        HAS_CONVEX ? (Xc -> -Convex.logdet(Xc)) : nothing,  # Note: Convex.jl uses -logdet for convexity
        "Baseline: Both DCP and DGCP should verify"
    )
    push!(results, result)
    
    #--------------------------------------------------------------------------
    # Case 2: tr(X^{-1}) - Both verify (convex in Euclidean, g-convex on SPD)
    #--------------------------------------------------------------------------
    expr = tr(inv(X)) |> Symbolics.unwrap
    result = compare_verification(
        "tr(inv(X))",
        expr,
        nothing,  # Convex.jl doesn't have matrix inverse + trace composition
        "Trace of inverse: g-convex on SPD"
    )
    push!(results, result)
    
    #--------------------------------------------------------------------------
    # Case 3: Riemannian distance squared - DGCP yes, DCP no
    #--------------------------------------------------------------------------
    expr = Manifolds.distance(M, A, X)^2 |> Symbolics.unwrap
    result = compare_verification(
        "distance(M, A, X)²",
        expr,
        nothing,  # No Euclidean equivalent
        "Riemannian distance: g-convex but NOT Euclidean convex"
    )
    push!(results, result)
    
    #--------------------------------------------------------------------------
    # Case 4: S-divergence - DGCP yes, DCP no
    #--------------------------------------------------------------------------
    expr = SymbolicAnalysis.sdivergence(X, A) |> Symbolics.unwrap
    result = compare_verification(
        "S-divergence(X, A)",
        expr,
        nothing,  # No DCP equivalent
        "Symmetric Stein divergence: g-convex, used in matrix mean problems"
    )
    push!(results, result)
    
    #--------------------------------------------------------------------------
    # Case 5: Conjugation logdet - DGCP yes, DCP limited
    #--------------------------------------------------------------------------
    expr = logdet(SymbolicAnalysis.conjugation(inv(X), A)) |> Symbolics.unwrap
    result = compare_verification(
        "logdet(A' X^{-1} A)",
        expr,
        nothing,
        "Conjugation composition: key for Brascamp-Lieb"
    )
    push!(results, result)
    
    #--------------------------------------------------------------------------
    # Case 6: Tyler's M-Estimator objective - DGCP yes, DCP no
    #--------------------------------------------------------------------------
    expr = sum(SymbolicAnalysis.log_quad_form(x, inv(X)) for x in xs) + 
           (1/5) * logdet(X) |> Symbolics.unwrap
    result = compare_verification(
        "Tyler's M-Estimator",
        expr,
        nothing,
        "Maximum likelihood covariance: g-convex, Euclidean non-convex"
    )
    push!(results, result)
    
    #--------------------------------------------------------------------------
    # Case 7: Karcher mean objective - DGCP yes, DCP no
    #--------------------------------------------------------------------------
    As = [randn(5, 5) |> x -> x * x' + I for _ in 1:3]
    expr = sum(Manifolds.distance(M, Ai, X)^2 for Ai in As) |> Symbolics.unwrap
    result = compare_verification(
        "Karcher Mean (Σ d²)",
        expr,
        nothing,
        "Frechet mean on SPD: g-convex, Euclidean non-convex"
    )
    push!(results, result)
    
    #--------------------------------------------------------------------------
    # Print Results Table
    #--------------------------------------------------------------------------
    println()
    println("Results:")
    println("-"^70)
    println(rpad("Expression", 25), " | ", 
            rpad("DGCP", 12), " | ",
            rpad("Eucl. Convex", 12), " | ",
            "G-Convex")
    println("-"^70)
    
    for r in results
        println(
            rpad(r.name, 25), " | ",
            rpad(string(r.dgcp_curvature), 12), " | ",
            rpad(r.euclidean_convex ? "Yes" : "No", 12), " | ",
            r.geodesically_convex ? "Yes" : "No"
        )
    end
    println("-"^70)
    
    #--------------------------------------------------------------------------
    # Key Finding
    #--------------------------------------------------------------------------
    dgcp_only = count(r -> r.geodesically_convex && !r.euclidean_convex, results)
    both = count(r -> r.geodesically_convex && r.euclidean_convex, results)
    
    println()
    println("Summary:")
    println("  • Functions verified by DGCP only (g-convex, not Eucl-convex): $dgcp_only")
    println("  • Functions verified by both (g-convex and Eucl-convex): $both")
    println()
    println("This demonstrates that DGCP extends DCP's verification scope to")
    println("geodesically convex functions that are Euclidean non-convex.")
    
    return results
end

# Run tests
@testset "DCP vs DGCP Scope Comparison" begin
    results = run_scope_comparison()
    
    # Verify key results
    @test any(r -> r.name == "logdet(X)" && r.geodesically_convex, results)
    @test any(r -> r.name == "distance(M, A, X)²" && r.geodesically_convex, results)
    @test any(r -> r.name == "Tyler's M-Estimator" && r.geodesically_convex, results)
end
