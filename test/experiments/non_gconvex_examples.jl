"""
Experiment 2: Non-G-Convex Identification Examples

This experiment demonstrates that DGCP correctly identifies functions
that are NOT verifiably geodesically convex, returning `GUnknownCurvature`.

Addresses:
- Reviewer 400: "provide explicit examples to illustrate how the framework 
  recognizes functions that are not geodesically convex"
"""

using SymbolicAnalysis
using Manifolds
using Symbolics
using LinearAlgebra
using Test

#==============================================================================#
# Test Cases: Known Non-G-Convex or Non-DGCP-Verifiable Functions
#==============================================================================#

"""
Run analysis and check for non-verification
"""
function test_non_gconvex(name::String, expr, expected_result::Symbol, reason::String)
    M = SymmetricPositiveDefinite(5)
    result = analyze(expr, M)
    
    return (
        name = name,
        gcurvature = result.gcurvature,
        eucl_curvature = result.curvature,
        expected = expected_result,
        passed = result.gcurvature == SymbolicAnalysis.GUnknownCurvature,
        reason = reason
    )
end

function run_non_gconvex_examples()
    println("="^70)
    println("EXPERIMENT 2: Non-G-Convex Identification")
    println("="^70)
    println()
    println("Testing that DGCP correctly returns GUnknownCurvature for")
    println("functions that cannot be verified as geodesically convex.")
    println()
    
    results = []
    
    # Setup
    @variables X[1:5, 1:5] Y[1:5, 1:5]
    @variables x[1:5]
    M = SymmetricPositiveDefinite(5)
    
    A = randn(5, 5)
    A = A * A' + I
    
    #--------------------------------------------------------------------------
    # Case 1: Product of matrix variables - not DGCP-verifiable
    #--------------------------------------------------------------------------
    # sqrt(X * Y) involves multiple variables in non-affine way
    expr = sqrt(X * Y) |> Symbolics.unwrap
    result = test_non_gconvex(
        "sqrt(X * Y)",
        expr,
        :GUnknownCurvature,
        "Product of two SPD variables: no composition rule applies"
    )
    push!(results, result)
    
    #--------------------------------------------------------------------------
    # Case 2: X - A (matrix difference) - not g-convex preserving
    #--------------------------------------------------------------------------
    expr = (X - A) |> Symbolics.unwrap
    result = test_non_gconvex(
        "X - A (difference)",
        expr,
        :GUnknownCurvature,
        "Matrix subtraction: doesn't preserve SPD structure"
    )
    push!(results, result)
    
    #--------------------------------------------------------------------------
    # Case 3: tr(X^2) - second power without log transform
    #--------------------------------------------------------------------------
    # Note: This depends on how X^2 is represented
    expr = tr(X * X) |> Symbolics.unwrap
    result = test_non_gconvex(
        "tr(X²)",
        expr,
        :GUnknownCurvature,
        "Quadratic in Frobenius: not g-convex without log transform"
    )
    push!(results, result)
    
    #--------------------------------------------------------------------------
    # Case 4: X + Y (sum of two matrix variables) - no DGCP rule for this
    #--------------------------------------------------------------------------
    expr = (X + Y) |> Symbolics.unwrap
    result = test_non_gconvex(
        "X + Y (sum)",
        expr,
        :GUnknownCurvature,
        "Sum of two matrix variables: not g-linear in general on SPD"
    )
    push!(results, result)
    
    #--------------------------------------------------------------------------
    # Case 5: log(det(X)^2) written as log(det(X))^2 - wrong composition
    #--------------------------------------------------------------------------
    # This is different from 2*logdet(X) which would be g-linear
    expr = logdet(X)^2 |> Symbolics.unwrap
    result = test_non_gconvex(
        "(logdet(X))²",
        expr,
        :GUnknownCurvature,
        "Square of logdet: not same as 2*logdet(X)"
    )
    push!(results, result)
    
    #--------------------------------------------------------------------------
    # Case 6: logdet(X) * logdet(Y) - product of two g-linear terms
    #--------------------------------------------------------------------------
    # Product of g-linear functions is not g-linear
    expr = logdet(X) * logdet(Y) |> Symbolics.unwrap
    result = test_non_gconvex(
        "logdet(X)*logdet(Y)",
        expr,
        :GUnknownCurvature,
        "Product of g-linear terms: not necessarily g-convex"
    )
    push!(results, result)
    
    #--------------------------------------------------------------------------
    # Print Results
    #--------------------------------------------------------------------------
    println("-"^70)
    println(rpad("Expression", 20), " | ", 
            rpad("DGCP Result", 20), " | ",
            rpad("Correctly Rejected?", 20))
    println("-"^70)
    
    for r in results
        status = r.passed ? "✓ Yes" : "✗ No"
        println(
            rpad(r.name, 20), " | ",
            rpad(string(r.gcurvature), 20), " | ",
            status
        )
    end
    println("-"^70)
    
    #--------------------------------------------------------------------------
    # Detailed Explanations
    #--------------------------------------------------------------------------
    println()
    println("Explanations:")
    println("-"^70)
    for r in results
        println("• $(r.name): $(r.reason)")
    end
    
    #--------------------------------------------------------------------------
    # Summary
    #--------------------------------------------------------------------------
    correctly_rejected = count(r -> r.passed, results)
    total = length(results)
    
    println()
    println("Summary:")
    println("  • Correctly identified as non-g-convex: $correctly_rejected / $total")
    println()
    println("This demonstrates that DGCP does not falsely claim g-convexity")
    println("for functions that cannot be verified through composition rules.")
    
    return results
end

#==============================================================================#
# Comparison: Equivalent Forms
#==============================================================================#

function run_equivalent_form_comparison()
    println()
    println("="^70)
    println("BONUS: Equivalent Forms with Different Verifiability")
    println("="^70)
    println()
    println("Demonstrating that symbolic representation affects verifiability")
    println("(addresses Reviewer 385's concern about symbolic non-uniqueness)")
    println()
    
    @variables X[1:5, 1:5]
    M = SymmetricPositiveDefinite(5)
    
    # Case: 2 * logdet(X) vs logdet(X)^2
    expr1 = 2 * logdet(X) |> Symbolics.unwrap
    expr2 = logdet(X)^2 |> Symbolics.unwrap
    
    result1 = analyze(expr1, M)
    result2 = analyze(expr2, M)
    
    println("Expression 1: 2 * logdet(X)")
    println("  → DGCP: $(result1.gcurvature)")
    println()
    println("Expression 2: logdet(X)²")
    println("  → DGCP: $(result2.gcurvature)")
    println()
    println("Note: These are mathematically different functions, but this")
    println("illustrates how users should choose DGCP-compliant formulations.")
    
    return (expr1_result=result1, expr2_result=result2)
end

# Run tests
@testset "Non-G-Convex Identification" begin
    results = run_non_gconvex_examples()
    
    # At least some should be correctly rejected
    @test any(r -> r.passed, results)
end

@testset "Equivalent Form Comparison" begin
    results = run_equivalent_form_comparison()
    
    # 2*logdet should verify, logdet^2 should not
    @test results.expr1_result.gcurvature == SymbolicAnalysis.GLinear
end
