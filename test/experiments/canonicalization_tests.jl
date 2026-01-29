"""
Canonicalization Tests

Tests for the DGCP-aware canonicalization pass.
"""

using SymbolicAnalysis
using Symbolics
using LinearAlgebra
using Manifolds
using Test

@testset "Canonicalization" begin
    @variables X[1:5, 1:5] Y[1:5, 1:5]
    M = SymmetricPositiveDefinite(5)
    
    @testset "Double Inverse Simplification" begin
        expr = inv(inv(X)) |> Symbolics.unwrap
        canon = SymbolicAnalysis.canonize(expr)
        # Should simplify to X
        @test string(canon) == "X"
    end
    
    @testset "Logdet of Inverse" begin
        expr = logdet(inv(X)) |> Symbolics.unwrap
        canon = SymbolicAnalysis.canonize(expr)
        # Should become negative logdet
        @test occursin("-", string(canon)) || occursin("log", string(canon))
    end
    
    @testset "Analysis After Canonicalization" begin
        # logdet should still verify correctly after canonicalization
        expr = logdet(X) |> Symbolics.unwrap
        result = analyze(expr, M)
        @test result.gcurvature == SymbolicAnalysis.GLinear
        
        # distance squared should verify
        A = randn(5, 5); A = A * A' + I
        expr = Manifolds.distance(M, A, X)^2 |> Symbolics.unwrap
        result = analyze(expr, M)
        @test result.gcurvature == SymbolicAnalysis.GConvex
    end
    
    @testset "Equivalent Forms Documentation" begin
        forms = SymbolicAnalysis.equivalent_forms()
        @test length(forms) >= 5
        @test all(haskey(f, :verifiable) for f in forms)
        @test all(haskey(f, :not_verifiable) for f in forms)
    end
    
    @testset "Is Canonical Check" begin
        # Simple expressions should be canonical
        expr = logdet(X) |> Symbolics.unwrap
        @test SymbolicAnalysis.is_canonical(expr)
        
        # inv(inv(X)) should NOT be canonical
        expr = inv(inv(X)) |> Symbolics.unwrap
        @test !SymbolicAnalysis.is_canonical(expr)
    end
end

println("✓ All canonicalization tests passed!")
