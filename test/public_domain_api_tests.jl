using SymbolicAnalysis
using Symbolics
using LinearAlgebra
using Test

@testset "Public domain declaration API" begin
    @testset "declared public" begin
        if isdefined(Base, :ispublic)
            @test Base.ispublic(SymbolicAnalysis, :definite_domain)
            @test Base.ispublic(SymbolicAnalysis, :semidefinite_domain)
            @test Base.ispublic(SymbolicAnalysis, :symmetric_domain)
        end
    end

    @testset "definite_domain unlocks the positive-definite matrix rules" begin
        @variables X[1:3, 1:3]
        Xd = setmetadata(X, SymbolicAnalysis.VarDomain, SymbolicAnalysis.definite_domain())
        @test analyze(logdet(Xd) - tr(Xd)).curvature == SymbolicAnalysis.Concave
    end

    @testset "semidefinite_domain unlocks the same rules as definite_domain" begin
        @variables X[1:3, 1:3]
        Xd = setmetadata(
            X, SymbolicAnalysis.VarDomain, SymbolicAnalysis.semidefinite_domain()
        )
        @test analyze(logdet(Xd) - tr(Xd)).curvature == SymbolicAnalysis.Concave
    end

    @testset "symmetric_domain unlocks the eigmax/eigmin rules" begin
        @variables X[1:3, 1:3]
        Xs = setmetadata(X, SymbolicAnalysis.VarDomain, SymbolicAnalysis.symmetric_domain())
        @test analyze(eigmax(Xs)).curvature == SymbolicAnalysis.Convex
        @test analyze(eigmin(Xs)).curvature == SymbolicAnalysis.Concave
    end
end
