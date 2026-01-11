using SymbolicAnalysis:
    propagate_curvature,
    propagate_sign,
    propagate_gcurvature,
    getcurvature,
    getsign,
    getgcurvature

using SafeTestsets, Test

@testset "DCP" begin
    include("test.jl")
end

@testset "DGCP - SPD Manifold" begin
    include("dgp.jl")
end

@testset "DGCP - Lorentz Manifold" begin
    include("lorentz.jl")
end

@testset "Interface Compatibility" begin
    include("interface_tests.jl")
end

# AllocCheck tests - run separately to avoid precompilation overhead
# These tests verify that key operations have minimal allocations
if get(ENV, "SYMBOLICANALYSIS_TEST_ALLOC", "true") == "true"
    @testset "Allocation Tests" begin
        include("alloc_tests.jl")
    end
end
