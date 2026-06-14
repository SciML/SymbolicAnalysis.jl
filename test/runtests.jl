using SymbolicAnalysis:
    propagate_curvature,
    propagate_sign,
    propagate_gcurvature,
    getcurvature,
    getsign,
    getgcurvature

using SafeTestsets, Test

@safetestset "DCP" begin
    include("test.jl")
end

@safetestset "DGCP - SPD Manifold" begin
    include("dgp.jl")
end

@safetestset "DGCP - Lorentz Manifold" begin
    include("lorentz.jl")
end

@safetestset "Interface Compatibility" begin
    include("interface_tests.jl")
end

# AllocCheck tests - run separately to avoid precompilation overhead
# These tests verify that key operations have minimal allocations
if get(ENV, "SYMBOLICANALYSIS_TEST_ALLOC", "true") == "true"
    @safetestset "Allocation Tests" begin
        include("alloc_tests.jl")
    end
end
