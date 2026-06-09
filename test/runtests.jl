using SafeTestsets, Test

const GROUP = get(ENV, "GROUP", "All")

if GROUP == "All" || GROUP == "Core"
    using SymbolicAnalysis:
        propagate_curvature,
        propagate_sign,
        propagate_gcurvature,
        getcurvature,
        getsign,
        getgcurvature

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
end

if GROUP == "QA"
    @testset "Quality Assurance" begin
        include("qa/qa.jl")
    end
end
