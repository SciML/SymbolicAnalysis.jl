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
