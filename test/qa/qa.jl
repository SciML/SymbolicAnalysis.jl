using SymbolicAnalysis, Aqua, JET, Test

@testset "Aqua" begin
    Aqua.test_all(SymbolicAnalysis)
end

@testset "JET" begin
    JET.test_package(SymbolicAnalysis; target_defined_modules = true)
end
