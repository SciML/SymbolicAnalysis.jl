# Interface compatibility tests for SymbolicAnalysis.jl
# Tests that helper functions work with BigFloat and AbstractVector types

using SymbolicAnalysis
using LinearAlgebra
using Test

@testset "BigFloat support" begin
    @testset "dotsort" begin
        x = BigFloat[3.0, 1.0, 2.0]
        y = BigFloat[6.0, 4.0, 5.0]
        result = SymbolicAnalysis.dotsort(x, y)
        @test result == BigFloat(32.0)
        @test eltype([result]) == BigFloat
    end

    @testset "invprod" begin
        x = BigFloat[2.0, 3.0, 4.0]
        result = SymbolicAnalysis.invprod(x)
        @test result isa BigFloat
        @test isapprox(result, BigFloat(1) / BigFloat(24))
    end

    @testset "matrix_frac" begin
        x = BigFloat[1.0, 2.0]
        P = BigFloat[2.0 0.5; 0.5 2.0]
        result = SymbolicAnalysis.matrix_frac(x, P)
        @test result isa BigFloat
    end

    @testset "quad_form" begin
        x = BigFloat[1.0, 2.0]
        P = BigFloat[2.0 0.5; 0.5 2.0]
        result = SymbolicAnalysis.quad_form(x, P)
        @test result isa BigFloat
        @test result == BigFloat(12.0)
    end

    @testset "quad_over_lin" begin
        x = BigFloat[1.0, 2.0, 3.0]
        y = BigFloat(2.0)
        result = SymbolicAnalysis.quad_over_lin(x, y)
        @test result isa BigFloat
        @test result == BigFloat(7.0)

        # Scalar version
        result_scalar = SymbolicAnalysis.quad_over_lin(BigFloat(3.0), BigFloat(2.0))
        @test result_scalar isa BigFloat
        @test result_scalar == BigFloat(4.5)
    end

    @testset "tv" begin
        x = BigFloat[1.0, 2.0, 5.0, 3.0]
        result = SymbolicAnalysis.tv(x)
        @test result isa BigFloat
        @test result == BigFloat(6.0)
    end

    @testset "sum_largest" begin
        x = BigFloat[1.0 2.0; 3.0 4.0]
        result = SymbolicAnalysis.sum_largest(x, 2)
        @test result isa BigFloat
        # sum_largest sums the k largest elements: sorted = [1,2,3,4], (end-k):end = 3:4+1 = 3 elements
        # Looking at the code: sort(vec(x))[(end - k):end] = [2, 3, 4] when k=2, so sum = 9
        @test result == BigFloat(9.0)
    end

    @testset "sum_smallest" begin
        x = BigFloat[1.0 2.0; 3.0 4.0]
        result = SymbolicAnalysis.sum_smallest(x, 2)
        @test result isa BigFloat
        @test result == BigFloat(3.0)
    end

    @testset "trinv" begin
        x = BigFloat[4.0 1.0; 1.0 3.0]
        result = SymbolicAnalysis.trinv(x)
        @test result isa BigFloat
    end

    @testset "huber" begin
        result = SymbolicAnalysis.huber(BigFloat(1.5), BigFloat(1.0))
        @test result isa BigFloat
        @test result == BigFloat(2.0)
    end

    @testset "rel_entr" begin
        result = SymbolicAnalysis.rel_entr(BigFloat(2.0), BigFloat(3.0))
        @test result isa BigFloat
    end

    @testset "perspective" begin
        result = SymbolicAnalysis.perspective(exp, BigFloat(2.0), BigFloat(3.0))
        @test result isa BigFloat
    end

    @testset "GDCP functions" begin
        @testset "conjugation" begin
            X = BigFloat[4.0 1.0; 1.0 3.0]
            B = BigFloat[1.0 0.0; 1.0 1.0]
            result = SymbolicAnalysis.conjugation(X, B)
            @test eltype(result) == BigFloat
        end

        @testset "scalar_mat" begin
            X = BigFloat[4.0 1.0; 1.0 3.0]
            result = SymbolicAnalysis.scalar_mat(X, 2)
            @test eltype(result) == BigFloat
        end

        @testset "sdivergence" begin
            X = BigFloat[4.0 1.0; 1.0 3.0]
            Y = BigFloat[3.0 0.5; 0.5 2.0]
            result = SymbolicAnalysis.sdivergence(X, Y)
            @test result isa BigFloat
        end

        @testset "log_quad_form" begin
            y = Vector(BigFloat[1.0, 2.0])
            X = Matrix(BigFloat[4.0 1.0; 1.0 3.0])
            result = SymbolicAnalysis.log_quad_form(y, X)
            @test result isa BigFloat
        end

        @testset "lorentz_homogeneous_diagonal" begin
            a = BigFloat[1.0, 2.0, -1.0]
            p = BigFloat[1.0, 0.5, sqrt(BigFloat(2))]
            result = SymbolicAnalysis.lorentz_homogeneous_diagonal(a, p)
            @test result isa BigFloat
        end

        @testset "lorentz_transform" begin
            O = BigFloat[1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0]
            p = BigFloat[0.0, 0.0, 1.0]
            result = SymbolicAnalysis.lorentz_transform(O, p)
            @test eltype(result) == BigFloat
        end
    end
end

@testset "AbstractVector support" begin
    @testset "dotsort with views" begin
        x = [3.0, 1.0, 2.0, 100.0]
        y = [6.0, 4.0, 5.0, 100.0]
        v_x = view(x, 1:3)
        v_y = view(y, 1:3)
        result = SymbolicAnalysis.dotsort(v_x, v_y)
        @test result == 32.0
    end

    @testset "quad_over_lin with views" begin
        x = [1.0, 2.0, 3.0, 100.0]
        v = view(x, 1:3)
        y = 2.0
        result = SymbolicAnalysis.quad_over_lin(v, y)
        @test result == 7.0
    end

    @testset "tv with views" begin
        x = [1.0, 2.0, 5.0, 3.0, 100.0]
        v = view(x, 1:4)
        result = SymbolicAnalysis.tv(v)
        @test result == 6.0
    end
end
