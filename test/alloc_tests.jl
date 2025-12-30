using SymbolicAnalysis
using SymbolicAnalysis: getsign, hassign, getcurvature, hascurvature, getgcurvature,
                        hasgcurvature,
                        add_sign, mul_sign, add_curvature, mul_curvature,
                        add_gcurvature, mul_gcurvature,
                        Sign, Positive, Negative, AnySign,
                        Curvature, Convex, Concave, Affine, UnknownCurvature,
                        GCurvature, GConvex, GConcave, GLinear, GUnknownCurvature
using Symbolics
using Test
using AllocCheck

# Test that sign/curvature accessor functions don't allocate
@testset "Allocation-free accessors" begin
    @variables x y z

    # Create expressions with metadata set
    ex_pos = SymbolicAnalysis.setsign(x, Positive)
    ex_neg = SymbolicAnalysis.setsign(y, Negative)

    # Test getsign doesn't allocate on numeric types
    @test (@allocations getsign(1.0)) == 0
    @test (@allocations getsign(-5)) == 0

    # Test hassign doesn't allocate on numeric types
    @test (@allocations hassign(1.0)) == 0
    @test (@allocations hassign(-5)) == 0

    # Test getcurvature doesn't allocate on numeric types
    @test (@allocations getcurvature(1.0)) == 0
    @test (@allocations getcurvature(-5)) == 0

    # Test hascurvature doesn't allocate on numeric types
    @test (@allocations hascurvature(1.0)) == 0
    @test (@allocations hascurvature(-5)) == 0

    # Test getgcurvature doesn't allocate on numeric types
    @test (@allocations getgcurvature(1.0)) == 0
    @test (@allocations getgcurvature(-5)) == 0

    # Test hasgcurvature doesn't allocate on numeric types
    @test (@allocations hasgcurvature(1.0)) == 0
    @test (@allocations hasgcurvature(-5)) == 0
end

@testset "Low allocation sign operations" begin
    # Test mul_sign with numeric arguments
    numeric_args = [1.0, -2.0, 3.0]
    # Warm up
    mul_sign(numeric_args)
    allocs = @allocations mul_sign(numeric_args)
    @test allocs <= 3  # Allow minimal allocations for iteration
end

@testset "Low allocation curvature operations" begin
    # Test with numeric arguments (should be Affine)
    numeric_args = [1.0, 2.0, 3.0]
    # Warm up
    add_curvature(numeric_args)
    mul_curvature(numeric_args)

    allocs_add = @allocations add_curvature(numeric_args)
    allocs_mul = @allocations mul_curvature(numeric_args)

    @test allocs_add <= 3
    @test allocs_mul <= 3
end

@testset "Low allocation gcurvature operations" begin
    # Test with numeric arguments (should be GLinear)
    numeric_args = [1.0, 2.0, 3.0]
    # Warm up
    add_gcurvature(numeric_args)
    mul_gcurvature(numeric_args)

    allocs_add = @allocations add_gcurvature(numeric_args)
    allocs_mul = @allocations mul_gcurvature(numeric_args)

    @test allocs_add <= 3
    @test allocs_mul <= 3
end
