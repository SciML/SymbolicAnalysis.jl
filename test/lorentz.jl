using Manifolds, Symbolics, SymbolicAnalysis, LinearAlgebra
using Test
using Symbolics: unwrap
using SymbolicAnalysis: propagate_sign, propagate_curvature, propagate_gcurvature

@testset "Lorentz Manifold" begin
    # Create a Lorentz manifold of dimension 2 (3D ambient space)
    M = Lorentz(2)

    # Define symbolic variables
    @variables p[1:3]

    # Test lorentz_distance
    q = [0.0, 0.0, 1.0]  # A point on the Lorentz model
    ex = Manifolds.distance(M, q, p) |> unwrap
    ex = propagate_sign(ex)
    ex = propagate_gcurvature(ex, M)
    @test SymbolicAnalysis.getgcurvature(ex) == SymbolicAnalysis.GConvex

    # Test analyze function
    analyze_res = analyze(ex, M)
    @test analyze_res.gcurvature == SymbolicAnalysis.GConvex

    # Test lorentz_homogeneous_quadratic
    A = [2.0 0.0 0.0; 0.0 2.0 0.0; 0.0 0.0 1.0]  # Positive definite matrix
    ex = SymbolicAnalysis.lorentz_homogeneous_quadratic(A, p) |> unwrap
    ex = propagate_sign(ex)
    ex = propagate_gcurvature(ex, M)
    @test SymbolicAnalysis.getgcurvature(ex) == SymbolicAnalysis.GConvex

    # Test lorentz_homogeneous_diagonal
    a = [2.0, 2.0, 1.0]  # min(a[1:2]) + a[3] = 2 + 1 = 3 ≥ 0
    ex = SymbolicAnalysis.lorentz_homogeneous_diagonal(a, p) |> unwrap
    ex = propagate_sign(ex)
    ex = propagate_gcurvature(ex, M)
    @test SymbolicAnalysis.getgcurvature(ex) == SymbolicAnalysis.GConvex

    # Test lorentz_least_squares
    # Create a valid X and y that satisfy the geodesic convexity conditions:
    # ∑^d_i=1(X'y)^2_i ≤ (X'y)^2_{d+1} and (X'y)_{d+1} ≤ 0
    X = [1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0]
    y = [0.0, 0.0, -1.0]  # X'y = [0.0, 0.0, -1.0], which satisfies both conditions
    ex = SymbolicAnalysis.lorentz_least_squares(X, y, p) |> unwrap
    ex = propagate_sign(ex)
    ex = propagate_gcurvature(ex, M)
    @test SymbolicAnalysis.getgcurvature(ex) == SymbolicAnalysis.GConvex

    @testset "Least Squares Problem" begin
        # Define variables for symbolic testing
        @variables p[1:3]
        M = Manifolds.Lorentz(2)  # 2D Lorentz model (3D ambient)

        # Create a valid test case with data that satisfies geodesic convexity conditions
        # We need X and y such that:
        # 1. ∑_{i=1}^d (X'y)²_i ≤ (X'y)²_{d+1}
        # 2. (X'y)_{d+1} ≤ 0
        X = [1.0 0.0 2.0; 0.0 1.0 3.0; 2.0 2.0 10.0]
        y = [1.0, 2.0, -5.0]

        # Verify conditions explicitly for the test
        Xty = X' * y
        @test sum(Xty[1:2] .^ 2) <= Xty[3]^2  # Condition 1
        @test Xty[3] <= 0                   # Condition 2

        # Compose the least squares problem from atoms
        A = X' * X  # Positive semidefinite, automatically ∂L-copositive
        b = -2 * X' * y  # Must be in Lorentz cone
        c = y' * y

        # Create expression using lorentz_nonhomogeneous_quadratic
        expr = SymbolicAnalysis.lorentz_nonhomogeneous_quadratic(A, b, c, p)

        # Verify geodesic convexity through DGCP framework
        expr = propagate_sign(expr)
        expr = propagate_gcurvature(expr, M)
        @test SymbolicAnalysis.getgcurvature(expr) == SymbolicAnalysis.GConvex

        # Verify that the composition matches the direct expansion
        # direct_expr = c - 2 * p' * X' * y + p' * X' * X * p
        # @test isequal(simplify(expr), simplify(direct_expr))
    end
    # Test composition of functions
    ex = 2.0 * Manifolds.distance(M, q, p) +
        Manifolds.distance(M, [1.0, 0.0, sqrt(2.0)], p) |> unwrap
    ex = propagate_sign(ex)
    ex = propagate_gcurvature(ex, M)
    @test SymbolicAnalysis.getgcurvature(ex) == SymbolicAnalysis.GConvex

    # Test lorentz_transform (should preserve geodesic convexity)
    # Create a Lorentz boost in the x-direction
    cosh_phi = 1.2
    sinh_phi = sqrt(cosh_phi^2 - 1)
    O = [1.0 0.0 0.0; 0.0 cosh_phi -sinh_phi; 0.0 -sinh_phi cosh_phi]

    # Create a compound expression with the transform
    q_transformed = O * q
    ex = Manifolds.distance(M, q_transformed, p) |> unwrap
    ex = propagate_sign(ex)
    ex = propagate_gcurvature(ex, M)
    @test SymbolicAnalysis.getgcurvature(ex) == SymbolicAnalysis.GConvex
end

# Sign audits: several Lorentz atoms were registered Positive while taking
# negative values on their stated geodesic-convexity domains.
@testset "Lorentz atom signs are not falsely Positive" begin
    @variables p[1:3]
    M = Manifolds.Lorentz(3)
    # -log(t-1) < 0 for t > 2 (t = p[end])
    @test analyze(SymbolicAnalysis.lorentz_log_barrier(p), M).sign !=
        SymbolicAnalysis.Positive
    A = Diagonal([1.0, 1.0, -0.5]) |> Matrix
    @test analyze(SymbolicAnalysis.lorentz_homogeneous_quadratic(A, p), M).sign !=
        SymbolicAnalysis.Positive
    a = [1.0, 1.0, -0.5]
    @test analyze(SymbolicAnalysis.lorentz_homogeneous_diagonal(a, p), M).sign !=
        SymbolicAnalysis.Positive
end

@testset "lorentz_log_barrier is not certified geodesically convex" begin
    @variables p[1:3]
    M = Manifolds.Lorentz(2)
    @test analyze(SymbolicAnalysis.lorentz_log_barrier(p), M).gcurvature ==
        SymbolicAnalysis.GUnknownCurvature
    ex = 2.0 * Manifolds.distance(M, [0.0, 0.0, 1.0], p) +
        SymbolicAnalysis.lorentz_log_barrier(p)
    @test analyze(ex, M).gcurvature == SymbolicAnalysis.GUnknownCurvature

    # Witness: along the geodesic through q in the direction v (v[end] = 0) the
    # barrier is strictly concave at q, with second derivative -q[end] / (q[end] - 1).
    H = Manifolds.Hyperbolic(2)
    q = [0.3, 0.4, sqrt(1.25)]
    v = [-0.8, 0.6, 0.0]
    a, b = exp(H, q, -v), exp(H, q, v)
    @test shortest_geodesic(H, a, b, 0.5) ≈ q
    f = SymbolicAnalysis.lorentz_log_barrier
    @test f(q) > (f(a) + f(b)) / 2
end
