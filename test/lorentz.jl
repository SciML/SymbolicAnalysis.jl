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
    @test analyze(Manifolds.distance(M, q, p), M).sign == SymbolicAnalysis.Positive
    @test analyze(Manifolds.distance(M, q, p)^2, M).gcurvature == SymbolicAnalysis.GConvex

    # Test lorentz_log_barrier
    ex = SymbolicAnalysis.lorentz_log_barrier(p) |> unwrap
    ex = propagate_sign(ex)
    ex = propagate_gcurvature(ex, M)
    @test SymbolicAnalysis.getgcurvature(ex) == SymbolicAnalysis.GConvex

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
    ex = 2.0 * Manifolds.distance(M, q, p) + SymbolicAnalysis.lorentz_log_barrier(p) |>
        unwrap
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

# SPD sign rules must not fire under Lorentz for a vector argument: p = e₃ lies
# on the hyperboloid and `quad_form(p, diag(1,1,-5)) = -5`.
@testset "SPD atoms do not certify Lorentz vectors" begin
    @variables p[1:3]
    Cind = Diagonal([1, 1, -5])
    for M in (Lorentz(2), Lorentz(3))
        @test analyze(SymbolicAnalysis.quad_form(p, Cind), M).sign !=
            SymbolicAnalysis.Positive
        @test analyze(tr(p), M).sign != SymbolicAnalysis.Positive
    end
    e3 = [0.0, 0.0, 1.0]
    @test SymbolicAnalysis.quad_form(e3, Matrix(Cind)) < 0
end

@testset "Lorentz parameters must match the point's dimension" begin
    @variables p[1:3] A3[1:3, 1:3] A4[1:4, 1:4]
    M = Lorentz(2)
    @test analyze(SymbolicAnalysis.lorentz_homogeneous_quadratic(A4, p), M).gcurvature ==
        SymbolicAnalysis.GUnknownCurvature
    @test analyze(SymbolicAnalysis.lorentz_homogeneous_quadratic(A3, p), M).gcurvature ==
        SymbolicAnalysis.GUnknownCurvature
end

@testset "Lorentz quadratic atoms need data satisfying their conditions" begin
    @variables p[1:3] A3[1:3, 1:3] a3[1:3] b3[1:3]
    M = Lorentz(2)
    gc(ex) = analyze(ex, M).gcurvature
    unknown = SymbolicAnalysis.GUnknownCurvature
    mink(u, v) = u[1] * v[1] + u[2] * v[2] - u[3] * v[3]
    function geodesic_mid(P, Q)
        d = acosh(-mink(P, Q))
        return (P .+ Q) ./ (2 * cosh(d / 2))
    end
    # Each `f` violates geodesic midpoint convexity between `P` and `Q`.
    s1, c1 = sinh(big(1)), cosh(big(1))
    pt(x1, x2) = BigFloat[x1, x2, sqrt(1 + big(x1)^2 + big(x2)^2)]
    p0, e2 = BigFloat[2.5, 0, sqrt(big(7.25))], BigFloat[0, 1, 0]
    cases = [
        (
            SymbolicAnalysis.lorentz_homogeneous_diagonal([1.0, 1.0, -3.0], p),
            x -> x[1]^2 + x[2]^2 - 3x[3]^2, pt(-1, 0), pt(1, 0),
        ),
        (
            SymbolicAnalysis.lorentz_homogeneous_quadratic([1.0 0 0; 0 1 0; 0 0 -3], p),
            x -> x[1]^2 + x[2]^2 - 3x[3]^2, pt(-1, 0), pt(1, 0),
        ),
        (
            SymbolicAnalysis.lorentz_nonhomogeneous_quadratic(
                [1.0 0 0; 0 1 0; 0 0 -1], [-1.0, 0, 0], 0.0, p
            ),
            x -> x[1]^2 + x[2]^2 - x[3]^2 - x[1], BigFloat[0, 0, 1],
            BigFloat[sinh(big(2)), 0, cosh(big(2))],
        ),
        (
            SymbolicAnalysis.lorentz_least_squares([1.0 0 0], [5.0], p),
            x -> (5 - x[1])^2, c1 .* p0 .- s1 .* e2, c1 .* p0 .+ s1 .* e2,
        ),
        (
            SymbolicAnalysis.lorentz_least_squares(Matrix(1.0I(3)), [0, 0, 5.0], p),
            x -> x[1]^2 + x[2]^2 + (5 - x[3])^2, BigFloat[s1, 0, c1],
            BigFloat[-s1, 0, c1],
        ),
    ]
    for (ex, f, P, Q) in cases
        @test f(geodesic_mid(P, Q)) > (f(P) + f(Q)) / 2
        @test gc(ex) == unknown
    end
    @test analyze(SymbolicAnalysis.lorentz_least_squares([1.0 0 0], [5.0], p), M).sign ==
        SymbolicAnalysis.Positive

    @test gc(SymbolicAnalysis.lorentz_homogeneous_quadratic(A3, p)) == unknown
    @test gc(SymbolicAnalysis.lorentz_homogeneous_diagonal(a3, p)) == unknown
    @test gc(SymbolicAnalysis.lorentz_nonhomogeneous_quadratic(Matrix(1.0I(3)), b3, 0.0, p)) ==
        unknown

    # Boundary data: `σ = -λ_min(Ā)` is admitted, one ulp below is not.
    @test gc(SymbolicAnalysis.lorentz_homogeneous_diagonal([1.0, 1.0, -1.0], p)) ==
        SymbolicAnalysis.GConvex
    @test gc(SymbolicAnalysis.lorentz_homogeneous_diagonal([1.0, 1.0, prevfloat(-1.0)], p)) ==
        unknown
    @test gc(SymbolicAnalysis.lorentz_homogeneous_quadratic([1.0 0 0; 0 1 0; 0 0 -1], p)) ==
        SymbolicAnalysis.GConvex
    @test gc(
        SymbolicAnalysis.lorentz_homogeneous_quadratic(
            [1.0 0 0; 0 1 0; 0 0 prevfloat(-1.0)], p
        )
    ) == unknown
    # `a ≠ 0` branch: σ + λ_min(Ā) > 2‖a‖ with `a = (0, 1)`, `Ā = I`.
    @test gc(SymbolicAnalysis.lorentz_homogeneous_quadratic([1.0 0 0; 0 1 1; 0 1 1.5], p)) ==
        SymbolicAnalysis.GConvex
    @test gc(SymbolicAnalysis.lorentz_homogeneous_quadratic([1.0 0 0; 0 1 1; 0 1 1], p)) ==
        unknown
end
