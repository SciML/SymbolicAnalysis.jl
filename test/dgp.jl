using Manifolds, Symbolics, SymbolicAnalysis, LinearAlgebra
using LinearAlgebra, PDMats
using Symbolics: unwrap
using Test, Zygote, ForwardDiff
using SymbolicAnalysis: propagate_sign, propagate_curvature, propagate_gcurvature

@variables X[1:5, 1:5]

M = Manifolds.SymmetricPositiveDefinite(5)

A = rand(5, 5)
A = A * A'

ex = SymbolicAnalysis.logdet(SymbolicAnalysis.conjugation(inv(X), A)) |> unwrap
ex = propagate_sign(ex)
ex = propagate_curvature(ex)
ex = propagate_gcurvature(ex, M)
SymbolicAnalysis.getcurvature(ex)
@test SymbolicAnalysis.getgcurvature(ex) == SymbolicAnalysis.GConvex

ex = SymbolicAnalysis.logdet(tr(inv(X))) |> unwrap
ex = propagate_sign(ex)
ex = propagate_curvature(ex)
ex = propagate_gcurvature(ex, M)
@test SymbolicAnalysis.getgcurvature(ex) == SymbolicAnalysis.GConvex
SymbolicAnalysis.getcurvature(ex)

@variables Sigma[1:5, 1:5]
xs = [rand(5) for i in 1:2]
ex = sum(SymbolicAnalysis.log_quad_form(x, inv(Sigma)) for x in xs) +
    1 / 5 * logdet(Sigma) |> Symbolics.unwrap
analyze_res = SymbolicAnalysis.analyze(ex, M)
@test analyze_res.gcurvature == SymbolicAnalysis.GConvex

##Brascamplieb Problem
M = SymmetricPositiveDefinite(5)
objective_expr = logdet(SymbolicAnalysis.conjugation(X, A)) - logdet(X) |> unwrap
objective_expr = SymbolicAnalysis.propagate_sign(objective_expr)
analyze_res = analyze(objective_expr, M)
@test analyze_res.gcurvature == SymbolicAnalysis.GConvex

objective_expr = SymbolicAnalysis.propagate_gcurvature(objective_expr, M)
@test SymbolicAnalysis.getgcurvature(objective_expr) == SymbolicAnalysis.GConvex

ex = SymbolicAnalysis.tr(SymbolicAnalysis.conjugation(X, A)) |> unwrap
ex = propagate_sign(ex)
ex = propagate_curvature(ex)
ex = propagate_gcurvature(ex, M)

@test analyze(ex, M).gcurvature == SymbolicAnalysis.GConvex

# using Convex

# X = Convex.Variable(5, 5)
# Y = Convex.Variable(5, 5)
# ex = sqrt(X*Y)
# vexity(ex)

## Karcher Mean
As = [rand(5, 5) for i in 1:5]
As = [As[i] * As[i]' for i in 1:5]

ex = SymbolicAnalysis.sdivergence(X, As[1]) |> unwrap
ex = SymbolicAnalysis.propagate_sign(ex)
ex = SymbolicAnalysis.propagate_gcurvature(ex, M)

@test SymbolicAnalysis.getgcurvature(ex) == SymbolicAnalysis.GConvex

ex = sum(SymbolicAnalysis.sdivergence(X, As[i]) for i in 1:5) |> Symbolics.unwrap
ex = SymbolicAnalysis.propagate_sign(ex)
ex = SymbolicAnalysis.propagate_gcurvature(ex, M)

@test SymbolicAnalysis.getgcurvature(ex) == SymbolicAnalysis.GConvex

ex = Manifolds.distance(M, As[1], X)^2 |> Symbolics.unwrap
ex = SymbolicAnalysis.propagate_sign(ex)
ex = SymbolicAnalysis.propagate_gcurvature(ex, M)

@test SymbolicAnalysis.getgcurvature(ex) == SymbolicAnalysis.GConvex

M = SymmetricPositiveDefinite(5)
objective_expr = sum(Manifolds.distance(M, As[i], X)^2 for i in 1:5) |> Symbolics.unwrap
analyze_res = analyze(objective_expr, M)
@test analyze_res.gcurvature == SymbolicAnalysis.GConvex

@variables Y[1:5, 1:5]
ex = sqrt(X * Y)
analyze_res = analyze(ex, M)
@test analyze_res.gcurvature == SymbolicAnalysis.GUnknownCurvature

# ex = exp(X*Y) |> unwrap
# ex = SymbolicAnalysis.propagate_sign(ex)
# @test_throws SymbolicUtils.RuleRewriteError SymbolicAnalysis.propagate_gcurvature(ex)

# using Manopt, Manifolds, Random, LinearAlgebra, ManifoldDiff
# using ManifoldDiff: grad_distance, prox_distance
# Random.seed!(42);

# m = 100
# σ = 0.005
# q = Matrix{Float64}(I, 5, 5) .+ 2.0
# data2 = [exp(M, q, σ * rand(M; vector_at=q)) for i in 1:m];

# f(M, x) = sum(distance(M, x, data2[i])^2 for i in 1:m)
# f(x) = sum(distance(M, x, data2[i])^2 for i in 1:m)

# using FiniteDifferences

# r_backend = ManifoldDiff.RiemannianProjectionBackend(
#     ManifoldDiff.FiniteDifferencesBackend()
# )
# gradf1_FD(M, p) = ManifoldDiff.gradient(M, f, p, r_backend)

# m1 = gradient_descent(M, f, gradf1_FD, data2[1]; maxiter=1000)

# ################################
using Optimization,
    OptimizationManopt, Symbolics, Manifolds, Random, LinearAlgebra, SymbolicAnalysis

M = SymmetricPositiveDefinite(5)
m = 100
σ = 0.005
q = Matrix{Float64}(LinearAlgebra.I(5)) .+ 2.0

data2 = [exp(M, q, σ * rand(M; vector_at = q)) for i in 1:m];

f(x, p = nothing) = sum(SymbolicAnalysis.distance(M, data2[i], x)^2 for i in 1:5)
optf = OptimizationFunction(f, Optimization.AutoZygote())
prob = OptimizationProblem(optf, data2[1]; manifold = M, structural_analysis = true)

opt = OptimizationManopt.GradientDescentOptimizer()
@time sol = solve(prob, opt, maxiters = 100)
@test sol.objective < 1.0e-2

M = SymmetricPositiveDefinite(5)
xs = [rand(5) for i in 1:5]

function f(S, p = nothing)
    return 1 / length(xs) * sum(SymbolicAnalysis.log_quad_form(x, S) for x in xs) +
        1 / 5 * logdet(inv(S))
end

optf = OptimizationFunction(f, Optimization.AutoZygote())
prob = OptimizationProblem(
    optf,
    Array{Float64}(LinearAlgebra.I(5));
    manifold = M,
    structural_analysis = true
)

opt = OptimizationManopt.GradientDescentOptimizer()
sol = solve(prob, opt, maxiters = 10)

A = randn(5, 5) #initialize random matrix
A = A * A' #make it a SPD matrix

function matsqrt(X, p = nothing) #setup objective function
    return SymbolicAnalysis.sdivergence(X, A) +
        SymbolicAnalysis.sdivergence(X, Matrix{Float64}(LinearAlgebra.I(5)))
end

optf = OptimizationFunction(matsqrt, Optimization.AutoZygote()) #setup oracles
prob = OptimizationProblem(optf, A / 2, manifold = M, structural_analysis = true) #setup problem with manifold and initial point

sol = solve(prob, GradientDescentOptimizer(), maxiters = 1000) #solve the problem
@test sqrt(A) ≈ sol.u rtol = 1.0e-3

ex = matsqrt(X) |> unwrap
ex = SymbolicAnalysis.propagate_sign(ex)
ex = SymbolicAnalysis.propagate_gcurvature(ex, M)

@test SymbolicAnalysis.getgcurvature(ex) == SymbolicAnalysis.GConvex

##Diagonal loading
@variables X[1:5, 1:5]

ex = tr(inv(X)) + logdet(X) |> unwrap
@test analyze(ex, M).gcurvature == SymbolicAnalysis.GConvex

γ = 1 / 2
ex = (tr(X + γ * I(5)))^(2) |> unwrap

@test analyze(ex, M).gcurvature == SymbolicAnalysis.GConvex

d = 10
n = 50
@variables X[1:d, 1:d]

@variables x[1:5] X[1:5, 1:5]
ex = SymbolicAnalysis.log_quad_form(x, inv(X)) |> unwrap
ex = SymbolicAnalysis.propagate_sign(ex)
ex = SymbolicAnalysis.propagate_gcurvature(ex, M)
@test SymbolicAnalysis.getgcurvature(ex) == SymbolicAnalysis.GConvex

ys = [rand(5) for i in 1:5]
ex = SymbolicAnalysis.log_quad_form(ys, X) |> unwrap
ex = SymbolicAnalysis.propagate_sign(ex)
ex = SymbolicAnalysis.propagate_gcurvature(ex, M)
@test SymbolicAnalysis.getgcurvature(ex) == SymbolicAnalysis.GConvex

ex = SymbolicAnalysis.log_quad_form(ys, inv(X)) |> unwrap
ex = SymbolicAnalysis.propagate_sign(ex)
ex = SymbolicAnalysis.propagate_gcurvature(ex, M)
@test SymbolicAnalysis.getgcurvature(ex) == SymbolicAnalysis.GConvex

ex = SymbolicAnalysis.log_quad_form(ys, X) |> unwrap
ex = SymbolicAnalysis.propagate_sign(ex)
ex = SymbolicAnalysis.propagate_gcurvature(ex, M)
@test SymbolicAnalysis.getgcurvature(ex) == SymbolicAnalysis.GConvex

ex = sum(SymbolicAnalysis.eigsummax(log(X), 2)) |> unwrap
anres = analyze(ex, M)
@test anres.gcurvature == SymbolicAnalysis.GConvex

ex = sum(SymbolicAnalysis.schatten_norm(log(X), 3)) |> unwrap
anres = analyze(ex, M)
@test anres.gcurvature == SymbolicAnalysis.GConvex

ex = exp(SymbolicAnalysis.eigsummax(log(X), 2)) |> unwrap
anres = analyze(ex, M)
@test anres.gcurvature == SymbolicAnalysis.GConvex

ex = SymbolicAnalysis.sum_log_eigmax(X, 2) |> unwrap
anres = analyze(ex, M)
@test anres.gcurvature == SymbolicAnalysis.GConvex

ex = SymbolicAnalysis.sum_log_eigmax(exp, X, 2) |> unwrap
anres = analyze(ex, M)
@test anres.gcurvature == SymbolicAnalysis.GConvex

B = rand(5, 5)
B = B * B'
Ys = [rand(5, 5) for i in 1:5]
Ys = [Y * Y' for Y in Ys]
ex = tr(SymbolicAnalysis.affine_map(SymbolicAnalysis.conjugation, X, B, Ys[1])) |> unwrap
anres = analyze(ex, M)
@test anres.gcurvature == SymbolicAnalysis.GConvex

ex = SymbolicAnalysis.hadamard_product(X, B) |> unwrap
anres = analyze(ex, M)
@test anres.gcurvature == SymbolicAnalysis.GConvex

A = rand(5, 5)
A = A * A'
ex = logdet(SymbolicAnalysis.affine_map(SymbolicAnalysis.hadamard_product, X, A, B)) |>
    unwrap
anres = analyze(ex, M)
@test anres.gcurvature == SymbolicAnalysis.GConvex

# The manifold-conditional signs in the geodesic rule table are gated on `M` being
# passed, so `analyze(ex, M)` must still see them (`tr` is Positive on the SPD cone).
@test analyze(tr(X), M).sign == SymbolicAnalysis.Positive
@test analyze(eigmax(X), M).sign == SymbolicAnalysis.Positive
@test analyze(tr(X)).sign == SymbolicAnalysis.AnySign

# A Euclidean rule does not imply a geodesic one. `find_gcurvature` used to fall
# back to the Euclidean rule table for any atom without a geodesic rule and treat
# that curvature as geodesic, which certified `eigmin(X)` GConcave on the SPD cone
# from its Euclidean Concave. Refuted: between A and B below, both SPD, the
# geodesic midpoint gives eigmin 2.3890 against a chord of 3.8712 — below the
# chord, so concavity fails.
@test analyze(unwrap(eigmin(X)), M).gcurvature == SymbolicAnalysis.GUnknownCurvature
@test analyze(unwrap(sqrt(X)), M).gcurvature == SymbolicAnalysis.GUnknownCurvature
@test analyze(unwrap(norm(X)), M).gcurvature == SymbolicAnalysis.GUnknownCurvature

# the numeric witness, so the assertion above is not just a recorded opinion
let A = [7.7517 1.132; 1.132 8.8903], B = [2.8936 0.3831; 0.3831 0.7551]
    M2 = Manifolds.SymmetricPositiveDefinite(2)
    γ = Manifolds.shortest_geodesic(M2, A, B, 0.5)
    @test eigmin(Symmetric(γ)) < (eigmin(Symmetric(A)) + eigmin(Symmetric(B))) / 2
end

# Atoms with a real geodesic rule keep it, and a Euclidean rule may still supply a
# composition over an argument that already carries one.
@test analyze(unwrap(tr(X)), M).gcurvature == SymbolicAnalysis.GConvex
@test analyze(unwrap(eigmax(X)), M).gcurvature == SymbolicAnalysis.GConvex
@test analyze(unwrap(logdet(X)), M).gcurvature == SymbolicAnalysis.GLinear
@test analyze(unwrap(log(tr(X))), M).gcurvature == SymbolicAnalysis.GConvex

# `logdet` on SPD is GLinear but not nonnegative (logdet(0.01I₂) ≈ -9.21). The
# false Positive sign certified `abs(logdet(X)+tr(X))` as GConvex; a geodesic
# midpoint can sit above the chord (convexity requires ≤).
@variables X2[1:2, 1:2]
M2 = SymmetricPositiveDefinite(2)
@test analyze(logdet(X2), M2).sign != SymbolicAnalysis.Positive
@test analyze(abs(logdet(X2) + tr(X2)), M2).gcurvature != SymbolicAnalysis.GConvex
let A = [
        0.19439548280694857 -0.027683931785232742;
        -0.027683931785232742 0.01892741465443008
    ],
        B = [
        0.9010733896285387 -1.365891310427168;
        -1.365891310427168 2.105439954050302
    ]
    f = P -> abs(logdet(P) + tr(P))
    γ = Manifolds.shortest_geodesic(M2, A, B, 0.5)
    @test f(γ) > (f(A) + f(B)) / 2
end
# Same false-Positive class: log(y'Xy) and sum(log λᵢ) go negative near 0.
@test analyze(SymbolicAnalysis.log_quad_form([1.0, 0.0], X2), M2).sign !=
    SymbolicAnalysis.Positive
@test analyze(SymbolicAnalysis.sum_log_eigmax(X2, 2), M2).sign !=
    SymbolicAnalysis.Positive

# Geodesic rules hold only at a point of the manifold: `tr(X .- 5)` can be
# negative on the SPD cone, and `tr(-X) = -tr(X)` is geodesically concave.
@testset "GDCP rules need a manifold-point argument" begin
    @variables X2[1:2, 1:2]
    xexpx = SymbolicAnalysis.LogExpFunctions.xexpx
    M2 = SymmetricPositiveDefinite(2)
    not_positive = [
        tr(X2 .- 5), tr(X2 - 5 * I(2)), tr(-X2), eigmax(-X2), tr(-X2 + I(2)),
        SymbolicAnalysis.eigsummax(log(X2), 1),
    ]
    for ex in not_positive
        @test analyze(ex, M2).sign != SymbolicAnalysis.Positive
    end
    not_convex = [
        xexpx(tr(X2 .- 5)), inv(tr(X2 .- 5)), xexpx(tr(X2 - 5 * I(2))), xexpx(tr(-X2)),
        abs(eigmax(-X2)), eigmax(-X2)^2,
    ]
    for ex in not_convex
        @test analyze(ex, M2).curvature != SymbolicAnalysis.Convex
    end
    not_gconvex = [
        tr(-X2), eigmax(-X2), tr(sin.(X2)), tr(-X2 + I(2)), tr(inv(-X2)),
        tr(inv(X2 .- 5)), tr(inv(X2 + I(2))), xexpx(tr(X2 .- 5)), tr(X2 - 5 * I(2))^2,
        abs(eigmax(-X2)), SymbolicAnalysis.eigsummax(log(X2), 1)^2,
        log(tr(X2 - I(2))),
    ]
    for ex in not_gconvex
        @test analyze(ex, M2).gcurvature != SymbolicAnalysis.GConvex
    end
    # A point's isometric images and PSD shifts keep the rule.
    for ex in [tr(X2), tr(inv(X2)), tr(2 * X2), tr(X2 + 0.5 * I(2)), eigmax(inv(X2))]
        @test analyze(ex, M2).sign == SymbolicAnalysis.Positive
    end
    for ex in [tr(inv(X2))^2, tr(2 * X2)^2, tr(X2 + 0.5 * I(2))^2, logdet(X2 + I(2))]
        @test analyze(ex, M2).gcurvature == SymbolicAnalysis.GConvex
    end
    Cns = [0.0 3.0; 0.0 0.0]
    for ex in [
            eigmax(X2 + Cns), eigmax(X2 .+ Cns), eigmax(X2 .- [0.0, 3.0]),
            SymbolicAnalysis.eigsummax(X2 + Cns, 1),
        ]
        @test analyze(ex, M2).gcurvature != SymbolicAnalysis.GConvex
    end
    @test analyze(eigmax(X2 .- 5), M2).gcurvature == SymbolicAnalysis.GConvex
    @test analyze(eigmax(X2 + [1.0 4.0; 4.0 -2.0]), M2).gcurvature ==
        SymbolicAnalysis.GConvex
    @test analyze(tr(X2), Lorentz(2)).sign != SymbolicAnalysis.Positive

    A = [0.5 0.1; 0.1 0.3]
    B = [4.0 -1.0; -1.0 3.0]
    @test tr(A .- 5) < 0
    γ = Manifolds.shortest_geodesic(M2, A, B, 0.5)
    @test tr(-γ) > (tr(-A) + tr(-B)) / 2
    xe = P -> (t = tr(P .- 5); t * exp(t))
    @test xe((A + B) / 2) > (xe(A) + xe(B)) / 2

    # det(Cind) = -1//10^18 but float.(Cind) rounds to the singular PSD ones(2, 2).
    Cind = [1 // 1 1 // 1; 1 // 1 (big(10)^18 - 1) // big(10)^18]
    @test analyze(tr(inv(X2 + Cind)), M2).sign != SymbolicAnalysis.Positive
    P = Matrix{Rational{BigInt}}(I, 2, 2) / (4 * big(10)^18)
    @test tr(inv(P + Cind)) < 0
    @test analyze(tr(inv(X2 + [1 // 1 1 // 2; 1 // 2 1 // 1])), M2).sign ==
        SymbolicAnalysis.Positive
    @test !SymbolicAnalysis.psd_const([1.0 1.0; 1.0 1 - eps() / 2])
    @test !SymbolicAnalysis.psd_const([Inf 0.0; 0.0 1.0])
    @test SymbolicAnalysis.psd_const([1 0; 0 0])

    # Curvature, like sign, needs each symbolic argument to be a point of `M`.
    @variables Y3[1:3, 1:3] p3[1:3]
    for (ex, M) in [
            (tr(Y3), M2), (tr(inv(Y3)), M2), (logdet(Y3 + I(3)), M2),
            (logdet(X2 + I(2)), Lorentz(2)), (log(tr(p3)), Lorentz(2)),
        ]
        @test analyze(ex, M).gcurvature == SymbolicAnalysis.GUnknownCurvature
    end
    @test analyze(tr(Y3), SymmetricPositiveDefinite(3)).gcurvature == SymbolicAnalysis.GConvex
    for ex in [SymbolicAnalysis.sdivergence(inv(X2), Y3), SymbolicAnalysis.sdivergence(Y3, inv(X2))]
        @test analyze(ex, M2).gcurvature == SymbolicAnalysis.GUnknownCurvature
    end

    # Non-point arguments of the SPD special cases must be constants.
    v = [sin(X2[1, 1]), 0]
    @test analyze(log(SymbolicAnalysis.quad_form(v, X2)), M2).gcurvature ==
        SymbolicAnalysis.GUnknownCurvature
    lq = P -> (w = [sin(P[1, 1]), 0]; log(dot(w, P * w)))
    I2 = Matrix(1.0I, 2, 2)
    @test lq(Manifolds.shortest_geodesic(M2, I2, 2I2, 0.5)) > (lq(I2) + lq(2I2)) / 2
    @test analyze(log(SymbolicAnalysis.quad_form([1.0, 2.0], X2)), M2).gcurvature ==
        SymbolicAnalysis.GConvex

    # A free parameter is admitted only if its known size fits its role.
    @variables v2[1:2] v3[1:3]
    for ex in [
            SymbolicAnalysis.log_quad_form(v3, X2), SymbolicAnalysis.quad_form(v3, X2),
            log(SymbolicAnalysis.quad_form(v3, X2)), SymbolicAnalysis.log_quad_form(v3, inv(X2)),
        ]
        @test analyze(ex, M2).gcurvature == SymbolicAnalysis.GUnknownCurvature
    end
    for ex in [
            SymbolicAnalysis.log_quad_form(v2, X2), SymbolicAnalysis.quad_form(v2, X2),
            log(SymbolicAnalysis.quad_form(v2, X2)),
        ]
        @test analyze(ex, M2).gcurvature == SymbolicAnalysis.GConvex
    end
    @variables z0[fill(1:1, 0)...]
    for y in [z0, X2, [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0]]]
        @test analyze(SymbolicAnalysis.log_quad_form(y, X2), M2).gcurvature ==
            SymbolicAnalysis.GUnknownCurvature
    end
    for ys in ([[1.0, 0.0]], [[1.0, 0.0], [0.0, 1.0]], [[1.0, 0.0], [0.0, 1.0], [1.0, 1.0]])
        @test analyze(SymbolicAnalysis.log_quad_form(ys, X2), M2).gcurvature ==
            SymbolicAnalysis.GConvex
    end
end

@testset "GDCP atom role and dimension validation" begin
    @variables X[1:2, 1:2] v[1:2] v3[1:3] p3[1:3] p4[1:4] A3[1:3, 1:3] A4[1:4, 1:4] s
    M = SymmetricPositiveDefinite(2)
    L2 = Lorentz(2)
    B2 = [2.0 0.5; 0.5 1.0]
    B3 = Matrix(1.0I, 3, 3)
    uv = unwrap(v)
    uX = unwrap(X)
    term = SymbolicAnalysis.SymbolicUtils.term

    bad_cases = [
        (tr(v), M),
        (term(sum, uv), M),
        (term(diag, uv), M),
        (eigmax(v), M),
        (logdet(v), M),
        (SymbolicAnalysis.eigsummax(v, 1), M),
        (term(SymbolicAnalysis.eigsummax, uX, uv), M),
        (term(SymbolicAnalysis.schatten_norm, uv, 2), M),
        (term(SymbolicAnalysis.schatten_norm, uX, uv), M),
        (term(SymbolicAnalysis.sum_log_eigmax, uv, 1), M),
        (term(SymbolicAnalysis.sum_log_eigmax, exp, uv, 1), M),
        (term(SymbolicAnalysis.sum_log_eigmax, uX, uv), M),
        (term(SymbolicAnalysis.scalar_mat, uv, 2), M),
        (term(SymbolicAnalysis.scalar_mat, uX, uv), M),
        (SymbolicAnalysis.sdivergence(X, rand(2, 3)), M),
        (SymbolicAnalysis.sdivergence(X, B3), M),
        (SymbolicAnalysis.sdivergence(X, v), M),
        (SymbolicAnalysis.sdivergence(v, X), M),
        (SymbolicAnalysis.conjugation(X, randn(3, 3)), M),
        (SymbolicAnalysis.conjugation(X, rand(2, 3)), M),
        (SymbolicAnalysis.conjugation(v, B2), M),
        (SymbolicAnalysis.hadamard_product(X, randn(3, 3)), M),
        (SymbolicAnalysis.hadamard_product(v, B2), M),
        (SymbolicAnalysis.affine_map(SymbolicAnalysis.conjugation, X, B3, B2), M),
        (SymbolicAnalysis.affine_map(SymbolicAnalysis.conjugation, X, B2, B3), M),
        (SymbolicAnalysis.affine_map(tr, X, B3), M),
        (SymbolicAnalysis.affine_map(diag, X, B3), M),
        (SymbolicAnalysis.affine_map(SymbolicAnalysis.hadamard_product, X, B3, B2), M),
        (SymbolicAnalysis.affine_map(tr, v, B2), M),
        (Manifolds.distance(M, I(3), X), M),
        (Manifolds.distance(M, I(2), v), M),
        (Manifolds.distance(L2, [0.0, 0.0, 1.0], p3), M),
        (SymbolicAnalysis.quad_form(v3, X), M),
        (SymbolicAnalysis.log_quad_form(v3, X), M),
        (SymbolicAnalysis.lorentz_log_barrier(p4), L2),
        (SymbolicAnalysis.lorentz_homogeneous_quadratic(A4, p3), L2),
        (SymbolicAnalysis.lorentz_homogeneous_diagonal(p4, p3), L2),
        (SymbolicAnalysis.lorentz_nonhomogeneous_quadratic(A4, p3, 1.0, p3), L2),
        (SymbolicAnalysis.lorentz_least_squares(randn(4, 3), randn(3), p3), L2),
        (Manifolds.distance(L2, [0.0, 0.0, 0.0, 1.0], p3), L2),
    ]
    for (ex, m) in bad_cases
        r = analyze(ex, m)
        @test r.gcurvature == SymbolicAnalysis.GUnknownCurvature
        @test r.sign == SymbolicAnalysis.AnySign
    end

    valid_controls = [
        (tr(X), M, SymbolicAnalysis.Positive, SymbolicAnalysis.GConvex),
        (term(sum, uX), M, SymbolicAnalysis.Positive, SymbolicAnalysis.GConvex),
        (term(diag, uX), M, SymbolicAnalysis.Positive, SymbolicAnalysis.GConvex),
        (eigmax(X), M, SymbolicAnalysis.Positive, SymbolicAnalysis.GConvex),
        (logdet(X), M, SymbolicAnalysis.AnySign, SymbolicAnalysis.GLinear),
        (SymbolicAnalysis.eigsummax(X, 1), M, SymbolicAnalysis.Positive, SymbolicAnalysis.GConvex),
        (SymbolicAnalysis.schatten_norm(X, 2), M, SymbolicAnalysis.Positive, SymbolicAnalysis.GConvex),
        (SymbolicAnalysis.sum_log_eigmax(X, 1), M, SymbolicAnalysis.AnySign, SymbolicAnalysis.GConvex),
        (SymbolicAnalysis.sum_log_eigmax(exp, X, 1), M, SymbolicAnalysis.AnySign, SymbolicAnalysis.GConvex),
        (SymbolicAnalysis.scalar_mat(X, 2), M, SymbolicAnalysis.Positive, SymbolicAnalysis.GConvex),
        (SymbolicAnalysis.sdivergence(X, B2), M, SymbolicAnalysis.Positive, SymbolicAnalysis.GConvex),
        (SymbolicAnalysis.conjugation(X, B2), M, SymbolicAnalysis.Positive, SymbolicAnalysis.GConvex),
        (SymbolicAnalysis.hadamard_product(X, B2), M, SymbolicAnalysis.Positive, SymbolicAnalysis.GConvex),
        (SymbolicAnalysis.affine_map(SymbolicAnalysis.conjugation, X, B2, B2), M, SymbolicAnalysis.Positive, SymbolicAnalysis.GConvex),
        (SymbolicAnalysis.affine_map(tr, X, B2), M, SymbolicAnalysis.Positive, SymbolicAnalysis.GConvex),
        (SymbolicAnalysis.affine_map(diag, X, B2), M, SymbolicAnalysis.Positive, SymbolicAnalysis.GConvex),
        (SymbolicAnalysis.affine_map(SymbolicAnalysis.hadamard_product, X, B2, B2), M, SymbolicAnalysis.Positive, SymbolicAnalysis.GConvex),
        (Manifolds.distance(M, I(2), X), M, SymbolicAnalysis.Positive, SymbolicAnalysis.GConvex),
        (SymbolicAnalysis.quad_form(v, X), M, SymbolicAnalysis.AnySign, SymbolicAnalysis.GConvex),
        (SymbolicAnalysis.log_quad_form(v, X), M, SymbolicAnalysis.AnySign, SymbolicAnalysis.GConvex),
        (SymbolicAnalysis.lorentz_log_barrier(p3), L2, SymbolicAnalysis.AnySign, SymbolicAnalysis.GConvex),
        (SymbolicAnalysis.lorentz_homogeneous_quadratic(A3, p3), L2, SymbolicAnalysis.AnySign, SymbolicAnalysis.GConvex),
        (SymbolicAnalysis.lorentz_homogeneous_diagonal([2.0, 2.0, 1.0], p3), L2, SymbolicAnalysis.AnySign, SymbolicAnalysis.GConvex),
        (SymbolicAnalysis.lorentz_nonhomogeneous_quadratic(A3, [0.0, 0.0, 1.0], s, p3), L2, SymbolicAnalysis.AnySign, SymbolicAnalysis.GConvex),
        (SymbolicAnalysis.lorentz_least_squares(B3, [0.0, 0.0, -1.0], p3), L2, SymbolicAnalysis.Positive, SymbolicAnalysis.GConvex),
        (Manifolds.distance(L2, [0.0, 0.0, 1.0], p3), L2, SymbolicAnalysis.Positive, SymbolicAnalysis.GConvex),
    ]
    for (ex, m, expected_sign, expected_gcurv) in valid_controls
        r = analyze(ex, m)
        @test r.sign == expected_sign
        @test r.gcurvature == expected_gcurv
    end
end
