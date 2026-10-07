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

# Hadamard and affine-map rules hold only for PSD constant factors (Schur product
# theorem; a PSD additive constant). Each claim below was refuted on random SPD pairs.
@testset "gDCP constant factors" begin
    had, aff, cj = SymbolicAnalysis.hadamard_product, SymbolicAnalysis.affine_map,
        SymbolicAnalysis.conjugation
    D = [1.0 0.0; 0.0 -1.0]
    Hind = [1.0 1.5; 1.5 1.0]
    Hpsd = [2.0 1.0; 1.0 1.0]
    I2 = Matrix(1.0I, 2, 2)
    A0 = [1.0 2.0; 0.5 3.0]
    @variables Bs[1:2, 1:2]
    gconvex(ex) = analyze(ex, M2).gcurvature == SymbolicAnalysis.GConvex
    positive(ex) = analyze(ex, M2).sign == SymbolicAnalysis.Positive

    @test !positive(had(X2, D))
    @test !gconvex(had(X2, D))
    @test !positive(tr(had(X2, D)))
    @test !gconvex(tr(had(X2, D)))
    @test !gconvex(sum(had(X2, D)))
    @test !gconvex(tr(aff(had, X2, D, I2)))
    @test !gconvex(logdet(had(X2, Hind)))
    @test !positive(aff(cj, X2, -0.5 * I2, I2))
    @test !gconvex(logdet(aff(cj, X2, -0.5 * I2, I2)))
    @test !positive(tr(aff(cj, X2, [1.0 4.0; 4.0 -2.0], A0)))
    @test !gconvex(log(tr(had(X2, D))))
    @test !gconvex(tr(inv(aff(cj, X2, -0.5 * I2, I2))))
    @test !positive(had(X2, Bs))
    @test !gconvex(logdet(aff(cj, X2, I2, Bs)))

    for ex in (
            had(X2, Hpsd), tr(had(X2, Hpsd)), aff(cj, X2, I2, A0),
            tr(aff(had, X2, Hpsd, I2)),
        )
        @test positive(ex)
        @test gconvex(ex)
    end
    @test gconvex(sum(had(X2, ones(2, 2))))
    @test gconvex(logdet(had(X2, Hpsd)))
    @test gconvex(log(tr(had(X2, Hpsd))))
    @test gconvex(logdet(aff(cj, X2, I2, A0)))

    # Matrices that are PSD only up to rounding are refused: Hadamard products
    # scale the negative eigenvalue by large entries of X.
    Beps = [1.0 0.0; 0.0 -2eps()]
    Bnear = [1.0 1.0; 1.0 1.0 - eps() / 2]
    for B in (Beps, Bnear)
        @test !positive(had(X2, B))
        @test !gconvex(had(X2, B))
        @test !gconvex(tr(had(X2, B)))
    end
    @test eigmin(Beps .* [1.0 0.0; 0.0 1.0e12]) < -1.0e-4
    @test det(Rational{BigInt}.(Bnear)) < 0
    for B in ([2 1; 1 1], [2 // 1 1; 1 1], big.(Hpsd))
        @test positive(had(X2, B))
        @test gconvex(had(X2, B))
    end

    Cind = [1.0 4.0; 4.0 -2.0]
    @test !gconvex(SymbolicAnalysis.sdivergence(X2, Cind))
    @test !positive(SymbolicAnalysis.sdivergence(X2, Cind))
    @test !gconvex(SymbolicAnalysis.sdivergence(X2, ones(2, 2)))
    @test gconvex(SymbolicAnalysis.sdivergence(X2, I2))

    let P = [2.0 0.0; 0.0 1.0], Q = [1.0 0.0; 0.0 4.0]
        @test eigmin(D .* P) < 0
        @test eigmin([0.1 0.0; 0.0 1.0] - 0.5I) < 0
        γ = Manifolds.shortest_geodesic(M2, P, Q, 0.5)
        @test tr(D .* γ) > (tr(D .* P) + tr(D .* Q)) / 2
    end
end
