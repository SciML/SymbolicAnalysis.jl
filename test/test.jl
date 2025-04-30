using SymbolicAnalysis
using Symbolics, SymbolicAnalysis.LogExpFunctions
using Symbolics: unwrap
using LinearAlgebra, Test

@variables x y
y = setmetadata(
    y,
    SymbolicAnalysis.VarDomain,
    Symbolics.DomainSets.HalfLine{Number,:open}(),
)
ex1 = exp(y) - log(y) |> unwrap
ex1 = propagate_curvature(propagate_sign(ex1))

@test getcurvature(ex1) == SymbolicAnalysis.Convex
@test getsign(ex1) == SymbolicAnalysis.AnySign

ex2 = -sqrt(x^2) |> unwrap
ex2 = propagate_curvature(propagate_sign(ex2))

@test getcurvature(ex2) == SymbolicAnalysis.UnknownCurvature
@test getsign(ex2) == SymbolicAnalysis.Negative

ex = -1 * LogExpFunctions.xlogx(x) |> unwrap
ex = propagate_curvature(propagate_sign(ex))
@test getcurvature(ex) == SymbolicAnalysis.Concave
@test getsign(ex) == SymbolicAnalysis.AnySign

ex = 2 * abs(x) - 1 |> unwrap
ex = propagate_curvature(propagate_sign(ex))
@test getcurvature(ex) == SymbolicAnalysis.Convex
@test getsign(ex) == SymbolicAnalysis.AnySign

# x = setmetadata(x, SymbolicAnalysis.Sign, SymbolicAnalysis.Positive)
ex = abs(x)^2 |> unwrap
ex = propagate_curvature(propagate_sign(ex))
@test getcurvature(ex) == SymbolicAnalysis.Convex
@test getsign(ex) == SymbolicAnalysis.Positive

ex = abs(x)^2 + abs(x)^3 |> unwrap
ex = propagate_curvature(propagate_sign(ex))
@test getcurvature(ex) == SymbolicAnalysis.Convex
@test getsign(ex) == SymbolicAnalysis.Positive

@variables x[1:3] y
ex = x .- y |> unwrap
ex = propagate_curvature(propagate_sign(ex))
@test_broken getcurvature(ex) == SymbolicAnalysis.Affine
@test getsign(ex) == SymbolicAnalysis.AnySign

ex = exp.(x) |> unwrap
ex = propagate_curvature(propagate_sign(ex))
@test getcurvature(ex) == SymbolicAnalysis.Convex
@test getsign(ex) == SymbolicAnalysis.Positive

##vector * scalar gets simplified

@variables x y z
obj = x^2 + y^2 + z^2 |> unwrap

ex = propagate_curvature(propagate_sign(obj))
@test_broken getcurvature(ex) == SymbolicAnalysis.Convex
@test getsign(ex) == SymbolicAnalysis.Positive

cons = [
    x + y + z ~ 10
    log1p(x)^2 - log1p(z) ≲ 0
]

ex = propagate_curvature(propagate_sign(cons[1].lhs |> unwrap))
@test getcurvature(ex) == SymbolicAnalysis.Affine

ex = propagate_curvature(propagate_sign(cons[2].lhs))
@test getcurvature(ex) == SymbolicAnalysis.Convex

@variables x y z

ex = SymbolicAnalysis.quad_over_lin(x - y, 1 - max(x, y)) |> unwrap
ex = propagate_curvature(propagate_sign(ex))
@test getcurvature(ex) == SymbolicAnalysis.Convex
