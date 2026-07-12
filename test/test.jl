using SymbolicAnalysis
using SymbolicAnalysis: propagate_curvature, propagate_sign, getcurvature, getsign
using Symbolics, SymbolicAnalysis.LogExpFunctions
using Symbolics: unwrap
using LinearAlgebra, Test

@variables x y
y = setmetadata(
    y,
    SymbolicAnalysis.VarDomain,
    Symbolics.DomainSets.HalfLine{Number, :open}()
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
@test getcurvature(ex) == SymbolicAnalysis.Affine
@test getsign(ex) == SymbolicAnalysis.AnySign

ex = exp.(x) |> unwrap
ex = propagate_curvature(propagate_sign(ex))
@test getcurvature(ex) == SymbolicAnalysis.Convex
@test getsign(ex) == SymbolicAnalysis.Positive

##vector * scalar gets simplified

@variables x y z
obj = x^2 + y^2 + z^2 |> unwrap

ex = propagate_curvature(propagate_sign(obj))
@test getcurvature(ex) == SymbolicAnalysis.Convex
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

# Composition-rule regressions: affine atoms over non-affine arguments,
# monotonicity tuples shorter than the argument list, and
# value-of-p-dependent curvature for norm.

@variables x y

# Single-entry monotonicity tuples apply to every argument; previously the
# second argument of max fell off the tuple (misclassified) and the concave
# branch indexed the tuple directly (min of two variables threw).
ex = max(y, x^2) |> unwrap
ex = propagate_curvature(propagate_sign(ex))
@test getcurvature(ex) == SymbolicAnalysis.Convex

ex = max(x^2, y) |> unwrap
ex = propagate_curvature(propagate_sign(ex))
@test getcurvature(ex) == SymbolicAnalysis.Convex

ex = min(x, y) |> unwrap
ex = propagate_curvature(propagate_sign(ex))
@test getcurvature(ex) == SymbolicAnalysis.Concave

ex = min(y, sqrt(x)) |> unwrap
ex = propagate_curvature(propagate_sign(ex))
@test getcurvature(ex) == SymbolicAnalysis.Concave

# Affine atoms compose as convex (concave) when their arguments do; previously
# any non-affine argument of an affine atom returned UnknownCurvature.
@variables X[1:3, 1:3]

ex = tr(exp.(X)) |> unwrap
ex = propagate_curvature(propagate_sign(ex))
@test getcurvature(ex) == SymbolicAnalysis.Convex

ex = tr(log.(X)) |> unwrap
ex = propagate_curvature(propagate_sign(ex))
@test getcurvature(ex) == SymbolicAnalysis.Concave

ex = tr(X) |> unwrap
ex = propagate_curvature(propagate_sign(ex))
@test getcurvature(ex) == SymbolicAnalysis.Affine

# norm(x, p) curvature depends on the value of p: convex for p >= 1, concave
# for 0 < p < 1 only on a nonnegative argument (unknown otherwise, previously
# claimed convex), and unknown for p <= 0.
@variables z[1:4]

ex = norm(z, 2) |> unwrap
ex = propagate_curvature(propagate_sign(ex))
@test getcurvature(ex) == SymbolicAnalysis.Convex

ex = norm(z, 0.5) |> unwrap
ex = propagate_curvature(propagate_sign(ex))
@test getcurvature(ex) == SymbolicAnalysis.UnknownCurvature

zpos = setmetadata(unwrap(z), SymbolicAnalysis.Sign, SymbolicAnalysis.Positive)
ex = norm(Symbolics.wrap(zpos), 0.5) |> unwrap
ex = propagate_curvature(ex)
@test getcurvature(ex) == SymbolicAnalysis.Concave

ex = norm(z, -1) |> unwrap
ex = propagate_curvature(propagate_sign(ex))
@test getcurvature(ex) == SymbolicAnalysis.UnknownCurvature

# sum/map over symbolic arrays trace to SymbolicUtils.Mapreducer/Mapper
# operations, not to `sum`/`map` themselves; previously they always analyzed
# as UnknownCurvature.
ex = sum(exp.(z)) |> unwrap
ex = propagate_curvature(propagate_sign(ex))
@test getcurvature(ex) == SymbolicAnalysis.Convex

ex = sum(log.(z)) |> unwrap
ex = propagate_curvature(propagate_sign(ex))
@test getcurvature(ex) == SymbolicAnalysis.Concave

ex = sum(z) |> unwrap
ex = propagate_curvature(propagate_sign(ex))
@test getcurvature(ex) == SymbolicAnalysis.Affine

ex = map(exp, z) |> unwrap
ex = propagate_curvature(propagate_sign(ex))
@test getcurvature(ex) == SymbolicAnalysis.Convex

# maximum/minimum over symbolic arrays reduce with max/min, tracing to
# SymbolicUtils.Mapreducer{identity, max}/{identity, min} rather than
# maximum/minimum; previously they analyzed as UnknownCurvature.
@variables z[1:4]

ex = maximum(z) |> unwrap
ex = propagate_curvature(propagate_sign(ex))
@test getcurvature(ex) == SymbolicAnalysis.Convex

ex = minimum(z) |> unwrap
ex = propagate_curvature(propagate_sign(ex))
@test getcurvature(ex) == SymbolicAnalysis.Concave

ex = maximum(exp.(z)) |> unwrap
ex = propagate_curvature(propagate_sign(ex))
@test getcurvature(ex) == SymbolicAnalysis.Convex
