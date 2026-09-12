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
# `log1p(x)^2` is not convex: d²/dx² log(1+x)² = 2(1 - log(1+x))/(1+x)², which is
# negative for x > e - 1. This certified Convex only because `log1p` declared a
# `Negative` sign, which made `increasing_if_positive` hand `^2` a `Decreasing`
# slot. With the sign corrected to AnySign the constraint is no longer DCP.
@test getcurvature(ex) == SymbolicAnalysis.UnknownCurvature

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

# logsumexp over a symbolic vector must stay an unevaluated atom — Symbolics'
# own method expands it to log(sum(exp, z)), erasing the atom and yielding
# UnknownCurvature; logsumexp itself is convex. (Also fixes the rule domain,
# which required ndims==2 and so never matched a vector.)
@variables z[1:4]

@test Symbolics.operation(LogExpFunctions.logsumexp(z) |> unwrap) === LogExpFunctions.logsumexp

ex = LogExpFunctions.logsumexp(z) |> unwrap
ex = propagate_curvature(propagate_sign(ex))
@test getcurvature(ex) == SymbolicAnalysis.Convex
# logistic (the sigmoid 1/(1+exp(-x))) is NOT globally convex: f'' is positive
# for x<0, zero at x=0, negative for x>0 (inflection at 0), so no single-curvature
# DCP rule is valid. The rule was dead today (logistic expands to `/`) but a latent
# trap: it must carry no rule and analyze as UnknownCurvature, not a false Convex.
@variables x
@test !SymbolicAnalysis.hasdcprule(logistic)
ex = logistic(x) |> unwrap
ex = propagate_curvature(propagate_sign(ex))
@test getcurvature(ex) == SymbolicAnalysis.UnknownCurvature
# eigmax/eigmin of a symbolic matrix must build symbolic atoms (Convex/Concave),
# not fall through to numeric eigmax and crash on eigvals!(::Matrix{Num}).
@variables X[1:3, 1:3]

ex = eigmax(X) |> unwrap
ex = propagate_curvature(propagate_sign(ex))
@test getcurvature(ex) == SymbolicAnalysis.Convex

ex = eigmin(X) |> unwrap
ex = propagate_curvature(propagate_sign(ex))
@test getcurvature(ex) == SymbolicAnalysis.Concave

# A matrix assembled from scalar variables traces to a SymbolicUtils.array_literal
# term. Without a rule for it the assembled matrix has no curvature, so every atom
# taking a matrix argument fails to compose: logdet/eigmax of such a matrix
# analyzed as UnknownCurvature even though each entry is affine.
@variables m[1:3]
ms = Symbolics.scalarize(m)
M = [ms[1] ms[2]; ms[2] ms[3]]

ex = SymbolicAnalysis.logdet(M) |> unwrap
ex = propagate_curvature(propagate_sign(ex))
@test getcurvature(ex) == SymbolicAnalysis.Concave

ex = -SymbolicAnalysis.logdet(M) |> unwrap        # the log-det barrier
ex = propagate_curvature(propagate_sign(ex))
@test getcurvature(ex) == SymbolicAnalysis.Convex

ex = eigmax(M) |> unwrap
ex = propagate_curvature(propagate_sign(ex))
@test getcurvature(ex) == SymbolicAnalysis.Convex

ex = tr(M) |> unwrap
ex = propagate_curvature(propagate_sign(ex))
@test getcurvature(ex) == SymbolicAnalysis.Affine

# The rule must not launder a non-affine entry into an affine matrix: `logdet` is
# concave with AnyMono monotonicity, so it composes only over affine arguments.
# Assembling from x^2 must therefore leave the result uncertified.
@variables q
Mq = [q^2 0.0; 0.0 q^2]
ex = SymbolicAnalysis.logdet(Mq) |> unwrap
ex = propagate_curvature(propagate_sign(ex))
@test getcurvature(ex) == SymbolicAnalysis.UnknownCurvature

# Unsound certificates (#156): these returned a curvature that was wrong in the
# permissive direction, which is worse than UnknownCurvature — a false Affine or
# Convex can route a non-convex problem to a conic solver.

# `find_curvature` fell through to Affine for a wrapped `Num` (a `Num` is not
# itself a call), so a container of non-affine expressions certified as Affine.
@variables w[1:3]
ws = Symbolics.scalarize(w)
@test SymbolicAnalysis.getcurvature(exp.(ws)) == SymbolicAnalysis.Convex
@test SymbolicAnalysis.getcurvature(log.(ws)) == SymbolicAnalysis.Concave
@test SymbolicAnalysis.find_curvature(exp(ws[1])) == SymbolicAnalysis.Convex
# a container mixing convex and concave entries has no single curvature
@test SymbolicAnalysis.getcurvature([exp(ws[1]), log(ws[2])]) ==
    SymbolicAnalysis.UnknownCurvature

# `dot` is bilinear: affine only when one side is constant. dot(x, x) is ‖x‖².
@test SymbolicAnalysis.analyze(unwrap(dot(w, w))).curvature ==
    SymbolicAnalysis.UnknownCurvature
@test SymbolicAnalysis.analyze(unwrap(dot(ones(3), ws))).curvature ==
    SymbolicAnalysis.Affine

# The scalar power laws must not be applied to a matrix base: tr(X^2) for an
# unconstrained X is indefinite (it contains cross terms X[i,j]*X[j,i]).
@variables Xm[1:3, 1:3]
@test SymbolicAnalysis.analyze(unwrap(tr(Xm * Xm))).curvature ==
    SymbolicAnalysis.UnknownCurvature
@test SymbolicAnalysis.analyze(unwrap(tr(Xm))).curvature == SymbolicAnalysis.Affine

# The elementwise `x .^ i` traces to `broadcast(^, x, i)` with an array base too,
# but applies the power pointwise, so the scalar laws do hold there and the
# matrix-power guard above must not reject it.
@test SymbolicAnalysis.analyze(unwrap(w .^ 2)).curvature == SymbolicAnalysis.Convex
@test SymbolicAnalysis.analyze(unwrap(sum(w .^ 2))).curvature == SymbolicAnalysis.Convex
@test SymbolicAnalysis.analyze(unwrap(sum(w .^ 3))).curvature == SymbolicAnalysis.Convex
@test SymbolicAnalysis.analyze(unwrap(sum(w .^ 0.5))).curvature ==
    SymbolicAnalysis.Concave
@test SymbolicAnalysis.analyze(unwrap(sum(Xm .^ 2))).curvature == SymbolicAnalysis.Convex
@test SymbolicAnalysis.analyze(unwrap(norm(w .^ 2, 1))).curvature ==
    SymbolicAnalysis.Convex

# A symbolic exponent has no fixed curvature law; it must degrade rather than throw
# from the `isinteger` comparisons.
@variables p
@test SymbolicAnalysis.analyze(unwrap(q^p)).curvature == SymbolicAnalysis.UnknownCurvature
@test SymbolicAnalysis.analyze(unwrap(sum(w .^ p))).curvature ==
    SymbolicAnalysis.UnknownCurvature

# A wrong sign becomes a wrong curvature one level up, because
# `increasing_if_positive` turns it into a monotonicity `abs` then composes over.

# A geodesic rule's sign holds only on its manifold — `tr` is Positive on the SPD
# cone — so it must not be consulted when no manifold was passed. It was, and
# `abs(tr(X) + a^2)` certified Convex while |−2 + t²| is concave near t = 0.
@test SymbolicAnalysis.getsign(propagate_sign(unwrap(tr(Xm)))) == SymbolicAnalysis.AnySign
@test SymbolicAnalysis.getsign(propagate_sign(unwrap(eigmax(Xm)))) ==
    SymbolicAnalysis.AnySign
@test SymbolicAnalysis.analyze(unwrap(abs(tr(Xm) + x^2))).curvature ==
    SymbolicAnalysis.UnknownCurvature
@test SymbolicAnalysis.analyze(unwrap(abs(eigmax(Xm) + x^2))).curvature ==
    SymbolicAnalysis.UnknownCurvature

# `log1p` is negative only on (-1, 0) and positive on all of (0, Inf), so the
# Negative it declared made `abs(log1p(a))` — strictly concave for a > 0 — Convex.
@test SymbolicAnalysis.getsign(propagate_sign(unwrap(log1p(x)))) ==
    SymbolicAnalysis.AnySign
@test SymbolicAnalysis.analyze(unwrap(abs(log1p(x)))).curvature ==
    SymbolicAnalysis.UnknownCurvature
# `lognormcdf` is the log of a CDF, so its Negative sign is genuine and
# `abs(lognormcdf(a)) = -lognormcdf(a)` stays Convex.
@test SymbolicAnalysis.analyze(unwrap(abs(SymbolicAnalysis.lognormcdf(x)))).curvature ==
    SymbolicAnalysis.Convex

# Directly: the squared concave atom above, numerically concave for x > e - 1
# (second difference at x = 3 is -0.0483 across h = 1e-2, 1e-3, 1e-4).
@test SymbolicAnalysis.analyze(unwrap(log1p(x)^2)).curvature ==
    SymbolicAnalysis.UnknownCurvature
