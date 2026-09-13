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

# `/` had no rule at all, so every division fell through to UnknownCurvature.
# It is bilinear like `*`: DCP only when one side is constant.
@variables dn dd dv[1:3]
# a declared positive domain — `analyze` overwrites Sign metadata during
# propagation, so a VarDomain declaration is how the precondition is stated
dpos = setmetadata(
    dn, SymbolicAnalysis.VarDomain, Symbolics.DomainSets.HalfLine{Number, :open}()
)

# a nonzero constant denominator is an affine rescaling
@test SymbolicAnalysis.analyze(unwrap(dn / 2)).curvature == SymbolicAnalysis.Affine
@test SymbolicAnalysis.analyze(unwrap(exp(dn) / 2)).curvature == SymbolicAnalysis.Convex
@test SymbolicAnalysis.analyze(unwrap(log(dn) / 2)).curvature == SymbolicAnalysis.Concave
@test SymbolicAnalysis.analyze(unwrap(exp(dn) / -2)).curvature == SymbolicAnalysis.Concave
@test SymbolicAnalysis.analyze(unwrap(exp(dn) / -2)).sign == SymbolicAnalysis.Negative
@test SymbolicAnalysis.analyze(unwrap(sum(dv) / 3)).curvature == SymbolicAnalysis.Affine

# Symbolics folds constant division into a rational coefficient, so `getsign`
# has to cover Rational — it threw a MethodError before.
@test SymbolicAnalysis.analyze(unwrap(dn / 3 + dn / 3)).curvature ==
    SymbolicAnalysis.Affine
@test SymbolicAnalysis.analyze(unwrap(-dn / 3 - dn / 3)).sign == SymbolicAnalysis.AnySign

# A constant numerator is the `inv` atom, and Symbolics rewrites `inv(x)` and
# `x^-n` through `/`, so this is also what gives negative powers their curvature.
# It fires only on an established-positive denominator: `1/x` is defined and
# concave below zero, so a sign-unknown argument gets no certificate.
@test SymbolicAnalysis.analyze(unwrap(1 / dpos)).curvature == SymbolicAnalysis.Convex
@test SymbolicAnalysis.analyze(unwrap(1 / dpos)).sign == SymbolicAnalysis.Positive
@test SymbolicAnalysis.analyze(unwrap(-1 / dpos)).curvature == SymbolicAnalysis.Concave
@test SymbolicAnalysis.analyze(unwrap(dpos^(-2))).curvature == SymbolicAnalysis.Convex
@test SymbolicAnalysis.analyze(unwrap(dpos^(-3))).curvature == SymbolicAnalysis.Convex
@test SymbolicAnalysis.analyze(unwrap(1 / dn)).curvature ==
    SymbolicAnalysis.UnknownCurvature
@test SymbolicAnalysis.analyze(unwrap(dn^(-2))).curvature ==
    SymbolicAnalysis.UnknownCurvature

# neither side constant
@test SymbolicAnalysis.analyze(unwrap(dn / dd)).curvature ==
    SymbolicAnalysis.UnknownCurvature
@test SymbolicAnalysis.analyze(unwrap(exp(dn) / dd)).curvature ==
    SymbolicAnalysis.UnknownCurvature
@test SymbolicAnalysis.analyze(unwrap(dn^2 / dd)).curvature ==
    SymbolicAnalysis.UnknownCurvature


# Bilinear atoms are convex only when the argument they are linear in is
# constant. The static rule table cannot express that, so each certified both
# arguments symbolic — the same defect `dot` was fixed for.
@variables bx[1:2] bv[1:2] bP[1:2, 1:2] bh
bxs = Symbolics.scalarize(bx)

# quad_form is quadratic in x but *linear* in P: `quad_form([a], [b;;])` is
# `a^2*b`, whose second difference along (1, -1) is -2.
@test SymbolicAnalysis.analyze(unwrap(SymbolicAnalysis.quad_form(bxs, collect(bP)))).curvature ==
    SymbolicAnalysis.UnknownCurvature
# A constant but indefinite P is `x[1]^2 - x[2]^2`.
@test SymbolicAnalysis.analyze(
    unwrap(SymbolicAnalysis.quad_form(bxs, [1.0 0.0; 0.0 -1.0]))
).curvature == SymbolicAnalysis.UnknownCurvature
@test SymbolicAnalysis.analyze(
    unwrap(SymbolicAnalysis.quad_form(bxs, [2.0 1.0; 1.0 2.0]))
).curvature == SymbolicAnalysis.Convex

# `x'Px` is nondecreasing in each x[i] over the nonnegative orthant only when P is
# *entrywise* nonnegative; positive definiteness alone does not give it. With
# P = [1 -0.9; -0.9 1] (positive definite) the composition certified
# `quad_form(exp.(v), P)` as Convex while the midpoint exceeds the chord by 1.92
# between v = [0.64, 2.0] and [1.64, 2.3]. Asserted on the rule directly: the
# end-to-end path currently degrades for an unrelated reason (an `array_literal`
# argument carries no sign, so `increasing_if_positive` yields AnyMono anyway).
@test SymbolicAnalysis.dcprule(
    SymbolicAnalysis.quad_form, bxs, [1.0 -0.9; -0.9 1.0]
)[1].monotonicity[1] === SymbolicAnalysis.AnyMono
@test SymbolicAnalysis.dcprule(
    SymbolicAnalysis.quad_form, bxs, [2.0 1.0; 1.0 2.0]
)[1].monotonicity[1] === SymbolicAnalysis.increasing_if_positive

# dotsort is a pointwise max of bilinear forms; for length-1 vectors it is x*y.
@test SymbolicAnalysis.analyze(
    unwrap(SymbolicAnalysis.dotsort(bxs, Symbolics.scalarize(bv)))
).curvature == SymbolicAnalysis.UnknownCurvature
@test SymbolicAnalysis.analyze(
    unwrap(SymbolicAnalysis.dotsort(bxs, [1.0, 2.0]))
).curvature == SymbolicAnalysis.Convex

# huber is concave in its threshold: for abs(x) > M it is 2M*abs(x) - M^2.
@test SymbolicAnalysis.analyze(unwrap(SymbolicAnalysis.huber(2.0, bh))).curvature ==
    SymbolicAnalysis.UnknownCurvature
@test SymbolicAnalysis.analyze(unwrap(SymbolicAnalysis.huber(bh, 1.0))).curvature ==
    SymbolicAnalysis.Convex


# `exp(X)` on a matrix is the matrix exponential, not the scalar law applied
# pointwise. `expm` is neither operator convex nor operator monotone, so
# `sum(exp(X))` is indefinite — it fails even on SPD X with symmetric directions.
@variables ex2[1:2, 1:2] ev[1:2] ea
@test SymbolicAnalysis.analyze(unwrap(sum(exp(ex2)))).curvature ==
    SymbolicAnalysis.UnknownCurvature
@test SymbolicAnalysis.analyze(unwrap(tr(exp(ex2)))).curvature ==
    SymbolicAnalysis.UnknownCurvature
# the elementwise form arrives through `broadcast` and keeps the scalar law
@test SymbolicAnalysis.analyze(unwrap(tr(exp.(ex2)))).curvature ==
    SymbolicAnalysis.Convex
@test SymbolicAnalysis.analyze(unwrap(sum(exp.(ev)))).curvature ==
    SymbolicAnalysis.Convex
@test SymbolicAnalysis.analyze(unwrap(exp(ea))).curvature == SymbolicAnalysis.Convex


# The matrix rules for `sqrt`, `log` and `inv` are Loewner-order statements.
# Operator concavity licenses PSD-weighted functionals — `tr(f(X))` and
# `sum(f(X)) = e'f(X)e` — but not a single entry, and the smallest entry of a
# positive definite matrix can be off-diagonal. `minimum(sqrt(X))` certified
# Concave while its second difference over the SPD cone takes both signs
# (-1.63 and +1.03 over a 400-sample probe with symmetric directions).
@variables LX[1:2, 1:2] Lv[1:2]
@test SymbolicAnalysis.analyze(unwrap(minimum(sqrt(LX)))).curvature ==
    SymbolicAnalysis.UnknownCurvature
@test SymbolicAnalysis.analyze(unwrap(minimum(log(LX)))).curvature ==
    SymbolicAnalysis.UnknownCurvature
# PSD-weighted consumers keep their certificates
@test SymbolicAnalysis.analyze(unwrap(tr(sqrt(LX)))).curvature == SymbolicAnalysis.Concave
@test SymbolicAnalysis.analyze(unwrap(sum(sqrt(LX)))).curvature == SymbolicAnalysis.Concave
@test SymbolicAnalysis.analyze(unwrap(tr(inv(LX)))).curvature == SymbolicAnalysis.Convex
# the elementwise and vector reductions are untouched
@test SymbolicAnalysis.analyze(unwrap(minimum(log.(LX)))).curvature ==
    SymbolicAnalysis.Concave
@test SymbolicAnalysis.analyze(unwrap(minimum(Lv))).curvature == SymbolicAnalysis.Concave
# `maximum` needs no guard: for positive definite M, M[i,j] <= sqrt(M[i,i]*M[j,j])
# <= max(M[i,i], M[j,j]), so the largest entry is on the diagonal.
@test SymbolicAnalysis.analyze(unwrap(maximum(inv(LX)))).curvature ==
    SymbolicAnalysis.Convex


# `find_curvature`'s `*` branch read `constval(args[1])` and required a `Number`,
# so a constant *matrix* coefficient was rejected: `A*x - b` did not certify while
# the broadcast `A*x .- b` did, and `norm(A*x - b)` — the least-squares objective —
# came back UnknownCurvature. It now delegates to `mul_curvature`, which handles a
# constant in any position and of any shape.
@variables lx[1:3] lx2[1:2] lXm[1:3, 1:3] la lb
lA = [1.0 2.0 3.0; 4.0 5.0 6.0]
lbv = [1.0, 2.0]
@test SymbolicAnalysis.analyze(unwrap(lA * lx - lbv)).curvature == SymbolicAnalysis.Affine
@test SymbolicAnalysis.analyze(unwrap(norm(lA * lx - lbv))).curvature ==
    SymbolicAnalysis.Convex
@test SymbolicAnalysis.analyze(unwrap(sum(lA * lx - lbv))).curvature ==
    SymbolicAnalysis.Affine
@test SymbolicAnalysis.analyze(unwrap(maximum(lA * lx - lbv))).curvature ==
    SymbolicAnalysis.Convex

# A linear map preserves a *curvature* only when its coefficients share a sign;
# a mixed-sign constant matrix sums convex and concave terms. An affine argument
# stays affine whatever the coefficients are.
lApos = [1.0 2.0; 3.0 4.0]
lAmix = [1.0 -2.0; 3.0 -4.0]
@test SymbolicAnalysis.analyze(unwrap(lApos * exp.(lx2))).curvature ==
    SymbolicAnalysis.Convex
@test SymbolicAnalysis.analyze(unwrap((-lApos) * exp.(lx2))).curvature ==
    SymbolicAnalysis.Concave
@test SymbolicAnalysis.analyze(unwrap(lAmix * exp.(lx2))).curvature ==
    SymbolicAnalysis.UnknownCurvature
@test SymbolicAnalysis.analyze(unwrap(lAmix * lx2)).curvature == SymbolicAnalysis.Affine

# two non-constant factors are not DCP in any position
@test SymbolicAnalysis.analyze(unwrap(la * lb)).curvature ==
    SymbolicAnalysis.UnknownCurvature
@test SymbolicAnalysis.analyze(unwrap(lXm * lx)).curvature ==
    SymbolicAnalysis.UnknownCurvature
@test SymbolicAnalysis.analyze(unwrap(sum(lXm * lx))).curvature ==
    SymbolicAnalysis.UnknownCurvature

# Binary `-` is increasing in its first argument; only the second is decreasing.
# A single declared monotonicity is broadcast to every slot, so `Decreasing`
# flipped the first one too and `exp.(v) .- w` — Hessian diag(exp(vᵢ)) ⪰ 0 —
# certified as Concave. Scalar `a - b` never reaches the rule (Symbolics rewrites
# it to `a + (-1)*b`); the broadcast form does.
@variables mv[1:2] mw[1:2] ma
@test SymbolicAnalysis.analyze(unwrap(mv .- mw)).curvature == SymbolicAnalysis.Affine
@test SymbolicAnalysis.analyze(unwrap(exp.(mv) .- mw)).curvature ==
    SymbolicAnalysis.Convex
@test SymbolicAnalysis.analyze(unwrap(sum(exp.(mv) .- mw))).curvature ==
    SymbolicAnalysis.Convex
@test SymbolicAnalysis.analyze(unwrap(sum(log.(mv) .- mw))).curvature ==
    SymbolicAnalysis.Concave
# convex - convex is indefinite, and unary `-` stays decreasing
@test SymbolicAnalysis.analyze(unwrap(sum(exp.(mv) .- exp.(mw)))).curvature ==
    SymbolicAnalysis.UnknownCurvature
@test SymbolicAnalysis.analyze(unwrap(sum(.-(exp.(mv))))).curvature ==
    SymbolicAnalysis.Concave
@test SymbolicAnalysis.analyze(unwrap(-exp(ma))).curvature == SymbolicAnalysis.Concave

# `analyze` must return a curvature for every well-formed symbolic expression
# (#156, "crash instead of degrade"): each block below threw before.

# `hasdcprule(broadcast)` is unconditional while the rule forwarded to the rule
# table, so any broadcasted function without a table entry was a KeyError.
@variables bx[1:3]
@test SymbolicAnalysis.analyze(unwrap(2.0 .* bx)).curvature == SymbolicAnalysis.Affine
@test SymbolicAnalysis.analyze(unwrap(bx .* 2.0)).curvature == SymbolicAnalysis.Affine
@test SymbolicAnalysis.analyze(unwrap(sum(2.0 .* bx))).curvature == SymbolicAnalysis.Affine
@test SymbolicAnalysis.analyze(unwrap([1.0, 2.0, 3.0] .* bx)).curvature ==
    SymbolicAnalysis.Affine
@test SymbolicAnalysis.analyze(unwrap(2.0 .* exp.(bx))).curvature == SymbolicAnalysis.Convex
@test SymbolicAnalysis.analyze(unwrap(-2.0 .* exp.(bx))).curvature ==
    SymbolicAnalysis.Concave
@test SymbolicAnalysis.analyze(unwrap(2.0 .* log.(bx))).curvature ==
    SymbolicAnalysis.Concave
# elementwise `*` of two non-constant factors is bilinear, and a broadcast of a
# function with no rule has no curvature: both must degrade, not crash.
@test SymbolicAnalysis.analyze(unwrap(bx .* bx)).curvature ==
    SymbolicAnalysis.UnknownCurvature
@test SymbolicAnalysis.analyze(unwrap(sin.(bx))).curvature ==
    SymbolicAnalysis.UnknownCurvature
@test SymbolicAnalysis.analyze(unwrap(abs2.(bx))).curvature ==
    SymbolicAnalysis.UnknownCurvature
# `./` by a constant forwards to the `/` rule and is an affine rescaling — this
# assertion exists to cover that forwarding path, not to claim a degradation.
@test SymbolicAnalysis.analyze(unwrap(bx ./ 2.0)).curvature == SymbolicAnalysis.Affine
# A mixed-sign constant factor flips the curvature of some elements and not
# others, so only an affine argument survives it.
@test SymbolicAnalysis.analyze(unwrap([1.0, -2.0, 3.0] .* bx)).curvature ==
    SymbolicAnalysis.Affine
@test SymbolicAnalysis.analyze(unwrap([1.0, -2.0, 3.0] .* exp.(bx))).curvature ==
    SymbolicAnalysis.UnknownCurvature
# A broadcast power is elementwise, so the scalar power laws apply to it even
# though the base is an array.
@test SymbolicAnalysis.analyze(unwrap(bx .^ 2)).curvature == SymbolicAnalysis.Convex
@test SymbolicAnalysis.analyze(unwrap(sum(bx .^ 2))).curvature == SymbolicAnalysis.Convex

# Symbolics folds `-c/2` into `(-1//2)*c`, so a coefficient can be a `Rational`.
@variables c d
@test SymbolicAnalysis.analyze(unwrap(-c / 2)).curvature == SymbolicAnalysis.Affine
@test SymbolicAnalysis.analyze(unwrap((1 // 2) * exp(c))).curvature ==
    SymbolicAnalysis.Convex
@test SymbolicAnalysis.analyze(unwrap(-exp(c) / 2)).curvature == SymbolicAnalysis.Concave

# A constant-folded expression is a `BasicSymbolic` that is neither `issym` nor
# `iscall`, so the propagation walk never annotates it.
@test SymbolicAnalysis.analyze(unwrap(c - c)).curvature == SymbolicAnalysis.Affine
@test SymbolicAnalysis.analyze(unwrap(0 * c)).curvature == SymbolicAnalysis.Affine
@test SymbolicAnalysis.analyze(unwrap(c^0)).curvature == SymbolicAnalysis.Affine
@test SymbolicAnalysis.analyze(unwrap(Num(3.0))).curvature == SymbolicAnalysis.Affine
@test SymbolicAnalysis.analyze(unwrap(exp(c) - exp(c))).curvature ==
    SymbolicAnalysis.Affine

# A symbolic exponent is not a number, so the power laws cannot branch on it.
# `c^g` with a constant base is `exp(g*log(c))`: convex, monotone with `log(c)`.
@test SymbolicAnalysis.analyze(unwrap(2^c)).curvature == SymbolicAnalysis.Convex
@test SymbolicAnalysis.analyze(unwrap(2^c)).sign == SymbolicAnalysis.Positive
@test SymbolicAnalysis.analyze(unwrap(2.0^(c + 1))).curvature == SymbolicAnalysis.Convex
@test SymbolicAnalysis.analyze(unwrap(2.0^exp(c))).curvature == SymbolicAnalysis.Convex
@test SymbolicAnalysis.analyze(unwrap(0.5^(-exp(c)))).curvature == SymbolicAnalysis.Convex
@test SymbolicAnalysis.analyze(unwrap(c^d)).curvature == SymbolicAnalysis.UnknownCurvature
@test SymbolicAnalysis.analyze(unwrap(exp(c)^d)).curvature ==
    SymbolicAnalysis.UnknownCurvature
# `2^g` composes only with an increasing slot and `0.5^g` only with a decreasing
# one; a negative base is not real-valued for every exponent.
@test SymbolicAnalysis.analyze(unwrap(2.0^(-exp(c)))).curvature ==
    SymbolicAnalysis.UnknownCurvature
@test SymbolicAnalysis.analyze(unwrap(0.5^exp(c))).curvature ==
    SymbolicAnalysis.UnknownCurvature
@test SymbolicAnalysis.analyze(unwrap((-2.0)^c)).curvature ==
    SymbolicAnalysis.UnknownCurvature

# `abs`/`conj`/`real`/`imag` are registered on ℂ, which DomainSets cannot compare
# against a `HalfLine`, so the package's own positivity annotation broke `abs`.
cpos = setmetadata(
    c, SymbolicAnalysis.VarDomain, Symbolics.DomainSets.HalfLine{Number, :open}()
)
creal = setmetadata(c, SymbolicAnalysis.VarDomain, Symbolics.DomainSets.RealLine())
@test SymbolicAnalysis.analyze(unwrap(abs(cpos))).curvature == SymbolicAnalysis.Convex
@test SymbolicAnalysis.analyze(unwrap(abs(cpos))).sign == SymbolicAnalysis.Positive
@test SymbolicAnalysis.analyze(unwrap(abs(creal))).curvature == SymbolicAnalysis.Convex
@test SymbolicAnalysis.analyze(unwrap(exp(abs(cpos)))).curvature == SymbolicAnalysis.Convex
@test SymbolicAnalysis.analyze(unwrap(log(cpos))).curvature == SymbolicAnalysis.Concave
# A genuine domain mismatch still selects no rule — `log` is not concave (or
# defined) on all of the reals — but it degrades instead of throwing.
@test SymbolicAnalysis.analyze(unwrap(log(creal))).curvature ==
    SymbolicAnalysis.UnknownCurvature

# The `logdet` registration asserted `symtype(X) <: Matrix{Num}`, which no
# derived matrix expression satisfies, so `logdet` of one threw at build time.
@variables Xl[1:3, 1:3]
Al = rand(3, 3)
@test SymbolicAnalysis.analyze(unwrap(logdet(2 * Xl))).curvature == SymbolicAnalysis.Concave
@test SymbolicAnalysis.analyze(unwrap(logdet(Al * Xl))).curvature ==
    SymbolicAnalysis.Concave
@test SymbolicAnalysis.analyze(unwrap(logdet(Al * Xl * Al'))).curvature ==
    SymbolicAnalysis.Concave
@test SymbolicAnalysis.analyze(unwrap(logdet(Xl + Xl'))).curvature ==
    SymbolicAnalysis.Concave
@test SymbolicAnalysis.analyze(unwrap(-logdet(2 * Xl))).curvature == SymbolicAnalysis.Convex
