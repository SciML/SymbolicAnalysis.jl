using SymbolicAnalysis, BenchmarkTools
using Symbolics, Manifolds, LinearAlgebra
using Symbolics: @variables

const SUITE = BenchmarkGroup()

# =============================================================================
# DCP curvature analysis (the core single-pass propagation engine)
# =============================================================================

SUITE["analyze"] = BenchmarkGroup()

# Small expression
@variables x
ex_small = Symbolics.unwrap(exp(x) + x^2 + log(1 + x^2))
SUITE["analyze"]["small"] = @benchmarkable SymbolicAnalysis.analyze($ex_small)

# Medium expression — same scaling workload the allocation regression test uses
vars = Symbolics.variables(:q, 1:200)
ex_large = Symbolics.unwrap(sum(exp(v) + v^2 for v in vars))
SUITE["analyze"]["large_200vars"] = @benchmarkable SymbolicAnalysis.analyze($ex_large)

# Geometric (gDCP) analysis on a manifold — Brascamp–Lieb-style objective
M = SymmetricPositiveDefinite(5)
@variables X[1:5, 1:5] A[1:5, 1:5]
ex_g = Symbolics.unwrap(
    logdet(SymbolicAnalysis.conjugation(X, A)) - logdet(X)
)
ex_g = SymbolicAnalysis.propagate_sign(ex_g)
SUITE["analyze"]["gdcp_manifold"] = @benchmarkable analyze($ex_g, $M)

# =============================================================================
# Sign propagation
# =============================================================================

SUITE["propagate"] = BenchmarkGroup()

SUITE["propagate"]["propagate_sign"] = @benchmarkable SymbolicAnalysis.propagate_sign(
    $ex_g
)
