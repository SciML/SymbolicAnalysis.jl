using SymbolicAnalysis, Aqua, JET, SciMLTesting

const SA = SymbolicAnalysis

# SymbolicAnalysis's entire purpose is to teach the symbolic-analysis machinery
# about existing functions: it registers DCP/gDCP curvature methods (and the
# Symbolics traversal/symtype/shape hooks they need) on functions owned by Base,
# LinearAlgebra, Symbolics/SymbolicUtils, Manifolds and LogExpFunctions. Those are
# intentional, by-design extensions of non-owned functions, so they are declared to
# Aqua's piracy check via `treat_as_own`. This is the only Aqua exception the
# package needs; ambiguities, stale-deps, etc. all pass cleanly.
const SYMBOLIC_OWN = Any[
    Base.:*, Base.log, Base.sqrt,
    SA.LinearAlgebra.inv, SA.LinearAlgebra.logdet,
    SA.Symbolics.arguments, SA.Symbolics.hasmetadata, SA.Symbolics.promote_symtype,
    SA.SymbolicUtils.promote_shape,
    SA.Manifolds.distance, SA.LogExpFunctions.xlogx, SA.LogExpFunctions.logsumexp,
]

run_qa(
    SymbolicAnalysis;
    Aqua = Aqua,
    JET = JET,
    jet = true,
    aqua_kwargs = (; piracies = (; treat_as_own = SYMBOLIC_OWN)),
    jet_kwargs = (; target_modules = (SymbolicAnalysis,), mode = :typo),
)
