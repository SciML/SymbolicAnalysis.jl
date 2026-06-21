using SymbolicAnalysis, Aqua, JET, SciMLTesting

const SA = SymbolicAnalysis

# Functions this package deliberately extends with symbolic-analysis methods.
# Registering DCP/GDCP methods on these non-owned functions (and the Symbolics
# traversal/symtype hooks) is the core purpose of the package, so they are
# intentional and declared to Aqua's piracy check via `treat_as_own`.
const SYMBOLIC_OWN = Any[
    Base.log, SA.LinearAlgebra.tr, SA.LinearAlgebra.inv, SA.LinearAlgebra.sqrt,
    SA.LinearAlgebra.logdet, SA.Manifolds.distance, SA.LogExpFunctions.xlogx,
    SA.Symbolics.arguments, SA.Symbolics.hasmetadata, SA.Symbolics.promote_symtype,
    # Symbolics v7 support: matrix `*` is re-wrapped to an `Arr` and `sqrt`/atom
    # registrations contribute `promote_shape` methods on non-owned functions.
    # These are intentional extensions of the symbolic machinery (the package's
    # core purpose), so declare them as own.
    Base.:*, SA.SymbolicUtils.promote_shape,
]

# The scalar/array `@register_symbolic`/`@register_array_symbolic` registrations
# for the package's own atoms generate overlapping signatures (e.g. a scalar and
# an array method of the same atom), which Aqua reports as internal ambiguities.
# These are all between SymbolicAnalysis's own methods (no cross-package
# ambiguity remains) and are benign; concrete calls dispatch unambiguously.
# Excluding only these atoms keeps the check live for everything else.
const ATOM_AMBIGUITIES = Any[
    SA.affine_map, SA.sdivergence, SA.lorentz_least_squares, SA.conjugation,
    SA.log_quad_form, SA.lorentz_homogeneous_quadratic, SA.hadamard_product,
    SA.lorentz_homogeneous_diagonal, SA.lorentz_transform, SA.quad_over_lin,
]

run_qa(
    SymbolicAnalysis;
    Aqua = Aqua,
    JET = JET,
    jet = true,
    aqua_kwargs = (;
        ambiguities = (; exclude = ATOM_AMBIGUITIES),
        piracies = (; treat_as_own = SYMBOLIC_OWN),
    ),
    jet_kwargs = (; target_modules = (SymbolicAnalysis,), mode = :typo),
)
