module SymbolicAnalysis

using DSP: conv
using Distributions: Normal, logcdf
import DomainSets
using DomainSets: HalfLine, RealLine, ℂ, ℤ
using IntervalSets: Domain, Interval
import LinearAlgebra
using LinearAlgebra: Diagonal, I, Symmetric, diag, diagm, dot, eigmax, eigmin, eigvals,
    isposdef, issymmetric, logdet, norm, opnorm, tr, tril, triu
import LogExpFunctions
using LogExpFunctions: xexpx, xlogx
import Manifolds
using Manifolds: Lorentz, SymmetricPositiveDefinite
using ManifoldsBase: AbstractManifold
using PrecompileTools: @compile_workload, @setup_workload
using SciMLPublic: @public
import StatsBase
using StatsBase: kldivergence
import Symbolics
using Symbolics: @variables, Num
import SymbolicUtils
using SymbolicUtils: @rule, BasicSymbolic, getmetadata, hasmetadata, iscall, issym,
    setmetadata, unwrap
using SymbolicUtils.Rewriters: Postwalk
using TermInterface: arguments, operation

# Symbolics v7 / SymbolicUtils v4 removed the `Symbolic` abstract type: every
# symbolic — scalar or array — is now a `BasicSymbolic{SymReal}`.
const Symbolic = BasicSymbolic

# The scalar-symbolic type Symbolics uses when dispatching `in(::symbolic, ::Domain)`.
# Matching it exactly lets the `in(::_, ::CustomDomain)` disambiguators below stay
# strictly more specific than Symbolics' `IntervalSets.Domain` method (which is
# keyed on `BasicSymbolic{SymReal}`).
const InDomainSymbolic = BasicSymbolic{SymbolicUtils.SymReal}

"""
    VarDomain

Metadata key used to attach a `DomainSets.Domain`
to a symbolic variable. DCP rule selection uses this metadata when a rule has
domain restrictions.

# Fields

`VarDomain` is a marker type and has no fields. Use the type itself as the
metadata key with `SymbolicUtils.setmetadata`.

# Examples

```jldoctest
julia> using SymbolicAnalysis, Symbolics, SymbolicUtils, DomainSets

julia> @variables x;

julia> x = setmetadata(x, SymbolicAnalysis.VarDomain, HalfLine{Number, :open}());

julia> hasmetadata(x, SymbolicAnalysis.VarDomain)
true
```
"""
struct VarDomain end

include("rules.jl")
include("atoms.jl")
include("gdcp/gdcp_rules.jl")
include("gdcp/spd.jl")
include("gdcp/lorentz.jl")
include("canon.jl")
include("fold.jl")

"""
    AnalysisResult

Result returned by [`analyze`](@ref). It contains the Euclidean sign and
curvature inferred from a symbolic expression and, when a supported manifold
is supplied, its geodesic curvature.

# Fields

- `curvature::Curvature`: Euclidean curvature of the expression.
- `sign::Sign`: inferred sign of the expression.
- `gcurvature::Union{GCurvature, Nothing}`: geodesic curvature when manifold
  analysis was requested, otherwise `nothing`.

The result is immutable. Read its fields directly rather than relying on the
internal metadata propagation functions.

# Examples

```jldoctest
julia> using SymbolicAnalysis, Symbolics

julia> @variables x;

julia> result = analyze(exp(x));

julia> (result.curvature, result.sign, result.gcurvature)
(SymbolicAnalysis.Convex, SymbolicAnalysis.Positive, nothing)
```
"""
struct AnalysisResult
    curvature::SymbolicAnalysis.Curvature
    sign::SymbolicAnalysis.Sign
    gcurvature::Union{SymbolicAnalysis.GCurvature, Nothing}
end

@public VarDomain, AnalysisResult
@public Sign, Curvature, Monotonicity, GCurvature, GMonotonicity
@public add_dcprule, add_gdcprule

"""
    analyze(ex, M = nothing) -> AnalysisResult

Analyze the symbolic expression `ex` and return its Euclidean curvature and sign.
When `M` is supplied, also determine the geodesic curvature on that manifold.

# Arguments

- `ex`: Symbolics expression to analyze.
- `M::Union{AbstractManifold, Nothing} = nothing`: optional manifold for geodesic
  curvature analysis. `SymmetricPositiveDefinite` and `Lorentz` manifolds are
  supported.

# Returns

An `AnalysisResult` with the fields:

- `curvature::SymbolicAnalysis.Curvature`: Euclidean curvature of `ex`.
- `sign::SymbolicAnalysis.Sign`: inferred sign of `ex`.
- `gcurvature::Union{SymbolicAnalysis.GCurvature, Nothing}`: geodesic curvature
  when `M` is supplied, or `nothing` otherwise.

# Throws

- `AssertionError`: if `M` is not a supported manifold.

# Examples

```jldoctest
julia> using SymbolicAnalysis, Symbolics

julia> @variables x;

julia> result = analyze(exp(x));

julia> result.curvature == SymbolicAnalysis.Convex
true

julia> result.gcurvature === nothing
true
```
"""
function analyze(ex, M::Union{AbstractManifold, Nothing} = nothing)
    ex = unwrap(ex)
    ex = canonize(ex)
    if !isnothing(M)
        @assert M isa SymmetricPositiveDefinite || M isa Lorentz "Only SymmetricPositiveDefinite and Lorentz manifolds are currently supported"
    end
    # One read-only memoized fold replaces the Postwalk rebuild passes. The
    # public `propagate_*` functions still annotate trees for callers that
    # read metadata; `analyze` only needs the root properties.
    props = analyze_fold(ex, M)
    # `props.sign` is ordinarily a `Sign`. Atoms that store callable
    # placeholders in the rule table (e.g. `perspective`) leave a `Function`
    # here; the `::Sign` assertion then fails exactly as `getsign` did on the
    # metadata-annotated tree under typed-rule-tables.
    return AnalysisResult(props.curvature, props.sign::Sign, props.gcurvature)
end

export analyze

@setup_workload begin
    @compile_workload begin
        @variables x y
        y_with_domain = setmetadata(
            y, VarDomain, DomainSets.HalfLine{Number, :open}()
        )

        ex1 = exp(y_with_domain) - log(y_with_domain) |> unwrap
        analyze(ex1)

        ex2 = abs(x)^2 + abs(x)^3 |> unwrap
        analyze(ex2)

        ex3 = 2 * abs(x) - 1 |> unwrap
        analyze(ex3)

        @variables z[1:3]
        ex4 = exp.(z) |> unwrap
        analyze(ex4)
    end
end

end
