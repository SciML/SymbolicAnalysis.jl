module SymbolicAnalysis

using DomainSets
using LinearAlgebra
using LogExpFunctions
using StatsBase
using Distributions
using DSP, DataStructures

using Symbolics
import Symbolics: Symbolic, issym, Term
using SymbolicUtils: iscall
using SymbolicUtils.Rewriters
SymbolicUtils.inspect_metadata[] = true

struct VarDomain end

include("rules.jl")
include("atoms.jl")
include("gdcp/gdcp_rules.jl")
include("gdcp/spd.jl")
include("gdcp/lorentz.jl")
include("canon.jl")

struct AnalysisResult
    curvature::SymbolicAnalysis.Curvature
    sign::SymbolicAnalysis.Sign
    gcurvature::Union{SymbolicAnalysis.GCurvature,Nothing}
end

"""
    analyze(ex)
    analyze(ex, M)

Analyze the expression `ex` and return the curvature and sign of the expression. If a manifold `M` from [Manifolds.jl](https://juliamanifolds.github.io/Manifolds.jl/stable/) is provided, also return the geodesic curvature of the expression.
Currently supports the `SymmetricPositiveDefinite` and `Lorentz` manifolds.

The returned `AnalysisResult` contains the following fields:
- `curvature::SymbolicAnalysis.Curvature`: The curvature of the expression.
- `sign::SymbolicAnalysis.Sign`: The sign of the expression.
- `gcurvature::Union{SymbolicAnalysis.GCurvature,Nothing}`: The geodesic curvature of the expression if `M` is provided. Otherwise, `nothing`.
"""
function analyze(ex, M::Union{AbstractManifold,Nothing} = nothing)
    ex = unwrap(ex)
    ex = canonize(ex)
    ex = propagate_sign(ex)
    ex = propagate_curvature(ex)
    if isnothing(M)
        return AnalysisResult(getcurvature(ex), getsign(ex), nothing)
    else
        @assert M isa SymmetricPositiveDefinite || M isa Lorentz "Only SymmetricPositiveDefinite and Lorentz manifolds are currently supported"
        ex = propagate_gcurvature(ex, M)
        return AnalysisResult(getcurvature(ex), getsign(ex), getgcurvature(ex))
    end
end

export analyze

end
