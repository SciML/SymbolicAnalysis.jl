"""
    Conic Form Generation

Transforms DCP-verified symbolic expressions into conic formulations
suitable for consumption by MathOptInterface (MOI) solvers.

The key idea: every DCP atom has a corresponding MOI cone. When we walk
the expression tree bottom-up, each atom call can be replaced by an
epigraph variable `t` plus a cone constraint linking `t` to the atom's arguments.
The result is a linear objective over epigraph variables subject to cone constraints.
"""

using MathOptInterface
const _MOI = MathOptInterface

"""
    ConeConstraint

A single conic constraint: `func ∈ cone`.

# Fields
- `func_vars::Vector{Symbol}` — variable names involved
- `func_coeffs::Vector{Float64}` — coefficients for each variable
- `func_constant::Float64` — constant offset
- `cone` — MOI cone type (e.g., `MOI.ExponentialCone`, `MOI.SecondOrderCone`)
- `description::String` — human-readable description
"""
struct ConeConstraint
    func_vars::Vector{Symbol}
    func_coeffs::Vector{Float64}
    func_constant::Float64
    cone::Any  # MOI.AbstractSet type
    description::String
end

"""
    ConicFormulation

The result of converting a DCP expression to conic form.

# Fields
- `objective_var::Symbol` — the top-level epigraph variable (minimize this for convex, maximize for concave)
- `objective_sense::Symbol` — `:minimize` or `:maximize`
- `constraints::Vector{ConeConstraint}` — cone constraints
- `epigraph_map::Dict{Symbol, Any}` — maps epigraph variable names to the expressions they represent
- `variables::Set{Symbol}` — all decision variables (original + epigraph)
- `original_variables::Set{Symbol}` — only the original (user) variables
"""
struct ConicFormulation
    objective_var::Symbol
    objective_sense::Symbol
    constraints::Vector{ConeConstraint}
    epigraph_map::Dict{Symbol, Any}
    variables::Set{Symbol}
    original_variables::Set{Symbol}
end

# Counter for generating unique epigraph variable names
const _epi_counter = Ref(0)

function _reset_epi_counter!()
    _epi_counter[] = 0
end

function _new_epi_var()
    _epi_counter[] += 1
    return Symbol("_t$(_epi_counter[])")
end

"""
    to_conic_form(ex)

Convert a DCP-verified symbolic expression to a `ConicFormulation`.

The expression `ex` should have already been analyzed via `analyze()` to confirm
DCP compliance. This function walks the expression tree bottom-up, introducing
epigraph variables and cone constraints for each atom.

# Returns
A `ConicFormulation` with:
- A linear objective over epigraph variables
- Cone constraints encoding each atom's epigraph
"""
function to_conic_form(ex)
    _reset_epi_counter!()
    ex = unwrap(ex)

    # First, analyze to get curvature
    analyzed = canonize(ex)
    analyzed = propagate_sign(analyzed)
    analyzed = propagate_curvature(analyzed)

    curv = getcurvature(analyzed)
    sense = if curv == Convex
        :minimize
    elseif curv == Concave
        :maximize
    else
        :minimize  # Affine can be either
    end

    original_vars = Set{Symbol}()
    _collect_variables!(analyzed, original_vars)

    constraints = ConeConstraint[]
    epigraph_map = Dict{Symbol, Any}()
    variables = copy(original_vars)

    obj_var = _process_node!(analyzed, constraints, epigraph_map, variables)

    return ConicFormulation(
        obj_var,
        sense,
        constraints,
        epigraph_map,
        variables,
        original_vars
    )
end

"""
    _collect_variables!(ex, vars)

Collect all symbolic variable names from an expression.
"""
function _collect_variables!(ex, vars::Set{Symbol})
    if issym(ex)
        push!(vars, Symbol(ex))
    elseif iscall(ex)
        for arg in arguments(ex)
            _collect_variables!(arg, vars)
        end
    end
end

"""
    _process_node!(ex, constraints, epigraph_map, variables)

Recursively process an expression node, emitting cone constraints and
returning the symbol for the epigraph variable that represents this node.
"""
function _process_node!(ex, constraints, epigraph_map, variables)
    # Base case: a symbolic variable
    if issym(ex)
        return Symbol(ex)
    end

    # Base case: a number
    if ex isa Number
        # Create an epigraph variable fixed to this constant
        t = _new_epi_var()
        push!(variables, t)
        epigraph_map[t] = ex
        # Add equality constraint: t == constant
        push!(constraints, ConeConstraint(
            [t], [1.0], -Float64(ex),
            _MOI.EqualTo(0.0),
            "constant: $t == $ex"
        ))
        return t
    end

    if !iscall(ex)
        # Wrapped Num or similar
        return _process_node!(unwrap(ex), constraints, epigraph_map, variables)
    end

    f = operation(ex)
    args = arguments(ex)

    # Handle addition: sum of subexpressions
    if Symbol(f) == :+
        child_vars = Symbol[]
        child_coeffs = Float64[]
        constant = 0.0
        for arg in args
            if arg isa Number
                constant += Float64(arg)
            else
                child = _process_node!(arg, constraints, epigraph_map, variables)
                push!(child_vars, child)
                push!(child_coeffs, 1.0)
            end
        end
        t = _new_epi_var()
        push!(variables, t)
        epigraph_map[t] = ex

        # t == sum of children + constant
        all_vars = vcat([t], child_vars)
        all_coeffs = vcat([1.0], [-c for c in child_coeffs])
        push!(constraints, ConeConstraint(
            all_vars, all_coeffs, -constant,
            _MOI.EqualTo(0.0),
            "sum: $t == $(join(child_vars, " + ")) + $constant"
        ))
        return t
    end

    # Handle multiplication by constant
    if Symbol(f) == :*
        # Find constant and non-constant parts
        constant = 1.0
        non_const = nothing
        for arg in args
            if arg isa Number
                constant *= Float64(arg)
            else
                non_const = arg
            end
        end

        if non_const !== nothing
            child = _process_node!(non_const, constraints, epigraph_map, variables)
            t = _new_epi_var()
            push!(variables, t)
            epigraph_map[t] = ex

            # t == constant * child
            push!(constraints, ConeConstraint(
                [t, child], [1.0, -constant], 0.0,
                _MOI.EqualTo(0.0),
                "scale: $t == $constant * $child"
            ))
            return t
        else
            # Pure constant
            t = _new_epi_var()
            push!(variables, t)
            epigraph_map[t] = constant
            push!(constraints, ConeConstraint(
                [t], [1.0], -constant,
                _MOI.EqualTo(0.0),
                "constant: $t == $constant"
            ))
            return t
        end
    end

    # Handle DCP atoms with cone annotations
    if hasdcprule(f)
        child_vars = Symbol[]
        for arg in args
            if arg isa Number
                child = _process_node!(arg, constraints, epigraph_map, variables)
                push!(child_vars, child)
            else
                child = _process_node!(arg, constraints, epigraph_map, variables)
                push!(child_vars, child)
            end
        end

        # Look up the cone for this atom
        rule, _ = dcprule(f, args...)
        cone = rule.cone
        curv = rule.curvature

        t = _new_epi_var()
        push!(variables, t)
        epigraph_map[t] = ex

        # Emit the cone constraint
        _emit_atom_constraint!(f, t, child_vars, cone, curv, constraints)

        return t
    end

    # Fallback: treat as opaque
    t = _new_epi_var()
    push!(variables, t)
    epigraph_map[t] = ex
    return t
end

"""
    _emit_atom_constraint!(f, t, child_vars, cone, curvature, constraints)

Emit the appropriate cone constraint for atom `f` with epigraph variable `t`
and argument variables `child_vars`.

For a convex atom f(x), the epigraph is: {(t, x) : f(x) ≤ t}
For a concave atom f(x), the hypograph is: {(t, x) : f(x) ≥ t}
"""
function _emit_atom_constraint!(f, t, child_vars, cone, curvature, constraints)
    fname = string(nameof(f))

    if cone === nothing || cone == _MOI.Reals
        # Linear or no specific cone — just record the relationship
        push!(constraints, ConeConstraint(
            vcat([t], child_vars),
            vcat([1.0], [-1.0 for _ in child_vars]),
            0.0,
            _MOI.Zeros(1 + length(child_vars)),
            "$fname: linear relationship"
        ))
        return
    end

    # Dispatch on specific atoms for proper conic reformulation
    if f === exp
        # exp(x) ≤ t  ⟺  (x, 1, t) ∈ ExponentialCone
        # MOI.ExponentialCone: (x, y, z) such that y * exp(x/y) ≤ z, y > 0
        @assert length(child_vars) == 1
        push!(constraints, ConeConstraint(
            [child_vars[1], t],  # x, t
            [1.0, 1.0],
            0.0,
            _MOI.ExponentialCone(),
            "$fname: ($(child_vars[1]), 1, $t) ∈ ExponentialCone"
        ))

    elseif f === log
        # log(x) ≥ t  ⟺  (t, 1, x) ∈ ExponentialCone
        @assert length(child_vars) == 1
        push!(constraints, ConeConstraint(
            [t, child_vars[1]],
            [1.0, 1.0],
            0.0,
            _MOI.ExponentialCone(),
            "$fname: ($t, 1, $(child_vars[1])) ∈ ExponentialCone"
        ))

    elseif f === abs
        # |x| ≤ t  ⟺  (t, x) ∈ NormOneCone(2)
        @assert length(child_vars) == 1
        push!(constraints, ConeConstraint(
            [t, child_vars[1]],
            [1.0, 1.0],
            0.0,
            _MOI.NormOneCone(2),
            "$fname: ($t, $(child_vars[1])) ∈ NormOneCone(2)"
        ))

    elseif f === norm
        # ‖x‖ ≤ t  ⟺  (t, x...) ∈ SecondOrderCone
        push!(constraints, ConeConstraint(
            vcat([t], child_vars),
            ones(1 + length(child_vars)),
            0.0,
            _MOI.SecondOrderCone(1 + length(child_vars)),
            "$fname: ($t, $(join(child_vars, ", "))) ∈ SOC"
        ))

    elseif f === sqrt
        # sqrt(x) ≥ t  ⟺  (t, 1, x) ∈ RotatedSecondOrderCone(3)
        # RSOC: t₁ * t₂ ≥ ‖x‖², t₁,t₂ ≥ 0
        @assert length(child_vars) == 1
        push!(constraints, ConeConstraint(
            [t, child_vars[1]],
            [1.0, 1.0],
            0.0,
            _MOI.RotatedSecondOrderCone(3),
            "$fname: ($t, 1, $(child_vars[1])) ∈ RSOC"
        ))

    elseif f === inv
        # inv(x) ≤ t, x > 0  ⟺  (t, x, 1) ∈ RotatedSecondOrderCone(3)
        @assert length(child_vars) == 1
        push!(constraints, ConeConstraint(
            [t, child_vars[1]],
            [1.0, 1.0],
            0.0,
            _MOI.RotatedSecondOrderCone(3),
            "$fname: ($t, $(child_vars[1]), 1) ∈ RSOC"
        ))

    elseif f === quad_over_lin
        # x²/y ≤ t  ⟺  (y, t, x) ∈ RotatedSecondOrderCone(3)
        @assert length(child_vars) == 2
        push!(constraints, ConeConstraint(
            [child_vars[2], t, child_vars[1]],
            [1.0, 1.0, 1.0],
            0.0,
            _MOI.RotatedSecondOrderCone(3),
            "$fname: ($(child_vars[2]), $t, $(child_vars[1])) ∈ RSOC"
        ))

    elseif f === rel_entr
        # x*log(x/y) ≤ t  ⟺  (-t, x, y) ∈ RelativeEntropyCone(3)
        @assert length(child_vars) == 2
        push!(constraints, ConeConstraint(
            [t, child_vars[1], child_vars[2]],
            [1.0, 1.0, 1.0],
            0.0,
            _MOI.RelativeEntropyCone(3),
            "$fname: ($t, $(child_vars[1]), $(child_vars[2])) ∈ RelativeEntropyCone"
        ))

    else
        # Generic: record the cone type without a specific reformulation
        sense_str = curvature == Convex ? "≤" : curvature == Concave ? "≥" : "=="
        push!(constraints, ConeConstraint(
            vcat([t], child_vars),
            ones(1 + length(child_vars)),
            0.0,
            cone isa DataType ? cone : typeof(cone),
            "$fname: $t $sense_str $fname($(join(child_vars, ", "))) via $(cone)"
        ))
    end
end

"""
    list_cone_annotations()

Return a list of all registered DCP and DGCP atoms with their MOI cone annotations.
"""
function list_cone_annotations()
    result = []
    for (f, rule) in dcprules_dict
        if rule isa Vector
            for r in rule
                push!(result, (atom = nameof(f), type = :DCP, cone = r.cone, curvature = r.curvature))
            end
        else
            push!(result, (atom = nameof(f), type = :DCP, cone = rule.cone, curvature = rule.curvature))
        end
    end
    for (f, rule) in gdcprules_dict
        push!(result, (atom = nameof(f), type = :GDCP, cone = rule.cone, gcurvature = rule.gcurvature))
    end
    return result
end

export to_conic_form, ConicFormulation, ConeConstraint, list_cone_annotations
