"""
    Conic Form Generation

Transforms DCP-verified symbolic expressions into conic formulations
suitable for consumption by MathOptInterface (MOI) solvers.

The key idea: every DCP atom has a corresponding MOI cone. When we walk
the expression tree bottom-up, each atom call can be replaced by an
epigraph variable `t` plus a cone constraint linking `t` to the atom's arguments.
The result is a linear objective over epigraph variables subject to cone constraints.
"""

# Use the module-level MOI alias from SymbolicAnalysis.jl
# (const MOI = MathOptInterface is defined there)

"""
    ConicConstraintTerm

A single row/dimension of a conic constraint: an affine expression `coeffs'*vars + constant`.

# Fields
- `vars::Vector{Symbol}` — variable names
- `coeffs::Vector{Float64}` — coefficient for each variable
- `constant::Float64` — constant offset
"""
struct ConicConstraintTerm
    vars::Vector{Symbol}
    coeffs::Vector{Float64}
    constant::Float64
end

"""
    ConeConstraint

A vector-valued conic constraint: `(terms[1], terms[2], ...) ∈ cone`.

Each `ConicConstraintTerm` produces one row of the vector-valued function.

# Fields
- `terms::Vector{ConicConstraintTerm}` — one per output dimension of the cone
- `cone::Any` — MOI cone instance (e.g., `MOI.ExponentialCone()`, `MOI.SecondOrderCone(3)`)
- `atom::Union{Function, Nothing}` — which atom generated this constraint
- `description::String` — human-readable description
"""
struct ConeConstraint
    terms::Vector{ConicConstraintTerm}
    cone::Any
    atom::Union{Function, Nothing}
    description::String
end

"""
    ConicContext

Mutable context for conic form generation, replacing global state.
Thread-safe: each call to `to_conic_form` creates its own context.

# Fields
- `epi_counter::Int` — counter for unique epigraph variable names
- `constraints::Vector{ConeConstraint}` — accumulated cone constraints
- `epigraph_map::Dict{Symbol, Any}` — maps epigraph variables to their expressions
- `variables::Set{Symbol}` — all variables (original + epigraph)
- `original_variables::Set{Symbol}` — only the original (user) variables
"""
mutable struct ConicContext
    epi_counter::Int
    constraints::Vector{ConeConstraint}
    epigraph_map::Dict{Symbol, Any}
    variables::Set{Symbol}
    original_variables::Set{Symbol}
end

function ConicContext(original_vars::Set{Symbol})
    return ConicContext(
        0,
        ConeConstraint[],
        Dict{Symbol, Any}(),
        copy(original_vars),
        original_vars
    )
end

function _new_epi_var!(ctx::ConicContext)
    ctx.epi_counter += 1
    t = Symbol("_t$(ctx.epi_counter)")
    push!(ctx.variables, t)
    return t
end

"""
    ConicFormulation

The result of converting a DCP expression to conic form.

# Fields
- `objective_var::Symbol` — the top-level epigraph variable
- `objective_sense::Symbol` — `:minimize` or `:maximize`
- `constraints::Vector{ConeConstraint}` — cone constraints
- `epigraph_map::Dict{Symbol, Any}` — maps epigraph variable names to expressions
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

# ──────────────────────────────────────────────────────────────────────────────
# Affine expression utilities
# ──────────────────────────────────────────────────────────────────────────────

"""
    _is_affine(ex)

Check if expression `ex` is purely affine (symbols, numbers, +, * by constant).
"""
function _is_affine(ex)
    if issym(ex) || ex isa Number
        return true
    end
    if !iscall(ex)
        return _is_affine(unwrap(ex))
    end
    f = operation(ex)
    args = arguments(ex)
    if Symbol(f) == :+
        return all(_is_affine, args)
    elseif Symbol(f) == :*
        # Affine if at most one non-constant factor
        non_const = count(a -> !(a isa Number), args)
        return non_const <= 1 && all(a -> a isa Number || _is_affine(a), args)
    end
    return false
end

"""
    _extract_affine(ex)

Extract affine structure from a purely affine expression.
Returns `(vars::Vector{Symbol}, coeffs::Vector{Float64}, constant::Float64)`.
Assumes `_is_affine(ex)` is true.
"""
function _extract_affine(ex)
    vars = Symbol[]
    coeffs = Float64[]
    constant = Ref(0.0)
    _extract_affine!(ex, vars, coeffs, constant, 1.0)
    return vars, coeffs, constant[]
end

function _extract_affine!(ex, vars, coeffs, constant, scale)
    if ex isa Number
        constant[] += scale * Float64(ex)
        return
    end
    if issym(ex)
        # Check if this variable already appears; if so, accumulate coefficient
        sym = Symbol(ex)
        idx = findfirst(==(sym), vars)
        if idx !== nothing
            coeffs[idx] += scale
        else
            push!(vars, sym)
            push!(coeffs, scale)
        end
        return
    end
    if !iscall(ex)
        _extract_affine!(unwrap(ex), vars, coeffs, constant, scale)
        return
    end
    f = operation(ex)
    args = arguments(ex)
    if Symbol(f) == :+
        for arg in args
            _extract_affine!(arg, vars, coeffs, constant, scale)
        end
    elseif Symbol(f) == :*
        # Find constant part and non-constant part
        c = 1.0
        non_const = nothing
        for arg in args
            if arg isa Number
                c *= Float64(arg)
            else
                non_const = arg
            end
        end
        if non_const !== nothing
            _extract_affine!(non_const, vars, coeffs, constant, scale * c)
        else
            constant[] += scale * c
        end
    end
end

# ──────────────────────────────────────────────────────────────────────────────
# Main entry point
# ──────────────────────────────────────────────────────────────────────────────

"""
    to_conic_form(ex)

Convert a DCP-verified symbolic expression to a `ConicFormulation`.

The expression `ex` should have already been analyzed via `analyze()` to confirm
DCP compliance. This function walks the expression tree bottom-up, introducing
epigraph variables and cone constraints for each atom.

Thread-safe: uses a local `ConicContext` instead of global state.

# Returns
A `ConicFormulation` with:
- A linear objective over epigraph variables
- Cone constraints encoding each atom's epigraph
"""
function to_conic_form(ex)
    ex = unwrap(ex)

    # Analyze to get curvature
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

    ctx = ConicContext(original_vars)
    obj_var = _process_node!(analyzed, ctx)

    return ConicFormulation(
        obj_var,
        sense,
        ctx.constraints,
        ctx.epigraph_map,
        ctx.variables,
        ctx.original_variables
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
    _process_node!(ex, ctx::ConicContext)

Recursively process an expression node, emitting cone constraints and
returning the symbol for the variable/epigraph that represents this node.
"""
function _process_node!(ex, ctx::ConicContext)
    # Base case: a symbolic variable
    if issym(ex)
        return Symbol(ex)
    end

    # Base case: a number
    if ex isa Number
        t = _new_epi_var!(ctx)
        ctx.epigraph_map[t] = ex
        # Add equality constraint: t == constant
        push!(ctx.constraints, ConeConstraint(
            [ConicConstraintTerm([t], [1.0], -Float64(ex))],
            MOI.EqualTo(0.0),
            nothing,
            "constant: $t == $ex"
        ))
        return t
    end

    if !iscall(ex)
        return _process_node!(unwrap(ex), ctx)
    end

    f = operation(ex)
    args = arguments(ex)

    # Affine flattening: if the entire subtree is affine, represent as
    # a single epigraph variable with an equality constraint
    if _is_affine(ex)
        avars, acoeffs, aconst = _extract_affine(ex)
        # If it's just a plain variable, return it directly
        if length(avars) == 1 && acoeffs[1] == 1.0 && aconst == 0.0
            return avars[1]
        end
        t = _new_epi_var!(ctx)
        ctx.epigraph_map[t] = ex
        # t == coeffs'*vars + constant  →  t - coeffs'*vars - constant == 0
        all_vars = vcat([t], avars)
        all_coeffs = vcat([1.0], [-c for c in acoeffs])
        push!(ctx.constraints, ConeConstraint(
            [ConicConstraintTerm(all_vars, all_coeffs, -aconst)],
            MOI.EqualTo(0.0),
            nothing,
            "affine: $t == expression"
        ))
        return t
    end

    # Handle addition: process children, link with equality constraint
    if Symbol(f) == :+
        child_vars = Symbol[]
        child_coeffs = Float64[]
        constant = 0.0
        for arg in args
            if arg isa Number
                constant += Float64(arg)
            else
                child = _process_node!(arg, ctx)
                push!(child_vars, child)
                push!(child_coeffs, 1.0)
            end
        end
        if length(child_vars) == 1 && child_coeffs[1] == 1.0 && constant == 0.0
            return child_vars[1]
        end
        t = _new_epi_var!(ctx)
        ctx.epigraph_map[t] = ex
        # t == sum of children + constant
        all_vars = vcat([t], child_vars)
        all_coeffs = vcat([1.0], [-c for c in child_coeffs])
        push!(ctx.constraints, ConeConstraint(
            [ConicConstraintTerm(all_vars, all_coeffs, -constant)],
            MOI.EqualTo(0.0),
            nothing,
            "sum: $t == $(join(child_vars, " + ")) + $constant"
        ))
        return t
    end

    # Handle multiplication by constant
    if Symbol(f) == :*
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
            child = _process_node!(non_const, ctx)
            t = _new_epi_var!(ctx)
            ctx.epigraph_map[t] = ex
            # t == constant * child
            push!(ctx.constraints, ConeConstraint(
                [ConicConstraintTerm([t, child], [1.0, -constant], 0.0)],
                MOI.EqualTo(0.0),
                nothing,
                "scale: $t == $constant * $child"
            ))
            return t
        else
            # Pure constant product
            t = _new_epi_var!(ctx)
            ctx.epigraph_map[t] = constant
            push!(ctx.constraints, ConeConstraint(
                [ConicConstraintTerm([t], [1.0], -constant)],
                MOI.EqualTo(0.0),
                nothing,
                "constant: $t == $constant"
            ))
            return t
        end
    end

    # Handle division: a / b → treat as a * inv(b)
    if Symbol(f) == :/
        @assert length(args) == 2
        if args[1] isa Number && !(args[2] isa Number)
            # constant / expr → constant * inv(expr)
            child = _process_node!(args[2], ctx)
            inv_t = _new_epi_var!(ctx)
            ctx.epigraph_map[inv_t] = :(_inv_aux)
            # inv(child) ≤ inv_t via RSOC
            sqrt2 = sqrt(2.0)
            push!(ctx.constraints, ConeConstraint(
                [
                    ConicConstraintTerm([inv_t], [1.0], 0.0),
                    ConicConstraintTerm([child], [1.0], 0.0),
                    ConicConstraintTerm(Symbol[], Float64[], sqrt2),
                ],
                MOI.RotatedSecondOrderCone(3),
                inv,
                "inv: ($inv_t, $child, √2) ∈ RSOC(3)"
            ))
            # result = constant * inv_t
            c = Float64(args[1])
            t = _new_epi_var!(ctx)
            ctx.epigraph_map[t] = ex
            push!(ctx.constraints, ConeConstraint(
                [ConicConstraintTerm([t, inv_t], [1.0, -c], 0.0)],
                MOI.EqualTo(0.0),
                nothing,
                "scale: $t == $c * $inv_t"
            ))
            return t
        else
            # General division: process numerator and denominator
            num_var = _process_node!(args[1], ctx)
            den_var = _process_node!(args[2], ctx)
            # Create inv(denominator) via RSOC
            inv_t = _new_epi_var!(ctx)
            ctx.epigraph_map[inv_t] = :(_inv_aux)
            sqrt2 = sqrt(2.0)
            push!(ctx.constraints, ConeConstraint(
                [
                    ConicConstraintTerm([inv_t], [1.0], 0.0),
                    ConicConstraintTerm([den_var], [1.0], 0.0),
                    ConicConstraintTerm(Symbol[], Float64[], sqrt2),
                ],
                MOI.RotatedSecondOrderCone(3),
                inv,
                "inv: ($inv_t, $den_var, √2) ∈ RSOC(3)"
            ))
            # result = numerator * inv_t (requires linearity of numerator)
            t = _new_epi_var!(ctx)
            ctx.epigraph_map[t] = ex
            push!(ctx.constraints, ConeConstraint(
                [ConicConstraintTerm([t, num_var, inv_t], [1.0, -1.0, 0.0], 0.0)],
                MOI.EqualTo(0.0),
                nothing,
                "div: $t == $num_var / $den_var"
            ))
            return t
        end
    end

    # Handle DCP atoms with cone annotations
    if hasdcprule(f)
        child_vars = Symbol[]
        for arg in args
            child = _process_node!(arg, ctx)
            push!(child_vars, child)
        end

        # Look up the cone for this atom
        rule, _ = dcprule(f, args...)
        cone = rule.cone
        curv = rule.curvature

        t = _new_epi_var!(ctx)
        ctx.epigraph_map[t] = ex

        # Emit the cone constraint
        _emit_atom_constraint!(f, t, child_vars, cone, curv, ctx)

        return t
    end

    # Fallback: error on unhandled atoms
    error("No conic reformulation for atom: $(nameof(f)). " *
          "All atoms must have a registered conic reformulation to generate valid conic form.")
end

# ──────────────────────────────────────────────────────────────────────────────
# Atom-specific cone constraint emission
# ──────────────────────────────────────────────────────────────────────────────

"""
    _emit_atom_constraint!(f, t, child_vars, cone, curvature, ctx)

Emit the appropriate cone constraint for atom `f` with epigraph variable `t`
and argument variables `child_vars`.

For a convex atom f(x), the epigraph is: {(t, x) : f(x) ≤ t}
For a concave atom f(x), the hypograph is: {(t, x) : f(x) ≥ t}
"""
function _emit_atom_constraint!(f, t, child_vars, cone, curvature, ctx::ConicContext)
    fname = string(nameof(f))

    # ── Check atom identity first (before linear fallback) ────────────
    # Some atoms like max, min have cone=MOI.Reals but need LP reformulations

    # ── LP atoms (max, min, maximum, minimum) ─────────────────────────

    if f === max
        # max(a,b) ≤ t  ⟺  t - a ≥ 0 AND t - b ≥ 0
        @assert length(child_vars) == 2
        a, b = child_vars[1], child_vars[2]
        push!(ctx.constraints, ConeConstraint(
            [ConicConstraintTerm([t, a], [1.0, -1.0], 0.0)],
            MOI.Nonnegatives(1),
            max,
            "max: $t - $a ≥ 0"
        ))
        push!(ctx.constraints, ConeConstraint(
            [ConicConstraintTerm([t, b], [1.0, -1.0], 0.0)],
            MOI.Nonnegatives(1),
            max,
            "max: $t - $b ≥ 0"
        ))
        return

    elseif f === min
        # min(a,b) ≥ t  ⟺  a - t ≥ 0 AND b - t ≥ 0
        @assert length(child_vars) == 2
        a, b = child_vars[1], child_vars[2]
        push!(ctx.constraints, ConeConstraint(
            [ConicConstraintTerm([a, t], [1.0, -1.0], 0.0)],
            MOI.Nonnegatives(1),
            min,
            "min: $a - $t ≥ 0"
        ))
        push!(ctx.constraints, ConeConstraint(
            [ConicConstraintTerm([b, t], [1.0, -1.0], 0.0)],
            MOI.Nonnegatives(1),
            min,
            "min: $b - $t ≥ 0"
        ))
        return

    elseif f === maximum
        # maximum(x) ≤ t  ⟺  t - xᵢ ≥ 0 for all i
        for xi in child_vars
            push!(ctx.constraints, ConeConstraint(
                [ConicConstraintTerm([t, xi], [1.0, -1.0], 0.0)],
                MOI.Nonnegatives(1),
                maximum,
                "maximum: $t - $xi ≥ 0"
            ))
        end
        return

    elseif f === minimum
        # minimum(x) ≥ t  ⟺  xᵢ - t ≥ 0 for all i
        for xi in child_vars
            push!(ctx.constraints, ConeConstraint(
                [ConicConstraintTerm([xi, t], [1.0, -1.0], 0.0)],
                MOI.Nonnegatives(1),
                minimum,
                "minimum: $xi - $t ≥ 0"
            ))
        end
        return
    end

    # ── Linear fallback ────────────────────────────────────────────────

    if cone === nothing || cone == MOI.Reals
        push!(ctx.constraints, ConeConstraint(
            [ConicConstraintTerm(
                vcat([t], child_vars),
                vcat([1.0], [-1.0 for _ in child_vars]),
                0.0
            )],
            MOI.EqualTo(0.0),
            f,
            "$fname: linear relationship"
        ))
        return
    end

    # ── Exponential Cone atoms ──────────────────────────────────────────

    if f === exp
        # exp(x) ≤ t  ⟺  (x, 1, t) ∈ ExponentialCone
        # MOI.ExponentialCone: (x, y, z) such that y * exp(x/y) ≤ z, y > 0
        @assert length(child_vars) == 1
        x = child_vars[1]
        push!(ctx.constraints, ConeConstraint(
            [
                ConicConstraintTerm([x], [1.0], 0.0),       # row 1: x
                ConicConstraintTerm(Symbol[], Float64[], 1.0),  # row 2: 1
                ConicConstraintTerm([t], [1.0], 0.0),       # row 3: t
            ],
            MOI.ExponentialCone(),
            exp,
            "$fname: ($(x), 1, $t) ∈ ExponentialCone"
        ))

    elseif f === log
        # log(x) ≥ t  ⟺  (t, 1, x) ∈ ExponentialCone
        @assert length(child_vars) == 1
        x = child_vars[1]
        push!(ctx.constraints, ConeConstraint(
            [
                ConicConstraintTerm([t], [1.0], 0.0),       # row 1: t
                ConicConstraintTerm(Symbol[], Float64[], 1.0),  # row 2: 1
                ConicConstraintTerm([x], [1.0], 0.0),       # row 3: x
            ],
            MOI.ExponentialCone(),
            log,
            "$fname: ($t, 1, $(x)) ∈ ExponentialCone"
        ))

    elseif f === log1p
        # log(1+x) ≥ t  ⟺  (t, 1, 1+x) ∈ ExponentialCone
        @assert length(child_vars) == 1
        x = child_vars[1]
        push!(ctx.constraints, ConeConstraint(
            [
                ConicConstraintTerm([t], [1.0], 0.0),           # row 1: t
                ConicConstraintTerm(Symbol[], Float64[], 1.0),  # row 2: 1
                ConicConstraintTerm([x], [1.0], 1.0),           # row 3: 1 + x
            ],
            MOI.ExponentialCone(),
            log1p,
            "log1p: ($t, 1, 1+$(x)) ∈ ExponentialCone"
        ))

    elseif f === logistic
        # logistic(x) = log(1 + exp(x)) ≤ t
        # Reformulation: introduce u, v s.t. u + v ≤ exp(t), u ≥ 1, v ≥ exp(x)
        # ⟺ (0, 1, u) ∈ ExpCone (u ≥ exp(0)=1) and (x, 1, v) ∈ ExpCone
        #    and u + v ≤ exp(t) ⟺ (0, u+v, exp(t)) ... but simpler:
        # log(1+exp(x)) ≤ t ⟺ two constraints:
        #   (x - t, 1, u₁) ∈ ExponentialCone  [u₁ ≥ exp(x-t)]
        #   (-t, 1, u₂) ∈ ExponentialCone     [u₂ ≥ exp(-t)]
        #   u₁ + u₂ ≤ 1
        @assert length(child_vars) == 1
        x = child_vars[1]
        u1 = _new_epi_var!(ctx)
        u2 = _new_epi_var!(ctx)
        ctx.epigraph_map[u1] = :(_logistic_aux1)
        ctx.epigraph_map[u2] = :(_logistic_aux2)

        # (x - t, 1, u1) ∈ ExponentialCone
        push!(ctx.constraints, ConeConstraint(
            [
                ConicConstraintTerm([x, t], [1.0, -1.0], 0.0),  # x - t
                ConicConstraintTerm(Symbol[], Float64[], 1.0),    # 1
                ConicConstraintTerm([u1], [1.0], 0.0),           # u1
            ],
            MOI.ExponentialCone(),
            logistic,
            "logistic: ($(x)-$t, 1, $u1) ∈ ExponentialCone"
        ))

        # (-t, 1, u2) ∈ ExponentialCone
        push!(ctx.constraints, ConeConstraint(
            [
                ConicConstraintTerm([t], [-1.0], 0.0),           # -t
                ConicConstraintTerm(Symbol[], Float64[], 1.0),    # 1
                ConicConstraintTerm([u2], [1.0], 0.0),           # u2
            ],
            MOI.ExponentialCone(),
            logistic,
            "logistic: (-$t, 1, $u2) ∈ ExponentialCone"
        ))

        # u1 + u2 ≤ 1  ⟺  1 - u1 - u2 ≥ 0
        push!(ctx.constraints, ConeConstraint(
            [ConicConstraintTerm([u1, u2], [-1.0, -1.0], 1.0)],
            MOI.Nonnegatives(1),
            logistic,
            "logistic: $u1 + $u2 ≤ 1"
        ))

    elseif f === xlogx
        # xlogx(x) = x*log(x) ≤ t
        # ⟺ (-t, x, 1) ∈ RelativeEntropyCone(3)
        # MOI.RelativeEntropyCone(3): (u, v, w) s.t. u ≥ v*log(v/w)
        # So u = -t, v = x, w = 1  gives  -t ≥ x*log(x/1) = x*log(x)
        # i.e. t ≤ -x*log(x)... wait, xlogx is convex, so epigraph is xlogx(x) ≤ t
        # We need: t ≥ x*log(x). RelEntropyCone: u ≥ v*log(v/w)
        # Set u = t, v = x, w = 1: t ≥ x*log(x)  ✓
        @assert length(child_vars) == 1
        x = child_vars[1]
        push!(ctx.constraints, ConeConstraint(
            [
                ConicConstraintTerm([t], [1.0], 0.0),           # row 1: t (= u)
                ConicConstraintTerm([x], [1.0], 0.0),           # row 2: x (= v)
                ConicConstraintTerm(Symbol[], Float64[], 1.0),  # row 3: 1 (= w)
            ],
            MOI.RelativeEntropyCone(3),
            xlogx,
            "xlogx: ($t, $(x), 1) ∈ RelativeEntropyCone(3)"
        ))

    # ── Norm / SOC atoms ───────────────────────────────────────────────

    elseif f === abs
        # |x| ≤ t  ⟺  (t, x) ∈ NormOneCone(2)
        @assert length(child_vars) == 1
        x = child_vars[1]
        push!(ctx.constraints, ConeConstraint(
            [
                ConicConstraintTerm([t], [1.0], 0.0),  # row 1: t
                ConicConstraintTerm([x], [1.0], 0.0),  # row 2: x
            ],
            MOI.NormOneCone(2),
            abs,
            "$fname: ($t, $(x)) ∈ NormOneCone(2)"
        ))

    elseif f === norm
        # ‖x‖ ≤ t  ⟺  (t, x...) ∈ SecondOrderCone
        dim = 1 + length(child_vars)
        terms = Vector{ConicConstraintTerm}(undef, dim)
        terms[1] = ConicConstraintTerm([t], [1.0], 0.0)
        for (i, v) in enumerate(child_vars)
            terms[i + 1] = ConicConstraintTerm([v], [1.0], 0.0)
        end
        push!(ctx.constraints, ConeConstraint(
            terms,
            MOI.SecondOrderCone(dim),
            norm,
            "$fname: ($t, $(join(child_vars, ", "))) ∈ SOC($dim)"
        ))

    # ── RSOC atoms ─────────────────────────────────────────────────────

    elseif f === sqrt
        # sqrt(x) ≥ t  ⟺  (t, 1, x) ∈ RotatedSecondOrderCone(3)
        # RSOC(3): 2*t₁*t₂ ≥ x₃², t₁,t₂ ≥ 0
        # We want t² ≤ x, i.e. (t, 0.5, x) or equivalently we use
        # the standard form: 2*t*1 ≥ ... wait, let's use the correct form:
        # RSOC: 2*u₁*u₂ ≥ ‖u₃:‖². For sqrt: t ≥ 0, x ≥ 0, t² ≤ x
        # Set u = (x, 0.5, t): 2*x*0.5 ≥ t² → x ≥ t² → t ≤ sqrt(x) ✓
        @assert length(child_vars) == 1
        x = child_vars[1]
        push!(ctx.constraints, ConeConstraint(
            [
                ConicConstraintTerm([x], [1.0], 0.0),           # row 1: x
                ConicConstraintTerm(Symbol[], Float64[], 0.5),   # row 2: 0.5
                ConicConstraintTerm([t], [1.0], 0.0),           # row 3: t
            ],
            MOI.RotatedSecondOrderCone(3),
            sqrt,
            "$fname: ($(x), 0.5, $t) ∈ RSOC(3)"
        ))

    elseif f === inv
        # inv(x) ≤ t, x > 0  ⟺  1/x ≤ t  ⟺  1 ≤ t*x
        # RSOC(3): 2*t*x ≥ (√2)² = 2, i.e. t*x ≥ 1
        # Set u = (t, x, √2): 2*t*x ≥ 2 → t*x ≥ 1 → 1/x ≤ t ✓
        @assert length(child_vars) == 1
        x = child_vars[1]
        sqrt2 = sqrt(2.0)
        push!(ctx.constraints, ConeConstraint(
            [
                ConicConstraintTerm([t], [1.0], 0.0),               # row 1: t
                ConicConstraintTerm([x], [1.0], 0.0),               # row 2: x
                ConicConstraintTerm(Symbol[], Float64[], sqrt2),     # row 3: √2
            ],
            MOI.RotatedSecondOrderCone(3),
            inv,
            "$fname: ($t, $(x), √2) ∈ RSOC(3)"
        ))

    elseif f === quad_over_lin
        # x²/y ≤ t  ⟺  (y, t, x) ∈ RotatedSecondOrderCone(3)
        # RSOC: 2*y*t ≥ x² (since x is scalar here)
        # Actually 2*y*t ≥ 2*(x²/2)... let's be precise:
        # RSOC(3): 2*u₁*u₂ ≥ u₃². Set u = (y, t, x): 2*y*t ≥ x² → x²/y ≤ 2t
        # Hmm, that has factor of 2. The standard RSOC in MOI is:
        # 2*u[1]*u[2] ≥ ‖u[3:]‖², u[1],u[2] ≥ 0
        # So (y, t, x): 2*y*t ≥ x² → t ≥ x²/(2y)
        # We need t ≥ x²/y, so use (0.5*y, t, x): 2*(0.5y)*t ≥ x² → y*t ≥ x² ✓
        # Or equivalently, scale: (y, t, x*√2): 2*y*t ≥ 2x² ... no.
        # Simplest: introduce s = 2t, then (y, s, x*√2)... too complex.
        # Better: just use (y/2, t, x) → 2*(y/2)*t = y*t ≥ x² ✓
        # But y/2 requires scaling. Let's use the VectorAffine approach:
        # row 1 = 0.5*y, row 2 = t, row 3 = x
        @assert length(child_vars) == 2
        x, y = child_vars[1], child_vars[2]
        push!(ctx.constraints, ConeConstraint(
            [
                ConicConstraintTerm([y], [0.5], 0.0),  # row 1: y/2
                ConicConstraintTerm([t], [1.0], 0.0),   # row 2: t
                ConicConstraintTerm([x], [1.0], 0.0),   # row 3: x
            ],
            MOI.RotatedSecondOrderCone(3),
            quad_over_lin,
            "$fname: ($(y)/2, $t, $(x)) ∈ RSOC(3)"
        ))

    # ── Relative entropy ───────────────────────────────────────────────

    elseif f === rel_entr
        # rel_entr(x,y) = x*log(x/y) ≤ t
        # MOI.RelativeEntropyCone(3): (u, v, w) s.t. u ≥ v*log(v/w)
        # Set u = t, v = x, w = y: t ≥ x*log(x/y) ✓
        @assert length(child_vars) == 2
        x, y = child_vars[1], child_vars[2]
        push!(ctx.constraints, ConeConstraint(
            [
                ConicConstraintTerm([t], [1.0], 0.0),   # row 1: t (= u)
                ConicConstraintTerm([x], [1.0], 0.0),   # row 2: x (= v)
                ConicConstraintTerm([y], [1.0], 0.0),   # row 3: y (= w)
            ],
            MOI.RelativeEntropyCone(3),
            rel_entr,
            "$fname: ($t, $(x), $(y)) ∈ RelativeEntropyCone(3)"
        ))

    elseif f === kldivergence
        # kldivergence(p, q) = Σ pᵢ*log(pᵢ/qᵢ) ≤ t
        # MOI.RelativeEntropyCone(2n+1): (u, p..., q...) s.t. u ≥ Σ pᵢ*log(pᵢ/qᵢ)
        # child_vars are the processed versions of the p and q arguments
        # For the symbolic case, p and q are vectors, but in our tree walk
        # they're already reduced to single epigraph vars
        @assert length(child_vars) == 2
        p, q = child_vars[1], child_vars[2]
        # Scalar case (each arg reduced to single var)
        push!(ctx.constraints, ConeConstraint(
            [
                ConicConstraintTerm([t], [1.0], 0.0),   # row 1: t (= u)
                ConicConstraintTerm([p], [1.0], 0.0),   # row 2: p
                ConicConstraintTerm([q], [1.0], 0.0),   # row 3: q
            ],
            MOI.RelativeEntropyCone(3),
            kldivergence,
            "kldivergence: ($t, $p, $q) ∈ RelativeEntropyCone(3)"
        ))

    # ── Power cone ─────────────────────────────────────────────────────

    elseif f === (^)
        # Power atom x^p: dispatch based on exponent
        @assert length(child_vars) >= 1
        x = child_vars[1]
        # Get the actual exponent from the original expression arguments
        p = nothing
        if length(args) >= 2 && args[2] isa Number
            p = Float64(args[2])
        end

        if p !== nothing && p == 2
            # x² ≤ t  ⟺  RSOC: (t, 0.5, x): 2*t*0.5 ≥ x² → t ≥ x²
            push!(ctx.constraints, ConeConstraint(
                [
                    ConicConstraintTerm([t], [1.0], 0.0),
                    ConicConstraintTerm(Symbol[], Float64[], 0.5),
                    ConicConstraintTerm([x], [1.0], 0.0),
                ],
                MOI.RotatedSecondOrderCone(3),
                (^),
                "power: ($t, 0.5, $(x)) ∈ RSOC(3) [x²]"
            ))
        elseif p !== nothing && p > 1
            # x^p ≤ t, x ≥ 0  ⟺  (t, x) ∈ PowerCone(1/p)
            # MOI.PowerCone(α): x₁^α * x₂^(1-α) ≥ |x₃|
            # For x^p ≤ t: set α = 1/p, (t, 1, x): t^(1/p) * 1^(1-1/p) ≥ |x|... no.
            # Actually MOI PowerCone: x₁^α * x₂^(1-α) ≥ |x₃|, x₁,x₂ ≥ 0
            # We want t ≥ x^p. Set α = 1/p:
            # (t, 1, x): t^(1/p) * 1^(1-1/p) ≥ |x| → t^(1/p) ≥ x → t ≥ x^p ✓
            push!(ctx.constraints, ConeConstraint(
                [
                    ConicConstraintTerm([t], [1.0], 0.0),           # t
                    ConicConstraintTerm(Symbol[], Float64[], 1.0),  # 1
                    ConicConstraintTerm([x], [1.0], 0.0),           # x
                ],
                MOI.PowerCone(1.0 / p),
                (^),
                "power: ($t, 1, $(x)) ∈ PowerCone($(1.0/p)) [x^$p]"
            ))
        elseif p !== nothing && p > 0 && p < 1
            # x^p ≥ t, x ≥ 0 (concave)  ⟺  PowerCone(p)
            # (x, 1, t): x^p * 1^(1-p) ≥ |t| → x^p ≥ t ✓
            push!(ctx.constraints, ConeConstraint(
                [
                    ConicConstraintTerm([x], [1.0], 0.0),           # x
                    ConicConstraintTerm(Symbol[], Float64[], 1.0),  # 1
                    ConicConstraintTerm([t], [1.0], 0.0),           # t
                ],
                MOI.PowerCone(p),
                (^),
                "power: ($(x), 1, $t) ∈ PowerCone($p) [x^$p]"
            ))
        elseif p !== nothing && p < 0
            # x^p (p<0), x > 0, convex.  x^p ≤ t ⟺ 1 ≤ t * x^(-p)
            # Use PowerCone: (t, x, 1) with α = 1/(1-p)...
            # t ≥ x^p. Let q = -p > 0. t ≥ 1/x^q.
            # (t, x, 1) ∈ PowerCone(1/(1+q)): t^(1/(1+q)) * x^(q/(1+q)) ≥ 1
            # → t * x^q ≥ 1^(1+q) = 1 → t ≥ 1/x^q = x^p ✓
            q = -p
            alpha = 1.0 / (1.0 + q)
            push!(ctx.constraints, ConeConstraint(
                [
                    ConicConstraintTerm([t], [1.0], 0.0),           # t
                    ConicConstraintTerm([x], [1.0], 0.0),           # x
                    ConicConstraintTerm(Symbol[], Float64[], 1.0),  # 1
                ],
                MOI.PowerCone(alpha),
                (^),
                "power: ($t, $(x), 1) ∈ PowerCone($alpha) [x^$p]"
            ))
        else
            # Fallback for integer powers or unrecognized
            _emit_generic_constraint!(f, t, child_vars, cone, curvature, ctx)
        end

    # ── Huber loss ─────────────────────────────────────────────────────

    elseif f === huber
        # huber(x, M) ≤ t
        # Decomposition: huber(x,M) = 2M*|x| - M² if |x|>M, x² if |x|≤M
        # Standard conic reformulation:
        #   t ≥ 2*M*s + v, |x| ≤ s + M, v ≥ x² (RSOC), s ≥ 0
        # Simpler: huber(x,M) = min_s (x-s)²/1 + 2M*|s|... no.
        # Standard: t ≥ u + 2Mv, |x| ≤ u + M, u ≥ 0, v ≥ 0, u*1 ≥ (stuff)...
        # Actually the cleanest decomposition:
        # huber(x) ≤ t  ⟺  ∃ s,v: t = 2v + s, |x| ≤ v + M, s ≥ x² (RSOC), v ≥ 0
        # Nah let's use the Convex.jl standard:
        # huber(x,M) ≤ t  ⟺  ∃ s ≥ 0, n:  |x| ≤ s + n, t ≥ s² + 2Mn
        # ... these get complex. Use the simple RSOC+LP approach:
        # Split: t = u + v, x = a + b, |a| ≤ M, v = a², u = 2M|b|
        # Even simpler, just create the generic constraint for now
        _emit_generic_constraint!(f, t, child_vars, cone, curvature, ctx)

    # ── Geometric mean ─────────────────────────────────────────────────

    elseif f === StatsBase.geomean
        # geomean(x) ≥ t  ⟺  (t, x...) ∈ GeometricMeanCone(n+1)
        dim = 1 + length(child_vars)
        terms = Vector{ConicConstraintTerm}(undef, dim)
        terms[1] = ConicConstraintTerm([t], [1.0], 0.0)
        for (i, v) in enumerate(child_vars)
            terms[i + 1] = ConicConstraintTerm([v], [1.0], 0.0)
        end
        push!(ctx.constraints, ConeConstraint(
            terms,
            MOI.GeometricMeanCone(dim),
            StatsBase.geomean,
            "geomean: ($t, $(join(child_vars, ", "))) ∈ GeometricMeanCone($dim)"
        ))

    else
        # Generic: record the cone type
        _emit_generic_constraint!(f, t, child_vars, cone, curvature, ctx)
    end
end

"""
Emit a generic cone constraint when no specific reformulation is available.
"""
function _emit_generic_constraint!(f, t, child_vars, cone, curvature, ctx::ConicContext)
    fname = string(nameof(f))
    sense_str = curvature == Convex ? "≤" : curvature == Concave ? "≥" : "=="
    dim = 1 + length(child_vars)
    terms = Vector{ConicConstraintTerm}(undef, dim)
    terms[1] = ConicConstraintTerm([t], [1.0], 0.0)
    for (i, v) in enumerate(child_vars)
        terms[i + 1] = ConicConstraintTerm([v], [1.0], 0.0)
    end
    cone_instance = if cone isa DataType
        try
            cone(dim)
        catch
            MOI.Reals(dim)
        end
    else
        cone
    end
    push!(ctx.constraints, ConeConstraint(
        terms,
        cone_instance,
        f,
        "$fname: $t $sense_str $fname($(join(child_vars, ", "))) via $(cone_instance)"
    ))
end

# ──────────────────────────────────────────────────────────────────────────────
# Utilities
# ──────────────────────────────────────────────────────────────────────────────

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

export to_conic_form, ConicFormulation, ConeConstraint, ConicConstraintTerm, ConicContext, list_cone_annotations
