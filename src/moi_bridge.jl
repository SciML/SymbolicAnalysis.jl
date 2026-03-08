"""
    MOI/JuMP Bridge

Converts a `ConicFormulation` (from `to_conic_form`) into an MOI model or JuMP model
that can be solved by any MOI-compatible solver.

Uses the new vector-valued `ConeConstraint` struct with `ConicConstraintTerm` rows,
enabling a single generic dispatch instead of per-cone-type if-elseif chains.
"""

import JuMP

# Use the module-level MOI alias from SymbolicAnalysis.jl

"""
    to_jump_model(cf::ConicFormulation; solver=nothing)

Convert a `ConicFormulation` to a JuMP `Model`.

# Arguments
- `cf::ConicFormulation` — the conic formulation from `to_conic_form`
- `solver` — optional solver (e.g., `SCS.Optimizer`). If `nothing`, creates model without solver.

# Returns
A JuMP `Model` with variables, objective, and cone constraints.
"""
function to_jump_model(cf::ConicFormulation; solver = nothing)
    model = solver === nothing ? JuMP.Model() : JuMP.Model(solver)

    # Create JuMP variables for all variables in the formulation
    jump_vars = Dict{Symbol, JuMP.VariableRef}()
    for v in cf.variables
        jump_vars[v] = JuMP.@variable(model, base_name = string(v))
    end

    # Set objective
    obj_var = jump_vars[cf.objective_var]
    if cf.objective_sense == :minimize
        JuMP.@objective(model, Min, obj_var)
    else
        JuMP.@objective(model, Max, obj_var)
    end

    # Add constraints
    for c in cf.constraints
        _add_jump_constraint!(model, c, jump_vars)
    end

    return model
end

"""
    _add_jump_constraint!(model, c::ConeConstraint, jump_vars)

Add a single ConeConstraint to a JuMP model using generic dispatch.
"""
function _add_jump_constraint!(model, c::ConeConstraint, jump_vars)
    return if c.cone isa MOI.AbstractScalarSet
        ct = only(c.terms)
        expr = JuMP.AffExpr(ct.constant)
        for (v, coeff) in zip(ct.vars, ct.coeffs)
            JuMP.add_to_expression!(expr, coeff, jump_vars[v])
        end
        JuMP.@constraint(model, expr in c.cone)
    else
        @assert c.cone isa MOI.AbstractVectorSet
        vec_expr = Vector{JuMP.AffExpr}(undef, length(c.terms))
        for (row, ct) in enumerate(c.terms)
            expr = JuMP.AffExpr(ct.constant)
            for (v, coeff) in zip(ct.vars, ct.coeffs)
                JuMP.add_to_expression!(expr, coeff, jump_vars[v])
            end
            vec_expr[row] = expr
        end
        JuMP.@constraint(model, vec_expr in c.cone)
    end
end

"""
    to_moi_model(cf::ConicFormulation)

Convert a `ConicFormulation` to a raw MOI model.

# Returns
A tuple `(model, variable_map)` where:
- `model` is an `MOI.Utilities.Model{Float64}`
- `variable_map` is a `Dict{Symbol, MOI.VariableIndex}`
"""
function to_moi_model(cf::ConicFormulation)
    model = MOI.Utilities.Model{Float64}()

    # Add variables
    var_map = Dict{Symbol, MOI.VariableIndex}()
    for v in cf.variables
        vi = MOI.add_variable(model)
        MOI.set(model, MOI.VariableName(), vi, string(v))
        var_map[v] = vi
    end

    # Set objective
    obj_vi = var_map[cf.objective_var]
    obj_func = MOI.ScalarAffineFunction([MOI.ScalarAffineTerm(1.0, obj_vi)], 0.0)
    sense = cf.objective_sense == :minimize ? MOI.MIN_SENSE : MOI.MAX_SENSE
    MOI.set(model, MOI.ObjectiveSense(), sense)
    MOI.set(model, MOI.ObjectiveFunction{typeof(obj_func)}(), obj_func)

    # Add constraints
    for c in cf.constraints
        _add_moi_constraint!(model, c, var_map)
    end

    return model, var_map
end

"""
    _add_moi_constraint!(model, c::ConeConstraint, var_map)

Add a single ConeConstraint to an MOI model using generic dispatch.
"""
function _add_moi_constraint!(model, c::ConeConstraint, var_map)
    return if c.cone isa MOI.AbstractScalarSet
        ct = only(c.terms)
        terms = [
            MOI.ScalarAffineTerm(coeff, var_map[v]) for
                (v, coeff) in zip(ct.vars, ct.coeffs)
        ]
        func = MOI.ScalarAffineFunction(terms, ct.constant)
        MOI.add_constraint(model, func, c.cone)
    else
        @assert c.cone isa MOI.AbstractVectorSet
        vat = MOI.VectorAffineTerm{Float64}[]
        for (row, ct) in enumerate(c.terms)
            for (v, coeff) in zip(ct.vars, ct.coeffs)
                push!(
                    vat,
                    MOI.VectorAffineTerm(row, MOI.ScalarAffineTerm(coeff, var_map[v])),
                )
            end
        end
        constants = [ct.constant for ct in c.terms]
        func = MOI.VectorAffineFunction(vat, constants)
        MOI.add_constraint(model, func, c.cone)
    end
end

"""
    extract_solution(cf::ConicFormulation, model, var_map)

Extract solution values from a solved MOI model back to the original variable names.

# Arguments
- `cf::ConicFormulation` — the conic formulation
- `model` — a solved MOI model
- `var_map::Dict{Symbol, MOI.VariableIndex}` — variable mapping from `to_moi_model`

# Returns
A `Dict{Symbol, Float64}` mapping original variable names to their optimal values.
"""
function extract_solution(cf::ConicFormulation, model, var_map)
    result = Dict{Symbol, Float64}()
    for v in cf.original_variables
        if haskey(var_map, v)
            val = MOI.get(model, MOI.VariablePrimal(), var_map[v])
            result[v] = val
        end
    end
    return result
end

"""
    print_conic_form(cf::ConicFormulation; io=stdout)

Pretty-print a conic formulation.
"""
function print_conic_form(cf::ConicFormulation; io = stdout)
    println(io, "Conic Formulation:")
    println(io, "  Objective: $(cf.objective_sense) $(cf.objective_var)")
    println(io, "  Original variables: $(join(sort(collect(cf.original_variables)), ", "))")
    println(
        io,
        "  Epigraph variables: $(join(sort(collect(setdiff(cf.variables, cf.original_variables))), ", "))",
    )
    println(io, "  Constraints ($(length(cf.constraints))):")
    for (i, c) in enumerate(cf.constraints)
        println(io, "    [$i] $(c.description)")
        for (j, term) in enumerate(c.terms)
            parts = String[]
            for (v, coeff) in zip(term.vars, term.coeffs)
                if coeff == 1.0
                    push!(parts, string(v))
                elseif coeff == -1.0
                    push!(parts, "-$(v)")
                else
                    push!(parts, "$(coeff)*$(v)")
                end
            end
            if term.constant != 0.0
                push!(parts, string(term.constant))
            end
            expr_str = isempty(parts) ? "0" : join(parts, " + ")
            println(io, "        row $j: $expr_str")
        end
    end
    return
end

export to_jump_model, to_moi_model, print_conic_form, extract_solution
