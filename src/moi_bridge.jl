"""
    MOI/JuMP Bridge

Converts a `ConicFormulation` (from `to_conic_form`) into an MOI model or JuMP model
that can be solved by any MOI-compatible solver.
"""

import JuMP
using MathOptInterface
const __MOI = MathOptInterface

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

Add a single ConeConstraint to a JuMP model.
"""
function _add_jump_constraint!(model, c::ConeConstraint, jump_vars)
    if c.cone isa __MOI.EqualTo
        expr = JuMP.AffExpr(c.func_constant)
        for (v, coeff) in zip(c.func_vars, c.func_coeffs)
            JuMP.add_to_expression!(expr, coeff, jump_vars[v])
        end
        JuMP.@constraint(model, expr == 0)

    elseif c.cone isa __MOI.Zeros
        expr = JuMP.AffExpr(c.func_constant)
        for (v, coeff) in zip(c.func_vars, c.func_coeffs)
            JuMP.add_to_expression!(expr, coeff, jump_vars[v])
        end
        JuMP.@constraint(model, expr == 0)

    elseif c.cone isa __MOI.ExponentialCone
        vars = [jump_vars[v] for v in c.func_vars]
        if length(vars) >= 2
            JuMP.@constraint(model, [vars[1], 1.0, vars[2]] in __MOI.ExponentialCone())
        end

    elseif c.cone isa __MOI.SecondOrderCone
        vars = [jump_vars[v] for v in c.func_vars]
        JuMP.@constraint(model, vars in JuMP.SecondOrderCone())

    elseif c.cone isa __MOI.RotatedSecondOrderCone
        vars = [jump_vars[v] for v in c.func_vars]
        if length(vars) == 2
            JuMP.@constraint(model, [vars[1], vars[2], 1.0] in JuMP.RotatedSecondOrderCone())
        else
            JuMP.@constraint(model, vars in JuMP.RotatedSecondOrderCone())
        end

    elseif c.cone isa __MOI.NormOneCone
        vars = [jump_vars[v] for v in c.func_vars]
        dim = c.cone.dimension
        JuMP.@constraint(model, vars in __MOI.NormOneCone(dim))

    elseif c.cone isa __MOI.RelativeEntropyCone
        vars = [jump_vars[v] for v in c.func_vars]
        dim = c.cone.dimension
        JuMP.@constraint(model, vars in __MOI.RelativeEntropyCone(dim))

    else
        # For other cone types, add a placeholder bound
        vars = [jump_vars[v] for v in c.func_vars]
        if length(vars) > 0
            JuMP.@constraint(model, sum(vars) >= 0)
        end
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
    model = __MOI.Utilities.Model{Float64}()

    # Add variables
    var_map = Dict{Symbol, __MOI.VariableIndex}()
    for v in cf.variables
        vi = __MOI.add_variable(model)
        __MOI.set(model, __MOI.VariableName(), vi, string(v))
        var_map[v] = vi
    end

    # Set objective
    obj_vi = var_map[cf.objective_var]
    obj_func = __MOI.ScalarAffineFunction(
        [__MOI.ScalarAffineTerm(1.0, obj_vi)],
        0.0
    )
    sense = cf.objective_sense == :minimize ? __MOI.MIN_SENSE : __MOI.MAX_SENSE
    __MOI.set(model, __MOI.ObjectiveSense(), sense)
    __MOI.set(model, __MOI.ObjectiveFunction{typeof(obj_func)}(), obj_func)

    # Add constraints
    for c in cf.constraints
        _add_moi_constraint!(model, c, var_map)
    end

    return model, var_map
end

"""
    _add_moi_constraint!(model, c::ConeConstraint, var_map)

Add a single ConeConstraint to an MOI model.
"""
function _add_moi_constraint!(model, c::ConeConstraint, var_map)
    if c.cone isa __MOI.EqualTo
        terms = [__MOI.ScalarAffineTerm(coeff, var_map[v])
                 for (v, coeff) in zip(c.func_vars, c.func_coeffs)]
        func = __MOI.ScalarAffineFunction(terms, c.func_constant)
        __MOI.add_constraint(model, func, __MOI.EqualTo(0.0))

    elseif c.cone isa __MOI.Zeros
        terms = [__MOI.VectorAffineTerm(1, __MOI.ScalarAffineTerm(coeff, var_map[v]))
                 for (v, coeff) in zip(c.func_vars, c.func_coeffs)]
        func = __MOI.VectorAffineFunction(terms, [c.func_constant])
        __MOI.add_constraint(model, func, __MOI.Zeros(1))

    elseif c.cone isa __MOI.ExponentialCone
        vars = [var_map[v] for v in c.func_vars]
        if length(vars) >= 2
            terms = [
                __MOI.VectorAffineTerm(1, __MOI.ScalarAffineTerm(1.0, vars[1])),
                __MOI.VectorAffineTerm(3, __MOI.ScalarAffineTerm(1.0, vars[2]))
            ]
            func = __MOI.VectorAffineFunction(terms, [0.0, 1.0, 0.0])
            __MOI.add_constraint(model, func, __MOI.ExponentialCone())
        end

    elseif c.cone isa __MOI.SecondOrderCone
        vars = [var_map[v] for v in c.func_vars]
        dim = length(vars)
        terms = [__MOI.VectorAffineTerm(i, __MOI.ScalarAffineTerm(1.0, vars[i]))
                 for i in 1:dim]
        func = __MOI.VectorAffineFunction(terms, zeros(dim))
        __MOI.add_constraint(model, func, __MOI.SecondOrderCone(dim))

    elseif c.cone isa __MOI.RotatedSecondOrderCone
        vars = [var_map[v] for v in c.func_vars]
        if length(vars) == 2
            terms = [
                __MOI.VectorAffineTerm(1, __MOI.ScalarAffineTerm(1.0, vars[1])),
                __MOI.VectorAffineTerm(2, __MOI.ScalarAffineTerm(1.0, vars[2]))
            ]
            func = __MOI.VectorAffineFunction(terms, [0.0, 0.0, 1.0])
            __MOI.add_constraint(model, func, __MOI.RotatedSecondOrderCone(3))
        else
            dim = length(vars)
            terms = [__MOI.VectorAffineTerm(i, __MOI.ScalarAffineTerm(1.0, vars[i]))
                     for i in 1:dim]
            func = __MOI.VectorAffineFunction(terms, zeros(dim))
            __MOI.add_constraint(model, func, __MOI.RotatedSecondOrderCone(dim))
        end

    elseif c.cone isa __MOI.NormOneCone
        vars = [var_map[v] for v in c.func_vars]
        dim = c.cone.dimension
        terms = [__MOI.VectorAffineTerm(i, __MOI.ScalarAffineTerm(1.0, vars[i]))
                 for i in 1:min(dim, length(vars))]
        func = __MOI.VectorAffineFunction(terms, zeros(dim))
        __MOI.add_constraint(model, func, __MOI.NormOneCone(dim))

    elseif c.cone isa __MOI.RelativeEntropyCone
        vars = [var_map[v] for v in c.func_vars]
        dim = c.cone.dimension
        terms = [__MOI.VectorAffineTerm(i, __MOI.ScalarAffineTerm(1.0, vars[i]))
                 for i in 1:min(dim, length(vars))]
        func = __MOI.VectorAffineFunction(terms, zeros(dim))
        __MOI.add_constraint(model, func, __MOI.RelativeEntropyCone(dim))
    end
end

"""
    print_conic_form(cf::ConicFormulation; io=stdout)

Pretty-print a conic formulation.
"""
function print_conic_form(cf::ConicFormulation; io = stdout)
    println(io, "Conic Formulation:")
    println(io, "  Objective: $(cf.objective_sense) $(cf.objective_var)")
    println(io, "  Original variables: $(join(sort(collect(cf.original_variables)), ", "))")
    println(io, "  Epigraph variables: $(join(sort(collect(setdiff(cf.variables, cf.original_variables))), ", "))")
    println(io, "  Constraints ($(length(cf.constraints))):")
    for (i, c) in enumerate(cf.constraints)
        println(io, "    [$i] $(c.description)")
    end
end

export to_jump_model, to_moi_model, print_conic_form
