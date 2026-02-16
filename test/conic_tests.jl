using SymbolicAnalysis
using SymbolicAnalysis:
    dcprules_dict,
    gdcprules_dict,
    propagate_curvature,
    propagate_sign,
    getcurvature,
    getsign,
    Convex,
    Concave,
    Affine,
    Positive,
    Negative

using Symbolics
using Symbolics: unwrap
using LinearAlgebra
using MathOptInterface
const MOI = MathOptInterface
import JuMP
using Test

@testset "Cone Annotations" begin
    @testset "DCP atoms have cone field" begin
        for (f, rule) in dcprules_dict
            if rule isa Vector
                for r in rule
                    @test hasproperty(r, :cone)
                end
            else
                @test hasproperty(rule, :cone)
            end
        end
    end

    @testset "GDCP atoms have cone field" begin
        for (f, rule) in gdcprules_dict
            @test hasproperty(rule, :cone)
        end
    end

    @testset "Specific DCP cone mappings" begin
        # Check that key atoms map to expected cones
        exp_rule = dcprules_dict[exp]
        if exp_rule isa Vector
            @test any(r -> r.cone == MOI.ExponentialCone, exp_rule)
        else
            @test exp_rule.cone == MOI.ExponentialCone
        end

        abs_rule = dcprules_dict[abs]
        if abs_rule isa Vector
            @test any(r -> r.cone == MOI.NormOneCone, abs_rule)
        else
            @test abs_rule.cone == MOI.NormOneCone
        end

        logdet_rule = dcprules_dict[LinearAlgebra.logdet]
        if logdet_rule isa Vector
            @test any(r -> r.cone == MOI.LogDetConeTriangle, logdet_rule)
        else
            @test logdet_rule.cone == MOI.LogDetConeTriangle
        end
    end

    @testset "Specific GDCP cone mappings" begin
        # SPD atoms
        @test gdcprules_dict[LinearAlgebra.logdet].cone == MOI.LogDetConeTriangle
        @test gdcprules_dict[LinearAlgebra.tr].cone == MOI.Reals
    end

    @testset "list_cone_annotations" begin
        annotations = list_cone_annotations()
        @test length(annotations) > 0
        # Check structure
        for a in annotations
            @test haskey(a, :atom)
            @test haskey(a, :type)
            @test haskey(a, :cone)
            @test a.type ∈ (:DCP, :GDCP)
        end
    end
end

@testset "Conic Form Generation" begin
    @variables x y

    @testset "exp(x) → ExponentialCone" begin
        cf = to_conic_form(exp(x) |> unwrap)
        @test cf isa ConicFormulation
        @test cf.objective_sense == :minimize
        @test :x ∈ cf.original_variables
        @test length(cf.constraints) >= 1

        # Should have an ExponentialCone constraint
        exp_cones = filter(c -> c.cone isa MOI.ExponentialCone, cf.constraints)
        @test length(exp_cones) >= 1
    end

    @testset "log(x) → ExponentialCone (concave)" begin
        cf = to_conic_form(log(x) |> unwrap)
        @test cf.objective_sense == :maximize
        exp_cones = filter(c -> c.cone isa MOI.ExponentialCone, cf.constraints)
        @test length(exp_cones) >= 1
    end

    @testset "abs(x) → NormOneCone" begin
        cf = to_conic_form(abs(x) |> unwrap)
        @test cf.objective_sense == :minimize
        norm_cones = filter(c -> c.cone isa MOI.NormOneCone, cf.constraints)
        @test length(norm_cones) >= 1
    end

    @testset "exp(x) + abs(x) → mixed cones" begin
        cf = to_conic_form((exp(x) + abs(x)) |> unwrap)
        @test cf.objective_sense == :minimize
        @test length(cf.constraints) >= 3  # abs + exp + sum

        exp_cones = filter(c -> c.cone isa MOI.ExponentialCone, cf.constraints)
        @test length(exp_cones) >= 1

        norm_cones = filter(c -> c.cone isa MOI.NormOneCone, cf.constraints)
        @test length(norm_cones) >= 1
    end

    @testset "2*abs(x) - 1 → scaled" begin
        cf = to_conic_form((2 * abs(x) - 1) |> unwrap)
        @test :x ∈ cf.original_variables
        @test length(cf.constraints) >= 2
    end

    @testset "Epigraph variables are distinct from original" begin
        cf = to_conic_form(exp(x) |> unwrap)
        epigraph = setdiff(cf.variables, cf.original_variables)
        @test length(epigraph) >= 1
        @test cf.objective_var ∈ epigraph
    end

    @testset "print_conic_form does not error" begin
        cf = to_conic_form(exp(x) |> unwrap)
        io = IOBuffer()
        print_conic_form(cf; io = io)
        output = String(take!(io))
        @test contains(output, "Conic Formulation")
        @test contains(output, "ExponentialCone")
    end
end

@testset "MOI Bridge" begin
    @variables x y

    @testset "to_jump_model creates valid model" begin
        cf = to_conic_form(exp(x) |> unwrap)
        model = to_jump_model(cf)
        @test model isa JuMP.Model
        # Should have variables
        @test JuMP.num_variables(model) >= 2  # x + at least 1 epigraph var
    end

    @testset "to_moi_model creates valid model" begin
        cf = to_conic_form(exp(x) |> unwrap)
        moi_model, var_map = to_moi_model(cf)
        @test length(var_map) >= 2
        # Should have an exponential cone constraint
        exp_ci = MOI.get(moi_model,
            MOI.ListOfConstraintIndices{
                MOI.VectorAffineFunction{Float64},
                MOI.ExponentialCone
            }())
        @test length(exp_ci) >= 1
    end

    @testset "abs(x) model has NormOneCone" begin
        cf = to_conic_form(abs(x) |> unwrap)
        moi_model, var_map = to_moi_model(cf)
        norm_ci = MOI.get(moi_model,
            MOI.ListOfConstraintIndices{
                MOI.VectorAffineFunction{Float64},
                MOI.NormOneCone
            }())
        @test length(norm_ci) >= 1
    end

    @testset "Composite expression model" begin
        cf = to_conic_form((exp(x) + abs(x)) |> unwrap)
        model = to_jump_model(cf)
        @test JuMP.num_variables(model) >= 3  # x + exp epi + abs epi + sum epi
    end
end

import JuMP

@testset "JuMP Model Structure" begin
    @variables x

    @testset "exp(x) JuMP model is minimization" begin
        cf = to_conic_form(exp(x) |> unwrap)
        model = to_jump_model(cf)
        @test JuMP.objective_sense(model) == MOI.MIN_SENSE
    end

    @testset "log(x) JuMP model is maximization" begin
        cf = to_conic_form(log(x) |> unwrap)
        model = to_jump_model(cf)
        @test JuMP.objective_sense(model) == MOI.MAX_SENSE
    end
end
