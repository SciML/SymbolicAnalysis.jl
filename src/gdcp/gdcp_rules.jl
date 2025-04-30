using Manifolds
using Symbolics: @register_symbolic, unwrap
using LinearAlgebra

# @enum GSign GPositive GNegative GAnySign
@enum GCurvature GConvex GConcave GLinear GUnknownCurvature
@enum GMonotonicity GIncreasing GDecreasing GAnyMono

const gdcprules_dict = Dict()

function add_gdcprule(f, manifold, sign, curvature, monotonicity)
    if !(monotonicity isa Tuple)
        monotonicity = (monotonicity,)
    end
    gdcprules_dict[f] = makegrule(manifold, sign, curvature, monotonicity)
end
makegrule(manifold, sign, curvature, monotonicity) =
    (manifold = manifold, sign = sign, gcurvature = curvature, gmonotonicity = monotonicity)

hasgdcprule(f::Function) = haskey(gdcprules_dict, f)
hasgdcprule(f) = false
gdcprule(f, args...) = gdcprules_dict[f], args

setgcurvature(ex::Union{Symbolic,Num}, curv) = setmetadata(ex, GCurvature, curv)
setgcurvature(ex, curv) = ex
getgcurvature(ex::Union{Symbolic,Num}) = getmetadata(ex, GCurvature)
getgcurvature(ex) = GLinear
hasgcurvature(ex::Union{Symbolic,Num}) = hasmetadata(ex, GCurvature)
hasgcurvature(ex) = ex isa Real

function mul_gcurvature(args)
    non_constants = findall(x -> issym(x) || iscall(x), args)
    constants = findall(x -> !issym(x) && !iscall(x), args)
    try
        @assert length(non_constants) <= 1
    catch
        @warn "DGCP does not support multiple non-constant arguments in multiplication"
        return UnknownGCurvature
    end
    if !isempty(non_constants)
        expr = args[first(non_constants)]
        curv = find_gcurvature(expr)
        return if prod(args[constants]) < 0
            # flip
            curv == GConvex ? GConcave : curv == GConcave ? GConvex : curv
        else
            curv
        end
    end
    return GLinear
end

function add_gcurvature(args)
    curvs = find_gcurvature.(args)
    all(==(GLinear), curvs) && return GLinear
    all(x -> x == GConvex || x == GLinear, curvs) && return GConvex
    all(x -> x == GConcave || x == GLinear, curvs) && return GConcave
    return GUnknownCurvature
end

function find_gcurvature(ex)
    if hasgcurvature(ex)
        return getgcurvature(ex)
    end
    if iscall(ex)
        f, args = operation(ex), arguments(ex)
        knowngcurv = false

        if hasgdcprule(f) && !any(iscall.(args))
            rule, args = gdcprule(f, args...)
            f_curvature = rule.gcurvature
            f_monotonicity = rule.gmonotonicity
            knowngcurv = true
        elseif f == LinearAlgebra.logdet
            if operation(args[1]) == conjugation ||
               operation(args[1]) == LinearAlgebra.diag ||
               Symbol(operation(args[1])) == :+ ||
               operation(args[1]) == affine_map ||
               operation(args[1]) == hadamard_product
                return GConvex
            end
        elseif f == log &&
               iscall(args[1]) &&
               (operation(args[1]) == LinearAlgebra.tr || operation(args[1]) == quad_form)
            return GConvex
        elseif (f == schatten_norm || f == eigsummax) && operation(args[1]) == log
            return GConvex
        elseif f == sum_log_eigmax && hasdcprule(args[1])
            if dcprule(operation(args[1])) == Convex
                return GConvex
            else
                return GUnknownCurvature
            end
        elseif f == affine_map
            if args[1] == tr || args[1] == conjugation || args[1] == diag
                return GConvex
            else
                return GUnknownCurvature
            end
        elseif hasgdcprule(f) && any(iscall.(args))
            for i in eachindex(args)
                if iscall(args[i])
                    if operation(args[i]) == inv
                        rule, args = gdcprule(f, args...)
                        f_curvature = rule.gcurvature
                        f_monotonicity = if rule.gmonotonicity == GIncreasing
                            GDecreasing
                        elseif rule.gmonotonicity == GDecreasing
                            GIncreasing
                        else
                            GAnyMono
                        end
                        knowngcurv = true
                    elseif operation(args[i]) == broadcast
                        rule, args = gdcprule(f, args...)
                        f_curvature = rule.gcurvature
                        f_monotonicity = rule.gmonotonicity
                        knowngcurv = true
                    elseif operation(args[i]) == affine_map
                        rule, args = gdcprule(f, args...)
                        f_curvature = rule.gcurvature
                        f_monotonicity = rule.gmonotonicity
                        knowngcurv = true
                    end
                end
            end
        elseif Symbol(f) == :*
            if args[1] isa Number && args[1] > 0
                return find_gcurvature(args[2])
            elseif args[1] isa Number && args[1] < 0
                argscurv = find_gcurvature(args[2])
                if argscurv == GConvex
                    return GConcave
                elseif argscurv == GConcave
                    return GConvex
                else
                    argscurv
                end
            else
                @warn "Disciplined Programming does not support multiple non-constant arguments in multiplication"
                return UnknownGCurvature
            end
        end

        if !(knowngcurv) && hasdcprule(f)
            rule, args = dcprule(f, args...)
            f_curvature = rule.curvature
            f_monotonicity = rule.monotonicity
        end

        if f_curvature == Convex || f_curvature == Affine
            if all(enumerate(args)) do (i, arg)
                arg_curv = find_gcurvature(arg)
                m = get_arg_property(f_monotonicity, i, args)
                # @show arg
                if arg_curv == GConvex
                    m == Increasing
                elseif arg_curv == GConcave
                    m == Decreasing
                else
                    arg_curv == GLinear
                end
            end
                return GConvex
            end
        elseif f_curvature == Concave || f_curvature == Affine
            if all(enumerate(args)) do (i, arg)
                arg_curv = find_gcurvature(arg)
                m = f_monotonicity[i]
                if arg_curv == GConcave
                    m == Increasing
                elseif arg_curv == GConvex
                    m == Decreasing
                else
                    arg_curv == GLinear
                end
            end
                return GConcave
            end
        elseif f_curvature == Affine
            if all(enumerate(args)) do (i, arg)
                arg_curv = find_gcurvature(arg)
                arg_curv == GLinear
            end
                return GLinear
            end
        elseif f_curvature isa GCurvature
            return f_curvature
        else
            return GUnknownCurvature
        end
    elseif hasfield(typeof(ex), :val) && operation(ex.val) in keys(gdcprules_dict)
        f, args = operation(ex.val), arguments(ex.val)
        rule, args = gdcprule(f, args...)
        return rule.gcurvature
    else
        return GLinear
    end
end

function propagate_gcurvature(ex, M::AbstractManifold)
    r = [
        @rule *(~~x) => setgcurvature(~MATCH, mul_gcurvature(~~x))
        @rule +(~~x) => setgcurvature(~MATCH, add_gcurvature(~~x))
        @rule ~x => setgcurvature(~x, find_gcurvature(~x))
    ]
    ex = Postwalk(Chain(r))(ex)
    ex = Prewalk(Chain(r))(ex)
    return ex
end
