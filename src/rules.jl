
@enum Sign Positive Negative AnySign
@enum Curvature Convex Concave Affine UnknownCurvature
@enum Monotonicity Increasing Decreasing AnyMono

struct CustomDomain{T} <: Domain{T}
    in::Function
end

Base.in(x, c::CustomDomain) = c.in(x)

function array_domain(element_domain)
    CustomDomain{AbstractArray}() do xs
        all(in(element_domain), xs)
    end
end

function array_domain(element_domain, N)
    CustomDomain{AbstractArray{<:Any, N}}() do xs
        ndims(xs) == N && all(in(element_domain), xs)
    end
end

function symmetric_domain()
    CustomDomain{AbstractArray{<:Any, 2}}(issymmetric)
end

function semidefinite_domain()
    CustomDomain{AbstractArray{<:Any, 2}}(isposdef) #not semi so needs to change
end

function negsemidefinite_domain()
    CustomDomain{AbstractArray{<:Any, 2}}(isposdef ∘ -) #not semi so needs to change
end

function definite_domain()
    CustomDomain{AbstractArray{<:Any, 2}}(isposdef)
end

function negdefinite_domain()
    CustomDomain{AbstractArray{<:Any, 2}}(isposdef ∘ -)
end

function function_domain()
    CustomDomain{Function}(x -> typeassert(x, Function))
end

function increasing_if_positive(x)
    sign = getsign(x)
    sign == AnySign ? AnyMono : sign == Positive ? Increasing : Decreasing
end

const dcprules_dict = Dict()

function add_dcprule(f, domain, sign, curvature, monotonicity)
    if !(monotonicity isa Tuple)
        monotonicity = (monotonicity,)
    end
    if f in keys(dcprules_dict)
        dcprules_dict[f] = vcat(
            dcprules_dict[f], makerule(domain, sign, curvature, monotonicity))
    else
        dcprules_dict[f] = makerule(domain, sign, curvature, monotonicity)
    end
end

function makerule(domain, sign, curvature, monotonicity)
    (; domain = domain, sign = sign, curvature = curvature, monotonicity = monotonicity)
end

hasdcprule(f::Function) = haskey(dcprules_dict, f)
hasdcprule(f) = false

Symbolics.hasmetadata(::Union{Real, AbstractArray{<:Real}}, args...) = false

function dcprule(f, args...)
    if all(hasmetadata.(args, Ref(VarDomain)))
        argsdomain = getmetadata.(args, Ref(VarDomain))
    else
        if dcprules_dict[f] isa Vector
            return dcprules_dict[f][1], args
        else
            return dcprules_dict[f], args
        end
    end

    if dcprules_dict[f] isa Vector
        for i in 1:length(dcprules_dict[f])
            if (dcprules_dict[f][i].domain isa Domain) &&
               all(issubset.(argsdomain, Ref(dcprules_dict[f][i].domain)))
                return dcprules_dict[f][i], args
            elseif !(dcprules_dict[f][i].domain isa Domain) &&
                   all(issubset.(argsdomain, dcprules_dict[f][i].domain))
                return dcprules_dict[f][i], args
            else
                throw(
                    ArgumentError(
                    "No DCP rule found for $f with arguments $args with domain $argsdomain",
                ),
                )
            end
        end
    elseif (dcprules_dict[f].domain isa Domain) &&
           all(issubset.(argsdomain, Ref(dcprules_dict[f].domain)))
        return dcprules_dict[f], args
    elseif dcprules_dict[f].domain isa Tuple &&
           all(issubset.(argsdomain, dcprules_dict[f].domain))
        return dcprules_dict[f], args
    else
        throw(ArgumentError("No DCP rule found for $f with arguments $args"))
    end
end

### Sign ###
setsign(ex::Union{Num, Symbolic}, sign) = setmetadata(ex, Sign, sign)
setsign(ex, sign) = ex

function getsign(ex::Union{Num, Symbolic})
    if hasmetadata(ex, Sign)
        return getmetadata(ex, Sign)
    end
    return AnySign
end

getsign(ex::Union{AbstractFloat, Integer}) = ex < 0 ? Negative : Positive

function getsign(ex::AbstractArray)
    if all(x -> getsign(x) == Negative, ex)
        return Negative
    elseif all(x -> getsign(x) == Positive, ex)
        return Positive
    else
        AnySign
    end
end

hassign(ex::Union{Num, Symbolic}) = hasmetadata(ex, Sign)
hassign(ex) = ex isa Real

hassign(ex::typeof(Base.broadcast)) = true
getsign(ex::typeof(Base.broadcast)) = Positive

Symbolics.arguments(x::Number) = x

function add_sign(args)
    if hassign(args)
        return getsign(args)
    end
    for i in eachindex(args)
        if iscall(args[i])
            args[i] = propagate_sign(args[i])
        end
    end
    signs = reduce(vcat, getsign.(args))
    if any(==(AnySign), signs)
        AnySign
    elseif all(==(Negative), signs)
        Negative
    elseif all(==(Positive), signs)
        Positive
    else
        AnySign
    end
end

function mul_sign(args)
    signs = getsign.(args)
    if any(==(AnySign), signs)
        AnySign
    elseif isodd(count(==(Negative), signs))
        Negative
    else
        Positive
    end
end

function propagate_sign(ex)
    # Step 1: set the sign of all variables to be AnySign
    rs = [@rule ~x::issym => setsign(~x, AnySign) where {hassign(~x)}
          @rule ~x::iscall => setsign(~x, AnySign) where {hassign(~x)}
          @rule ~x::issym => setsign(~x, (dcprule(~x))[1].sign) where {hasdcprule(~x)}
          @rule ~x::issym => setsign(~x, (gdcprule(~x))[1].sign) where {hasgdcprule(~x)}
          @rule ~x::iscall => setsign(
              ~x,
              (dcprule(operation(~x), arguments(~x)...)[1].sign)
          ) where {hasdcprule(operation(~x))}
          @rule ~x::iscall => setsign(
              ~x,
              (gdcprule(operation(~x), arguments(~x)...)[1].sign)
          ) where {hasgdcprule(operation(~x))}
          @rule *(~~x) => setsign(~MATCH, mul_sign(~~x))
          @rule +(~~x) => setsign(~MATCH, add_sign(~~x))]
    rc = Chain(rs)
    ex = Postwalk(rc)(ex)
    ex = Prewalk(rc)(ex)
    return ex
end

### Curvature ###

setcurvature(ex::Union{Num, Symbolic}, curv) = setmetadata(ex, Curvature, curv)
setcurvature(ex, curv) = ex
getcurvature(ex::Union{Num, Symbolic}) = getmetadata(ex, Curvature)
getcurvature(ex) = Affine
hascurvature(ex::Union{Num, Symbolic}) = hasmetadata(ex, Curvature)
hascurvature(ex) = ex isa Real

function mul_curvature(args)
    # all but one arg is constant
    non_constants = findall(x -> issym(x) || iscall(x), args)
    constants = findall(x -> !issym(x) && !iscall(x), args)
    try
        @assert length(non_constants) <= 1
    catch
        @warn "DCP does not support multiple non-constant arguments in multiplication"
        return UnknownCurvature
    end

    if !isempty(non_constants)
        expr = args[first(non_constants)]
        curv = find_curvature(expr)
        return if getsign(prod(args[constants])) == Negative
            # flip
            curv == Convex ? Concave : curv == Concave ? Convex : curv
        else
            curv
        end
    end
    return Affine
end

function add_curvature(args)
    curvs = find_curvature.(args)
    all(==(Affine), curvs) && return Affine
    all(x -> x == Convex || x == Affine, curvs) && return Convex
    all(x -> x == Concave || x == Affine, curvs) && return Concave
    return UnknownCurvature
end

function propagate_curvature(ex)
    rs = [@rule *(~~x) => setcurvature(~MATCH, mul_curvature(~~x))
          @rule +(~~x) => setcurvature(~MATCH, add_curvature(~~x))
          @rule ~x => setcurvature(~x, find_curvature(~x))]
    rc = Chain(rs)
    ex = Postwalk(rc)(ex)
    ex = Prewalk(rc)(ex)
    # SymbolicUtils.inspect(ex, metadata = true)
    return ex
end

function get_arg_property(monotonicity, i, args)
    @label start
    if monotonicity isa Function
        monotonicity(args[i])
    elseif monotonicity isa Tuple && i <= length(monotonicity)
        monotonicity = monotonicity[i]
        @goto start
    else
        monotonicity
    end
end

function find_curvature(ex)
    if hascurvature(ex)
        return getcurvature(ex)
    end

    if iscall(ex)
        f, args = operation(ex), arguments(ex)
        # @show f
        if hasdcprule(f)
            rule, args = dcprule(f, args...)
        elseif Symbol(f) == :*
            if args[1] isa Number && args[1] > 0
                return find_curvature(args[2])
            elseif args[1] isa Number && args[1] < 0
                argscurv = find_curvature(args[2])
                if argscurv == Convex
                    return Concave
                elseif argscurv == Concave
                    return Convex
                else
                    argscurv
                end
            else
                @warn "DCP does not support multiple non-constant arguments in multiplication"
                return UnknownCurvature
            end
        else
            return UnknownCurvature
        end
        f_curvature = rule.curvature
        f_monotonicity = rule.monotonicity

        if f_curvature == Affine
            if all(enumerate(args)) do (i, arg)
                arg_curv = find_curvature(arg)
                arg_curv == Affine
            end
                return Affine
            end
        elseif f_curvature == Convex || f_curvature == Affine
            if all(enumerate(args)) do (i, arg)
                arg_curv = find_curvature(arg)
                m = get_arg_property(f_monotonicity, i, args)
                # @show f_monotonicity
                # @show arg
                # @show m
                if arg_curv == Convex
                    m == Increasing
                elseif arg_curv == Concave
                    m == Decreasing
                else
                    arg_curv == Affine
                end
            end
                return Convex
            end
        elseif f_curvature == Concave || f_curvature == Affine
            if all(enumerate(args)) do (i, arg)
                arg_curv = find_curvature(arg)
                m = f_monotonicity[i]
                if arg_curv == Concave
                    m == Increasing
                elseif arg_curv == Convex
                    m == Decreasing
                else
                    arg_curv == Affine
                end
            end
                return Concave
            end
        end
        return UnknownCurvature
    else
        return Affine
    end
end
