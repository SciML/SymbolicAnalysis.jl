@enum Sign Positive Negative AnySign
@enum Curvature Convex Concave Affine UnknownCurvature
@enum Monotonicity Increasing Decreasing AnyMono

# Symbolics v7 wraps numeric literals (multiplication coefficients, exponents,
# broadcasted functions, ...) as constant `BasicSymbolic`s, so a plain
# `arg isa Number` test no longer recognises them. `Symbolics.value` collapses a
# wrapped constant back to its underlying value (identity on already-unwrapped
# values), so use it before any `isa Number` dispatch.
constval(x) = Symbolics.value(x)

struct CustomDomain{T} <: Domain{T}
    in::Function
end

Base.in(x, c::CustomDomain) = c.in(x)
# Disambiguate against Symbolics' `in(::Num/::Symbolic, ::Domain)`, since
# `CustomDomain <: Domain` makes the symbolic-variable calls match both methods.
# `InDomainSymbolic` matches Symbolics' own dispatch type so these remain strictly
# more specific (no ambiguity).
Base.in(x::Union{Symbolics.Num, InDomainSymbolic}, c::CustomDomain) = c.in(x)
Base.in(x::NTuple{N, Union{Symbolics.Num, InDomainSymbolic}}, c::CustomDomain) where {N} = c.in(x)

function array_domain(element_domain)
    return CustomDomain{AbstractArray}() do xs
        all(in(element_domain), xs)
    end
end

function array_domain(element_domain, N)
    return CustomDomain{AbstractArray{<:Any, N}}() do xs
        ndims(xs) == N && all(in(element_domain), xs)
    end
end

function symmetric_domain()
    return CustomDomain{AbstractArray{<:Any, 2}}(issymmetric)
end

function semidefinite_domain()
    return CustomDomain{AbstractArray{<:Any, 2}}(isposdef) #not semi so needs to change
end

function negsemidefinite_domain()
    return CustomDomain{AbstractArray{<:Any, 2}}(isposdef ∘ -) #not semi so needs to change
end

function definite_domain()
    return CustomDomain{AbstractArray{<:Any, 2}}(isposdef)
end

function negdefinite_domain()
    return CustomDomain{AbstractArray{<:Any, 2}}(isposdef ∘ -)
end

function function_domain()
    return CustomDomain{Function}(x -> typeassert(x, Function))
end

function increasing_if_positive(x)
    sign = getsign(x)
    return sign == AnySign ? AnyMono : sign == Positive ? Increasing : Decreasing
end

const dcprules_dict = Dict()

function add_dcprule(f, domain, sign, curvature, monotonicity)
    if !(monotonicity isa Tuple)
        monotonicity = (monotonicity,)
    end
    return if f in keys(dcprules_dict)
        dcprules_dict[f] = vcat(dcprules_dict[f], makerule(domain, sign, curvature, monotonicity))
    else
        dcprules_dict[f] = makerule(domain, sign, curvature, monotonicity)
    end
end

function makerule(domain, sign, curvature, monotonicity)
    return (; domain = domain, sign = sign, curvature = curvature, monotonicity = monotonicity)
end

hasdcprule(f::Function) = haskey(dcprules_dict, f)
hasdcprule(f) = false

# Only symbolic leaves can carry a `VarDomain`, and `arguments` always yields
# `BasicSymbolic`s, so a plain type guard replaces the former
# `Symbolics.hasmetadata(::Union{Real, AbstractArray{<:Real}}, args...) = false`
# pirate and is strictly more robust (false for any non-symbolic argument).
_has_vardomain(x) = (x isa Union{Num, Symbolic}) && hasmetadata(x, VarDomain)

function dcprule(f, args...)
    if all(_has_vardomain, args)
        argsdomain = getmetadata.(args, Ref(VarDomain))
    else
        if dcprules_dict[f] isa Vector
            return dcprules_dict[f][1], args
        else
            return dcprules_dict[f], args
        end
    end

    if dcprules_dict[f] isa Vector
        for rule in dcprules_dict[f]
            if (rule.domain isa Domain) &&
                    all(issubset.(argsdomain, Ref(rule.domain)))
                return rule, args
            elseif !(rule.domain isa Domain) &&
                    all(issubset.(argsdomain, rule.domain))
                return rule, args
            end
        end
        throw(
            ArgumentError(
                "No DCP rule found for $f with arguments $args with domain $argsdomain",
            ),
        )
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

function add_sign(args)
    if hassign(args)
        return getsign(args)
    end
    # Check signs without allocating intermediate arrays
    has_anysign = false
    all_negative = true
    all_positive = true
    for i in eachindex(args)
        arg = args[i]
        # The bottom-up pass annotates children before their parent, so each
        # argument already carries its sign metadata here.
        s = getsign(arg)
        if s == AnySign
            has_anysign = true
            break
        elseif s == Negative
            all_positive = false
        elseif s == Positive
            all_negative = false
        else
            all_negative = false
            all_positive = false
        end
    end
    return if has_anysign
        AnySign
    elseif all_negative
        Negative
    elseif all_positive
        Positive
    else
        AnySign
    end
end

function mul_sign(args)
    # Avoid allocating intermediate arrays
    neg_count = 0
    for arg in args
        v = constval(arg)
        s = v isa Number ? getsign(v) : getsign(arg)
        if s == AnySign
            return AnySign
        elseif s == Negative
            neg_count += 1
        end
    end
    return isodd(neg_count) ? Negative : Positive
end

# Sign of a single node whose children have already been annotated by the
# bottom-up walk. Mirrors the priority the old rule chain established through
# successive overwrites: `*`/`+` aggregation wins over a GDCP rule, which wins
# over a DCP rule, which wins over the `AnySign` default (also overwriting any
# stale metadata from a previous analysis).
function node_sign(ex)
    if iscall(ex)
        f = operation(ex)
        if Symbol(f) == :*
            return mul_sign(arguments(ex))
        elseif Symbol(f) == :+
            return add_sign(arguments(ex))
        elseif hasgdcprule(f)
            return gdcprule(f, arguments(ex)...)[1].sign
        elseif hasdcprule(f)
            return dcprule(f, arguments(ex)...)[1].sign
        end
    elseif issym(ex)
        if hasgdcprule(ex)
            return gdcprule(ex)[1].sign
        elseif hasdcprule(ex)
            return dcprule(ex)[1].sign
        end
    end
    return AnySign
end

function propagate_sign(ex)
    # Work on the raw symbolic so the sign metadata survives the walk on
    # Symbolics v7 (a `Num`/`Arr` wrapper round-trips through wrap/unwrap and drops
    # it). `analyze` already unwraps; mirror that for direct callers.
    ex = Symbolics.unwrap(ex)
    # A single bottom-up walk: `Postwalk` rebuilds each node from its annotated
    # children before the annotation function sees it, so every node's sign is
    # computed exactly once. Only symbols and calls are annotated — wrapped
    # constants can't carry metadata on older SymbolicUtils v4 releases, and
    # every getter already falls back correctly for them.
    return Postwalk(x -> issym(x) || iscall(x) ? setsign(x, node_sign(x)) : x)(ex)
end

### Curvature ###

setcurvature(ex::Union{Num, Symbolic}, curv) = setmetadata(ex, Curvature, curv)
setcurvature(ex, curv) = ex
getcurvature(ex::Union{Num, Symbolic}) = getmetadata(ex, Curvature)
getcurvature(ex) = Affine
hascurvature(ex::Union{Num, Symbolic}) = hasmetadata(ex, Curvature)
hascurvature(ex) = ex isa Real

function mul_curvature(args)
    # Avoid allocations by not using findall
    # All but one arg must be constant
    non_constant_expr = nothing
    non_constant_count = 0
    constant_prod = one(Float64)
    for arg in args
        if issym(arg) || iscall(arg)
            non_constant_count += 1
            non_constant_expr = arg
            if non_constant_count > 1
                @warn "DCP does not support multiple non-constant arguments in multiplication"
                return UnknownCurvature
            end
        else
            constant_prod *= constval(arg)
        end
    end

    if non_constant_expr !== nothing
        curv = find_curvature(non_constant_expr)
        return if getsign(constant_prod) == Negative
            # flip
            curv == Convex ? Concave : curv == Concave ? Convex : curv
        else
            curv
        end
    end
    return Affine
end

function add_curvature(args)
    # Avoid allocating intermediate arrays - check curvatures in one pass
    has_convex = false
    has_concave = false
    for arg in args
        curv = find_curvature(arg)
        if curv == Affine
            continue
        elseif curv == Convex
            has_convex = true
            if has_concave
                return UnknownCurvature
            end
        elseif curv == Concave
            has_concave = true
            if has_convex
                return UnknownCurvature
            end
        else
            return UnknownCurvature
        end
    end
    if has_convex
        return Convex
    elseif has_concave
        return Concave
    else
        return Affine
    end
end

# Curvature of a single node whose children have already been annotated by the
# bottom-up walk (so the `find_curvature` recursion terminates at the children's
# cached metadata).
function node_curvature(ex)
    if iscall(ex)
        f = operation(ex)
        if Symbol(f) == :*
            return mul_curvature(arguments(ex))
        elseif Symbol(f) == :+
            return add_curvature(arguments(ex))
        end
    end
    return find_curvature(ex)
end

function propagate_curvature(ex)
    # See `propagate_sign`: unwrap so the curvature metadata survives the v7 walk,
    # and annotate symbols and calls in a single bottom-up pass.
    ex = Symbolics.unwrap(ex)
    return Postwalk(x -> issym(x) || iscall(x) ? setcurvature(x, node_curvature(x)) : x)(ex)
end

function get_arg_property(monotonicity, i, args)
    @label start
    return if monotonicity isa Function
        monotonicity(args[i])
    elseif monotonicity isa Tuple
        # A rule may declare fewer monotonicities than the call has arguments —
        # in particular `add_dcprule` wraps a scalar monotonicity into a
        # single-entry tuple that applies to every argument — so the last
        # declared entry covers the remaining arguments.
        monotonicity = monotonicity[min(i, length(monotonicity))]
        @goto start
    else
        monotonicity
    end
end

# The DCP composition rule: an atom of curvature `target` keeps that curvature
# when every argument is affine, matches `target` under an increasing slot, or
# matches the flipped curvature under a decreasing slot.
function composes_as(target::Curvature, f_monotonicity, args)
    flipped = target == Convex ? Concave : Convex
    return all(enumerate(args)) do (i, arg)
        arg_curv = find_curvature(arg)
        arg_curv == Affine && return true
        m = get_arg_property(f_monotonicity, i, args)
        if arg_curv == target
            m == Increasing
        elseif arg_curv == flipped
            m == Decreasing
        else
            false
        end
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
            a1 = constval(args[1])
            if a1 isa Number && a1 > 0
                return find_curvature(args[2])
            elseif a1 isa Number && a1 < 0
                argscurv = find_curvature(args[2])
                if argscurv == Convex
                    return Concave
                elseif argscurv == Concave
                    return Convex
                else
                    return argscurv
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

        if f_curvature == Convex
            return composes_as(Convex, f_monotonicity, args) ? Convex : UnknownCurvature
        elseif f_curvature == Concave
            return composes_as(Concave, f_monotonicity, args) ? Concave : UnknownCurvature
        elseif f_curvature == Affine
            # An affine atom composes both as a convex and as a concave one; it
            # stays affine only when both hold, i.e. every argument is affine.
            cvx = composes_as(Convex, f_monotonicity, args)
            ccv = composes_as(Concave, f_monotonicity, args)
            return cvx && ccv ? Affine : cvx ? Convex : ccv ? Concave : UnknownCurvature
        end
        return UnknownCurvature
    else
        return Affine
    end
end
