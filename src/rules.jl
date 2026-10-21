"""
    Sign

Sign classification propagated by the DCP analyzer.

# Values

- `Positive`: the expression is known to be positive.
- `Negative`: the expression is known to be negative.
- `AnySign`: no sign restriction is known.

These values are accepted by [`add_dcprule`](@ref) and are returned by
`AnalysisResult.sign`.
"""
@enum Sign Positive Negative AnySign

"""
    Curvature

Euclidean curvature classification propagated by the DCP analyzer.

# Values

- `Convex`: the expression is convex.
- `Concave`: the expression is concave.
- `Affine`: the expression is both convex and concave.
- `UnknownCurvature`: no supported curvature rule applies.

These values are accepted by [`add_dcprule`](@ref) and are returned by
`AnalysisResult.curvature`.
"""
@enum Curvature Convex Concave Affine UnknownCurvature

"""
    Monotonicity

Monotonicity classification for an argument of a DCP atom.

# Values

- `Increasing`: the atom is nondecreasing in the argument.
- `Decreasing`: the atom is nonincreasing in the argument.
- `AnyMono`: the atom's monotonicity is unrestricted for the rule.

These values are accepted by [`add_dcprule`](@ref).
"""
@enum Monotonicity Increasing Decreasing AnyMono

# Symbolics v7 wraps numeric literals (multiplication coefficients, exponents,
# broadcasted functions, ...) as constant `BasicSymbolic`s, so a plain
# `arg isa Number` test no longer recognises them. `Symbolics.value` collapses a
# wrapped constant back to its underlying value (identity on already-unwrapped
# values), so use it before any `isa Number` dispatch.
constval(x) = Symbolics.value(x)

# Whether an argument is a constant rather than something to optimize over.
# Bilinear atoms (`dot`, `dotsort`, `quad_form`, `huber`'s threshold) are convex
# only in that case, and the static rule table has no way to say so.
# A constant array arrives as a wrapped `Matrix`/`Vector` — neither `issym` nor
# `iscall` — while a symbolic one arrives as an `array_literal` call, so the test
# has to look through `constval` rather than at the node kind.
function isconstarg(x)
    v = constval(x)
    v isa Number && return true
    v isa AbstractArray && return all(y -> constval(y) isa Number, v)
    return false
end

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

"""
    DCPRule

Immutable descriptor for a registered Euclidean DCP atom rule. Fields are
concretely typed so rule-table lookups do not box through `Any`.

`sign` / `curvature` are ordinarily [`Sign`](@ref) / [`Curvature`](@ref); a
small number of atoms (e.g. `perspective`) store callable placeholders instead,
matching the historical NamedTuple storage.
"""
struct DCPRule
    domain::Any
    sign::Union{Sign, Function}
    curvature::Union{Curvature, Function}
    monotonicity::Any
end

# Identity-keyed so function objects used as atom keys resolve without hashing
# through `Any`. Every entry is a `Vector{DCPRule}` — one rule for most atoms,
# several when distinct domains are registered for the same operation.
const dcprules_dict = IdDict{Any, Vector{DCPRule}}()

"""
    add_dcprule(f, domain, sign, curvature, monotonicity)

Register a Euclidean disciplined-convex-programming (DCP) rule for the
symbolic atom `f`. This is the extension point for packages that add a new
symbolic atom whose sign, curvature, and argument monotonicity are known.

Register the symbolic operation with `Symbolics.@register_symbolic` before
registering its DCP rule. The function object passed to `add_dcprule` must be
the same operation stored in the symbolic term.

# Arguments

- `f`: function used as the operation of the symbolic atom.
- `domain`: a `DomainSets.Domain` for a single argument, or a tuple of domains
  for a multi-argument atom. Domain checks use `VarDomain` metadata when it is
  attached to the symbolic arguments.
- `sign::Sign`: sign guaranteed by the atom on its declared domain.
- `curvature::Curvature`: Euclidean curvature of the atom.
- `monotonicity`: one [`Monotonicity`](@ref) value, a tuple with one value per
  argument, or a function that computes the value from an argument. A single
  value applies to every argument.

# Returns

The registered rule descriptor. Registering another rule for the same
function appends a descriptor so that distinct domains can be supported.

# Examples

```jldoctest
julia> using SymbolicAnalysis, Symbolics, DomainSets

julia> test_atom(x) = x;

julia> @register_symbolic test_atom(x)

julia> SymbolicAnalysis.add_dcprule(
           test_atom, RealLine(), SymbolicAnalysis.Positive,
           SymbolicAnalysis.Convex, SymbolicAnalysis.Increasing
       );

julia> @variables x;

julia> analyze(test_atom(x)).curvature == SymbolicAnalysis.Convex
true
```

This registration hook is public for packages extending SymbolicAnalysis. The
rule dictionaries and the `dcprule` dispatch function are implementation
details and should not be accessed directly.
"""
function add_dcprule(f, domain, sign, curvature, monotonicity)
    if !(monotonicity isa Tuple)
        monotonicity = (monotonicity,)
    end
    rule = makerule(domain, sign, curvature, monotonicity)
    if haskey(dcprules_dict, f)
        push!(dcprules_dict[f], rule)
        return dcprules_dict[f]
    else
        return dcprules_dict[f] = DCPRule[rule]
    end
end

function makerule(domain, sign, curvature, monotonicity)
    return DCPRule(domain, sign, curvature, monotonicity)
end

hasdcprule(f::Function) = haskey(dcprules_dict, f)
hasdcprule(f) = false

# Only symbolic leaves can carry a `VarDomain`, and `arguments` always yields
# `BasicSymbolic`s, so a plain type guard replaces the former
# `Symbolics.hasmetadata(::Union{Real, AbstractArray{<:Real}}, args...) = false`
# pirate and is strictly more robust (false for any non-symbolic argument).
_has_vardomain(x) = (x isa Union{Num, Symbolic}) && hasmetadata(x, VarDomain)

# Domain matching for rule selection. `ℂ` is how an atom that is unrestricted on
# the reals is spelled (`abs`, `conj`, `real`, `imag`), but DomainSets cannot
# decide `issubset(::HalfLine, ::ComplexNumbers)` — it answers `false` for some
# element types and throws a `MethodError` for others — so a real argument domain
# always matches it. Any other incomparable pair degrades to "no match".
function subdomain(argdomain, ruledomain)
    ruledomain === ℂ && return true
    return try
        issubset(argdomain, ruledomain)
    catch e
        e isa MethodError ? false : rethrow()
    end
end

# Selected when no registered domain covers the declared argument domains. A
# missing rule is not an error: `analyze` must return a curvature for every
# well-formed expression, and `UnknownCurvature` is the conservative answer.
function no_rule(argsdomain, args)
    return makerule(argsdomain, AnySign, UnknownCurvature, AnyMono), args
end

function dcprule(f, args...)
    rules = dcprules_dict[f]
    if all(_has_vardomain, args)
        argsdomain = getmetadata.(args, Ref(VarDomain))
        for rule in rules
            if (rule.domain isa Domain) &&
                    all(subdomain.(argsdomain, Ref(rule.domain)))
                return rule, args
            elseif !(rule.domain isa Domain) &&
                    all(subdomain.(argsdomain, rule.domain))
                return rule, args
            end
        end
        return no_rule(argsdomain, args)
    else
        return rules[1], args
    end
end

### Sign ###
setsign(ex::Union{Num, Symbolic}, sign) = setmetadata(ex, Sign, sign)
setsign(ex, sign) = ex

function symbolic_sign(ex)::Sign
    # An assembled array aggregates its elements' signs, as `getsign(::AbstractArray)`
    # does for a plain container. This must win over both the fold memo and any
    # Sign metadata: `node_sign` writes the rule-table `AnySign` onto an
    # `array_literal`, but callers (e.g. `norm`'s `increasing_if_positive`) need
    # the elementwise aggregation. Argument 1 is the size tuple.
    if iscall(ex) && operation(ex) === SymbolicUtils.array_literal
        return add_sign(@view arguments(ex)[2:end])
    end
    # During the read-only fold, child results live in the memo rather than in
    # Sign metadata. Prefer the memo so rule logic that calls `getsign` on an
    # already-folded child sees the same value the Postwalk engine wrote.
    hit = fold_lookup(ex)
    hit !== nothing && return hit.sign::Sign
    if hasmetadata(ex, Sign)
        return getmetadata(ex, Sign)::Sign
    end
    return AnySign
end

# Split over `Num` and `Symbolic` instead of their union: `Num <: Real` while
# `BasicSymbolic` is not, so a union method would be ambiguous with the `Real`
# method below rather than more specific than it.
getsign(ex::Num)::Sign = symbolic_sign(unwrap(ex))
getsign(ex::Symbolic)::Sign = symbolic_sign(ex)

# `Real`, not `Union{AbstractFloat, Integer}`: Symbolics folds `-a/2` into
# `(-1//2)*a`, so a coefficient reaching here can be any real subtype.
getsign(ex::Real)::Sign = ex < 0 ? Negative : Positive

function getsign(ex::AbstractArray)::Sign
    if all(x -> getsign(x) == Negative, ex)
        return Negative
    elseif all(x -> getsign(x) == Positive, ex)
        return Positive
    else
        return AnySign
    end
end

function hassign(ex::Union{Num, Symbolic})
    hit = fold_lookup(ex isa Num ? unwrap(ex) : ex)
    hit !== nothing && return true
    return hasmetadata(ex, Sign)
end
hassign(ex) = ex isa Real

hassign(ex::typeof(Base.broadcast)) = true
getsign(ex::typeof(Base.broadcast))::Sign = Positive

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

function _is_scalar_reduction_op(f)
    return f === tr || f === logdet || f === norm || f === opnorm ||
        f === sum || f === prod || f === maximum || f === minimum ||
        f === eigmax || f === eigmin || f === dot || f === dotsort ||
        f === quad_form || f === quad_over_lin || f === matrix_frac ||
        f === trinv || f === tv || f === kldivergence ||
        f === LogExpFunctions.logsumexp || f === invprod ||
        f === StatsBase.geomean || f === StatsBase.harmmean ||
        f === eigsummax || f === eigsummin || f === sum_largest || f === sum_smallest ||
        f === log_quad_form || f === sum_log_eigmax || f === Base.getindex ||
        f === Manifolds.distance || f isa SymbolicUtils.Mapreducer
end

function is_nonscalar_array(x)
    x isa Real && !(x isa Num) && return false
    x isa AbstractArray && return ndims(x) > 0
    u = unwrap(x)
    v = constval(u)
    v isa Number && return false
    v isa AbstractArray && return ndims(v) > 0
    SymbolicUtils.symtype(u) <: AbstractArray && return true
    sh = SymbolicUtils.shape(u)
    sh isa AbstractVector && return !isempty(sh)
    if iscall(u)
        f = operation(u)
        if f === conjugation || f === hadamard_product || f === affine_map ||
                f === scalar_mat || f === lorentz_transform ||
                f === SymbolicUtils.array_literal || f === diag || f === diagm ||
                f === cumsum || f === diff || f === hcat || f === vcat ||
                f === kron || f === reshape || f === tril || f === triu ||
                f === vec || f isa SymbolicUtils.Mapper
            return true
        elseif f === (+) || f === (-) || f === (*) || f === (/) ||
                f === broadcast || f === inv || f === sqrt || f === log ||
                f === exp || f === adjoint || f === transpose
            return any(is_nonscalar_array, arguments(u))
        elseif SymbolicUtils.symtype(u) === Any && !_is_scalar_reduction_op(f)
            return any(is_nonscalar_array, arguments(u))
        end
    end
    return false
end

function _entrywise_domain_sign(dom)::Sign
    if dom isa Domain && !(dom isa CustomDomain)
        subdomain(dom, HalfLine()) && return Positive
        subdomain(dom, Interval{:open, :closed}(-Inf, 0)) && return Negative
        return AnySign
    elseif dom isa CustomDomain && isdefined(dom.in, :element_domain)
        return _entrywise_domain_sign(getfield(dom.in, :element_domain))
    end
    return AnySign
end

function _entrywise_domain_nonnegative(dom)
    return _entrywise_domain_sign(dom) == Positive
end

@inline function _ensure_entrywise_memo(
        memo::Union{Nothing, IdDict{Any, Sign}}, args, start_idx::Int
    )
    memo !== nothing && return memo
    for i in start_idx:length(args)
        if iscall(unwrap(args[i]))
            return IdDict{Any, Sign}()
        end
    end
    return nothing
end

@inline function _safe_even_integer(p::Real)::Bool
    return try
        isfinite(p) && isinteger(p) && iseven(p)
    catch e
        (e isa MethodError || e isa InexactError || e isa DomainError) ? false : rethrow()
    end
end

@inline function _safe_odd_integer(p::Real)::Bool
    return try
        isfinite(p) && isinteger(p) && isodd(p)
    catch e
        (e isa MethodError || e isa InexactError || e isa DomainError) ? false : rethrow()
    end
end

function _combine_entrywise_add(args, start_idx::Int, memo)::Sign
    start_idx > length(args) && return AnySign
    all_pos = true
    all_neg = true
    for i in start_idx:length(args)
        s = _entrywise_sign(args[i], memo)
        if s == Positive
            all_neg = false
        elseif s == Negative
            all_pos = false
        else
            return AnySign
        end
    end
    return all_pos ? Positive : all_neg ? Negative : AnySign
end

function _combine_entrywise_mul(args, start_idx::Int, memo)::Sign
    start_idx > length(args) && return AnySign
    neg_count = 0
    for i in start_idx:length(args)
        s = _entrywise_sign(args[i], memo)
        if s == AnySign
            return AnySign
        elseif s == Negative
            neg_count += 1
        end
    end
    return isodd(neg_count) ? Negative : Positive
end

function _elementwise_call_sign(@nospecialize(bf), args, start_idx::Int, memo)::Sign
    n = length(args) - start_idx + 1
    n <= 0 && return AnySign
    if bf === exp || bf === exp2 || bf === exp10 || bf === abs || bf === abs2 ||
            bf === cosh || bf === hypot || bf === sqrt || bf === LogExpFunctions.log1pexp
        return Positive
    elseif bf === lognormcdf
        return Negative
    elseif bf === (+)
        return _combine_entrywise_add(args, start_idx, memo)
    elseif bf === (-)
        if n == 1
            s = _entrywise_sign(args[start_idx], memo)
            return s == Positive ? Negative : s == Negative ? Positive : AnySign
        elseif n == 2
            s1 = _entrywise_sign(args[start_idx], memo)
            s2 = _entrywise_sign(args[start_idx + 1], memo)
            return (s1 == Positive && s2 == Negative) ? Positive :
                (s1 == Negative && s2 == Positive) ? Negative : AnySign
        end
        return AnySign
    elseif bf === (*)
        return _combine_entrywise_mul(args, start_idx, memo)
    elseif bf === (/) && n == 2
        num = args[start_idx]
        den = args[start_idx + 1]
        dv = constval(den)
        dv isa Number && iszero(dv) && return AnySign
        nv = constval(num)
        nv isa Number && iszero(nv) && !(dv isa Number) && return AnySign
        s1 = _entrywise_sign(num, memo)
        s1 == AnySign && return AnySign
        s2 = _entrywise_sign(den, memo)
        s2 == AnySign && return AnySign
        return s1 == s2 ? Positive : Negative
    elseif bf === inv && n == 1
        return _entrywise_sign(args[start_idx], memo)
    elseif bf === (^) && n == 2
        base = args[start_idx]
        p = constval(args[start_idx + 1])
        if p isa Real
            try
                isfinite(p) || return AnySign
                isone(p) && return _entrywise_sign(base, memo)
                if _safe_even_integer(p) && p > 0
                    return Positive
                end
                sb = _entrywise_sign(base, memo)
                sb == Positive && return Positive
                if _safe_odd_integer(p) && sb == Negative
                    return Negative
                end
                return AnySign
            catch e
                (e isa MethodError || e isa InexactError || e isa DomainError) &&
                    return AnySign
                rethrow()
            end
        else
            bval = constval(base)
            return (bval isa Real && isfinite(bval) && bval > 0) ? Positive : AnySign
        end
    elseif bf === LogExpFunctions.xexpx && n == 1
        return _entrywise_sign(args[start_idx], memo) == Positive ? Positive : AnySign
    elseif bf === huber && n == 2
        m = args[start_idx + 1]
        return (isconstarg(m) && _entrywise_sign(m, memo) == Positive) ? Positive : AnySign
    elseif haskey(dcprules_dict, bf) && !hasmethod(dcprule, Tuple{typeof(bf), Vararg{Any}})
        r = dcprules_dict[bf][1]
        if (r.domain === RealLine() || r.domain === ℂ) && r.sign isa Sign
            return r.sign::Sign
        end
    end
    return AnySign
end

function _entrywise_sign(x, memo::Union{Nothing, IdDict{Any, Sign}} = nothing)::Sign
    x isa Real && !(x isa Num) && return getsign(x)
    x isa AbstractArray && return getsign(x)
    u = unwrap(x)
    v = constval(u)
    v isa Number && return getsign(v)
    v isa AbstractArray && return getsign(v)
    if issym(u)
        if _has_vardomain(u)
            s = _entrywise_domain_sign(getmetadata(u, VarDomain))
            s != AnySign && return s
        end
        return is_nonscalar_array(u) ? AnySign : getsign(u)
    elseif iscall(u)
        if memo !== nothing
            cached = get(memo, u, nothing)
            cached !== nothing && return cached
        end
        f = operation(u)
        args = arguments(u)
        if f === broadcast && !isempty(args)
            bf = constval(args[1])
            if bf === (+) || bf === (-) || bf === (*) || bf === (/)
                memo = _ensure_entrywise_memo(memo, args, 2)
            end
            res = _elementwise_call_sign(bf, args, 2, memo)
        elseif f === (+) || f === hcat || f === vcat
            memo = _ensure_entrywise_memo(memo, args, 1)
            res = _combine_entrywise_add(args, 1, memo)
        elseif f === (-)
            memo = _ensure_entrywise_memo(memo, args, 1)
            res = _elementwise_call_sign(-, args, 1, memo)
        elseif f === vec || f === reshape || f === tril || f === triu ||
                f === diag || f === diagm || f === cumsum ||
                f === adjoint || f === transpose
            res = isempty(args) ? AnySign : _entrywise_sign(args[1], memo)
        elseif f === (*) || f === kron
            memo = _ensure_entrywise_memo(memo, args, 1)
            res = _combine_entrywise_mul(args, 1, memo)
        elseif f === (/) && length(args) == 2
            memo = _ensure_entrywise_memo(memo, args, 1)
            res = !is_nonscalar_array(args[2]) ? _elementwise_call_sign(/, args, 1, memo) : AnySign
        elseif f === SymbolicUtils.array_literal
            memo = _ensure_entrywise_memo(memo, args, 2)
            res = _combine_entrywise_add(args, 2, memo)
        elseif f === Base.getindex && !isempty(args)
            res = _entrywise_sign(args[1], memo)
            if res == AnySign && !is_nonscalar_array(u)
                res = getsign(u)
            end
        elseif f isa SymbolicUtils.Mapper
            res = _elementwise_call_sign(f.f, args, 1, memo)
        else
            res = !is_nonscalar_array(u) ? getsign(u) : AnySign
        end
        if memo !== nothing
            memo[u] = res
        end
        return res
    end
    return AnySign
end

entrywise_sign(x)::Sign = _entrywise_sign(x, nothing)
is_entrywise_nonnegative(x) = entrywise_sign(x) == Positive

function _scalar_factor_sign(arg)::Sign
    arg isa Real && !(arg isa Num) && return getsign(arg)
    v = constval(arg)
    v isa Number && return getsign(v)
    s = getsign(arg)
    s != AnySign && return s
    u = unwrap(arg)
    if issym(u) && _has_vardomain(u)
        return _entrywise_domain_sign(getmetadata(u, VarDomain))
    end
    return AnySign
end

function mul_sign(args)
    # A product with two or more non-scalar array factors is matrix/array
    # multiplication: Loewner-positive (PSD) factors do not have a PSD or
    # entrywise-positive product in general, so only propagate a sign when every
    # array factor is entrywise nonnegative.
    array_count = 0
    neg_count = 0
    first_array_arg = nothing
    for arg in args
        if is_nonscalar_array(arg)
            array_count += 1
            if array_count == 1
                first_array_arg = arg
            elseif array_count == 2
                _entrywise_sign(first_array_arg, nothing) == Positive || return AnySign
                (arg === first_array_arg || _entrywise_sign(arg, nothing) == Positive) ||
                    return AnySign
            else
                (arg === first_array_arg || _entrywise_sign(arg, nothing) == Positive) ||
                    return AnySign
            end
        else
            s = _scalar_factor_sign(arg)
            if s == AnySign
                return AnySign
            elseif s == Negative
                neg_count += 1
            end
        end
    end
    if array_count == 1
        v = constval(first_array_arg)
        if v isa Number || v isa AbstractArray
            s = getsign(v)
        else
            s = getsign(first_array_arg)
            if s == AnySign
                s = entrywise_sign(first_array_arg)
            end
        end
        if s == AnySign
            return AnySign
        elseif s == Negative
            neg_count += 1
        end
    end
    return isodd(neg_count) ? Negative : Positive
end

function broadcast_mul_sign(args)
    array_count = 0
    scalar_neg_count = 0
    first_array_arg = nothing
    for arg in args
        if is_nonscalar_array(arg)
            array_count += 1
            if array_count == 1
                first_array_arg = arg
            end
        else
            s = _scalar_factor_sign(arg)
            if s == AnySign
                return AnySign
            elseif s == Negative
                scalar_neg_count += 1
            end
        end
    end
    if array_count == 0
        return isodd(scalar_neg_count) ? Negative : Positive
    elseif array_count == 1
        s = _entrywise_sign(first_array_arg, nothing)
        if s == AnySign
            s = _loewner_sign(first_array_arg)
        end
        if s == AnySign
            s = getsign(first_array_arg)
        end
        s == AnySign && return AnySign
        s == Negative && (scalar_neg_count += 1)
        return isodd(scalar_neg_count) ? Negative : Positive
    end

    entrywise_neg_count = scalar_neg_count
    all_entrywise = true
    for arg in args
        is_nonscalar_array(arg) || continue
        s = _entrywise_sign(arg, nothing)
        if s == AnySign
            all_entrywise = false
            break
        elseif s == Negative
            entrywise_neg_count += 1
        end
    end
    all_entrywise && return isodd(entrywise_neg_count) ? Negative : Positive

    loewner_neg_count = scalar_neg_count
    for arg in args
        is_nonscalar_array(arg) || continue
        s = _loewner_sign(arg)
        if s == AnySign
            return AnySign
        elseif s == Negative
            loewner_neg_count += 1
        end
    end
    return isodd(loewner_neg_count) ? Negative : Positive
end

function _loewner_domain_sign(dom)::Sign
    if dom isa CustomDomain
        dom.in === isposdef && return Positive
        dom.in == (isposdef ∘ -) && return Negative
    end
    return AnySign
end

function _loewner_sign(x)::Sign
    x isa Real && !(x isa Num) && return AnySign
    if x isa AbstractArray
        if x isa AbstractMatrix{<:Real} && size(x, 1) == size(x, 2) && issymmetric(x)
            isposdef(x) && return Positive
            isposdef(-x) && return Negative
        end
        return AnySign
    end
    u = unwrap(x)
    v = constval(u)
    v isa Number && return AnySign
    if v isa AbstractArray
        if v isa AbstractMatrix{<:Real} && size(v, 1) == size(v, 2) && issymmetric(v)
            isposdef(v) && return Positive
            isposdef(-v) && return Negative
        end
        return AnySign
    end
    is_nonscalar_array(u) || return AnySign
    if issym(u)
        if _has_vardomain(u)
            return _loewner_domain_sign(getmetadata(u, VarDomain))
        end
        return AnySign
    elseif iscall(u)
        f = operation(u)
        args = arguments(u)
        if f === (+)
            isempty(args) && return AnySign
            s = _loewner_sign(args[1])
            s == AnySign && return AnySign
            for i in 2:length(args)
                _loewner_sign(args[i]) == s || return AnySign
            end
            return s
        elseif f === (-)
            if length(args) == 1
                s = _loewner_sign(args[1])
                return s == Positive ? Negative : s == Negative ? Positive : AnySign
            elseif length(args) == 2
                s1 = _loewner_sign(args[1])
                s1 == AnySign && return AnySign
                s2 = _loewner_sign(args[2])
                return s2 == (s1 == Positive ? Negative : Positive) ? s1 : AnySign
            end
            return AnySign
        elseif f === adjoint || f === transpose
            return isempty(args) ? AnySign : _loewner_sign(args[1])
        elseif f === (*) || f === kron
            neg_count = 0
            array_count = 0
            for arg in args
                if is_nonscalar_array(arg)
                    array_count += 1
                    if f === (*) && array_count > 1
                        return AnySign
                    end
                    s = _loewner_sign(arg)
                    s == AnySign && return AnySign
                    s == Negative && (neg_count += 1)
                else
                    s = _scalar_factor_sign(arg)
                    s == AnySign && return AnySign
                    s == Negative && (neg_count += 1)
                end
            end
            array_count == 0 && return AnySign
            return isodd(neg_count) ? Negative : Positive
        elseif f === (/) && length(args) == 2
            if is_nonscalar_array(args[1]) && !is_nonscalar_array(args[2])
                s1 = _loewner_sign(args[1])
                s1 == AnySign && return AnySign
                s2 = _scalar_factor_sign(args[2])
                s2 == AnySign && return AnySign
                return s1 == s2 ? Positive : Negative
            end
            return AnySign
        elseif f === broadcast && !isempty(args)
            bf = constval(args[1])
            if bf === (+) && length(args) >= 2
                s = _loewner_sign(args[2])
                s == AnySign && return AnySign
                for i in 3:length(args)
                    _loewner_sign(args[i]) == s || return AnySign
                end
                return s
            elseif bf === (-)
                if length(args) == 2
                    s = _loewner_sign(args[2])
                    return s == Positive ? Negative : s == Negative ? Positive : AnySign
                elseif length(args) == 3
                    s1 = _loewner_sign(args[2])
                    s1 == AnySign && return AnySign
                    s2 = _loewner_sign(args[3])
                    return s2 == (s1 == Positive ? Negative : Positive) ? s1 : AnySign
                end
                return AnySign
            elseif bf === (*) && length(args) >= 2
                neg_count = 0
                array_count = 0
                for i in 2:length(args)
                    arg = args[i]
                    if is_nonscalar_array(arg)
                        array_count += 1
                        s = _loewner_sign(arg)
                        s == AnySign && return AnySign
                        s == Negative && (neg_count += 1)
                    else
                        s = _scalar_factor_sign(arg)
                        s == AnySign && return AnySign
                        s == Negative && (neg_count += 1)
                    end
                end
                array_count == 0 && return AnySign
                return isodd(neg_count) ? Negative : Positive
            elseif bf === (/) && length(args) == 3
                if is_nonscalar_array(args[2]) && !is_nonscalar_array(args[3])
                    s1 = _loewner_sign(args[2])
                    s1 == AnySign && return AnySign
                    s2 = _scalar_factor_sign(args[3])
                    s2 == AnySign && return AnySign
                    return s1 == s2 ? Positive : Negative
                end
                return AnySign
            end
            return AnySign
        elseif f === vec || f === reshape || f === tril || f === triu ||
                f === diag || f === diagm || f === cumsum || f === diff ||
                f === hcat || f === vcat || f === SymbolicUtils.array_literal ||
                f === Base.getindex || f isa SymbolicUtils.Mapper
            return AnySign
        else
            return hasgdcprule(f) ? getsign(u) : AnySign
        end
    end
    return AnySign
end

# Sign of a single node whose children have already been annotated by the
# bottom-up walk. Mirrors the priority the old rule chain established through
# successive overwrites: `*`/`+` aggregation wins over a GDCP rule, which wins
# over a DCP rule, which wins over the `AnySign` default (also overwriting any
# stale metadata from a previous analysis).
#
# A GDCP rule's sign only holds on its manifold — `eigmax` is registered
# `Positive` because every eigenvalue of an SPD matrix is — so it is consulted
# only when manifold analysis was actually requested. Without the gate it leaked
# into Euclidean analysis, where `eigmax(X)` for an unconstrained symmetric `X`
# has no sign; that made `abs`'s `increasing_if_positive` resolve to `Increasing`
# and certified `abs(eigmax(X))` as `Convex`, which is false.
function node_sign(ex, M = nothing)
    usegdcp = !isnothing(M)
    if iscall(ex)
        f = operation(ex)
        if f === (*)
            return mul_sign(arguments(ex))
        elseif f === (+)
            return add_sign(arguments(ex))
        elseif f === broadcast && !isempty(arguments(ex)) &&
                constval(arguments(ex)[1]) === (*)
            return broadcast_mul_sign(@view arguments(ex)[2:end])
        elseif usegdcp && hasgdcprule(f)
            return gdcprule(f, arguments(ex)...)[1].sign
        elseif hasdcprule(f)
            return dcprule(f, arguments(ex)...)[1].sign
        end
    elseif issym(ex)
        if usegdcp && hasgdcprule(ex)
            return gdcprule(ex)[1].sign
        elseif hasdcprule(ex)
            return dcprule(ex)[1].sign
        end
    end
    return AnySign
end

function _propagate_sign_walk(x, M)
    if iscall(x)
        orig_type = SymbolicUtils.symtype(x)
        orig_shape = SymbolicUtils.shape(x)
        args = arguments(x)
        new_args = nothing
        for i in eachindex(args)
            arg = args[i]
            newarg = _propagate_sign_walk(arg, M)
            if new_args === nothing
                arg === newarg && continue
                new_args = Any[args[j] for j in eachindex(args)]
            end
            new_args[i] = newarg
        end
        if new_args !== nothing
            f = operation(x)
            m = metadata(x)
            x = maketerm(typeof(x), f, new_args, m; type = orig_type)
            if SymbolicUtils.symtype(x) !== orig_type || SymbolicUtils.shape(x) != orig_shape
                x = SymbolicUtils.term(
                    f, new_args...; type = orig_type, shape = orig_shape
                )
                if m isa AbstractDict
                    for (k, v) in m
                        x = setmetadata(x, k, v)
                    end
                end
            end
        end
        return setsign(x, node_sign(x, M))
    elseif issym(x)
        return setsign(x, node_sign(x, M))
    else
        return x
    end
end

function propagate_sign(ex, M = nothing)
    # Work on the raw symbolic so the sign metadata survives the walk on
    # Symbolics v7 (a `Num`/`Arr` wrapper round-trips through wrap/unwrap and drops
    # it). `analyze` already unwraps; mirror that for direct callers.
    ex = SymbolicUtils.unwrap(ex)
    return _propagate_sign_walk(ex, M)
end

### Curvature ###

setcurvature(ex::Union{Num, Symbolic}, curv) = setmetadata(ex, Curvature, curv)
setcurvature(ex, curv) = ex
function getcurvature(ex::Union{Num, Symbolic})::Curvature
    v = unwrap(ex)
    hit = fold_lookup(v)
    hit !== nothing && return hit.curvature
    hasmetadata(ex, Curvature) && return getmetadata(ex, Curvature)::Curvature
    # A constant-folded expression (`a - a`, `0*a`, `Num(3.0)`) is a
    # `BasicSymbolic` that is neither `issym` nor `iscall`, so the propagation
    # walk never annotates it; it is a constant, hence affine. Anything else
    # arriving here was never analyzed, and carries no certificate.
    return (issym(v) || iscall(v)) ? UnknownCurvature : Affine
end
getcurvature(ex)::Curvature = Affine

# A container is as curved as its elements, so aggregate rather than falling
# through to the `Affine` default for non-symbolic values: a plain
# `Vector{Num}`/`Matrix{Num}` of non-affine entries is not affine. Mirrors
# `getsign(::AbstractArray)`; elementwise curvatures combine as in a sum, since
# every element must hold the claimed curvature for the container to.
function getcurvature(ex::AbstractArray)::Curvature
    has_convex = false
    has_concave = false
    for e in ex
        curv = find_curvature(e)
        if curv == Affine
            continue
        elseif curv == Convex
            has_concave && return UnknownCurvature
            has_convex = true
        elseif curv == Concave
            has_convex && return UnknownCurvature
            has_concave = true
        else
            return UnknownCurvature
        end
    end
    return has_convex ? Convex : has_concave ? Concave : Affine
end
function hascurvature(ex::Union{Num, Symbolic})
    v = unwrap(ex)
    hit = fold_lookup(v)
    hit !== nothing && return true
    return hasmetadata(ex, Curvature)
end
hascurvature(ex) = ex isa Real

# Sign of the product of the constant factors of a multiplication, i.e. the
# direction in which the product is monotone in its one non-constant factor.
# `mul_sign` cannot answer this: it reads sign *metadata*, which a wrapped
# constant array never carries.
function const_factor_sign(args)
    constant_prod = one(Float64)
    for arg in args
        (issym(arg) || iscall(arg)) && continue
        constant_prod = constant_prod .* constval(arg)
    end
    return getsign(constant_prod)
end

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
        # A linear map keeps an affine argument affine whatever the coefficients
        # are, but it preserves a curvature only when every coefficient shares a
        # sign: a mixed-sign constant matrix sums convex and concave terms.
        curv == Affine && return Affine
        constant_sign = getsign(constant_prod)
        return if constant_sign == Negative
            # flip
            curv == Convex ? Concave : curv == Concave ? Convex : curv
        elseif constant_sign == Positive
            curv
        else
            UnknownCurvature
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
        if f === (*)
            return mul_curvature(arguments(ex))
        elseif f === (+)
            return add_curvature(arguments(ex))
        end
    end
    return find_curvature(ex)
end

function propagate_curvature(ex)
    # See `propagate_sign`: unwrap so the curvature metadata survives the v7 walk,
    # and annotate symbols and calls in a single bottom-up pass.
    ex = SymbolicUtils.unwrap(ex)
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

function find_curvature(ex)::Curvature
    if hascurvature(ex)
        return getcurvature(ex)
    end

    # A `Num` is not itself a call, so without unwrapping it would fall through to
    # the `Affine` default below and certify e.g. `exp(x)` as affine.
    if ex isa Num
        return find_curvature(Symbolics.unwrap(ex))
    end

    if iscall(ex)
        f, args = operation(ex), arguments(ex)
        # @show f
        if hasdcprule(f)
            rule, args = dcprule(f, args...)
        elseif f === (*)
            # `mul_curvature` handles a constant in any position and of any shape;
            # the scalar-only version this replaced read `constval(args[1])` and
            # required it to be a `Number`, so a constant *matrix* coefficient was
            # rejected and `A*x - b` did not certify while `A*x .- b` did.
            return mul_curvature(args)
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
