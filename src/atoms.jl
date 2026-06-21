### DCP atom rules

add_dcprule(+, RealLine(), AnySign, Affine, Increasing)
add_dcprule(-, RealLine(), AnySign, Affine, Decreasing)

add_dcprule(Base.Ref, RealLine(), AnySign, Affine, AnyMono)

add_dcprule(
    dot,
    (array_domain(RealLine()), array_domain(RealLine())),
    AnySign,
    Affine,
    Increasing
)

"""
    dotsort(x, y)

Sorts `x` and `y` and returns the dot product of the sorted vectors.

# Arguments

    - `x::AbstractVector`: A vector.
    - `y::AbstractVector`: A vector.
"""
function dotsort(x::AbstractVector, y::AbstractVector)
    if length(x) != length(y)
        throw(DimensionMismatch("AbstractVectors must have same length"))
    end
    return dot(sort(x), sort(y))
end
Symbolics.@register_symbolic dotsort(x::AbstractVector, y::AbstractVector)
add_dcprule(
    dotsort,
    (array_domain(RealLine(), 1), array_domain(RealLine(), 1)),
    AnySign,
    Convex,
    (AnyMono, increasing_if_positive ∘ minimum)
)

add_dcprule(
    StatsBase.geomean,
    array_domain(HalfLine{Real, :open}(), 1),
    Positive,
    Concave,
    Increasing
)
add_dcprule(
    StatsBase.harmmean,
    array_domain(HalfLine{Real, :open}(), 1),
    Positive,
    Concave,
    Increasing
)

"""
    invprod(x::AbstractVector)

Returns the inverse of the product of the elements of `x`.

# Arguments

    - `x::AbstractVector`: A vector.
"""
function invprod(x::AbstractVector)
    if any(iszero(x))
        throw(DivideError())
    end
    return inv(prod(x))
end
Symbolics.@register_symbolic invprod(x::AbstractVector)

add_dcprule(invprod, array_domain(HalfLine{Real, :open}()), Positive, Convex, Decreasing)

add_dcprule(eigmax, symmetric_domain(), AnySign, Convex, AnyMono)

add_dcprule(eigmin, symmetric_domain(), AnySign, Concave, AnyMono)

"""
    eigsummax(m::Symmetric, k)

Returns the sum of the `k` largest eigenvalues of `m`.

# Arguments

    - `m::Symmetric`: A symmetric matrix.
    - `k::Int`: The number of largest eigenvalues to sum.
"""
function eigsummax(m::Symmetric, k::Int)
    if k < 1 || k > size(m, 1)
        throw(DomainError(k, "k must be between 1 and size(m, 1)"))
    end
    nrows = size(m, 1)
    return sum(eigvals(m, (nrows - k + 1):nrows))
end
Symbolics.@register_symbolic eigsummax(m::Matrix, k::Int)
add_dcprule(eigsummax, (array_domain(RealLine(), 2), RealLine()), AnySign, Convex, AnyMono)

"""
    eigsummin(m::Symmetric, k)

Returns the sum of the `k` smallest eigenvalues of `m`.

# Arguments

    - `m::Symmetric`: A symmetric matrix.
    - `k::Int`: The number of smallest eigenvalues to sum.
"""
function eigsummin(m::Symmetric, k::Int)
    if k < 1 || k > size(m, 1)
        throw(DomainError(k, "k must be between 1 and size(m, 1)"))
    end
    return sum(eigvals(m, 1:k))
end
Symbolics.@register_symbolic eigsummin(m::Matrix, k::Int)
add_dcprule(eigsummin, (array_domain(RealLine(), 2), RealLine()), AnySign, Concave, AnyMono)

add_dcprule(logdet, semidefinite_domain(), AnySign, Concave, AnyMono)

add_dcprule(
    LogExpFunctions.logsumexp,
    array_domain(RealLine(), 2),
    AnySign,
    Convex,
    Increasing
)

"""
    matrix_frac(x::AbstractVector, P::AbstractMatrix)

Returns the quadratic form `x' * P^{-1} * x`.

# Arguments

    - `x::AbstractVector`: A vector.
    - `P::AbstractMatrix`: A matrix.
"""
function matrix_frac(x::AbstractVector, P::AbstractMatrix)
    if length(x) != size(P, 1)
        throw(DimensionMismatch("x and P must have same length"))
    end
    return x' * inv(P) * x
end
Symbolics.@register_symbolic AbstractMatrix_frac(x::AbstractVector, P::AbstractMatrix)
add_dcprule(
    matrix_frac,
    (array_domain(RealLine(), 1), definite_domain()),
    AnySign,
    Convex,
    AnyMono
)

add_dcprule(maximum, array_domain(RealLine()), AnySign, Convex, Increasing)

add_dcprule(minimum, array_domain(RealLine()), AnySign, Concave, Increasing)

#incorrect for p<1
add_dcprule(
    norm,
    (array_domain(RealLine()), Interval{:closed, :open}(1, Inf)),
    Positive,
    Convex,
    increasing_if_positive
)
add_dcprule(
    norm,
    (array_domain(RealLine()), Interval{:closed, :open}(0, 1)),
    Positive,
    Convex,
    increasing_if_positive
)

"""
    perspective(f::Function, x, s::Real)

Returns the perspective function `s * f(x / s)`.

# Arguments

    - `f::Function`: A function.
    - `x`: A Real.
    - `s::Real`: A positive Real.
"""
function perspective(f::Function, x, s::Real)
    if s < 0
        throw(DomainError(s, "s must be positive"))
    end
    if s == 0
        return zero(typeof(f(x)))
    end
    return s * f(x / s)
end
Symbolics.@register_symbolic perspective(f::Function, x, s::Real)
add_dcprule(
    perspective,
    (function_domain(), RealLine(), Positive),
    getsign,
    getcurvature,
    AnyMono
)

"""
    quad_form(x::AbstractVector, P::AbstractMatrix)

Returns the quadratic form `x' * P * x`.

# Arguments

    - `x::AbstractVector`: A vector.
    - `P::AbstractMatrix`: A matrix.
"""
function quad_form(x::AbstractVector, P::AbstractMatrix)
    if length(x) != size(P, 1)
        throw(DimensionMismatch("x and P must have same length"))
    end
    return x' * P * x
end
Symbolics.@register_symbolic quad_form(x::AbstractVector, P::AbstractMatrix)
add_dcprule(
    quad_form,
    (array_domain(RealLine(), 1), semidefinite_domain()),
    Positive,
    Convex,
    (increasing_if_positive, Increasing)
)

function quad_over_lin(x::AbstractVector{<:Real}, y::Real)
    if getsign(y) == Negative
        throw(DomainError(y, "y must be positive"))
    end
    return sum(x .^ 2) / y
end

# On Symbolics v7 both scalar- and array-valued registered symbolics share the
# concrete type `BasicSymbolic{SymReal}`, so this vector registration would
# generate the same `quad_over_lin(::BasicSymbolic{SymReal}, ::Real)` method as
# the scalar one below and collide. Only register it on v6, where vector and
# scalar symbolics are distinct `BasicSymbolic` types.
if pkgversion(Symbolics) < v"7"
    Symbolics.@register_symbolic quad_over_lin(x::AbstractVector, y::Real) false
end

"""
    quad_over_lin(x::Real, y::Real)

Returns the quadratic over linear form `x^2 / y`.

# Arguments

    - `x`: A Real or a vector.
    - `y::Real`: A positive Real.
"""
function quad_over_lin(x::Real, y::Real)
    if getsign(y) == Negative
        throw(DomainError(y, "y must be positive"))
    end
    return x^2 / y
end

Symbolics.@register_symbolic quad_over_lin(x::Real, y::Real)

add_dcprule(
    quad_over_lin,
    (array_domain(RealLine()), HalfLine{Real, :open}()),
    Positive,
    Convex,
    (increasing_if_positive, Decreasing)
)

add_dcprule(
    quad_over_lin,
    (RealLine(), HalfLine{Real, :open}()),
    Positive,
    Convex,
    (increasing_if_positive, Decreasing)
)

add_dcprule(sum, array_domain(RealLine(), 2), AnySign, Affine, Increasing)

"""
    sum_largest(x::AbstractMatrix, k)

Returns the sum of the `k` largest elements of `x`.

# Arguments

    - `x::AbstractMatrix`: A matrix.
    - `k::Int`: The number of largest elements to sum.
"""
function sum_largest(x::AbstractMatrix, k::Integer)
    return sum(sort(vec(x))[(end - k):end])
end
Symbolics.@register_symbolic sum_largest(x::AbstractMatrix, k::Integer)
add_dcprule(sum_largest, (array_domain(RealLine(), 2), ℤ), AnySign, Convex, Increasing)

"""
    sum_smallest(x::AbstractMatrix, k)

Returns the sum of the `k` smallest elements of `x`.

# Arguments

    - `x::AbstractMatrix`: A matrix.
    - `k::Int`: The number of smallest elements to sum.
"""
function sum_smallest(x::AbstractMatrix, k::Integer)
    return sum(sort(vec(x))[1:k])
end

Symbolics.@register_symbolic sum_smallest(x::AbstractArray, k::Integer)
add_dcprule(sum_smallest, (array_domain(RealLine(), 2), ℤ), AnySign, Concave, Increasing)

add_dcprule(tr, array_domain(RealLine(), 2), AnySign, Affine, Increasing)

"""
    trinv(x::AbstractMatrix)

Returns the trace of the inverse of `x`.

# Arguments

    - `x::AbstractMatrix`: A matrix.
"""
function trinv(x::AbstractMatrix)
    return tr(inv(x))
end
Symbolics.@register_symbolic trinv(x::AbstractMatrix)
add_dcprule(trinv, definite_domain(), Positive, Convex, AnyMono)

"""
    tv(x::AbstractVector{<:Real})

Returns the total variation of `x`, defined as `sum_i |x_{i+1} - x_i|`.

# Arguments

    - `x::AbstractVector`: A vector.
"""
function tv(x::AbstractVector{<:Real})
    return sum(abs.(x[2:end] - x[1:(end - 1)]))
end
Symbolics.@register_symbolic tv(x::AbstractVector) false
add_dcprule(tv, array_domain(RealLine(), 1), Positive, Convex, AnyMono)

"""
    tv(x::AbstractVector{<:AbstractMatrix})

Returns the total variation of `x`, defined as `sum_{i,j} |x_{k+1}[i,j] - x_k[i,j]|`.

# Arguments

    - `x::AbstractVector`: A vector of matrices.
"""
function tv(x::AbstractVector{<:AbstractMatrix})
    return sum(
        map(1:(size(x, 1) - 1)) do i
            map(1:(size(x, 2) - 1)) do j
                norm([x[k][i + 1, j] - x[k][i, j] for k in eachindex(x)])
            end
        end
    )
end
add_dcprule(tv, array_domain(array_domain(RealLine(), 2), 1), Positive, Convex, AnyMono)

add_dcprule(abs, ℂ, Positive, Convex, increasing_if_positive)

add_dcprule(conj, ℂ, AnySign, Affine, AnyMono)

add_dcprule(exp, RealLine(), Positive, Convex, Increasing)

Symbolics.@register_symbolic LogExpFunctions.xlogx(x::Real)
add_dcprule(xlogx, RealLine(), AnySign, Convex, AnyMono)

"""
    huber(x, M=1)

Returns the Huber loss function of `x` with threshold `M`.

# Arguments

    - `x::Real`: A Real.
    - `M::Real`: The threshold.
"""
function huber(x::Real, M::Real = 1)
    if M < 0
        throw(DomainError(M, "M must be positive"))
    end

    if abs(x) <= M
        return x^2
    else
        return 2 * M * abs(x) - M^2
    end
end
Symbolics.@register_symbolic huber(x::Real, M::Real)
add_dcprule(huber, (RealLine(), HalfLine()), Positive, Convex, increasing_if_positive)

add_dcprule(imag, ℂ, AnySign, Affine, AnyMono)

add_dcprule(inv, HalfLine{Real, :open}(), Positive, Convex, Decreasing)
add_dcprule(log, HalfLine{Real, :open}(), AnySign, Concave, Increasing)

# Matrix-valued atoms (`log`, `inv`, `sqrt` of a symbolic matrix). On Symbolics
# v7 / SymbolicUtils v4 a symbolic matrix unwraps to `BasicSymbolic{SymReal}` —
# the same concrete type as a symbolic scalar — so the `@register_symbolic`
# macros would emit a `f(::BasicSymbolic{SymReal})` method that overwrites
# SymbolicUtils' own scalar `f` and aborts precompilation. Build the matrix term
# directly off the `Arr` wrapper instead; `SymbolicUtils.term(f, x; type)` is the
# stable constructor on both v6 and v7 and leaves the scalar methods untouched.
# On Symbolics v7 some scalar atoms (e.g. `sqrt`) register a `promote_shape` that
# rejects matrix shapes, which `term` would invoke unless the shape is supplied up
# front. v6's `term` does not take a `shape` kwarg (and doesn't need it), so only
# pass it on v7.
@static if pkgversion(Symbolics) < v"7"
    function matrix_atom(f, A::Symbolics.Arr)
        a = Symbolics.unwrap(A)
        return Symbolics.wrap(SymbolicUtils.term(f, a; type = SymbolicUtils.symtype(a)))
    end
    function matrix_atom(f, A::AbstractMatrix{<:Num})
        a = Symbolics.unwrap.(A)
        return Symbolics.wrap(SymbolicUtils.term(f, a; type = Matrix{Real}))
    end
else
    function matrix_atom(f, A::Symbolics.Arr)
        a = Symbolics.unwrap(A)
        return Symbolics.wrap(
            SymbolicUtils.term(
                f, a; type = SymbolicUtils.symtype(a), shape = Symbolics.shape(a)
            )
        )
    end
    function matrix_atom(f, A::AbstractMatrix{<:Num})
        a = Symbolics.unwrap.(A)
        return Symbolics.wrap(
            SymbolicUtils.term(f, a; type = Matrix{Real}, shape = map(Base.OneTo, size(A)))
        )
    end
end

@static if pkgversion(Symbolics) >= v"7"
    # Symbolics v6 returns a wrapped `Arr` from a symbolic matrix–matrix product;
    # v7 returns a bare `BasicSymbolic`. Re-wrap the 2-argument `Arr` product so it
    # reaches the matrix-atom methods below (which dispatch on `Arr`), restoring the
    # v6 shape. This 2-arg method is more specific than Symbolics' variadic `*`, so
    # it does not overwrite it.
    function Base.:*(x::Symbolics.Arr{<:Any, 2}, y::Symbolics.Arr{<:Any, 2})
        return Symbolics.wrap(Symbolics.unwrap(x) * Symbolics.unwrap(y))
    end

    # SymbolicUtils v4 gives `log` a matrix-permissive `promote_shape` but leaves
    # `sqrt` scalar-only, so rewriting a matrix `sqrt` term (e.g. during the
    # analysis walk's `maketerm`) throws "Invalid shapes for sqrt". Add a matrix
    # rule for `sqrt`; dispatching on the concrete `ShapeVecT` (rather than the
    # `ShapeT` union SymbolicUtils uses) makes this strictly more specific, so it
    # extends rather than overwrites the existing method.
    function SymbolicUtils.promote_shape(::typeof(sqrt), sh::SymbolicUtils.ShapeVecT)
        (length(sh) == 0 || length(sh) == 2) && return sh
        return SymbolicUtils._throw_array(sqrt, sh)
    end
end

Base.log(A::Symbolics.Arr) = matrix_atom(log, A)
Base.log(A::Matrix{Num}) = matrix_atom(log, A)
add_dcprule(log, array_domain(RealLine(), 2), Positive, Concave, Increasing)

LinearAlgebra.inv(A::Symbolics.Arr) = matrix_atom(inv, A)
add_dcprule(inv, semidefinite_domain(), AnySign, Convex, Decreasing)

LinearAlgebra.sqrt(A::Symbolics.Arr) = matrix_atom(sqrt, A)
add_dcprule(sqrt, semidefinite_domain(), Positive, Concave, Increasing)

add_dcprule(
    kldivergence,
    (array_domain(HalfLine{Real, :open}, 1), array_domain(HalfLine{Real, :open}, 1)),
    Positive,
    Convex,
    AnyMono
)

"""
    lognormcdf(x::Real)

Returns the log of the normal cumulative distribution function of `x`.

# Arguments

    - `x::Real`: A Real.
"""
function lognormcdf(x::Real)
    return logcdf(Normal(), x)
end
Symbolics.@register_symbolic lognormcdf(x::Real)
add_dcprule(lognormcdf, RealLine(), Negative, Concave, Increasing)

add_dcprule(log1p, Interval{:open, :open}(-1, Inf), Negative, Concave, Increasing)

add_dcprule(logistic, RealLine(), Positive, Convex, Increasing)

add_dcprule(max, (RealLine(), RealLine()), AnySign, Convex, Increasing)
add_dcprule(min, (RealLine(), RealLine()), AnySign, Concave, Increasing)

# special cases which depend on arguments:
function dcprule(::typeof(^), x::Symbolic, i)
    # On Symbolics v7 a literal exponent is wrapped as a constant `BasicSymbolic`,
    # so `isinteger`/`isone`/comparisons below would operate on a symbolic and
    # error; `Symbolics.value` unwraps it to the underlying number (identity on v6
    # and for already-numeric exponents).
    i = Symbolics.value(i)
    args = (x, i)
    if isone(i)
        return makerule(RealLine(), AnySign, Affine, Increasing), args
    elseif isinteger(i) && iseven(i)
        return makerule(RealLine(), Positive, Convex, increasing_if_positive), args
    elseif isinteger(i) && isodd(i)
        return makerule(HalfLine(), Positive, Convex, Increasing), args
    elseif i >= 1
        return makerule(HalfLine(), Positive, Convex, Increasing), args
    elseif i > 0 && i < 1
        return makerule(HalfLine(), Positive, Concave, Increasing), args
    elseif i < 0
        return makerule(HalfLine{Float64, :closed}(), Positive, Convex, Increasing), args
    end
end
dcprule(::typeof(Base.literal_pow), f, x...) = dcprule(^, x...)

hasdcprule(::typeof(^)) = true

add_dcprule(real, ℂ, AnySign, Affine, Increasing)

function rel_entr(x::Real, y::Real)
    if x < 0 || y < 0
        throw(DomainError((x, y), "x and y must be positive"))
    end
    if x == 0
        return 0
    end
    return x * log(x / y)
end
Symbolics.@register_symbolic rel_entr(x::Real, y::Real)
add_dcprule(
    rel_entr,
    (HalfLine{Real, :open}(), HalfLine{Real, :open}()),
    AnySign,
    Convex,
    (AnyMono, Decreasing)
)

add_dcprule(sqrt, HalfLine(), Positive, Concave, Increasing)

add_dcprule(xexpx, HalfLine, Positive, Convex, Increasing)

add_dcprule(
    conv,
    (array_domain(RealLine(), 1), array_domain(RealLine(), 1)),
    AnySign,
    Affine,
    AnyMono
)

add_dcprule(cumsum, array_domain(RealLine()), AnySign, Affine, Increasing)

add_dcprule(diagm, array_domain(RealLine(), 1), AnySign, Affine, Increasing)

add_dcprule(diag, array_domain(RealLine(), 2), AnySign, Affine, Increasing)

add_dcprule(diff, array_domain(RealLine()), AnySign, Affine, Increasing)

add_dcprule(hcat, array_domain(array_domain(RealLine(), 1), 1), AnySign, Affine, Increasing)

add_dcprule(
    kron,
    (array_domain(RealLine(), 2), array_domain(RealLine(), 2)),
    AnySign,
    Affine,
    Increasing
)

add_dcprule(reshape, array_domain(RealLine(), 2), AnySign, Affine, Increasing)

add_dcprule(triu, array_domain(RealLine(), 2), AnySign, Affine, Increasing)

add_dcprule(vec, array_domain(RealLine(), 2), AnySign, Affine, Increasing)

add_dcprule(vcat, array_domain(array_domain(RealLine(), 1), 1), AnySign, Affine, Increasing)

function dcprule(::typeof(broadcast), f, x...)
    # On Symbolics v7 the broadcasted function is wrapped as a constant symbolic
    # (e.g. `broadcast(exp, z)` carries a symbolic `exp`); `Symbolics.value`
    # recovers the underlying function (identity on v6 / for a plain function).
    return dcprule(Symbolics.value(f), x...)
end
hasdcprule(::typeof(broadcast)) = true

# add_dcprule(broadcast, (function_domain, array_domain(RealLine())), AnySign, Affine, (AnyMono, AnyMono))

add_dcprule(LinearAlgebra.adjoint, array_domain(RealLine(), 1), AnySign, Affine, Increasing)
add_dcprule(Base.getindex, array_domain(RealLine(), 1), AnySign, Affine, AnyMono)
