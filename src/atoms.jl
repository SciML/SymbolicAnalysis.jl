### DCP atom rules

# Linear atoms — no cone needed (MOI.Reals / linear constraints)
add_dcprule(+, RealLine(), AnySign, Affine, Increasing; cone = MOI.Reals)
add_dcprule(-, RealLine(), AnySign, Affine, Decreasing; cone = MOI.Reals)

add_dcprule(Base.Ref, RealLine(), AnySign, Affine, AnyMono; cone = MOI.Reals)

add_dcprule(
    dot,
    (array_domain(RealLine()), array_domain(RealLine())),
    AnySign,
    Affine,
    Increasing;
    cone = MOI.Reals
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
    (AnyMono, increasing_if_positive ∘ minimum);
    cone = MOI.Reals  # LP reformulation
)

add_dcprule(
    StatsBase.geomean,
    array_domain(HalfLine{Real, :open}(), 1),
    Positive,
    Concave,
    Increasing;
    cone = MOI.GeometricMeanCone
)
add_dcprule(
    StatsBase.harmmean,
    array_domain(HalfLine{Real, :open}(), 1),
    Positive,
    Concave,
    Increasing;
    cone = MOI.RotatedSecondOrderCone
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

add_dcprule(invprod, array_domain(HalfLine{Real, :open}()), Positive, Convex, Decreasing;
    cone = MOI.RotatedSecondOrderCone)

add_dcprule(eigmax, symmetric_domain(), AnySign, Convex, AnyMono;
    cone = MOI.PositiveSemidefiniteConeTriangle)

add_dcprule(eigmin, symmetric_domain(), AnySign, Concave, AnyMono;
    cone = MOI.PositiveSemidefiniteConeTriangle)

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
add_dcprule(eigsummax, (array_domain(RealLine(), 2), RealLine()), AnySign, Convex, AnyMono;
    cone = MOI.PositiveSemidefiniteConeTriangle)

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
add_dcprule(eigsummin, (array_domain(RealLine(), 2), RealLine()), AnySign, Concave, AnyMono;
    cone = MOI.PositiveSemidefiniteConeTriangle)

add_dcprule(logdet, semidefinite_domain(), AnySign, Concave, AnyMono;
    cone = MOI.LogDetConeTriangle)

add_dcprule(
    LogExpFunctions.logsumexp,
    array_domain(RealLine(), 2),
    AnySign,
    Convex,
    Increasing;
    cone = MOI.ExponentialCone
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
Symbolics.@register_symbolic matrix_frac(x::AbstractVector, P::AbstractMatrix)
add_dcprule(
    matrix_frac,
    (array_domain(RealLine(), 1), definite_domain()),
    AnySign,
    Convex,
    AnyMono;
    cone = MOI.PositiveSemidefiniteConeTriangle
)

add_dcprule(maximum, array_domain(RealLine()), AnySign, Convex, Increasing;
    cone = MOI.Reals)  # LP reformulation

add_dcprule(minimum, array_domain(RealLine()), AnySign, Concave, Increasing;
    cone = MOI.Reals)  # LP reformulation

# Note: p-norms for p < 1 are not convex (they are not even norms).
# Only p >= 1 is registered as convex.
add_dcprule(
    norm,
    (array_domain(RealLine()), Interval{:closed, :open}(1, Inf)),
    Positive,
    Convex,
    increasing_if_positive;
    cone = nothing
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
    (increasing_if_positive, Increasing);
    cone = MOI.PositiveSemidefiniteConeTriangle
)

function quad_over_lin(x::AbstractVector{<:Real}, y::Real)
    if getsign(y) == Negative
        throw(DomainError(y, "y must be positive"))
    end
    return sum(x .^ 2) / y
end

Symbolics.@register_symbolic quad_over_lin(x::AbstractVector, y::Real) false

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
    (increasing_if_positive, Decreasing);
    cone = MOI.RotatedSecondOrderCone
)

add_dcprule(
    quad_over_lin,
    (RealLine(), HalfLine{Real, :open}()),
    Positive,
    Convex,
    (increasing_if_positive, Decreasing);
    cone = MOI.RotatedSecondOrderCone
)

add_dcprule(sum, array_domain(RealLine(), 2), AnySign, Affine, Increasing;
    cone = MOI.Reals)

"""
    sum_largest(x::AbstractMatrix, k)

Returns the sum of the `k` largest elements of `x`.

# Arguments

    - `x::AbstractMatrix`: A matrix.
    - `k::Int`: The number of largest elements to sum.
"""
function sum_largest(x::AbstractMatrix, k::Integer)
    return sum(sort(vec(x))[(end - k + 1):end])
end
Symbolics.@register_symbolic sum_largest(x::AbstractMatrix, k::Integer)
add_dcprule(sum_largest, (array_domain(RealLine(), 2), ℤ), AnySign, Convex, Increasing;
    cone = MOI.Reals)  # LP reformulation

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
add_dcprule(sum_smallest, (array_domain(RealLine(), 2), ℤ), AnySign, Concave, Increasing;
    cone = MOI.Reals)  # LP reformulation

add_dcprule(tr, array_domain(RealLine(), 2), AnySign, Affine, Increasing;
    cone = MOI.Reals)

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
add_dcprule(trinv, definite_domain(), Positive, Convex, AnyMono;
    cone = MOI.PositiveSemidefiniteConeTriangle)

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
add_dcprule(tv, array_domain(RealLine(), 1), Positive, Convex, AnyMono;
    cone = MOI.NormOneCone)

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
add_dcprule(tv, array_domain(array_domain(RealLine(), 2), 1), Positive, Convex, AnyMono;
    cone = MOI.SecondOrderCone)

add_dcprule(abs, ℂ, Positive, Convex, increasing_if_positive;
    cone = MOI.NormOneCone)

add_dcprule(conj, ℂ, AnySign, Affine, AnyMono;
    cone = MOI.Reals)

add_dcprule(exp, RealLine(), Positive, Convex, Increasing;
    cone = MOI.ExponentialCone)

Symbolics.@register_symbolic LogExpFunctions.xlogx(x::Real)
add_dcprule(xlogx, RealLine(), AnySign, Convex, AnyMono;
    cone = MOI.RelativeEntropyCone)

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
add_dcprule(huber, (RealLine(), HalfLine()), Positive, Convex, increasing_if_positive;
    cone = MOI.SecondOrderCone)

add_dcprule(imag, ℂ, AnySign, Affine, AnyMono;
    cone = MOI.Reals)

add_dcprule(inv, HalfLine{Real, :open}(), Positive, Convex, Decreasing;
    cone = MOI.RotatedSecondOrderCone)
add_dcprule(log, HalfLine{Real, :open}(), AnySign, Concave, Increasing;
    cone = MOI.ExponentialCone)

@register_symbolic Base.log(A::Symbolics.Arr)
add_dcprule(log, array_domain(RealLine(), 2), Positive, Concave, Increasing;
    cone = MOI.ExponentialCone)

@register_symbolic LinearAlgebra.inv(A::Symbolics.Arr)
add_dcprule(inv, semidefinite_domain(), AnySign, Convex, Decreasing;
    cone = MOI.PositiveSemidefiniteConeTriangle)

@register_symbolic LinearAlgebra.sqrt(A::Symbolics.Arr)
add_dcprule(sqrt, semidefinite_domain(), Positive, Concave, Increasing;
    cone = MOI.PositiveSemidefiniteConeTriangle)

add_dcprule(
    kldivergence,
    (array_domain(HalfLine{Real, :open}, 1), array_domain(HalfLine{Real, :open}, 1)),
    Positive,
    Convex,
    AnyMono;
    cone = MOI.RelativeEntropyCone
)

"""
    lognormcdf(x::Real)

Returns the log of the normal cumulative distribution function of `x`.

# Arguments

    - `x::Real`: A Real.
"""
function lognormcdf(x::Real)
    return logcdf(Normal, x)
end
Symbolics.@register_symbolic lognormcdf(x::Real)
add_dcprule(lognormcdf, RealLine(), Negative, Concave, Increasing)

add_dcprule(log1p, Interval{:open, :open}(-1, Inf), Negative, Concave, Increasing;
    cone = MOI.ExponentialCone)

add_dcprule(logistic, RealLine(), Positive, Convex, Increasing;
    cone = MOI.ExponentialCone)

add_dcprule(max, (RealLine(), RealLine()), AnySign, Convex, Increasing;
    cone = MOI.Reals)  # LP reformulation
add_dcprule(min, (RealLine(), RealLine()), AnySign, Concave, Increasing;
    cone = MOI.Reals)  # LP reformulation

# special cases which depend on arguments:
function dcprule(::typeof(^), x::Symbolic, i)
    args = (x, i)
    if isone(i)
        return makerule(RealLine(), AnySign, Affine, Increasing; cone = MOI.Reals), args
    elseif i == 2
        return makerule(RealLine(), Positive, Convex, increasing_if_positive;
            cone = MOI.RotatedSecondOrderCone), args
    elseif isinteger(i) && iseven(i)
        return makerule(RealLine(), Positive, Convex, increasing_if_positive;
            cone = nothing), args
    elseif isinteger(i) && isodd(i)
        return makerule(HalfLine(), Positive, Convex, Increasing;
            cone = MOI.PowerCone), args
    elseif i >= 1
        return makerule(HalfLine(), Positive, Convex, Increasing;
            cone = MOI.PowerCone), args
    elseif i > 0 && i < 1
        return makerule(HalfLine(), Positive, Concave, Increasing;
            cone = MOI.PowerCone), args
    elseif i < 0
        return makerule(HalfLine{Float64, :closed}(), Positive, Convex, Increasing;
            cone = MOI.PowerCone), args
    end
end
dcprule(::typeof(Base.literal_pow), f, x...) = dcprule(^, x...)

hasdcprule(::typeof(^)) = true

add_dcprule(real, ℂ, AnySign, Affine, Increasing;
    cone = MOI.Reals)

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
    (AnyMono, Decreasing);
    cone = MOI.RelativeEntropyCone
)

add_dcprule(sqrt, HalfLine(), Positive, Concave, Increasing;
    cone = MOI.RotatedSecondOrderCone)

add_dcprule(xexpx, HalfLine, Positive, Convex, Increasing;
    cone = MOI.ExponentialCone)

add_dcprule(
    conv,
    (array_domain(RealLine(), 1), array_domain(RealLine(), 1)),
    AnySign,
    Affine,
    AnyMono;
    cone = MOI.Reals
)

add_dcprule(cumsum, array_domain(RealLine()), AnySign, Affine, Increasing;
    cone = MOI.Reals)

add_dcprule(diagm, array_domain(RealLine(), 1), AnySign, Affine, Increasing;
    cone = MOI.Reals)

add_dcprule(diag, array_domain(RealLine(), 2), AnySign, Affine, Increasing;
    cone = MOI.Reals)

add_dcprule(diff, array_domain(RealLine()), AnySign, Affine, Increasing;
    cone = MOI.Reals)

add_dcprule(hcat, array_domain(array_domain(RealLine(), 1), 1), AnySign, Affine, Increasing;
    cone = MOI.Reals)

add_dcprule(
    kron,
    (array_domain(RealLine(), 2), array_domain(RealLine(), 2)),
    AnySign,
    Affine,
    Increasing;
    cone = MOI.Reals
)

add_dcprule(reshape, array_domain(RealLine(), 2), AnySign, Affine, Increasing;
    cone = MOI.Reals)

add_dcprule(triu, array_domain(RealLine(), 2), AnySign, Affine, Increasing;
    cone = MOI.Reals)

add_dcprule(vec, array_domain(RealLine(), 2), AnySign, Affine, Increasing;
    cone = MOI.Reals)

add_dcprule(vcat, array_domain(array_domain(RealLine(), 1), 1), AnySign, Affine, Increasing;
    cone = MOI.Reals)

function dcprule(::typeof(broadcast), f, x...)
    return dcprule(f, x...)
end
hasdcprule(::typeof(broadcast)) = true

# add_dcprule(broadcast, (function_domain, array_domain(RealLine())), AnySign, Affine, (AnyMono, AnyMono))

add_dcprule(LinearAlgebra.adjoint, array_domain(RealLine(), 1), AnySign, Affine, Increasing;
    cone = MOI.Reals)
add_dcprule(Base.getindex, array_domain(RealLine(), 1), AnySign, Affine, AnyMono;
    cone = MOI.Reals)
