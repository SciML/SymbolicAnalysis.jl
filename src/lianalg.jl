using LinearAlgebra
using Symbolics: Num, simplify

function LinearAlgebra.ishermitian(A::AbstractMatrix{Num}; kwargs...)
    indsm, indsn = axes(A)
    if indsm != indsn
        return false
    end
    for i in indsn, j in i:last(indsn)
        d = simplify(A[i, j] - adjoint(A[j, i]))
        if !isapprox(d, 0.0; kwargs...)
            return false
        end
    end
    return true
end

## Numbers
function LinearAlgebra._chol!(x::Num, uplo)
    rx = real(x)
    rxr = sqrt(abs(rx))
    rval = convert(promote_type(typeof(x), typeof(rxr)), rxr)
    d = rx - abs(x) |> simplify
    println(d)
    isapprox(d, 0.0) ? (rval, convert(BlasInt, 0)) : (rval, convert(BlasInt, 1))
end
