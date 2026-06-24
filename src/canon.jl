function canonize(ex)
    # Symbolics v7 flattens `x' * (Y * x)` / `(B' * X) * B` to a single
    # scalar/matrix `*` term, so match the flattened `*` directly.
    rs = [
        @rule adjoint(~x) * ~Y * ~x => quad_form(~x, ~Y)
        @rule adjoint(~B) * ~X * ~B => conjugation(~X, ~B)
    ]
    try
        rc = SymbolicUtils.Chain(rs)
        ex = SymbolicUtils.Postwalk(rc)(ex)
        ex = SymbolicUtils.Prewalk(rc)(ex)
        return ex
    catch
        return ex
    end
end
