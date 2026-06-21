function canonize(ex)
    # Symbolics v6 keeps `x' * (Y * x)` / `(B' * X) * B` as (1×1 or full) arrays,
    # so the canonicalization rules index into them. Symbolics v7 flattens these
    # to a single scalar/matrix `*` term and rejects `getindex` on a slot pattern
    # of unknown shape (it errors at rule-construction time), so v7 matches the
    # flattened `*` directly without indexing.
    rs = if pkgversion(Symbolics) < v"7"
        [
            @rule (adjoint(~x) * (~Y * ~x))[1] => quad_form(~x, ~Y)
            @rule ((adjoint(~B) * ~X) * ~B)[
                Base.OneTo(size(~B, 2)), Base.OneTo(
                    size(
                        ~B, 1
                    )
                ),
            ] => conjugation(~X, ~B)
        ]
    else
        [
            @rule adjoint(~x) * ~Y * ~x => quad_form(~x, ~Y)
            @rule adjoint(~B) * ~X * ~B => conjugation(~X, ~B)
        ]
    end
    try
        rc = SymbolicUtils.Chain(rs)
        ex = SymbolicUtils.Postwalk(rc)(ex)
        ex = SymbolicUtils.Prewalk(rc)(ex)
        return ex
    catch
        return ex
    end
end
