function canonize(ex)
    rs = [@rule (adjoint(~x) * (~Y * ~x))[1] => quad_form(~x, ~Y)
          @rule ((adjoint(~B) * ~X) * ~B)[Base.OneTo(size(~B, 2)), Base.OneTo(size(
              ~B, 1))] => conjugation(~X, ~B)]
    try
        rc = SymbolicUtils.Chain(rs)
        ex = SymbolicUtils.Postwalk(rc)(ex)
        ex = SymbolicUtils.Prewalk(rc)(ex)
        return ex
    catch
        return ex
    end
end
