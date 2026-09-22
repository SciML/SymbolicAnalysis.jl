---
---

# Special functions {#Special-functions}

Since some atoms are not available in the base language or other packages we have implemented them here.
<details class='jldocstring custom-block' open>
<summary><a id='SymbolicAnalysis.dcprule-Tuple{typeof-, Any}' href='#SymbolicAnalysis.dcprule-Tuple{typeof-, Any}'><span class="jlbinding">SymbolicAnalysis.dcprule</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
dcprule(::typeof(-), x[, y])
```


Unary `-x` is decreasing, but binary `x - y` is _increasing_ in `x` and decreasing only in `y`. A single declared monotonicity is broadcast to every argument by `get_arg_property`, so registering `Decreasing` flipped the first slot too: `exp.(v) .- w`, whose Hessian is `diag(exp(vᵢ)) ⪰ 0`, certified as `Concave`. Scalar `a - b` never reaches here — Symbolics rewrites it to `a + (-1)*b` — but the broadcast form does.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/SciML/SymbolicAnalysis.jl" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='SymbolicAnalysis.dcprule-Tuple{typeof/, Any, Any}' href='#SymbolicAnalysis.dcprule-Tuple{typeof/, Any, Any}'><span class="jlbinding">SymbolicAnalysis.dcprule</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
dcprule(::typeof(/), num, den)
```


`/` is bilinear in the same sense as `*`, so it is DCP only when one side is constant. A nonzero constant denominator is an affine rescaling, flipping curvature and sign when it is negative.

A constant numerator is the `inv` atom: `c/x` is convex-decreasing for `c > 0` and concave-increasing for `c < 0` — but **only on `x > 0`**, and the denominator must be _known_ positive for the rule to fire. That is deliberately stricter than the unchecked-precondition convention `log`, `sqrt` and `geomean` use, and the asymmetry is the point: those are only wrong where the function does not exist, whereas `1/x` exists for every nonzero `x` and is _concave_ below zero. Certifying it `Convex` for a sign-unknown argument would be a false certificate over the function's own real domain — `f(-1) = -1`, `f(2) = 0.5`, and the midpoint `f(0.5) = 2` is far above the chord `-0.25`.

Symbolics rewrites `inv(x)` to `/(1, x)` and `x^-n` to `/(1, x)^n`, so this rule is also what gives those their curvature.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/SciML/SymbolicAnalysis.jl" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='SymbolicAnalysis.dcprule-Tuple{typeofLinearAlgebra.dot, Any, Any}' href='#SymbolicAnalysis.dcprule-Tuple{typeofLinearAlgebra.dot, Any, Any}'><span class="jlbinding">SymbolicAnalysis.dcprule</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
dcprule(::typeof(dot), x, y)
```


`dot` is bilinear, so it is affine only when one side is constant: `dot(c, x)` is affine in `x`, but `dot(x, x)` is the quadratic `‖x‖²`. Registering it as unconditionally affine certified `dot(x, x)` as `Affine`, which is not a valid certificate, so the curvature depends on the arguments like `^` and `norm`.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/SciML/SymbolicAnalysis.jl" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='SymbolicAnalysis.dcprule-Tuple{typeofLogExpFunctions.xexpx, Any}' href='#SymbolicAnalysis.dcprule-Tuple{typeofLogExpFunctions.xexpx, Any}'><span class="jlbinding">SymbolicAnalysis.dcprule</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
dcprule(::typeof(xexpx), x)
```


`x*exp(x)` is convex only on `[-2, Inf)`: at `x = -4` the second difference is `-0.0366`, at `x = -3` it is `-0.0498`, and the inflection sits exactly at `-2`. So unlike `log` or `sqrt`, this atom **exists** outside its declared domain and is concave there — the same shape that made `/` establish its precondition rather than assume it, rather than the `log`/`sqrt` case where the declared domain is only where the function is defined at all.

Registering it unguarded would have been a false certificate the moment the atom became reachable: until this rule existed `xexpx(s)` traced to `*` and degraded, which is the only reason the unenforced domain never bit.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/SciML/SymbolicAnalysis.jl" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='SymbolicAnalysis.dcprule-Tuple{typeofSymbolicAnalysis.dotsort, Any, Any}' href='#SymbolicAnalysis.dcprule-Tuple{typeofSymbolicAnalysis.dotsort, Any, Any}'><span class="jlbinding">SymbolicAnalysis.dcprule</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
dcprule(::typeof(dotsort), x, y)
```


`dotsort` is a pointwise maximum of bilinear forms, so — like `dot` — it is a convex atom only when one side is constant. With both sides symbolic it is indefinite: for length-1 vectors it is literally `x[1]*y[1]`.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/SciML/SymbolicAnalysis.jl" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='SymbolicAnalysis.dcprule-Tuple{typeofSymbolicAnalysis.huber, Any, Any}' href='#SymbolicAnalysis.dcprule-Tuple{typeofSymbolicAnalysis.huber, Any, Any}'><span class="jlbinding">SymbolicAnalysis.dcprule</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
dcprule(::typeof(huber), x, M)
```


Huber is convex in `x`, but for `abs(x) > M` it equals `2M*abs(x) - M^2`, which is **concave** in `M` (`d²/dM² = -2`). The threshold must therefore be a constant for the `Convex` certificate to hold.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/SciML/SymbolicAnalysis.jl" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='SymbolicAnalysis.dcprule-Tuple{typeofSymbolicAnalysis.quad_form, Any, Any}' href='#SymbolicAnalysis.dcprule-Tuple{typeofSymbolicAnalysis.quad_form, Any, Any}'><span class="jlbinding">SymbolicAnalysis.dcprule</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
dcprule(::typeof(quad_form), x, P)
```


`x'Px` is quadratic in `x` but **linear** in `P`, so it is convex only for a constant positive definite `P`. Registering it as unconditionally `Convex` certified both the both-symbolic form (indefinite: `quad_form([a], [b;;])` is `a^2*b`, whose second difference along `(1, -1)` is `-2`) and a constant indefinite `P` (`quad_form(x, [1 0; 0 -1])` is `x[1]^2 - x[2]^2`).

`P` is checked with `isposdef`, matching the `semidefinite_domain()` the rule declares; a singular positive semidefinite `P` therefore gets no certificate.

Slot 1 is only `increasing_if_positive` when `P` is **entrywise** nonnegative. `x'Px = Σ P[i,j]·x[i]·x[j]` is nondecreasing in each `x[i]` over the nonnegative orthant only then; positive definiteness alone does not give it. With `P = [1 -0.9; -0.9 1]` (which is positive definite) the composition certified `quad_form(exp.(v), P)` as `Convex` while it is not — the midpoint exceeds the chord by 1.92 between `v = [0.64, 2.0]` and `[1.64, 2.3]`.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/SciML/SymbolicAnalysis.jl" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='SymbolicAnalysis.dcprule-Tuple{typeofbroadcast, Any, Vararg{Any}}' href='#SymbolicAnalysis.dcprule-Tuple{typeofbroadcast, Any, Vararg{Any}}'><span class="jlbinding">SymbolicAnalysis.dcprule</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
dcprule(::typeof(broadcast), f, x...)
```


A broadcast carries its function as the first argument, so the rule is the one `f` has when applied elementwise. Two cases are not a rule lookup: `*` has no table entry (it is special-cased by the multiplication helpers, which pick out the constant factor and flip the curvature on a negative one), and a broadcasted function with no rule at all must degrade to `UnknownCurvature` rather than fail the lookup.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/SciML/SymbolicAnalysis.jl" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='SymbolicAnalysis.dcprule-Tuple{typeofexp, Any}' href='#SymbolicAnalysis.dcprule-Tuple{typeofexp, Any}'><span class="jlbinding">SymbolicAnalysis.dcprule</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
dcprule(::typeof(exp), x)
```


`exp(X)` on a matrix is the matrix exponential, not the scalar law applied pointwise: `expm` is neither operator convex nor operator monotone, so `sum(exp(X))` and `tr(exp(X))` are indefinite even on the SPD cone. The elementwise `exp.(X)` arrives through `broadcast` and keeps the scalar law via `elementwise_dcprule`.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/SciML/SymbolicAnalysis.jl" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='SymbolicAnalysis.dcprule-Tuple{typeofminimum, Any}' href='#SymbolicAnalysis.dcprule-Tuple{typeofminimum, Any}'><span class="jlbinding">SymbolicAnalysis.dcprule</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
dcprule(::typeof(minimum), x)
```


An entrywise `minimum` must not consume a Loewner-order matrix rule. The smallest entry of a positive definite matrix can be an off-diagonal one, so `minimum(sqrt(X))` and `minimum(log(X))` are indefinite over the SPD cone (second differences of −0.274 and −0.753 at one base point, +0.462 and +1.305 at another) while both certified as `Concave`.

`maximum` needs no such guard: for a positive definite `M`, `M[i,j] ≤ √(M[i,i]·M[j,j]) ≤ max(M[i,i], M[j,j])`, so the largest entry is always on the diagonal and `maximum(inv(X))` is a max of convex functions.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/SciML/SymbolicAnalysis.jl" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='SymbolicAnalysis.dotsort-Tuple{AbstractVector, AbstractVector}' href='#SymbolicAnalysis.dotsort-Tuple{AbstractVector, AbstractVector}'><span class="jlbinding">SymbolicAnalysis.dotsort</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
dotsort(x, y)
```


Sorts `x` and `y` and returns the dot product of the sorted vectors.

**Arguments**

```julia
- `x::AbstractVector`: A vector.
- `y::AbstractVector`: A vector.
```



<Badge type="info" class="source-link" text="source"><a href="https://github.com/SciML/SymbolicAnalysis.jl" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='SymbolicAnalysis.eigsummax-Tuple{LinearAlgebra.Symmetric, Int64}' href='#SymbolicAnalysis.eigsummax-Tuple{LinearAlgebra.Symmetric, Int64}'><span class="jlbinding">SymbolicAnalysis.eigsummax</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
eigsummax(m::Symmetric, k)
```


Returns the sum of the `k` largest eigenvalues of `m`.

**Arguments**

```julia
- `m::Symmetric`: A symmetric matrix.
- `k::Int`: The number of largest eigenvalues to sum.
```



<Badge type="info" class="source-link" text="source"><a href="https://github.com/SciML/SymbolicAnalysis.jl" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='SymbolicAnalysis.eigsummin-Tuple{LinearAlgebra.Symmetric, Int64}' href='#SymbolicAnalysis.eigsummin-Tuple{LinearAlgebra.Symmetric, Int64}'><span class="jlbinding">SymbolicAnalysis.eigsummin</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
eigsummin(m::Symmetric, k)
```


Returns the sum of the `k` smallest eigenvalues of `m`.

**Arguments**

```julia
- `m::Symmetric`: A symmetric matrix.
- `k::Int`: The number of smallest eigenvalues to sum.
```



<Badge type="info" class="source-link" text="source"><a href="https://github.com/SciML/SymbolicAnalysis.jl" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='SymbolicAnalysis.elementwise_dcprule-Tuple{Any, Vararg{Any}}' href='#SymbolicAnalysis.elementwise_dcprule-Tuple{Any, Vararg{Any}}'><span class="jlbinding">SymbolicAnalysis.elementwise_dcprule</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
elementwise_dcprule(f, args...)
```


The rule for `f` applied elementwise, as under a broadcast. Defaults to `f`'s own rule. Atoms that mean something different on a matrix than pointwise — `^` is a matrix power, `exp` is the matrix exponential — guard against the matrix meaning in `dcprule` and restore the scalar law here.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/SciML/SymbolicAnalysis.jl" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='SymbolicAnalysis.huber' href='#SymbolicAnalysis.huber'><span class="jlbinding">SymbolicAnalysis.huber</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
huber(x, M=1)
```


Returns the Huber loss function of `x` with threshold `M`.

**Arguments**

```julia
- `x::Real`: A Real.
- `M::Real`: The threshold.
```



<Badge type="info" class="source-link" text="source"><a href="https://github.com/SciML/SymbolicAnalysis.jl" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='SymbolicAnalysis.invprod-Tuple{AbstractVector}' href='#SymbolicAnalysis.invprod-Tuple{AbstractVector}'><span class="jlbinding">SymbolicAnalysis.invprod</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
invprod(x::AbstractVector)
```


Returns the inverse of the product of the elements of `x`.

**Arguments**

```julia
- `x::AbstractVector`: A vector.
```



<Badge type="info" class="source-link" text="source"><a href="https://github.com/SciML/SymbolicAnalysis.jl" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='SymbolicAnalysis.lognormcdf-Tuple{Real}' href='#SymbolicAnalysis.lognormcdf-Tuple{Real}'><span class="jlbinding">SymbolicAnalysis.lognormcdf</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
lognormcdf(x::Real)
```


Returns the log of the normal cumulative distribution function of `x`.

**Arguments**

```julia
- `x::Real`: A Real.
```



<Badge type="info" class="source-link" text="source"><a href="https://github.com/SciML/SymbolicAnalysis.jl" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='SymbolicAnalysis.matrix_frac-Tuple{AbstractVector, AbstractMatrix}' href='#SymbolicAnalysis.matrix_frac-Tuple{AbstractVector, AbstractMatrix}'><span class="jlbinding">SymbolicAnalysis.matrix_frac</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
matrix_frac(x::AbstractVector, P::AbstractMatrix)
```


Returns the quadratic form `x' * P^{-1} * x`.

**Arguments**

```julia
- `x::AbstractVector`: A vector.
- `P::AbstractMatrix`: A matrix.
```



<Badge type="info" class="source-link" text="source"><a href="https://github.com/SciML/SymbolicAnalysis.jl" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='SymbolicAnalysis.multilinear_rule-Tuple{Any, Any}' href='#SymbolicAnalysis.multilinear_rule-Tuple{Any, Any}'><span class="jlbinding">SymbolicAnalysis.multilinear_rule</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
multilinear_rule(domain, args)
```


Shared by `dot`, `conv` and `kron`: each is linear in every argument separately but not jointly, so each is affine only when all but one argument is constant. `conv(x, x)` and `kron(X, X)` are quadratic in the same way `dot(x, x)` is.

The affine case declares `AnyMono` rather than `Increasing` because the constant side holds the coefficients of the linear map and they may have either sign. That only refuses a _curved_ argument such as `kron(C, exp.(X))`; an affine one composes without consulting monotonicity at all.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/SciML/SymbolicAnalysis.jl" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='SymbolicAnalysis.perspective-Tuple{Function, Any, Real}' href='#SymbolicAnalysis.perspective-Tuple{Function, Any, Real}'><span class="jlbinding">SymbolicAnalysis.perspective</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
perspective(f::Function, x, s::Real)
```


Returns the perspective function `s * f(x / s)`.

**Arguments**

```julia
- `f::Function`: A function.
- `x`: A Real.
- `s::Real`: A positive Real.
```



<Badge type="info" class="source-link" text="source"><a href="https://github.com/SciML/SymbolicAnalysis.jl" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='SymbolicAnalysis.quad_form-Tuple{AbstractVector, AbstractMatrix}' href='#SymbolicAnalysis.quad_form-Tuple{AbstractVector, AbstractMatrix}'><span class="jlbinding">SymbolicAnalysis.quad_form</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
quad_form(x::AbstractVector, P::AbstractMatrix)
```


Returns the quadratic form `x' * P * x`.

**Arguments**

```julia
- `x::AbstractVector`: A vector.
- `P::AbstractMatrix`: A matrix.
```



<Badge type="info" class="source-link" text="source"><a href="https://github.com/SciML/SymbolicAnalysis.jl" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='SymbolicAnalysis.quad_over_lin-Tuple{Real, Real}' href='#SymbolicAnalysis.quad_over_lin-Tuple{Real, Real}'><span class="jlbinding">SymbolicAnalysis.quad_over_lin</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
quad_over_lin(x::Real, y::Real)
```


Returns the quadratic over linear form `x^2 / y`.

**Arguments**

```julia
- `x`: A Real or a vector.
- `y::Real`: A positive Real.
```



<Badge type="info" class="source-link" text="source"><a href="https://github.com/SciML/SymbolicAnalysis.jl" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='SymbolicAnalysis.sum_largest-Tuple{AbstractMatrix, Integer}' href='#SymbolicAnalysis.sum_largest-Tuple{AbstractMatrix, Integer}'><span class="jlbinding">SymbolicAnalysis.sum_largest</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
sum_largest(x::AbstractMatrix, k)
```


Returns the sum of the `k` largest elements of `x`.

**Arguments**

```julia
- `x::AbstractMatrix`: A matrix.
- `k::Int`: The number of largest elements to sum.
```



<Badge type="info" class="source-link" text="source"><a href="https://github.com/SciML/SymbolicAnalysis.jl" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='SymbolicAnalysis.sum_smallest-Tuple{AbstractMatrix, Integer}' href='#SymbolicAnalysis.sum_smallest-Tuple{AbstractMatrix, Integer}'><span class="jlbinding">SymbolicAnalysis.sum_smallest</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
sum_smallest(x::AbstractMatrix, k)
```


Returns the sum of the `k` smallest elements of `x`.

**Arguments**

```julia
- `x::AbstractMatrix`: A matrix.
- `k::Int`: The number of smallest elements to sum.
```



<Badge type="info" class="source-link" text="source"><a href="https://github.com/SciML/SymbolicAnalysis.jl" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='SymbolicAnalysis.trinv-Tuple{AbstractMatrix}' href='#SymbolicAnalysis.trinv-Tuple{AbstractMatrix}'><span class="jlbinding">SymbolicAnalysis.trinv</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
trinv(x::AbstractMatrix)
```


Returns the trace of the inverse of `x`.

**Arguments**

```julia
- `x::AbstractMatrix`: A matrix.
```



<Badge type="info" class="source-link" text="source"><a href="https://github.com/SciML/SymbolicAnalysis.jl" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='SymbolicAnalysis.tv-Tuple{AbstractVector{<:AbstractMatrix}}' href='#SymbolicAnalysis.tv-Tuple{AbstractVector{<:AbstractMatrix}}'><span class="jlbinding">SymbolicAnalysis.tv</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
tv(x::AbstractVector{<:AbstractMatrix})
```


Returns the total variation of `x`, defined as `sum_{i,j} |x_{k+1}[i,j] - x_k[i,j]|`.

**Arguments**

```julia
- `x::AbstractVector`: A vector of matrices.
```



<Badge type="info" class="source-link" text="source"><a href="https://github.com/SciML/SymbolicAnalysis.jl" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='SymbolicAnalysis.tv-Tuple{AbstractVector{<:Real}}' href='#SymbolicAnalysis.tv-Tuple{AbstractVector{<:Real}}'><span class="jlbinding">SymbolicAnalysis.tv</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
tv(x::AbstractVector{<:Real})
```


Returns the total variation of `x`, defined as `sum_i |x_{i+1} - x_i|`.

**Arguments**

```julia
- `x::AbstractVector`: A vector.
```



<Badge type="info" class="source-link" text="source"><a href="https://github.com/SciML/SymbolicAnalysis.jl" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='SymbolicAnalysis.affine_map-Tuple{typeofSymbolicAnalysis.conjugation, Matrix, Matrix, Matrix}' href='#SymbolicAnalysis.affine_map-Tuple{typeofSymbolicAnalysis.conjugation, Matrix, Matrix, Matrix}'><span class="jlbinding">SymbolicAnalysis.affine_map</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
affine_map(f, X, B, Y)
affine_map(f, X, B, Ys)
```


Affine map, i.e., `B + f(X, Y)` or `B + sum(f(X, Y) for Y in Ys)` for a function `f` where `f` is a positive linear operator.

**Arguments**

```julia
- `f::Function`: One of the following functions: `conjugation`, `diag`, `tr` and `hadamard_product`.
- `X::Matrix`: A symmetric positive definite matrix.
- `B::Matrix`: A matrix.
- `Y::Matrix`: A matrix.
- `Ys::Vector{<:Matrix}`: A vector of matrices.
```



<Badge type="info" class="source-link" text="source"><a href="https://github.com/SciML/SymbolicAnalysis.jl" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='SymbolicAnalysis.conjugation-Tuple{Any, Any}' href='#SymbolicAnalysis.conjugation-Tuple{Any, Any}'><span class="jlbinding">SymbolicAnalysis.conjugation</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
conjugation(X, B)
```


Conjugation of a matrix `X` by a matrix `B` is defined as `B'X*B`.

**Arguments**

```julia
- `X::Matrix`: A symmetric positive definite matrix.
- `B::Matrix`: A matrix.
```



<Badge type="info" class="source-link" text="source"><a href="https://github.com/SciML/SymbolicAnalysis.jl" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='SymbolicAnalysis.hadamard_product-Tuple{AbstractMatrix, AbstractMatrix}' href='#SymbolicAnalysis.hadamard_product-Tuple{AbstractMatrix, AbstractMatrix}'><span class="jlbinding">SymbolicAnalysis.hadamard_product</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
hadamard_product(X, B)
```


Hadamard product or element-wise multiplication of a symmetric positive definite matrix `X` by a positive semi-definite matrix `B`.

**Arguments**

```julia
- `X::Matrix`: A symmetric positive definite matrix.
- `B::Matrix`: A positive semi-definite matrix.
```



<Badge type="info" class="source-link" text="source"><a href="https://github.com/SciML/SymbolicAnalysis.jl" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='SymbolicAnalysis.log_quad_form-Tuple{Vector{<:Number}, Matrix}' href='#SymbolicAnalysis.log_quad_form-Tuple{Vector{<:Number}, Matrix}'><span class="jlbinding">SymbolicAnalysis.log_quad_form</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
log_quad_form(y, X)
log_quad_form(ys, X)
```


Log of the quadratic form of a symmetric positive definite matrix `X` and a vector `y` is defined as `log(y'*X*y)` or for a vector of vectors `ys` as `log(sum(y'*X*y for y in ys))`.

**Arguments**

```julia
- `y::Vector`: A vector of `Number`s or a `Vector` of `Vector`s.
- `X::Matrix`: A symmetric positive definite matrix.
```



<Badge type="info" class="source-link" text="source"><a href="https://github.com/SciML/SymbolicAnalysis.jl" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='SymbolicAnalysis.scalar_mat' href='#SymbolicAnalysis.scalar_mat'><span class="jlbinding">SymbolicAnalysis.scalar_mat</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
scalar_mat(X, k=size(X, 1))
```


Scalar matrix of a symmetric positive definite matrix `X` is defined as `tr(X)*I(k)`.

**Arguments**

```julia
- `X::Matrix`: A symmetric positive definite matrix.
- `k::Int`: The size of the identity matrix.
```



<Badge type="info" class="source-link" text="source"><a href="https://github.com/SciML/SymbolicAnalysis.jl" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='SymbolicAnalysis.schatten_norm' href='#SymbolicAnalysis.schatten_norm'><span class="jlbinding">SymbolicAnalysis.schatten_norm</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
schatten_norm(X, p=2)
```


Schatten norm of a symmetric positive definite matrix `X`.

**Arguments**

```julia
- `X::Matrix`: A symmetric positive definite matrix.
- `p::Int`: The p-norm.
```



<Badge type="info" class="source-link" text="source"><a href="https://github.com/SciML/SymbolicAnalysis.jl" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='SymbolicAnalysis.sdivergence-Tuple{Any, Any}' href='#SymbolicAnalysis.sdivergence-Tuple{Any, Any}'><span class="jlbinding">SymbolicAnalysis.sdivergence</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
sdivergence(X, Y)
```


Symmetric divergence of two symmetric positive definite matrices `X` and `Y` is defined as `logdet((X+Y)/2) - 1/2*logdet(X*Y)`.

**Arguments**

```julia
- `X::Matrix`: A symmetric positive definite matrix.
- `Y::Matrix`: A symmetric positive definite matrix.
```



<Badge type="info" class="source-link" text="source"><a href="https://github.com/SciML/SymbolicAnalysis.jl" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='SymbolicAnalysis.sum_log_eigmax-Tuple{Function, AbstractMatrix, Int64}' href='#SymbolicAnalysis.sum_log_eigmax-Tuple{Function, AbstractMatrix, Int64}'><span class="jlbinding">SymbolicAnalysis.sum_log_eigmax</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
sum_log_eigmax(X, k)
sum_log_eigmax(f, X, k)
```


Sum of the log of the maximum eigenvalues of a symmetric positive definite matrix `X`. If a function `f` is provided, the sum is over `f` applied to the log of the eigenvalues.

**Arguments**

```julia
- `f::Function`: A function.
- `X::Matrix`: A symmetric positive definite matrix.
- `k::Int`: The number of eigenvalues to consider.
```



<Badge type="info" class="source-link" text="source"><a href="https://github.com/SciML/SymbolicAnalysis.jl" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='SymbolicAnalysis.lorentz_homogeneous_diagonal-Tuple{AbstractVector, AbstractVector}' href='#SymbolicAnalysis.lorentz_homogeneous_diagonal-Tuple{AbstractVector, AbstractVector}'><span class="jlbinding">SymbolicAnalysis.lorentz_homogeneous_diagonal</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
lorentz_homogeneous_diagonal(a, p) -> Real
```


Evaluate the diagonal quadratic `sum(a[i] * p[i]^2)` on the Lorentz model. Geodesic convexity requires `minimum(a[1:end-1]) + a[end] >= 0`.

**Arguments**
- `a::AbstractVector`: `(d + 1)` diagonal coefficients.
    
- `p::AbstractVector`: point on the Lorentz manifold.
    

**Returns**
- A scalar containing the diagonal quadratic value.
    

**Throws**
- `DimensionMismatch`: if `a` and `p` have different lengths.
    
- `ArgumentError`: if `a` does not satisfy the geodesic-convexity condition.
    

**Examples**

```julia
julia> using SymbolicAnalysis

julia> lorentz_homogeneous_diagonal([1.0, 2.0, 0.0], [0.0, 0.0, 1.0])
0.0
```



<Badge type="info" class="source-link" text="source"><a href="https://github.com/SciML/SymbolicAnalysis.jl" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='SymbolicAnalysis.lorentz_homogeneous_quadratic-Tuple{AbstractMatrix, AbstractVector}' href='#SymbolicAnalysis.lorentz_homogeneous_quadratic-Tuple{AbstractMatrix, AbstractVector}'><span class="jlbinding">SymbolicAnalysis.lorentz_homogeneous_quadratic</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
lorentz_homogeneous_quadratic(A, p) -> Real
```


Evaluate the homogeneous quadratic `transpose(p) * A * p` on the Lorentz model. The matrix must satisfy one of the implemented geodesic-convexity conditions.

**Arguments**
- `A::AbstractMatrix`: symmetric `(d + 1)` by `(d + 1)` coefficient matrix.
    
- `p::AbstractVector`: point on the Lorentz manifold.
    

**Returns**
- A scalar containing the quadratic value.
    

**Throws**
- `ArgumentError`: if `A` does not satisfy the geodesic-convexity conditions.
    

**Examples**

```julia
julia> using SymbolicAnalysis, LinearAlgebra

julia> A = Matrix{Float64}(I, 3, 3);

julia> lorentz_homogeneous_quadratic(A, [0.0, 0.0, 1.0])
1.0
```



<Badge type="info" class="source-link" text="source"><a href="https://github.com/SciML/SymbolicAnalysis.jl" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='SymbolicAnalysis.lorentz_least_squares-Tuple{AbstractMatrix, AbstractVector, AbstractVector}' href='#SymbolicAnalysis.lorentz_least_squares-Tuple{AbstractMatrix, AbstractVector, AbstractVector}'><span class="jlbinding">SymbolicAnalysis.lorentz_least_squares</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
lorentz_least_squares(X, y, p) -> Real
```


Evaluate the squared residual norm `sum(abs2, y - X * p)` on the Lorentz model. The derived quadratic and linear terms must satisfy the implemented geodesic-convexity conditions.

**Arguments**
- `X::AbstractMatrix`: design matrix with `d + 1` columns.
    
- `y::AbstractVector`: response vector with one entry per row of `X`.
    
- `p::AbstractVector`: point on the Lorentz manifold.
    

**Returns**
- A scalar containing the squared residual norm.
    

**Throws**
- `ArgumentError`: if the derived homogeneous or linear term does not satisfy the geodesic-convexity conditions.
    
- `DimensionMismatch`: if the dimensions of `X`, `y`, and `p` are incompatible.
    

**Examples**

```julia
julia> using SymbolicAnalysis, LinearAlgebra

julia> X = Matrix{Float64}(I, 3, 3);

julia> lorentz_least_squares(X, [0.0, 0.0, -1.0], [0.0, 0.0, 1.0])
4.0
```



<Badge type="info" class="source-link" text="source"><a href="https://github.com/SciML/SymbolicAnalysis.jl" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='SymbolicAnalysis.lorentz_log_barrier-Tuple{AbstractVector}' href='#SymbolicAnalysis.lorentz_log_barrier-Tuple{AbstractVector}'><span class="jlbinding">SymbolicAnalysis.lorentz_log_barrier</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
lorentz_log_barrier(p) -> Real
```


Evaluate the log-barrier `-log(-1 - <a, p>_L)` for the Lorentz model, where `a = (0, ..., 0, 1)` and `<., .>_L` is the Lorentzian inner product.

**Arguments**
- `p::AbstractVector`: point on the Lorentz manifold. Its last coordinate must be greater than one for the barrier to be finite and real.
    

**Returns**
- A real scalar containing the barrier value.
    

**Throws**
- `DomainError`: if the last coordinate of `p` is outside the real logarithm's domain.
    

**Examples**

```julia
julia> using SymbolicAnalysis

julia> lorentz_log_barrier([0.0, 2.0]) == 0.0
true
```



<Badge type="info" class="source-link" text="source"><a href="https://github.com/SciML/SymbolicAnalysis.jl" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='SymbolicAnalysis.lorentz_nonhomogeneous_quadratic-Tuple{AbstractMatrix, AbstractVector, Real, AbstractVector}' href='#SymbolicAnalysis.lorentz_nonhomogeneous_quadratic-Tuple{AbstractMatrix, AbstractVector, Real, AbstractVector}'><span class="jlbinding">SymbolicAnalysis.lorentz_nonhomogeneous_quadratic</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
lorentz_nonhomogeneous_quadratic(A::AbstractMatrix, b::AbstractVector, c::Real, p::AbstractVector)
```


Computes the non-homogeneous quadratic function f(p) = p'Ap + b'p + c on the Lorentz model. For geodesic convexity, p'Ap must be geodesically convex and b must be in the Lorentz cone L.

**Arguments**

```julia
- `A::AbstractMatrix`: A symmetric matrix in R^((d+1)×(d+1)).
- `b::AbstractVector`: A vector in R^(d+1) which must be in the Lorentz cone.
- `c::Real`: A constant term.
- `p::AbstractVector`: A point on the Lorentz manifold.
```



<Badge type="info" class="source-link" text="source"><a href="https://github.com/SciML/SymbolicAnalysis.jl" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='SymbolicAnalysis.lorentz_transform-Tuple{AbstractMatrix, AbstractVector}' href='#SymbolicAnalysis.lorentz_transform-Tuple{AbstractMatrix, AbstractVector}'><span class="jlbinding">SymbolicAnalysis.lorentz_transform</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
lorentz_transform(O, p) -> AbstractVector
```


Apply the orthochronous Lorentz transformation `O` to the point `p`.

**Arguments**
- `O::AbstractMatrix`: element of the orthochronous Lorentz group.
    
- `p::AbstractVector`: point on the Lorentz manifold.
    

**Returns**
- The transformed point `O * p`.
    

**Throws**
- `ArgumentError`: if `O` does not preserve the Lorentz metric or the positive time direction.
    

**Examples**

```julia
julia> using SymbolicAnalysis, LinearAlgebra

julia> O = Matrix{Float64}(I, 3, 3);

julia> lorentz_transform(O, [0.0, 0.0, 1.0])
3-element Vector{Float64}:
 0.0
 0.0
 1.0
```



<Badge type="info" class="source-link" text="source"><a href="https://github.com/SciML/SymbolicAnalysis.jl" target="_blank" rel="noreferrer">source</a></Badge>

</details>

