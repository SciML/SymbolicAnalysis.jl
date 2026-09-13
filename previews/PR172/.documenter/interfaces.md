---
---

# Interfaces {#Interfaces}

The interfaces in this page are the supported extension points for packages that add symbolic atoms to SymbolicAnalysis. They describe the contract at the same level as the generic analyzer: register the symbolic operation, register its property rule, and call [`analyze`](/index#SymbolicAnalysis.analyze) to inspect the result.

`dcprule`, `gdcprule`, the rule dictionaries, and the propagation passes are implementation details. They are used to implement the interfaces below and are not a supported basis for downstream code.

## Analysis data {#Analysis-data}
<details class='jldocstring custom-block' open>
<summary><a id='SymbolicAnalysis.Sign' href='#SymbolicAnalysis.Sign'><span class="jlbinding">SymbolicAnalysis.Sign</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



```julia
Sign
```


Sign classification propagated by the DCP analyzer.

**Values**
- `Positive`: the expression is known to be positive.
  
- `Negative`: the expression is known to be negative.
  
- `AnySign`: no sign restriction is known.
  

These values are accepted by [`add_dcprule`](/interfaces#SymbolicAnalysis.add_dcprule) and are returned by `AnalysisResult.sign`.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/SciML/SymbolicAnalysis.jl" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='SymbolicAnalysis.Curvature' href='#SymbolicAnalysis.Curvature'><span class="jlbinding">SymbolicAnalysis.Curvature</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



```julia
Curvature
```


Euclidean curvature classification propagated by the DCP analyzer.

**Values**
- `Convex`: the expression is convex.
  
- `Concave`: the expression is concave.
  
- `Affine`: the expression is both convex and concave.
  
- `UnknownCurvature`: no supported curvature rule applies.
  

These values are accepted by [`add_dcprule`](/interfaces#SymbolicAnalysis.add_dcprule) and are returned by `AnalysisResult.curvature`.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/SciML/SymbolicAnalysis.jl" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='SymbolicAnalysis.Monotonicity' href='#SymbolicAnalysis.Monotonicity'><span class="jlbinding">SymbolicAnalysis.Monotonicity</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



```julia
Monotonicity
```


Monotonicity classification for an argument of a DCP atom.

**Values**
- `Increasing`: the atom is nondecreasing in the argument.
  
- `Decreasing`: the atom is nonincreasing in the argument.
  
- `AnyMono`: the atom's monotonicity is unrestricted for the rule.
  

These values are accepted by [`add_dcprule`](/interfaces#SymbolicAnalysis.add_dcprule).


<Badge type="info" class="source-link" text="source"><a href="https://github.com/SciML/SymbolicAnalysis.jl" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='SymbolicAnalysis.GCurvature' href='#SymbolicAnalysis.GCurvature'><span class="jlbinding">SymbolicAnalysis.GCurvature</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



```julia
GCurvature
```


Geodesic-curvature classification propagated by the gDCP analyzer.

**Values**
- `GConvex`: geodesically convex.
  
- `GConcave`: geodesically concave.
  
- `GLinear`: geodesically affine.
  
- `GUnknownCurvature`: no supported geodesic-curvature rule applies.
  

These values are accepted by [`add_gdcprule`](/interfaces#SymbolicAnalysis.add_gdcprule) and are returned in `AnalysisResult.gcurvature` when manifold analysis is requested.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/SciML/SymbolicAnalysis.jl" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='SymbolicAnalysis.GMonotonicity' href='#SymbolicAnalysis.GMonotonicity'><span class="jlbinding">SymbolicAnalysis.GMonotonicity</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



```julia
GMonotonicity
```


Monotonicity classification for an argument of a geodesic DCP atom.

**Values**
- `GIncreasing`: the atom is nondecreasing along the manifold argument.
  
- `GDecreasing`: the atom is nonincreasing along the manifold argument.
  
- `GAnyMono`: the atom's monotonicity is unrestricted for the rule.
  

These values are accepted by [`add_gdcprule`](/interfaces#SymbolicAnalysis.add_gdcprule).


<Badge type="info" class="source-link" text="source"><a href="https://github.com/SciML/SymbolicAnalysis.jl" target="_blank" rel="noreferrer">source</a></Badge>

</details>


## DCP extension {#DCP-extension}

Use [`add_dcprule`](/interfaces#SymbolicAnalysis.add_dcprule) to register the Euclidean sign, curvature, and monotonicity of a symbolic atom. The `domain` argument describes where the rule is valid. Symbolic arguments may carry a `VarDomain` metadata entry so that rules with different domain restrictions can be selected.
<details class='jldocstring custom-block' open>
<summary><a id='SymbolicAnalysis.add_dcprule' href='#SymbolicAnalysis.add_dcprule'><span class="jlbinding">SymbolicAnalysis.add_dcprule</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
add_dcprule(f, domain, sign, curvature, monotonicity)
```


Register a Euclidean disciplined-convex-programming (DCP) rule for the symbolic atom `f`. This is the extension point for packages that add a new symbolic atom whose sign, curvature, and argument monotonicity are known.

Register the symbolic operation with `Symbolics.@register_symbolic` before registering its DCP rule. The function object passed to `add_dcprule` must be the same operation stored in the symbolic term.

**Arguments**
- `f`: function used as the operation of the symbolic atom.
  
- `domain`: a `DomainSets.Domain` for a single argument, or a tuple of domains for a multi-argument atom. Domain checks use `VarDomain` metadata when it is attached to the symbolic arguments.
  
- `sign::Sign`: sign guaranteed by the atom on its declared domain.
  
- `curvature::Curvature`: Euclidean curvature of the atom.
  
- `monotonicity`: one [`Monotonicity`](/interfaces#SymbolicAnalysis.Monotonicity) value, a tuple with one value per argument, or a function that computes the value from an argument. A single value applies to every argument.
  

**Returns**

The registered rule descriptor. Registering another rule for the same function appends a descriptor so that distinct domains can be supported.

**Examples**

```julia
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


This registration hook is public for packages extending SymbolicAnalysis. The rule dictionaries and the `dcprule` dispatch function are implementation details and should not be accessed directly.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/SciML/SymbolicAnalysis.jl" target="_blank" rel="noreferrer">source</a></Badge>

</details>


## Geodesic DCP extension {#Geodesic-DCP-extension}

Use [`add_gdcprule`](/interfaces#SymbolicAnalysis.add_gdcprule) to register the geodesic curvature of an atom on a particular manifold type. A gDCP rule is selected when `analyze` is called with that manifold.
<details class='jldocstring custom-block' open>
<summary><a id='SymbolicAnalysis.add_gdcprule' href='#SymbolicAnalysis.add_gdcprule'><span class="jlbinding">SymbolicAnalysis.add_gdcprule</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
add_gdcprule(f, manifold, sign, curvature, monotonicity)
```


Register a geodesic disciplined-convex-programming (gDCP) rule for the symbolic atom `f` on a manifold type. This is the extension point for packages that add atoms with known geodesic curvature.

Register the symbolic operation with `Symbolics.@register_symbolic` before registering its rule. The function object passed here must be the operation stored in the symbolic term.

**Arguments**
- `f`: function used as the operation of the symbolic atom.
  
- `manifold`: manifold type on which the rule applies, such as `Manifolds.SymmetricPositiveDefinite` or `Manifolds.Lorentz`.
  
- `sign::Sign`: sign guaranteed by the atom.
  
- `curvature::GCurvature`: geodesic curvature of the atom.
  
- `monotonicity::GMonotonicity`: geodesic monotonicity of the atom. A tuple may be supplied for a multi-argument atom.
  

**Returns**

The registered geodesic-rule descriptor.

**Examples**

```julia
julia> using SymbolicAnalysis, Symbolics, Manifolds

julia> test_gdcp_atom(p) = p[1];

julia> @register_symbolic test_gdcp_atom(p::Vector{Num})

julia> SymbolicAnalysis.add_gdcprule(
           test_gdcp_atom, Manifolds.Lorentz, SymbolicAnalysis.Positive,
           SymbolicAnalysis.GConvex, SymbolicAnalysis.GIncreasing
       );

julia> @variables p[1:3];

julia> result = analyze(test_gdcp_atom(p), Manifolds.Lorentz(2));

julia> result.gcurvature == SymbolicAnalysis.GConvex
true
```


This registration hook is public for packages extending SymbolicAnalysis. The rule dictionary and `gdcprule` lookup function are implementation details and should not be accessed directly.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/SciML/SymbolicAnalysis.jl" target="_blank" rel="noreferrer">source</a></Badge>

</details>


## Extension rules {#Extension-rules}

An extension package should:
2. define the ordinary numerical operation;
  
3. register its symbolic operation with `Symbolics.@register_symbolic`;
  
4. register its DCP or gDCP rule with the corresponding function above; and
  
5. verify the generic path by constructing a symbolic expression and calling `analyze`.
  

The operation supplied to the rule registration function must be identical to the operation stored in the symbolic term. Rule registration is global to the loaded Julia process, so extension packages should register each atom once at module load time and use distinct function objects for distinct atoms.
