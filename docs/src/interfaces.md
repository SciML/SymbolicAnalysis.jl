# Interfaces

The interfaces in this page are the supported extension points for packages
that add symbolic atoms to SymbolicAnalysis. They describe the contract at the
same level as the generic analyzer: register the symbolic operation, register
its property rule, and call [`analyze`](@ref) to inspect the result.

`dcprule`, `gdcprule`, the rule dictionaries, and the propagation passes are
implementation details. They are used to implement the interfaces below and
are not a supported basis for downstream code.

## Analysis data

```@docs
SymbolicAnalysis.Sign
SymbolicAnalysis.Curvature
SymbolicAnalysis.Monotonicity
SymbolicAnalysis.GCurvature
SymbolicAnalysis.GMonotonicity
```

## DCP extension

Use [`add_dcprule`](@ref SymbolicAnalysis.add_dcprule) to register the
Euclidean sign, curvature, and monotonicity of a symbolic atom. The `domain`
argument describes where the rule is valid. Symbolic arguments may carry a
`VarDomain` metadata entry so that rules with different domain restrictions
can be selected.

```@docs
SymbolicAnalysis.add_dcprule
```

## Geodesic DCP extension

Use [`add_gdcprule`](@ref SymbolicAnalysis.add_gdcprule) to register the
geodesic curvature of an atom on a particular manifold type. A gDCP rule is
selected when `analyze` is called with that manifold.

```@docs
SymbolicAnalysis.add_gdcprule
```

## Extension rules

An extension package should:

1. define the ordinary numerical operation;
2. register its symbolic operation with `Symbolics.@register_symbolic`;
3. register its DCP or gDCP rule with the corresponding function above; and
4. verify the generic path by constructing a symbolic expression and calling
   `analyze`.

The operation supplied to the rule registration function must be identical to
the operation stored in the symbolic term. Rule registration is global to the
loaded Julia process, so extension packages should register each atom once at
module load time and use distinct function objects for distinct atoms.
