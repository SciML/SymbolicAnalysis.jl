# Atoms for DCP and DGCP

This page is intended to be a reference for the atoms that are currently implemented in this package with their respective properties. As much as possible atoms are created with functions from base, standard libraries and popular packages, but we also inherit a few functions from the CVX family of packages such as `quad_form`, `quad_over_lin` etc. and also introduce some new functions in this package. Description of all such special functions implemented in this package is available in the [Special functions](@ref) section of the documentation.

Use the **search box** to filter by name, or the **dropdown menus** to filter by curvature, sign, or monotonicity. Click any **column header** to sort.

## DCP Atoms

```@raw html
<AtomTable src="dcp" />
```

### Special Cases for `^(x, i)`

```@raw html
<AtomTable src="dcp-power" />
```

## DGCP Atoms (SPD)

```@raw html
<AtomTable src="dgcp-spd" />
```

## DGCP Atoms (Lorentz)

```@raw html
<AtomTable src="dgcp-lorentz" />
```

!!! note
    `lorentz_transform` does not have specific geodesic curvature properties by itself, 
    but it preserves geodesic convexity when applied to geodesically convex functions.
