# Porting DGCP to Python (via CVXPY) or Matlab

This guide provides practical instructions for adding Disciplined Geodesically Convex Programming (DGCP) verification to CVXPY and Matlab. Rather than building a symbolic system from scratch, the approach extends CVXPY's existing DCP infrastructure -- expression trees, atom library, sign/curvature propagation, and composition rules -- with a geodesic curvature layer.

The SymbolicAnalysis.jl implementation serves as the reference architecture.

## Why Extend CVXPY Instead of Building from Scratch

CVXPY already implements three of the four pipeline stages DGCP needs:

```
CVXPY today (3 phases):
  Expression tree  -->  Sign propagation  -->  Curvature propagation (DCP)
                              |                         |
                        atom.sign()              atom.func_curvature()
                                                 + composition rules

DGCP extension (adds 4th phase):
  Expression tree  -->  Sign propagation  -->  Curvature propagation  -->  G-Curvature propagation
                              |                         |                           |
                        atom.sign()              atom.func_curvature()       atom.g_curvature()
                                                 + DCP composition           + DGCP composition
                                                                             + fallback to DCP
```

CVXPY provides:
- **Expression trees** (`Expression` base class with `args`, recursive structure)
- **Atom base class** with `sign_from_args()`, `func_curvature()`, `monotonicity()`, `is_atom_convex()`, etc.
- **Curvature class** (`AFFINE`, `CONVEX`, `CONCAVE`, `UNKNOWN`) with arithmetic (`+` combines curvatures)
- **Monotonicity class** (`INCREASING`, `DECREASING`, `NONMONOTONIC`)
- **Composition rules** in `dcp_attr()` that check `f(g(x))` validity

DGCP adds `GCurvature`, `GMonotonicity`, a fourth propagation pass, and a registry of g-convex atoms.

---

## Step 1: Add GCurvature and GMonotonicity Enums

CVXPY defines `Curvature` and `Monotonicity` in `cvxpy/utilities/`. Add parallel geodesic types alongside them.

In SymbolicAnalysis.jl these are defined in `src/gdcp/gdcp_rules.jl`:

```julia
# Julia reference
@enum GCurvature GConvex GConcave GLinear GUnknownCurvature
@enum GMonotonicity GIncreasing GDecreasing GAnyMono
```

The Python equivalent, placed in `cvxpy/utilities/gcurvature.py`:

```python
# cvxpy/utilities/gcurvature.py
class GCurvature:
    """Geodesic curvature for manifold-valued expressions."""
    G_CONVEX = "G_CONVEX"
    G_CONCAVE = "G_CONCAVE"
    G_LINEAR = "G_LINEAR"       # Analogous to Affine in DCP
    G_UNKNOWN = "G_UNKNOWN"

    @staticmethod
    def combine(gcurvatures):
        """Combine g-curvatures under addition (same logic as DCP).

        Maps to add_gcurvature() in gdcp_rules.jl:
          - All GLinear -> GLinear
          - Mix of GConvex/GLinear -> GConvex
          - Mix of GConcave/GLinear -> GConcave
          - Any conflict -> GUnknown
        """
        has_gconvex = False
        has_gconcave = False
        for gc in gcurvatures:
            if gc == GCurvature.G_LINEAR:
                continue
            elif gc == GCurvature.G_CONVEX:
                has_gconvex = True
                if has_gconcave:
                    return GCurvature.G_UNKNOWN
            elif gc == GCurvature.G_CONCAVE:
                has_gconcave = True
                if has_gconvex:
                    return GCurvature.G_UNKNOWN
            else:
                return GCurvature.G_UNKNOWN
        if has_gconvex:
            return GCurvature.G_CONVEX
        elif has_gconcave:
            return GCurvature.G_CONCAVE
        return GCurvature.G_LINEAR

    @staticmethod
    def negate(gcurv):
        """Negate g-curvature (multiplication by negative scalar).

        Maps to mul_gcurvature() in gdcp_rules.jl.
        """
        if gcurv == GCurvature.G_CONVEX:
            return GCurvature.G_CONCAVE
        elif gcurv == GCurvature.G_CONCAVE:
            return GCurvature.G_CONVEX
        return gcurv

    @staticmethod
    def from_dcp(curvature):
        """Fall back from DCP curvature to g-curvature.

        Maps to the fallback logic in find_gcurvature() in gdcp_rules.jl:
          Convex -> GConvex, Concave -> GConcave, Affine -> GLinear
        """
        from cvxpy.utilities.curvature import Curvature
        mapping = {
            Curvature.AFFINE: GCurvature.G_LINEAR,
            Curvature.CONVEX: GCurvature.G_CONVEX,
            Curvature.CONCAVE: GCurvature.G_CONCAVE,
        }
        return mapping.get(curvature, GCurvature.G_UNKNOWN)


class GMonotonicity:
    """Geodesic monotonicity for manifold-valued expressions."""
    G_INCREASING = "G_INCREASING"
    G_DECREASING = "G_DECREASING"
    G_ANY_MONO = "G_ANY_MONO"
```

---

## Step 2: Extend the Atom Base Class

CVXPY atoms inherit from `cvxpy.atoms.atom.Atom`. Each atom defines `func_curvature()`, `sign_from_args()`, and `monotonicity()`. DGCP adds two new methods.

```python
# Add to cvxpy/atoms/atom.py (or a DGCP mixin)
class Atom(Expression):
    # ... existing methods ...

    def g_curvature(self):
        """Return the geodesic curvature of this atom.

        Default: fall back to DCP curvature via GCurvature.from_dcp().
        DGCP atoms override this to return their specific g-curvature.
        """
        return GCurvature.from_dcp(self.func_curvature())

    def g_monotonicity(self):
        """Return geodesic monotonicity (list, one per argument).

        Default: convert from DCP monotonicity.
        DGCP atoms override this.
        """
        return [GMonotonicity.G_ANY_MONO] * len(self.args)

    @property
    def manifold(self):
        """Return the manifold this atom operates on, or None for Euclidean."""
        return None
```

The default `g_curvature()` returns `GCurvature.from_dcp(self.func_curvature())`, which means every existing CVXPY atom automatically gets a valid g-curvature without modification. This mirrors the fallback in `find_gcurvature()` from `gdcp_rules.jl` (lines 181-185):

```julia
# Julia reference: when no GDCP rule exists, fall back to DCP
if !(knowngcurv) && hasdcprule(f)
    rule, args = dcprule(f, args...)
    f_curvature = rule.curvature
    f_monotonicity = rule.monotonicity
end
```

---

## Step 3: Register DGCP Atoms

Each DGCP atom is a CVXPY `Atom` subclass that overrides `g_curvature()`, `g_monotonicity()`, and `manifold`. The properties come directly from `add_gdcprule()` calls in `src/gdcp/spd.jl` and `src/gdcp/lorentz.jl`.

### SPD Manifold Atoms

```python
# cvxpy/atoms/dgcp/logdet_spd.py
import numpy as np
from cvxpy.atoms.atom import Atom
from cvxpy.utilities.gcurvature import GCurvature, GMonotonicity

class LogDetSPD(Atom):
    """log(det(X)) on the SPD manifold.

    Maps to: add_gdcprule(logdet, SymmetricPositiveDefinite,
                          AnySign, GLinear, GIncreasing)
    """
    def func_curvature(self):
        return Curvature.CONCAVE       # Standard DCP property

    def g_curvature(self):
        return GCurvature.G_LINEAR     # Key DGCP property: geodesically linear

    def g_monotonicity(self):
        return [GMonotonicity.G_INCREASING]

    @property
    def manifold(self):
        return "SymmetricPositiveDefinite"

    def sign_from_args(self):
        return (True, True)  # Can be positive or negative (AnySign)

    def numeric(self, values):
        return np.log(np.linalg.det(values[0]))


class Conjugation(Atom):
    """B' @ X @ B on the SPD manifold.

    Maps to: add_gdcprule(conjugation, SymmetricPositiveDefinite,
                          Positive, GConvex, GIncreasing)
    """
    def func_curvature(self):
        return Curvature.CONVEX

    def g_curvature(self):
        return GCurvature.G_CONVEX

    def g_monotonicity(self):
        return [GMonotonicity.G_INCREASING]

    @property
    def manifold(self):
        return "SymmetricPositiveDefinite"

    def sign_from_args(self):
        return (True, False)  # Positive

    def numeric(self, values):
        X, B = values
        return B.T @ X @ B


class TraceSPD(Atom):
    """tr(X) on the SPD manifold.

    Maps to: add_gdcprule(tr, SymmetricPositiveDefinite,
                          Positive, GConvex, GIncreasing)
    """
    def func_curvature(self):
        return Curvature.AFFINE        # Affine in Euclidean DCP

    def g_curvature(self):
        return GCurvature.G_CONVEX     # But g-convex on SPD!

    def g_monotonicity(self):
        return [GMonotonicity.G_INCREASING]

    @property
    def manifold(self):
        return "SymmetricPositiveDefinite"

    def sign_from_args(self):
        return (True, False)  # Positive on SPD

    def numeric(self, values):
        return np.trace(values[0])


class InvSPD(Atom):
    """inv(X) on the SPD manifold.

    Maps to: add_gdcprule(inv, SymmetricPositiveDefinite,
                          Positive, GConvex, GDecreasing)
    """
    def func_curvature(self):
        return Curvature.CONVEX

    def g_curvature(self):
        return GCurvature.G_CONVEX

    def g_monotonicity(self):
        return [GMonotonicity.G_DECREASING]  # Note: decreasing

    @property
    def manifold(self):
        return "SymmetricPositiveDefinite"

    def sign_from_args(self):
        return (True, False)

    def numeric(self, values):
        return np.linalg.inv(values[0])


class SDivergence(Atom):
    """S-divergence: logdet((X+Y)/2) - 0.5*logdet(X*Y).

    Maps to: add_gdcprule(sdivergence, SymmetricPositiveDefinite,
                          Positive, GConvex, GIncreasing)
    """
    def g_curvature(self):
        return GCurvature.G_CONVEX

    def g_monotonicity(self):
        return [GMonotonicity.G_INCREASING, GMonotonicity.G_INCREASING]

    @property
    def manifold(self):
        return "SymmetricPositiveDefinite"

    def sign_from_args(self):
        return (True, False)

    def numeric(self, values):
        X, Y = values
        return (np.log(np.linalg.det((X + Y) / 2))
                - 0.5 * np.log(np.linalg.det(X @ Y)))


class DistanceSPD(Atom):
    """Riemannian distance on SPD manifold.

    Maps to: add_gdcprule(distance, SymmetricPositiveDefinite,
                          Positive, GConvex, GAnyMono)
    """
    def g_curvature(self):
        return GCurvature.G_CONVEX

    def g_monotonicity(self):
        return [GMonotonicity.G_ANY_MONO, GMonotonicity.G_ANY_MONO]

    @property
    def manifold(self):
        return "SymmetricPositiveDefinite"

    def sign_from_args(self):
        return (True, False)
```

### Lorentz Manifold Atoms

```python
# cvxpy/atoms/dgcp/lorentz.py

class LorentzDistance(Atom):
    """Riemannian distance on Lorentz (hyperbolic) manifold.

    Maps to: add_gdcprule(distance, Lorentz, Positive, GConvex, GAnyMono)
    """
    def g_curvature(self):
        return GCurvature.G_CONVEX

    def g_monotonicity(self):
        return [GMonotonicity.G_ANY_MONO, GMonotonicity.G_ANY_MONO]

    @property
    def manifold(self):
        return "Lorentz"

    def sign_from_args(self):
        return (True, False)


class LorentzLogBarrier(Atom):
    """Log-barrier for Lorentz cone: -log(-1 - <a, p>_L).

    Maps to: add_gdcprule(lorentz_log_barrier, Lorentz,
                          Positive, GConvex, GIncreasing)
    """
    def g_curvature(self):
        return GCurvature.G_CONVEX

    def g_monotonicity(self):
        return [GMonotonicity.G_INCREASING]

    @property
    def manifold(self):
        return "Lorentz"

    def numeric(self, values):
        p = values[0]
        return -np.log(-1 + p[-1])


class LorentzHomogeneousQuadratic(Atom):
    """p'Ap on the Lorentz manifold (with convexity conditions on A).

    Maps to: add_gdcprule(lorentz_homogeneous_quadratic, Lorentz,
                          Positive, GConvex, GAnyMono)
    """
    def g_curvature(self):
        return GCurvature.G_CONVEX

    def g_monotonicity(self):
        return [GMonotonicity.G_ANY_MONO]

    @property
    def manifold(self):
        return "Lorentz"


class LorentzLeastSquares(Atom):
    """||y - Xp||^2 on the Lorentz manifold.

    Maps to: add_gdcprule(lorentz_least_squares, Lorentz,
                          Positive, GConvex, AnyMono)
    """
    def g_curvature(self):
        return GCurvature.G_CONVEX

    def g_monotonicity(self):
        return [GMonotonicity.G_ANY_MONO]

    @property
    def manifold(self):
        return "Lorentz"
```

The full set of atoms to register is listed in the reference table at the end of this document.

---

## Step 4: Add the G-Curvature Propagation Pass

CVXPY's DCP check lives in `cvxpy/problems/problem.py` and `cvxpy/reductions/dcp2cone/`. It walks the expression tree bottom-up and applies composition rules via `dcp_attr()`. DGCP adds a parallel `gdcp_attr()` pass.

This maps to `find_gcurvature()` in `src/gdcp/gdcp_rules.jl`, which:
1. Looks up a GDCP rule for the atom
2. If none exists, falls back to the DCP rule
3. Applies composition rules using g-curvature of arguments

```python
# cvxpy/reductions/dgcp/dgcp_attr.py
from cvxpy.utilities.gcurvature import GCurvature, GMonotonicity


def is_dgcp(problem, manifold):
    """Check if a problem is DGCP-compliant on the given manifold.

    Mirrors SymbolicAnalysis.jl's propagate_gcurvature(ex, M).
    """
    obj_gcurv = expr_gcurvature(problem.objective.expr, manifold)

    if problem.objective.NAME == "minimize":
        if obj_gcurv not in (GCurvature.G_CONVEX, GCurvature.G_LINEAR):
            return False
    elif problem.objective.NAME == "maximize":
        if obj_gcurv not in (GCurvature.G_CONCAVE, GCurvature.G_LINEAR):
            return False

    # Constraints: each must be g-convex (for <= 0) or g-linear (for == 0)
    for constr in problem.constraints:
        c_gcurv = expr_gcurvature(constr.expr, manifold)
        if isinstance(constr, ZeroConstraint):
            if c_gcurv != GCurvature.G_LINEAR:
                return False
        else:
            if c_gcurv not in (GCurvature.G_CONVEX, GCurvature.G_LINEAR):
                return False
    return True


def expr_gcurvature(expr, manifold):
    """Determine the g-curvature of an expression tree.

    Applies DGCP composition rules bottom-up, mirroring
    find_gcurvature() in gdcp_rules.jl.
    """
    from cvxpy.atoms.atom import Atom
    from cvxpy.atoms.affine.add_expr import AddExpression
    from cvxpy.atoms.affine.multiply import multiply
    from cvxpy.atoms.constants import Constant

    # Base cases: variables and constants are g-linear
    if expr.is_constant():
        return GCurvature.G_LINEAR
    if expr.is_var():
        return GCurvature.G_LINEAR

    # Addition: combine child g-curvatures
    if isinstance(expr, AddExpression):
        child_gcurvs = [expr_gcurvature(arg, manifold) for arg in expr.args]
        return GCurvature.combine(child_gcurvs)

    # Scalar multiplication: negate if constant is negative
    if isinstance(expr, multiply):
        if expr.args[0].is_constant():
            child_gcurv = expr_gcurvature(expr.args[1], manifold)
            if expr.args[0].is_nonneg():
                return child_gcurv
            elif expr.args[0].is_nonpos():
                return GCurvature.negate(child_gcurv)
            return GCurvature.G_UNKNOWN
        if expr.args[1].is_constant():
            child_gcurv = expr_gcurvature(expr.args[0], manifold)
            if expr.args[1].is_nonneg():
                return child_gcurv
            elif expr.args[1].is_nonpos():
                return GCurvature.negate(child_gcurv)
            return GCurvature.G_UNKNOWN
        return GCurvature.G_UNKNOWN

    # Atom: apply DGCP or DCP composition rules
    if isinstance(expr, Atom):
        # Check if this atom has manifold-specific g-curvature
        if expr.manifold is not None and expr.manifold != manifold:
            return GCurvature.G_UNKNOWN

        f_gcurv = expr.g_curvature()
        f_gmono = expr.g_monotonicity()

        return _apply_composition(f_gcurv, f_gmono, expr.args, manifold)

    return GCurvature.G_UNKNOWN


def _apply_composition(f_gcurv, f_gmono, args, manifold):
    """Apply DGCP composition rules for f(g1(x), g2(x), ...).

    Same rules as DCP but using g-curvature values:
      - f g-convex, g g-convex, f g-increasing -> g-convex
      - f g-convex, g g-concave, f g-decreasing -> g-convex
      - f g-concave, g g-concave, f g-increasing -> g-concave
      - f g-concave, g g-convex, f g-decreasing -> g-concave
      - f g-linear: result = g's curvature

    This mirrors the composition logic in find_gcurvature() in gdcp_rules.jl
    (lines 191-232).
    """
    if f_gcurv == GCurvature.G_LINEAR:
        # G-linear composed with anything preserves the argument's g-curvature
        if len(args) == 0:
            return GCurvature.G_LINEAR
        child_gcurvs = [expr_gcurvature(a, manifold) for a in args]
        return GCurvature.combine(child_gcurvs)

    if f_gcurv == GCurvature.G_CONVEX:
        for i, arg in enumerate(args):
            arg_gcurv = expr_gcurvature(arg, manifold)
            mono = f_gmono[i] if i < len(f_gmono) else f_gmono[-1]

            if arg_gcurv == GCurvature.G_CONVEX:
                if mono not in (GMonotonicity.G_INCREASING, "INCREASING"):
                    return GCurvature.G_UNKNOWN
            elif arg_gcurv == GCurvature.G_CONCAVE:
                if mono not in (GMonotonicity.G_DECREASING, "DECREASING"):
                    return GCurvature.G_UNKNOWN
            elif arg_gcurv == GCurvature.G_LINEAR:
                continue  # G-linear arguments are always OK
            else:
                return GCurvature.G_UNKNOWN
        return GCurvature.G_CONVEX

    if f_gcurv == GCurvature.G_CONCAVE:
        for i, arg in enumerate(args):
            arg_gcurv = expr_gcurvature(arg, manifold)
            mono = f_gmono[i] if i < len(f_gmono) else f_gmono[-1]

            if arg_gcurv == GCurvature.G_CONCAVE:
                if mono not in (GMonotonicity.G_INCREASING, "INCREASING"):
                    return GCurvature.G_UNKNOWN
            elif arg_gcurv == GCurvature.G_CONVEX:
                if mono not in (GMonotonicity.G_DECREASING, "DECREASING"):
                    return GCurvature.G_UNKNOWN
            elif arg_gcurv == GCurvature.G_LINEAR:
                continue
            else:
                return GCurvature.G_UNKNOWN
        return GCurvature.G_CONCAVE

    return GCurvature.G_UNKNOWN
```

### Special Composition Rules

The Julia `find_gcurvature()` also handles special compound expressions. For example, `logdet(conjugation(X, B))` is recognized as g-convex even though `logdet` alone is g-linear. These are hardcoded pattern matches in `gdcp_rules.jl` (lines 110-136):

```python
# Additional patterns to handle in expr_gcurvature():
def _check_special_compositions(expr, manifold):
    """Handle compound expressions from gdcp_rules.jl lines 110-136."""
    if manifold != "SymmetricPositiveDefinite":
        return None

    # logdet(conjugation(...)) -> GConvex
    # logdet(diag(...)) -> GConvex
    # logdet(affine_map(...)) -> GConvex
    # logdet(hadamard_product(...)) -> GConvex
    # logdet(X + Y) -> GConvex
    if isinstance(expr, LogDetSPD) and len(expr.args) == 1:
        inner = expr.args[0]
        if isinstance(inner, (Conjugation, DiagSPD, AffineMap,
                              HadamardProduct, AddExpression)):
            return GCurvature.G_CONVEX

    # log(tr(X)) -> GConvex
    # log(quad_form(y, X)) -> GConvex
    if isinstance(expr, log) and len(expr.args) == 1:
        inner = expr.args[0]
        if isinstance(inner, (TraceSPD, QuadFormSPD)):
            return GCurvature.G_CONVEX

    return None
```

---

## Step 5: Wire into CVXPY's Problem Interface

Add a `is_dgcp()` method to `Problem`, parallel to the existing `is_dcp()`:

```python
# Add to cvxpy/problems/problem.py
class Problem:
    # ... existing methods ...

    def is_dgcp(self, manifold="SymmetricPositiveDefinite"):
        """Check if this problem satisfies DGCP rules on the given manifold.

        This extends CVXPY's is_dcp() with a 4th verification phase:
        sign propagation -> curvature propagation -> g-curvature propagation.
        """
        from cvxpy.reductions.dgcp.dgcp_attr import is_dgcp
        return is_dgcp(self, manifold)
```

Usage:

```python
import cvxpy as cp
import numpy as np

n = 3
X = cp.Variable((n, n), symmetric=True)
Y_data = np.eye(n)

# A problem that is DGCP but not DCP
prob = cp.Problem(
    cp.Minimize(logdet_spd(conjugation(X, B)) + tr_spd(X)),
    [X >> 0]
)

print(prob.is_dcp())   # False -- logdet(B'XB) + tr(X) is not DCP
print(prob.is_dgcp())  # True  -- g-convex on SPD manifold
```

---

## Porting DGCP to Matlab

### Using Symbolic Math Toolbox

Matlab's Symbolic Math Toolbox provides `sym` objects and expression manipulation.

### Step 1: Define Curvature Types

```matlab
% dgcp_types.m
classdef CurvatureType
    enumeration
        Convex, Concave, Affine, Unknown
    end
end

classdef GCurvatureType
    enumeration
        GConvex, GConcave, GLinear, GUnknown
    end
end

classdef SignType
    enumeration
        Positive, Negative, AnySign
    end
end

classdef MonotonicityType
    enumeration
        Increasing, Decreasing, AnyMono
    end
end
```

### Step 2: Create Atom Registry

```matlab
% DCPAtomRegistry.m
classdef DCPAtomRegistry < handle
    properties
        rules containers.Map
        gdcp_rules containers.Map
    end

    methods
        function obj = DCPAtomRegistry()
            obj.rules = containers.Map('KeyType', 'char', 'ValueType', 'any');
            obj.gdcp_rules = containers.Map('KeyType', 'char', 'ValueType', 'any');
            obj.registerDefaultAtoms();
        end

        function addRule(obj, funcName, sign, curvature, monotonicity)
            rule = struct('sign', sign, ...
                         'curvature', curvature, ...
                         'monotonicity', monotonicity);
            obj.rules(funcName) = rule;
        end

        function addGDCPRule(obj, funcName, manifold, sign, gcurvature, gmonotonicity)
            rule = struct('manifold', manifold, ...
                         'sign', sign, ...
                         'gcurvature', gcurvature, ...
                         'gmonotonicity', gmonotonicity);
            obj.gdcp_rules(funcName) = rule;
        end

        function registerDefaultAtoms(obj)
            % Standard DCP atoms
            obj.addRule('exp', SignType.Positive, CurvatureType.Convex, MonotonicityType.Increasing);
            obj.addRule('log', SignType.AnySign, CurvatureType.Concave, MonotonicityType.Increasing);
            obj.addRule('abs', SignType.Positive, CurvatureType.Convex, MonotonicityType.AnyMono);
            obj.addRule('sqrt', SignType.Positive, CurvatureType.Concave, MonotonicityType.Increasing);
            obj.addRule('norm', SignType.Positive, CurvatureType.Convex, MonotonicityType.AnyMono);

            % DGCP atoms for SPD manifold
            obj.addGDCPRule('logdet', 'SPD', SignType.AnySign, ...
                           GCurvatureType.GLinear, MonotonicityType.Increasing);
            obj.addGDCPRule('trace', 'SPD', SignType.Positive, ...
                           GCurvatureType.GConvex, MonotonicityType.Increasing);
            obj.addGDCPRule('conjugation', 'SPD', SignType.Positive, ...
                           GCurvatureType.GConvex, MonotonicityType.Increasing);
            obj.addGDCPRule('inv', 'SPD', SignType.Positive, ...
                           GCurvatureType.GConvex, MonotonicityType.Decreasing);
            obj.addGDCPRule('sdivergence', 'SPD', SignType.Positive, ...
                           GCurvatureType.GConvex, MonotonicityType.Increasing);
            obj.addGDCPRule('distance', 'SPD', SignType.Positive, ...
                           GCurvatureType.GConvex, MonotonicityType.AnyMono);
        end

        function rule = getRule(obj, funcName)
            if obj.rules.isKey(funcName)
                rule = obj.rules(funcName);
            else
                rule = [];
            end
        end

        function rule = getGDCPRule(obj, funcName)
            if obj.gdcp_rules.isKey(funcName)
                rule = obj.gdcp_rules(funcName);
            else
                rule = [];
            end
        end
    end
end
```

### Step 3: G-Curvature Propagation

```matlab
% findGCurvature.m
function gcurvature = findGCurvature(expr, manifold, registry)
    % Base case
    if isnumeric(expr) || isempty(symvar(expr))
        gcurvature = GCurvatureType.GLinear;
        return;
    end

    [op, args] = getOpAndArgs(expr);

    % Handle addition
    if strcmp(op, 'plus')
        gcurvatures = arrayfun(@(a) findGCurvature(a, manifold, registry), args);
        gcurvature = combineGCurvatures(gcurvatures);
        return;
    end

    % Handle scalar multiplication
    if strcmp(op, 'times') || strcmp(op, 'mtimes')
        gcurvature = handleGMultiplication(args, manifold, registry);
        return;
    end

    % Look up GDCP rule
    rule = registry.getGDCPRule(op);
    if isempty(rule) || ~strcmp(rule.manifold, manifold)
        % Fall back to DCP curvature
        gcurvature = dcpToGdcp(findCurvature(expr, registry));
        return;
    end

    gcurvature = rule.gcurvature;
end

function gcurvature = combineGCurvatures(gcurvatures)
    if all(gcurvatures == GCurvatureType.GLinear)
        gcurvature = GCurvatureType.GLinear;
    elseif all(gcurvatures == GCurvatureType.GConvex | gcurvatures == GCurvatureType.GLinear)
        gcurvature = GCurvatureType.GConvex;
    elseif all(gcurvatures == GCurvatureType.GConcave | gcurvatures == GCurvatureType.GLinear)
        gcurvature = GCurvatureType.GConcave;
    else
        gcurvature = GCurvatureType.GUnknown;
    end
end

function gcurvature = dcpToGdcp(curvature)
    switch curvature
        case CurvatureType.Convex
            gcurvature = GCurvatureType.GConvex;
        case CurvatureType.Concave
            gcurvature = GCurvatureType.GConcave;
        case CurvatureType.Affine
            gcurvature = GCurvatureType.GLinear;
        otherwise
            gcurvature = GCurvatureType.GUnknown;
    end
end
```

### Step 4: Main Analysis Function

```matlab
% analyze.m
function result = analyze(expr, manifold)
    arguments
        expr sym
        manifold string = ""
    end

    registry = DCPAtomRegistry();
    curvature = findCurvature(expr, registry);
    sign = propagateSign(expr, registry);

    if manifold ~= ""
        gcurvature = findGCurvature(expr, manifold, registry);
    else
        gcurvature = [];
    end

    result = struct('curvature', curvature, ...
                    'sign', sign, ...
                    'gcurvature', gcurvature);
end
```

### Example Usage

```matlab
syms x y positive

registry = DCPAtomRegistry();

% DCP analysis
expr1 = exp(x) + exp(y);
result1 = analyze(expr1);
fprintf('exp(x) + exp(y): %s\n', string(result1.curvature));

% DGCP analysis
syms X [3 3] matrix
result2 = analyze(trace(X), 'SPD');
fprintf('trace(X) on SPD: %s\n', string(result2.gcurvature));
```

---

## DGCP Atoms Reference

### Symmetric Positive Definite (SPD) Manifold

| Atom | Sign | G-Curvature | G-Monotonicity | Julia Function |
|------|------|-------------|----------------|----------------|
| logdet(X) | AnySign | GLinear | GIncreasing | `LinearAlgebra.logdet` |
| tr(X) | Positive | GConvex | GIncreasing | `LinearAlgebra.tr` |
| conjugation(X, B) = B'XB | Positive | GConvex | GIncreasing | `conjugation` |
| diag(X) | Positive | GConvex | GIncreasing | `LinearAlgebra.diag` |
| inv(X) | Positive | GConvex | GDecreasing | `inv` |
| quad_form(y, X) = y'Xy | Positive | GConvex | GIncreasing | `quad_form` |
| log_quad_form(y, X) | Positive | GConvex | GIncreasing | `log_quad_form` |
| eigmax(X) | Positive | GConvex | GIncreasing | `LinearAlgebra.eigmax` |
| schatten_norm(X, p) | Positive | GConvex | GIncreasing | `schatten_norm` |
| sum_log_eigmax(X, k) | Positive | GConvex | GIncreasing | `sum_log_eigmax` |
| affine_map(f, X, B, Y) | Positive | GConvex | GIncreasing | `affine_map` |
| hadamard_product(X, B) | Positive | GConvex | GIncreasing | `hadamard_product` |
| sdivergence(X, Y) | Positive | GConvex | GIncreasing | `sdivergence` |
| distance(M, X, Y) | Positive | GConvex | GAnyMono | `Manifolds.distance` |

### Lorentz Manifold (Hyperbolic Space)

| Atom | Sign | G-Curvature | G-Monotonicity | Julia Function |
|------|------|-------------|----------------|----------------|
| distance(M, p, q) | Positive | GConvex | GAnyMono | `Manifolds.distance` |
| lorentz_log_barrier(p) | Positive | GConvex | GIncreasing | `lorentz_log_barrier` |
| lorentz_homogeneous_quadratic(A, p) | Positive | GConvex | GAnyMono | `lorentz_homogeneous_quadratic` |
| lorentz_homogeneous_diagonal(a, p) | Positive | GConvex | GAnyMono | `lorentz_homogeneous_diagonal` |
| lorentz_least_squares(X, y, p) | Positive | GConvex | GAnyMono | `lorentz_least_squares` |

---

## Implementation Checklist

When extending CVXPY with DGCP:

- [ ] **Add `GCurvature` and `GMonotonicity`** in `cvxpy/utilities/` alongside existing `Curvature`/`Monotonicity`
- [ ] **Extend `Atom` base class** with `g_curvature()`, `g_monotonicity()`, and `manifold` (defaults fall back to DCP)
- [ ] **Implement DGCP atom subclasses** for SPD atoms (logdet, conjugation, tr, inv, distance, sdivergence, etc.)
- [ ] **Implement DGCP atom subclasses** for Lorentz atoms (distance, log_barrier, homogeneous_quadratic, etc.)
- [ ] **Add `expr_gcurvature()` propagation** that walks the tree bottom-up applying DGCP composition rules
- [ ] **Handle special compositions** (logdet of conjugation, log of trace, etc.) as pattern matches
- [ ] **Add `Problem.is_dgcp(manifold)`** as the public API
- [ ] **Test with known expressions** from the paper's experiments (SPD matrix means, Lorentz regression)

## Key Design Decisions

1. **Fallback to DCP**: Every existing CVXPY atom gets automatic DGCP support via `GCurvature.from_dcp()`. Only atoms with manifold-specific properties need overrides.
2. **Manifold as parameter**: The manifold tag on atoms prevents mixing SPD and Lorentz atoms in the same expression.
3. **Separate pass, not interleaved**: G-curvature propagation runs as a distinct 4th phase after DCP, not interleaved with it. This keeps CVXPY's existing DCP logic untouched.
4. **Composition rules are identical**: The DCP and DGCP composition rule tables have the same structure -- only the enum values differ (`Convex`/`GConvex`, `Increasing`/`GIncreasing`). This is by design in the theory.
5. **Return `G_UNKNOWN` rather than throwing**: Mirrors Julia's approach of returning `GUnknownCurvature` when rules do not apply, rather than raising exceptions.
