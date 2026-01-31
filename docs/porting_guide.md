# Porting DGCP to Python or Matlab

> **Verification Note**: This document was verified against the source code on 2026-01-30.
> The architecture description, enumerations, and composition rules have been confirmed to match:
> - `src/gdcp/gdcp_rules.jl` (GCurvature, GMonotonicity enums, find_gcurvature, propagate_gcurvature)
> - `src/rules.jl` (Sign, Curvature, Monotonicity enums, find_curvature, propagate_curvature)
> - `src/gdcp/spd.jl` and `src/gdcp/lorentz.jl` (atom registrations)

This guide provides practical instructions for implementing Disciplined Geodesically Convex Programming (DGCP) in Python or Matlab. The SymbolicAnalysis.jl implementation serves as the reference architecture.

## Architecture Overview

DGCP verification follows a four-stage pipeline:

```
Expression → Canonize → Sign Propagation → Curvature Propagation → G-Curvature Propagation
                ↓              ↓                    ↓                        ↓
          Pattern rewrite  Metadata attach    DCP rules apply         DGCP rules apply
```

### Core Components

1. **Expression Tree Representation**: Symbolic expressions as trees with operations and arguments
2. **Metadata System**: Attach curvature/sign/monotonicity properties to expression nodes
3. **Atom Registry**: Dictionary mapping functions to their DCP/DGCP properties
4. **Rule-Based Rewriting**: Tree traversal applying composition rules
5. **Curvature Propagation**: Bottom-up inference following DCP composition rules

### Key Enumerations

```
Sign:        Positive | Negative | AnySign
Curvature:   Convex | Concave | Affine | UnknownCurvature
GCurvature:  GConvex | GConcave | GLinear | GUnknownCurvature
Monotonicity: Increasing | Decreasing | AnyMono
GMonotonicity: GIncreasing | GDecreasing | GAnyMono
```

### Composition Rules (DCP/DGCP)

For a composite function `f(g(x))`:

| f curvature | g curvature | f monotonicity | Result |
|-------------|-------------|----------------|--------|
| Convex | Convex | Increasing | Convex |
| Convex | Concave | Decreasing | Convex |
| Concave | Concave | Increasing | Concave |
| Concave | Convex | Decreasing | Concave |
| Affine | Any | Any | Same as g |

The same rules apply for geodesic curvature (GConvex, GConcave, GLinear).

---

## Porting DGCP to Python

### Recommended Library: SymPy

SymPy provides expression trees, pattern matching, and a metadata system that maps well to the Julia implementation.

### Step 1: Define Enumerations

```python
from enum import Enum, auto

class Sign(Enum):
    POSITIVE = auto()
    NEGATIVE = auto()
    ANY_SIGN = auto()

class Curvature(Enum):
    CONVEX = auto()
    CONCAVE = auto()
    AFFINE = auto()
    UNKNOWN = auto()

class GCurvature(Enum):
    G_CONVEX = auto()
    G_CONCAVE = auto()
    G_LINEAR = auto()
    G_UNKNOWN = auto()

class Monotonicity(Enum):
    INCREASING = auto()
    DECREASING = auto()
    ANY_MONO = auto()

class GMonotonicity(Enum):
    G_INCREASING = auto()
    G_DECREASING = auto()
    G_ANY_MONO = auto()
```

### Step 2: Create the Atom Registry

```python
from dataclasses import dataclass
from typing import Dict, Tuple, Callable, Any, Union
import sympy as sp

@dataclass
class DCPRule:
    """DCP rule for a function atom."""
    sign: Sign
    curvature: Curvature
    monotonicity: Tuple[Monotonicity, ...]  # One per argument

@dataclass
class GDCPRule:
    """GDCP rule for a geodesically convex atom."""
    manifold: str  # e.g., "SymmetricPositiveDefinite", "Lorentz"
    sign: Sign
    gcurvature: GCurvature
    gmonotonicity: Tuple[GMonotonicity, ...]

# DCP atom registry
dcp_rules: Dict[Callable, DCPRule] = {}

# GDCP atom registry
gdcp_rules: Dict[Callable, GDCPRule] = {}

def add_dcp_rule(func: Callable, sign: Sign, curvature: Curvature,
                 monotonicity: Union[Monotonicity, Tuple[Monotonicity, ...]]):
    """Register a DCP rule for a function."""
    if not isinstance(monotonicity, tuple):
        monotonicity = (monotonicity,)
    dcp_rules[func] = DCPRule(sign, curvature, monotonicity)

def add_gdcp_rule(func: Callable, manifold: str, sign: Sign,
                  gcurvature: GCurvature,
                  gmonotonicity: Union[GMonotonicity, Tuple[GMonotonicity, ...]]):
    """Register a GDCP rule for a function."""
    if not isinstance(gmonotonicity, tuple):
        gmonotonicity = (gmonotonicity,)
    gdcp_rules[func] = GDCPRule(manifold, sign, gcurvature, gmonotonicity)
```

### Step 3: Register Atom Rules

```python
import numpy as np

# Standard DCP atoms
add_dcp_rule(sp.exp, Sign.POSITIVE, Curvature.CONVEX, Monotonicity.INCREASING)
add_dcp_rule(sp.log, Sign.ANY_SIGN, Curvature.CONCAVE, Monotonicity.INCREASING)
add_dcp_rule(sp.Abs, Sign.POSITIVE, Curvature.CONVEX, Monotonicity.ANY_MONO)
add_dcp_rule(sp.sqrt, Sign.POSITIVE, Curvature.CONCAVE, Monotonicity.INCREASING)

# DGCP atoms for SPD manifold
def logdet(X):
    """Log-determinant of a matrix."""
    return np.log(np.linalg.det(X))

def conjugation(X, B):
    """Conjugation B' @ X @ B."""
    return B.T @ X @ B

def trace(X):
    """Matrix trace."""
    return np.trace(X)

add_gdcp_rule(logdet, "SymmetricPositiveDefinite", Sign.ANY_SIGN,
              GCurvature.G_LINEAR, GMonotonicity.G_INCREASING)
add_gdcp_rule(conjugation, "SymmetricPositiveDefinite", Sign.POSITIVE,
              GCurvature.G_CONVEX, GMonotonicity.G_INCREASING)
add_gdcp_rule(trace, "SymmetricPositiveDefinite", Sign.POSITIVE,
              GCurvature.G_CONVEX, GMonotonicity.G_INCREASING)
```

### Step 4: Expression Tree Traversal

```python
from typing import Optional

class ExpressionNode:
    """Wrapper for SymPy expressions with curvature metadata."""

    def __init__(self, expr: sp.Expr):
        self.expr = expr
        self.sign: Optional[Sign] = None
        self.curvature: Optional[Curvature] = None
        self.gcurvature: Optional[GCurvature] = None

    @property
    def is_atom(self) -> bool:
        """Check if this is a leaf node (symbol or number)."""
        return self.expr.is_Symbol or self.expr.is_Number

    @property
    def operation(self) -> Optional[Callable]:
        """Get the operation (function) of this node."""
        if self.is_atom:
            return None
        return type(self.expr)

    @property
    def arguments(self) -> list:
        """Get the arguments of this node."""
        if self.is_atom:
            return []
        return [ExpressionNode(arg) for arg in self.expr.args]


def find_curvature(node: ExpressionNode) -> Curvature:
    """
    Recursively determine the curvature of an expression.
    Implements DCP composition rules.
    """
    # Base case: symbols and numbers are affine
    if node.is_atom:
        return Curvature.AFFINE

    op = node.operation
    args = node.arguments

    # Handle addition: preserves curvature if all same type
    if op == sp.Add:
        curvs = [find_curvature(arg) for arg in args]
        if all(c == Curvature.AFFINE for c in curvs):
            return Curvature.AFFINE
        if all(c in (Curvature.CONVEX, Curvature.AFFINE) for c in curvs):
            return Curvature.CONVEX
        if all(c in (Curvature.CONCAVE, Curvature.AFFINE) for c in curvs):
            return Curvature.CONCAVE
        return Curvature.UNKNOWN

    # Handle multiplication: only valid if at most one non-constant
    if op == sp.Mul:
        non_constants = [arg for arg in args if not arg.expr.is_Number]
        if len(non_constants) > 1:
            return Curvature.UNKNOWN
        if len(non_constants) == 0:
            return Curvature.AFFINE

        # Get the non-constant's curvature
        nc = non_constants[0]
        curv = find_curvature(nc)

        # Check sign of constant multiplier
        constants = [arg.expr for arg in args if arg.expr.is_Number]
        const_prod = sp.prod(constants) if constants else 1

        if const_prod < 0:
            # Flip curvature
            if curv == Curvature.CONVEX:
                return Curvature.CONCAVE
            elif curv == Curvature.CONCAVE:
                return Curvature.CONVEX
        return curv

    # Look up DCP rule for this operation
    # Note: Need to map SymPy type to registered function
    func = _get_registered_function(op)
    if func is None or func not in dcp_rules:
        return Curvature.UNKNOWN

    rule = dcp_rules[func]
    f_curv = rule.curvature
    f_mono = rule.monotonicity

    # Apply composition rules
    if f_curv == Curvature.AFFINE:
        # Affine composed with anything preserves inner curvature
        return find_curvature(args[0]) if args else Curvature.AFFINE

    if f_curv == Curvature.CONVEX:
        # Check all arguments satisfy composition rule
        for i, arg in enumerate(args):
            arg_curv = find_curvature(arg)
            mono = f_mono[i] if i < len(f_mono) else f_mono[-1]

            if arg_curv == Curvature.CONVEX and mono != Monotonicity.INCREASING:
                return Curvature.UNKNOWN
            if arg_curv == Curvature.CONCAVE and mono != Monotonicity.DECREASING:
                return Curvature.UNKNOWN
            if arg_curv == Curvature.UNKNOWN:
                return Curvature.UNKNOWN
        return Curvature.CONVEX

    if f_curv == Curvature.CONCAVE:
        for i, arg in enumerate(args):
            arg_curv = find_curvature(arg)
            mono = f_mono[i] if i < len(f_mono) else f_mono[-1]

            if arg_curv == Curvature.CONCAVE and mono != Monotonicity.INCREASING:
                return Curvature.UNKNOWN
            if arg_curv == Curvature.CONVEX and mono != Monotonicity.DECREASING:
                return Curvature.UNKNOWN
            if arg_curv == Curvature.UNKNOWN:
                return Curvature.UNKNOWN
        return Curvature.CONCAVE

    return Curvature.UNKNOWN


def _get_registered_function(sympy_type):
    """Map SymPy expression type to registered function."""
    type_map = {
        sp.exp: sp.exp,
        sp.log: sp.log,
        sp.Abs: sp.Abs,
        sp.sqrt: sp.sqrt,
    }
    return type_map.get(sympy_type)
```

### Step 5: GDCP Analysis (Geodesic Curvature)

```python
def find_gcurvature(node: ExpressionNode, manifold: str) -> GCurvature:
    """
    Determine geodesic curvature for manifold-valued expressions.
    """
    if node.is_atom:
        return GCurvature.G_LINEAR

    op = node.operation
    args = node.arguments

    # Handle addition
    if op == sp.Add:
        gcurvs = [find_gcurvature(arg, manifold) for arg in args]
        if all(gc == GCurvature.G_LINEAR for gc in gcurvs):
            return GCurvature.G_LINEAR
        if all(gc in (GCurvature.G_CONVEX, GCurvature.G_LINEAR) for gc in gcurvs):
            return GCurvature.G_CONVEX
        if all(gc in (GCurvature.G_CONCAVE, GCurvature.G_LINEAR) for gc in gcurvs):
            return GCurvature.G_CONCAVE
        return GCurvature.G_UNKNOWN

    # Handle multiplication (scalar * expression)
    if op == sp.Mul:
        non_constants = [arg for arg in args if not arg.expr.is_Number]
        if len(non_constants) > 1:
            return GCurvature.G_UNKNOWN
        if len(non_constants) == 0:
            return GCurvature.G_LINEAR

        nc = non_constants[0]
        gcurv = find_gcurvature(nc, manifold)

        constants = [arg.expr for arg in args if arg.expr.is_Number]
        const_prod = sp.prod(constants) if constants else 1

        if const_prod < 0:
            if gcurv == GCurvature.G_CONVEX:
                return GCurvature.G_CONCAVE
            elif gcurv == GCurvature.G_CONCAVE:
                return GCurvature.G_CONVEX
        return gcurv

    # Look up GDCP rule
    func = _get_registered_gdcp_function(op)
    if func is None or func not in gdcp_rules:
        # Fall back to DCP rule if available
        return _fallback_to_dcp(node, manifold)

    rule = gdcp_rules[func]
    if rule.manifold != manifold:
        return GCurvature.G_UNKNOWN

    return rule.gcurvature


def _get_registered_gdcp_function(sympy_type):
    """Map to registered GDCP function."""
    # Custom mapping for matrix operations
    return None  # Implement based on your registered functions


def _fallback_to_dcp(node: ExpressionNode, manifold: str) -> GCurvature:
    """Use standard DCP curvature when no GDCP rule exists."""
    curv = find_curvature(node)
    curv_map = {
        Curvature.CONVEX: GCurvature.G_CONVEX,
        Curvature.CONCAVE: GCurvature.G_CONCAVE,
        Curvature.AFFINE: GCurvature.G_LINEAR,
        Curvature.UNKNOWN: GCurvature.G_UNKNOWN,
    }
    return curv_map.get(curv, GCurvature.G_UNKNOWN)
```

### Step 6: Main Analysis Function

```python
@dataclass
class AnalysisResult:
    """Result of DGCP analysis."""
    curvature: Curvature
    sign: Sign
    gcurvature: Optional[GCurvature] = None

def analyze(expr: sp.Expr, manifold: Optional[str] = None) -> AnalysisResult:
    """
    Analyze a symbolic expression for DCP/DGCP compliance.

    Args:
        expr: A SymPy expression
        manifold: Optional manifold name for GDCP analysis
                  ("SymmetricPositiveDefinite" or "Lorentz")

    Returns:
        AnalysisResult with curvature, sign, and optionally gcurvature
    """
    node = ExpressionNode(expr)

    # Step 1: Propagate sign
    sign = propagate_sign(node)

    # Step 2: Determine curvature
    curvature = find_curvature(node)

    # Step 3: Determine geodesic curvature if manifold specified
    gcurvature = None
    if manifold is not None:
        gcurvature = find_gcurvature(node, manifold)

    return AnalysisResult(curvature, sign, gcurvature)


def propagate_sign(node: ExpressionNode) -> Sign:
    """Propagate sign through the expression tree."""
    if node.expr.is_Number:
        return Sign.POSITIVE if node.expr > 0 else Sign.NEGATIVE
    if node.is_atom:
        return Sign.ANY_SIGN

    op = node.operation
    args = node.arguments

    if op == sp.Add:
        signs = [propagate_sign(arg) for arg in args]
        if all(s == Sign.POSITIVE for s in signs):
            return Sign.POSITIVE
        if all(s == Sign.NEGATIVE for s in signs):
            return Sign.NEGATIVE
        return Sign.ANY_SIGN

    if op == sp.Mul:
        signs = [propagate_sign(arg) for arg in args]
        neg_count = sum(1 for s in signs if s == Sign.NEGATIVE)
        if any(s == Sign.ANY_SIGN for s in signs):
            return Sign.ANY_SIGN
        return Sign.NEGATIVE if neg_count % 2 == 1 else Sign.POSITIVE

    # Look up rule for sign
    func = _get_registered_function(op)
    if func and func in dcp_rules:
        return dcp_rules[func].sign

    return Sign.ANY_SIGN
```

### Complete Example Usage

```python
import sympy as sp

# Define symbolic variables
x = sp.Symbol('x', positive=True)
y = sp.Symbol('y', positive=True)

# Example 1: DCP analysis
expr1 = sp.exp(x) + sp.log(y)
result1 = analyze(expr1)
print(f"exp(x) + log(y): curvature={result1.curvature}")
# Output: curvature=Curvature.UNKNOWN (convex + concave)

expr2 = sp.exp(x) + sp.exp(y)
result2 = analyze(expr2)
print(f"exp(x) + exp(y): curvature={result2.curvature}")
# Output: curvature=Curvature.CONVEX

# Example 2: DGCP analysis for SPD manifold
# For matrix expressions, you would extend with numpy/scipy
result3 = analyze(expr2, manifold="SymmetricPositiveDefinite")
print(f"DGCP result: gcurvature={result3.gcurvature}")
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
            % Add a DCP rule for a function
            rule = struct('sign', sign, ...
                         'curvature', curvature, ...
                         'monotonicity', monotonicity);
            obj.rules(funcName) = rule;
        end

        function addGDCPRule(obj, funcName, manifold, sign, gcurvature, gmonotonicity)
            % Add a GDCP rule for a function
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

### Step 3: Expression Tree Analysis

```matlab
% findCurvature.m
function curvature = findCurvature(expr, registry)
    % Find the curvature of a symbolic expression
    %
    % Args:
    %   expr: A symbolic expression (sym)
    %   registry: DCPAtomRegistry instance
    %
    % Returns:
    %   curvature: CurvatureType enumeration value

    % Base case: numbers are affine
    if isnumeric(expr) || isempty(symvar(expr))
        curvature = CurvatureType.Affine;
        return;
    end

    % Get the operation and arguments
    [op, args] = getOpAndArgs(expr);

    % Handle addition
    if strcmp(op, 'plus')
        curvatures = arrayfun(@(a) findCurvature(a, registry), args);
        curvature = combineCurvatures(curvatures);
        return;
    end

    % Handle multiplication
    if strcmp(op, 'times') || strcmp(op, 'mtimes')
        curvature = handleMultiplication(args, registry);
        return;
    end

    % Look up rule for this operation
    rule = registry.getRule(op);
    if isempty(rule)
        curvature = CurvatureType.Unknown;
        return;
    end

    % Apply composition rules
    curvature = applyCompositionRules(rule, args, registry);
end

function [op, args] = getOpAndArgs(expr)
    % Extract operation and arguments from symbolic expression
    str = char(expr);

    % Try to identify the outermost operation
    if contains(str, '+')
        op = 'plus';
        args = children(expr);
    elseif contains(str, '*')
        op = 'times';
        args = children(expr);
    else
        % Function call
        op = func2str(symFunType(expr));
        args = argnames(expr);
    end
end

function curvature = combineCurvatures(curvatures)
    % Combine curvatures for addition
    if all(curvatures == CurvatureType.Affine)
        curvature = CurvatureType.Affine;
    elseif all(curvatures == CurvatureType.Convex | curvatures == CurvatureType.Affine)
        curvature = CurvatureType.Convex;
    elseif all(curvatures == CurvatureType.Concave | curvatures == CurvatureType.Affine)
        curvature = CurvatureType.Concave;
    else
        curvature = CurvatureType.Unknown;
    end
end

function curvature = handleMultiplication(args, registry)
    % Handle multiplication - at most one non-constant allowed
    nonConstants = [];
    constProd = 1;

    for i = 1:length(args)
        if isnumeric(args(i)) || isempty(symvar(args(i)))
            constProd = constProd * double(args(i));
        else
            nonConstants = [nonConstants, args(i)];
        end
    end

    if length(nonConstants) > 1
        curvature = CurvatureType.Unknown;
        return;
    end

    if isempty(nonConstants)
        curvature = CurvatureType.Affine;
        return;
    end

    curv = findCurvature(nonConstants(1), registry);

    % Flip if multiplied by negative
    if constProd < 0
        if curv == CurvatureType.Convex
            curvature = CurvatureType.Concave;
        elseif curv == CurvatureType.Concave
            curvature = CurvatureType.Convex;
        else
            curvature = curv;
        end
    else
        curvature = curv;
    end
end

function curvature = applyCompositionRules(rule, args, registry)
    % Apply DCP composition rules
    fCurv = rule.curvature;
    fMono = rule.monotonicity;

    if fCurv == CurvatureType.Affine
        if isempty(args)
            curvature = CurvatureType.Affine;
        else
            curvature = findCurvature(args(1), registry);
        end
        return;
    end

    if fCurv == CurvatureType.Convex
        for i = 1:length(args)
            argCurv = findCurvature(args(i), registry);
            mono = getMono(fMono, i);

            if argCurv == CurvatureType.Convex && mono ~= MonotonicityType.Increasing
                curvature = CurvatureType.Unknown;
                return;
            end
            if argCurv == CurvatureType.Concave && mono ~= MonotonicityType.Decreasing
                curvature = CurvatureType.Unknown;
                return;
            end
            if argCurv == CurvatureType.Unknown
                curvature = CurvatureType.Unknown;
                return;
            end
        end
        curvature = CurvatureType.Convex;
        return;
    end

    if fCurv == CurvatureType.Concave
        for i = 1:length(args)
            argCurv = findCurvature(args(i), registry);
            mono = getMono(fMono, i);

            if argCurv == CurvatureType.Concave && mono ~= MonotonicityType.Increasing
                curvature = CurvatureType.Unknown;
                return;
            end
            if argCurv == CurvatureType.Convex && mono ~= MonotonicityType.Decreasing
                curvature = CurvatureType.Unknown;
                return;
            end
            if argCurv == CurvatureType.Unknown
                curvature = CurvatureType.Unknown;
                return;
            end
        end
        curvature = CurvatureType.Concave;
        return;
    end

    curvature = CurvatureType.Unknown;
end

function mono = getMono(fMono, idx)
    if iscell(fMono)
        if idx <= length(fMono)
            mono = fMono{idx};
        else
            mono = fMono{end};
        end
    else
        mono = fMono;
    end
end
```

### Step 4: GDCP Analysis for Manifolds

```matlab
% findGCurvature.m
function gcurvature = findGCurvature(expr, manifold, registry)
    % Find the geodesic curvature of a symbolic expression
    %
    % Args:
    %   expr: A symbolic expression (sym)
    %   manifold: Manifold name ('SPD' or 'Lorentz')
    %   registry: DCPAtomRegistry instance
    %
    % Returns:
    %   gcurvature: GCurvatureType enumeration value

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

    % Handle multiplication
    if strcmp(op, 'times') || strcmp(op, 'mtimes')
        gcurvature = handleGMultiplication(args, manifold, registry);
        return;
    end

    % Look up GDCP rule
    rule = registry.getGDCPRule(op);
    if isempty(rule) || ~strcmp(rule.manifold, manifold)
        % Fall back to DCP
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

### Step 5: Main Analysis Function

```matlab
% analyze.m
function result = analyze(expr, manifold)
    % Analyze a symbolic expression for DCP/DGCP compliance
    %
    % Args:
    %   expr: A symbolic expression
    %   manifold: Optional manifold name ('SPD' or 'Lorentz')
    %
    % Returns:
    %   result: struct with curvature, sign, and gcurvature fields

    arguments
        expr sym
        manifold string = ""
    end

    registry = DCPAtomRegistry();

    % Determine curvature
    curvature = findCurvature(expr, registry);

    % Determine sign
    sign = propagateSign(expr, registry);

    % Determine geodesic curvature if manifold specified
    if manifold ~= ""
        gcurvature = findGCurvature(expr, manifold, registry);
    else
        gcurvature = [];
    end

    result = struct('curvature', curvature, ...
                    'sign', sign, ...
                    'gcurvature', gcurvature);
end

function sign = propagateSign(expr, registry)
    % Propagate sign through expression tree
    if isnumeric(expr)
        if expr > 0
            sign = SignType.Positive;
        else
            sign = SignType.Negative;
        end
        return;
    end

    if isempty(symvar(expr))
        val = double(expr);
        if val > 0
            sign = SignType.Positive;
        else
            sign = SignType.Negative;
        end
        return;
    end

    sign = SignType.AnySign;
end
```

### Complete Matlab Example

```matlab
% Example usage
syms x y positive

% Create registry
registry = DCPAtomRegistry();

% Analyze expressions
expr1 = exp(x) + exp(y);
result1 = analyze(expr1);
fprintf('exp(x) + exp(y): %s\n', string(result1.curvature));

expr2 = log(x) + log(y);
result2 = analyze(expr2);
fprintf('log(x) + log(y): %s\n', string(result2.curvature));

% DGCP analysis
syms X [3 3] matrix
result3 = analyze(trace(X), 'SPD');
fprintf('trace(X) on SPD: %s\n', string(result3.gcurvature));
```

---

## DGCP Atoms Reference

### Symmetric Positive Definite (SPD) Manifold

| Atom | Sign | G-Curvature | Monotonicity | Julia Function |
|------|------|-------------|--------------|----------------|
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

| Atom | Sign | G-Curvature | Monotonicity | Julia Function |
|------|------|-------------|--------------|----------------|
| distance(M, p, q) | Positive | GConvex | GAnyMono | `Manifolds.distance` |
| lorentz_log_barrier(p) | Positive | GConvex | GIncreasing | `lorentz_log_barrier` |
| lorentz_homogeneous_quadratic(A, p) | Positive | GConvex | GAnyMono | `lorentz_homogeneous_quadratic` |
| lorentz_homogeneous_diagonal(a, p) | Positive | GConvex | GAnyMono | `lorentz_homogeneous_diagonal` |
| lorentz_least_squares(X, y, p) | Positive | GConvex | AnyMono | `lorentz_least_squares` |

---

## Implementation Checklist

When porting DGCP to a new language:

- [ ] **Define enumerations** for Sign, Curvature, GCurvature, Monotonicity
- [ ] **Create atom registry** as a dictionary mapping functions to properties
- [ ] **Implement expression tree traversal** (bottom-up for curvature propagation)
- [ ] **Handle special cases**: addition (combine curvatures), multiplication (flip if negative constant)
- [ ] **Implement composition rules** for convex/concave functions with monotonicity checks
- [ ] **Add canonization pass** to normalize expressions (e.g., x'Ax -> quad_form)
- [ ] **Register atoms** with their DCP and DGCP properties
- [ ] **Extend for manifolds** by adding manifold-specific GDCP rules
- [ ] **Test with known expressions** from the paper's experiments

## Key Design Decisions

1. **Metadata attachment**: Use language-specific metadata/attribute systems (Python `__dict__`, Matlab `properties`)
2. **Pattern matching**: Use symbolic library's rewriting capabilities or implement manual tree traversal
3. **Extensibility**: Make atom registry a global or singleton that users can extend
4. **Error handling**: Return `Unknown` curvature rather than throwing when rules don't apply
5. **Manifold support**: Keep manifold as a parameter to allow extension to new Riemannian geometries
