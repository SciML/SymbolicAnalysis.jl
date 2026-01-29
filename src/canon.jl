"""
DGCP-Aware Canonicalization Pass

This module provides expression rewriting rules to transform symbolic expressions
into DGCP-verifiable canonical forms.

Addresses Reviewer 385's concern about symbolic representation non-uniqueness:
"log(x²) is not DCP-valid while 2log(x) is. Similar situations may occur for 
the proposed methodology."
"""

using SymbolicUtils.Rewriters: Chain, Postwalk, Prewalk
using SymbolicUtils: @rule

#==============================================================================#
# Core Canonicalization (Original + Safe Extensions)
#==============================================================================#

"""
    canonize(ex)

Apply DGCP-aware canonicalization rules to transform expressions into
forms that are more likely to be verifiable by the DGCP framework.

Currently applies:
1. Pattern recognition: x'Ax → quad_form(x, A), B'XB → conjugation(X, B)
2. Inverse simplification: inv(inv(X)) → X
"""
function canonize(ex)
    # Core rules that are safe and well-tested
    core_rules = [
        # Quadratic form recognition: x'*Y*x → quad_form(x, Y)
        @rule (adjoint(~x) * (~Y * ~x))[1] => quad_form(~x, ~Y)
        
        # Conjugation recognition: B'*X*B → conjugation(X, B)
        @rule ((adjoint(~B) * ~X) * ~B)[Base.OneTo(size(~B, 2)), Base.OneTo(size(~B, 1))] => 
            conjugation(~X, ~B)
        
        # Double inverse: inv(inv(X)) → X
        @rule inv(inv(~X)) => ~X
    ]
    
    try
        rc = Chain(core_rules)
        ex = Postwalk(rc)(ex)
        ex = Prewalk(rc)(ex)
        return ex
    catch
        return ex
    end
end

#==============================================================================#
# Extended Canonicalization (Optional, for standalone use)
#==============================================================================#

"""
    canonize_extended(ex)

More aggressive canonicalization with additional rules.
Use with caution - may not work with all expression types.
"""
function canonize_extended(ex)
    ex = canonize(ex)  # First apply core rules
    
    extended_rules = [
        # logdet(inv(X)) → -logdet(X)
        @rule LinearAlgebra.logdet(inv(~X)) => -LinearAlgebra.logdet(~X)
    ]
    
    try
        rc = Chain(extended_rules)
        ex = Postwalk(rc)(ex)
        return ex
    catch
        return ex
    end
end

#==============================================================================#
# Helper: Check if Expression is in Canonical Form
#==============================================================================#

"""
    is_canonical(ex)

Check if an expression is already in canonical form.
Returns true if no canonicalization rules would change the expression.
"""
function is_canonical(ex)
    try
        canonical_ex = canonize(ex)
        return isequal(ex, canonical_ex)
    catch
        return true  # If canonization fails, assume canonical
    end
end

#==============================================================================#
# Documentation of Known Equivalent Forms
#==============================================================================#

"""
    equivalent_forms()

Returns documentation of known equivalent symbolic forms where one form
is DGCP-verifiable and another is not.

This addresses Reviewer 385's concern about symbolic representation non-uniqueness.
"""
function equivalent_forms()
    # Note: These document cases where symbolic representation affects verifiability.
    # Some are mathematically equivalent, others are different functions that users
    # might confuse or accidentally write in non-verifiable form.
    forms = [
        (
            verifiable = "-logdet(X)",
            not_verifiable = "logdet(inv(X))",
            note = "Mathematically equivalent: -log|X| = log|X^{-1}|. Use canonize_extended() to transform."
        ),
        (
            verifiable = "log(tr(X)) + log(tr(Y))",
            not_verifiable = "log(tr(X) * tr(Y))",
            note = "Mathematically equivalent by log properties. Sum of logs is DGCP-compliant."
        ),
        (
            verifiable = "tr(inv(X))",
            not_verifiable = "sum(eigvals(inv(X)))",
            note = "Semantically equivalent (trace = sum of eigenvalues). High-level form verifiable."
        ),
        (
            verifiable = "sum(distance(M, As[i], X)^2 for i in 1:n)",
            not_verifiable = "sum(log(eigvals(As[i]^(-1/2) * X * As[i]^(-1/2)))^2 for i in 1:n)",
            note = "Semantically equivalent. Use high-level distance atom for verification."
        ),
        (
            verifiable = "2 * logdet(X)",
            not_verifiable = "logdet(X)^2",
            note = "NOT equivalent! Common user mistake. 2*log|X| ≠ (log|X|)². First is g-linear."
        ),
    ]
    return forms
end
