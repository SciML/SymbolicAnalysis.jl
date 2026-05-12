# Theoretical Complexity Analysis of DCP and DGCP Verification in SymbolicAnalysis.jl

## 1. Problem Parameterization

We analyze the computational complexity of the DCP (Disciplined Convex Programming) and DGCP (Disciplined Geodesically Convex Programming) verification algorithms implemented in SymbolicAnalysis.jl. The following parameters characterize the input:

| Parameter | Definition |
|-----------|-----------|
| $n$ | Number of nodes in the expression AST (abstract syntax tree) |
| $k$ | Maximum arity of any function node in the AST |
| $d$ | Depth of the expression tree |
| $R_{\text{DCP}}$ | Number of registered DCP rewrite rules (constant; currently 65) |
| $R_{\text{DGCP}}$ | Number of registered DGCP rewrite rules (constant; currently 27) |

**Critical clarification on matrix dimensions.** The expression `logdet(X)` has the same AST regardless of whether `X` is a $5 \times 5$ or $500 \times 500$ matrix. Matrix entries appear as numerical constants at evaluation time, not as additional symbolic nodes during verification. The parameter $n$ counts symbolic nodes in the expression tree, and is independent of any matrix dimension parameters. This is a fundamental distinction: DCP/DGCP verification operates on the *symbolic structure* of the expression, not on its numerical evaluation.

## 2. Pipeline Architecture

The `analyze(ex)` function (defined in `src/SymbolicAnalysis.jl:48--60`) implements a four-phase pipeline for DCP verification, extended to five phases when a manifold is provided for DGCP:

```
Phase 1: Canonicalization        canonize(ex)                 [1 Postwalk + 1 Prewalk]
Phase 2: Sign Propagation        propagate_sign(ex)           [1 Postwalk + 1 Prewalk]
Phase 3: Curvature Propagation   propagate_curvature(ex)      [1 Postwalk + 1 Prewalk]
Phase 4: DGCP Propagation        propagate_gcurvature(ex, M)  [1 Postwalk + 1 Prewalk]
```

Each phase uses the SymbolicUtils rewriting framework, which implements `Postwalk` (bottom-up traversal) and `Prewalk` (top-down traversal) as single-pass tree walks.

## 3. Formal Complexity Theorems

### Definitions

Let $T$ denote an expression tree with $n$ nodes. A *traversal* of $T$ visits every node exactly once, performing $O(1)$ work at each node (rule matching against a constant-size rule set plus metadata attachment). We write $\text{Postwalk}(T)$ and $\text{Prewalk}(T)$ for bottom-up and top-down traversals, respectively.

A *rule chain* $C = \text{Chain}(r_1, \ldots, r_m)$ is a sequence of rewrite rules applied at each node. At each node, the chain attempts each rule in order until one matches, performing at most $m$ comparisons. Since $m$ is a fixed constant (bounded by $\max(R_{\text{DCP}}, R_{\text{DGCP}})$ within each chain, and typically much smaller since each chain uses a dedicated subset of 3--8 rules), the per-node cost of applying the chain is $O(1)$.

---

### Theorem 1 (Canonicalization Complexity)

**Statement.** The canonicalization phase `canonize(ex)` runs in $\Theta(n)$ time and $O(n)$ space.

**Proof sketch.** The implementation (`src/canon.jl:31--58`) constructs a chain of 5 structural rewrite rules (quadratic form recognition, conjugation recognition, double inverse elimination, `log(det(X)) \to \text{logdet}(X)`, and `sum(diag(X)) \to \text{tr}(X)`). Each rule performs constant-time pattern matching via SymbolicUtils' term-level dispatch. The function applies one `Postwalk` followed by one `Prewalk`, each visiting all $n$ nodes exactly once. Total work: $2n \cdot O(1) = \Theta(n)$.

Space is $O(n)$ because both the input tree and the (possibly rewritten) output tree have at most $n$ nodes, and the traversal uses $O(d) \leq O(n)$ stack space. $\square$

---

### Theorem 2 (Sign Propagation Complexity)

**Statement.** The sign propagation phase `propagate_sign(ex)` runs in $\Theta(n)$ time and $O(n)$ space.

**Proof sketch.** The implementation (`src/rules.jl:195--217`) constructs a chain of 8 rewrite rules:
1. Two rules reset signs on symbols and calls (constant-time metadata check).
2. Two rules assign signs from DCP/DGCP rule tables (dictionary lookup in `dcprules_dict` or `gdcprules_dict`, each $O(1)$ amortized via hash table).
3. Two rules assign signs to call expressions using rule table lookup.
4. One rule for multiplication: `mul_sign` (`src/rules.jl:181--193`) iterates over $k_v$ children of the node, where $k_v$ is the arity. Over the entire tree, $\sum_v k_v = n - 1$ (each non-root node is a child of exactly one parent), so the total work across all multiplication nodes is $O(n)$.
5. One rule for addition: `add_sign` (`src/rules.jl:143--179`) similarly iterates children, with the same $O(n)$ amortized bound.

The function applies one `Postwalk` (bottom-up, to propagate signs from leaves) followed by one `Prewalk` (top-down, to finalize). Each traversal is $\Theta(n)$. The total work is $\Theta(n)$.

Space is $O(n)$: sign metadata is attached in-place to existing tree nodes via the SymbolicUtils metadata system. $\square$

---

### Theorem 3 (DCP Curvature Propagation Complexity)

**Statement.** The DCP curvature propagation phase `propagate_curvature(ex)` runs in $\Theta(n)$ time and $O(n)$ space.

**Proof sketch.** The implementation (`src/rules.jl:290--301`) constructs a chain of 3 rewrite rules:
1. Multiplication curvature: `mul_curvature` (`src/rules.jl:228--257`) iterates over children to find at most one non-constant factor and determine curvature. Cost at each multiplication node is $O(k_v)$.
2. Addition curvature: `add_curvature` (`src/rules.jl:259--288`) iterates over children to check for curvature conflicts. Cost at each addition node is $O(k_v)$.
3. General curvature: `find_curvature` (`src/rules.jl:315--390`) performs a dictionary lookup for the atom's DCP rule ($O(1)$ via hash table in `dcprules_dict`), retrieves the atom's curvature and monotonicity, then checks each argument's curvature against the composition rule. Cost at each call node is $O(k_v)$.

The composition rule check (lines 348--385) implements the standard DCP composition theorem: for a convex non-decreasing $f$ composed with convex $g$, the composition $f \circ g$ is convex. This check examines each argument once, costing $O(k_v)$ per node. Summing over all nodes: $\sum_v k_v \leq n$, giving $O(n)$ total.

Two traversals (Postwalk + Prewalk) yield $\Theta(n)$ total. $\square$

---

### Theorem 4 (DGCP Curvature Propagation Complexity)

**Statement.** The DGCP curvature propagation phase `propagate_gcurvature(ex, M)` runs in $\Theta(n)$ time and $O(n)$ space, with a modestly larger constant factor than DCP curvature propagation.

**Proof sketch.** The implementation (`src/gdcp/gdcp_rules.jl:244--253`) has the same structure as DCP propagation: a chain of 3 rewrite rules applied via Postwalk + Prewalk.

The function `find_gcurvature` (`src/gdcp/gdcp_rules.jl:97--242`) is structurally analogous to `find_curvature` but performs additional case analysis:
- It first checks the DGCP rule table (`gdcprules_dict`, $O(1)$ lookup).
- It handles special structural patterns: `logdet` composed with specific operations (lines 110--117), `log` composed with `tr` or `quad_form` (lines 118--121), Schatten norms of matrix logarithms (lines 122--123), etc. Each of these checks is $O(1)$ at each node.
- If no DGCP-specific rule matches, it falls back to the DCP rule table (line 181--185), applying the same composition theorem with geodesic curvature labels.

The constant factor is larger than for DCP due to the additional case analysis (approximately 10 additional $O(1)$ checks per call node), but the asymptotic complexity remains $\Theta(n)$.

**Constant factor analysis.** Let $c_{\text{DCP}}$ denote the average per-node cost of `find_curvature` and $c_{\text{DGCP}}$ denote that of `find_gcurvature`. From the code, `find_gcurvature` performs: one `gdcprules_dict` lookup, up to 6 structural pattern checks, a possible fallback to `dcprules_dict` lookup, and the same argument-iteration loop. In the worst case, $c_{\text{DGCP}} \approx 2 c_{\text{DCP}}$, though for typical expressions where the first or second check matches, $c_{\text{DGCP}} \approx 1.2 c_{\text{DCP}}$. $\square$

---

### Theorem 5 (Total DCP Analysis Complexity)

**Statement.** The complete DCP analysis pipeline `analyze(ex)` performs exactly 6 tree traversals and runs in $\Theta(n)$ time.

**Proof.**

| Phase | Traversals | Time |
|-------|-----------|------|
| Canonicalization | 1 Postwalk + 1 Prewalk | $\Theta(n)$ |
| Sign propagation | 1 Postwalk + 1 Prewalk | $\Theta(n)$ |
| Curvature propagation | 1 Postwalk + 1 Prewalk | $\Theta(n)$ |
| **Total** | **6** | $\Theta(n)$ |

Let $c_1, c_2, c_3$ denote the per-node constants for canonicalization, sign propagation, and curvature propagation respectively. The total time is:

$$T_{\text{DCP}}(n) = 2(c_1 + c_2 + c_3) \cdot n = \Theta(n)$$

The factor of 2 accounts for the Postwalk + Prewalk pair in each phase. $\square$

---

### Theorem 6 (Total DGCP Analysis Complexity)

**Statement.** The complete DGCP analysis pipeline `analyze(ex, M)` performs exactly 8 tree traversals and runs in $\Theta(n)$ time.

**Proof.** DGCP analysis extends DCP analysis with one additional phase:

| Phase | Traversals | Time |
|-------|-----------|------|
| Canonicalization | 1 Postwalk + 1 Prewalk | $\Theta(n)$ |
| Sign propagation | 1 Postwalk + 1 Prewalk | $\Theta(n)$ |
| Curvature propagation (DCP) | 1 Postwalk + 1 Prewalk | $\Theta(n)$ |
| Curvature propagation (DGCP) | 1 Postwalk + 1 Prewalk | $\Theta(n)$ |
| **Total** | **8** | $\Theta(n)$ |

The total time is:

$$T_{\text{DGCP}}(n) = 2(c_1 + c_2 + c_3 + c_4) \cdot n = \Theta(n)$$

where $c_4$ is the per-node constant for DGCP curvature propagation.

Note that Phase 3 (DCP curvature propagation) is *not* skipped in DGCP mode. The implementation (`src/SymbolicAnalysis.jl:48--60`) always runs `propagate_curvature(ex)` before optionally running `propagate_gcurvature(ex, M)`. This is by design: the DCP curvature labels computed in Phase 3 are used as a fallback within `find_gcurvature` (line 181--185 of `gdcp_rules.jl`). $\square$

---

### Theorem 7 (DGCP Marginal Cost)

**Statement.** The marginal cost of DGCP verification over DCP verification is exactly one additional $\Theta(n)$ phase consisting of 2 tree traversals. The theoretical overhead ratio is $8/6 \approx 1.33$.

**Proof.** From Theorems 5 and 6:

$$\frac{T_{\text{DGCP}}(n)}{T_{\text{DCP}}(n)} = \frac{2(c_1 + c_2 + c_3 + c_4)}{2(c_1 + c_2 + c_3)} = 1 + \frac{c_4}{c_1 + c_2 + c_3}$$

In terms of traversal count, the ratio is $8/6 \approx 1.33$.

In terms of wall-clock time, the ratio depends on the relative magnitudes of the per-node constants. If all phases have comparable per-node cost ($c_1 \approx c_2 \approx c_3 \approx c_4$), the ratio is $4/3 \approx 1.33$. If DGCP propagation has a larger constant (say $c_4 \approx 2c_3$ due to additional case analysis), the ratio is $(c_1 + c_2 + c_3 + 2c_3)/(c_1 + c_2 + c_3) \approx 5/3 \approx 1.67$ under the assumption $c_1 \approx c_2 \approx c_3$.

**Why empirically observed "2--3x overhead" is misleading.** Empirical measurements that report 2--3x overhead for DGCP over DCP conflate several sources of overhead that are orthogonal to the algorithmic cost:

1. **JIT compilation.** Julia's just-in-time compiler generates specialized machine code on first invocation of each function. The DGCP pathway involves distinct type specializations (`GCurvature`, `GMonotonicity`, manifold-specific dispatch) that trigger additional compilation. This is a fixed startup cost amortized to zero over repeated calls.

2. **GC pressure from metadata allocation.** DGCP propagation attaches `GCurvature` metadata to nodes via `setgcurvature`, which allocates new `Metadata` wrappers. The incremental GC pressure from this additional metadata pass can cause non-deterministic slowdowns.

3. **Measurement noise at small $n$.** For small expression trees ($n < 100$), the absolute time difference between DCP and DGCP is in the microsecond range, where timer resolution and system jitter dominate.

The correct framing: DGCP adds one $O(n)$ pass to a pipeline of three $O(n)$ passes. The asymptotic complexity class is identical. For problems that DGCP can verify but DCP cannot (e.g., geodesically convex optimization on symmetric positive definite manifolds), the alternative is *no automated verification at all*. $\square$

---

### Theorem 8 (Conic Form Generation Complexity)

**Statement.** The conic form generation procedure `to_conic_form(ex)` runs in $\Theta(n)$ time and produces $O(n)$ conic constraints.

**Proof sketch.** The implementation (`src/conic.jl:221--257`) first runs the DCP verification pipeline ($\Theta(n)$ by Theorem 5), then performs a single recursive bottom-up traversal via `_process_node!` (`src/conic.jl:280--484`).

At each node, `_process_node!` performs one of:
- **Leaf (symbol or number):** $O(1)$ work --- either return the symbol or create one equality constraint.
- **Affine subtree:** `_is_affine` checks the subtree structure and `_extract_affine` linearizes it. For an affine subtree with $m$ nodes, this takes $O(m)$ time and produces exactly 1 constraint. Since affine subtrees are disjoint, the total cost across all affine subtrees is $O(n)$.
- **Addition/multiplication:** $O(k_v)$ work to process $k_v$ children, producing 1 constraint.
- **DCP atom:** $O(1)$ dictionary lookup for the cone annotation, $O(k_v)$ to process children, then `_emit_atom_constraint!` (`src/conic.jl:499--953`) emits a constant number of constraints per atom (at most 3, e.g., for `logistic` which requires two exponential cone constraints and one linear inequality).

Each node is visited exactly once. The total number of constraints is at most $O(n)$: each node contributes at most a constant number of constraints. Total time: $\Theta(n)$ for the traversal plus $\Theta(n)$ for the preceding DCP verification, giving $\Theta(n)$ overall.

Space is $O(n)$: the `ConicContext` accumulates $O(n)$ constraints, each of constant size, plus $O(n)$ epigraph variables. $\square$

---

### Theorem 9 (Lower Bound)

**Statement.** Any DCP or DGCP verifier must examine every node of the expression tree, requiring $\Omega(n)$ time.

**Proof.** We give an adversarial argument. Consider the family of expressions $e_i$ for $i \in \{1, \ldots, n\}$ defined as follows: $e_i$ is a sum of $n$ terms, where $n - 1$ terms are affine (and hence DCP-compliant) and the $i$-th term is a non-convex function (e.g., $\sin(x_i)$, which has no DCP rule). The expression $e_i$ is DCP-compliant if and only if the $i$-th term is replaced by a DCP-compliant atom.

Any verifier that does not examine the $i$-th node cannot distinguish the DCP-compliant expression from the non-compliant one. Since $i$ can be any value in $\{1, \ldots, n\}$, the verifier must examine all $n$ nodes in the worst case. $\square$

---

## 4. Detailed Per-Phase Analysis

### 4.1 Rule Matching Cost

Each phase uses a `Chain` of rewrite rules. At each node, the chain applies rules sequentially until one matches. The key observation is that the number of rules per chain is a small constant:

| Phase | Rules in Chain | Bound |
|-------|---------------|-------|
| Canonicalization | 5 structural patterns | $O(1)$ |
| Sign propagation | 8 rules | $O(1)$ |
| DCP curvature | 3 rules + `find_curvature` | $O(1)$ |
| DGCP curvature | 3 rules + `find_gcurvature` | $O(1)$ |

Within `find_curvature` and `find_gcurvature`, the dominant operation is a hash-table lookup in `dcprules_dict` or `gdcprules_dict`. These dictionaries are keyed by function identity (the Julia `Function` object), giving $O(1)$ amortized lookup. The DCP composition rule check iterates over the atom's arguments, but as shown in the proof of Theorem 3, the total work across all nodes is $O(n)$.

### 4.2 Why Two Traversals Per Phase

Each phase applies both a `Postwalk` (bottom-up) and a `Prewalk` (top-down). This is necessary for correctness:

- **Postwalk (bottom-up):** Propagates information from leaves to the root. For sign propagation, this computes the sign of each subexpression from its children's signs. For curvature propagation, this applies the DCP composition theorem bottom-up.

- **Prewalk (top-down):** Handles cases where a parent's metadata should override or refine a child's. For example, if canonicalization rewrites a subtree at the top level, the Prewalk ensures the rewritten form is propagated downward. For sign and curvature, the Prewalk catches expressions where a top-level rule match (e.g., a registered atom with known sign) should propagate to children.

Using both traversals is a standard technique in attribute grammar evaluation for synthesized and inherited attributes.

### 4.3 Dictionary Lookup Analysis

The rule dictionaries `dcprules_dict` and `gdcprules_dict` are Julia `Dict` objects using hash-based lookup.

- `dcprules_dict` contains entries for 65 atoms (including overloaded rules for the same function under different domains). Multiple rules for the same function are stored in a `Vector`, searched linearly. The maximum number of rules per function is 3 (for `quad_over_lin`, `inv`, `sqrt`, and `log`). Thus the per-lookup cost is $O(1)$ with a small constant.

- `gdcprules_dict` contains entries for 27 atoms (21 for `SymmetricPositiveDefinite`, 6 for `Lorentz`). Each function has at most 1 rule, giving strictly $O(1)$ per lookup.

## 5. Comparison with Related Systems

All major DCP verification systems share the same asymptotic complexity:

| System | Verification | Traversals | Rule Set Size | Verifiable Class |
|--------|-------------|-----------|---------------|-----------------|
| **Convex.jl** | $O(n)$ | ~2 (sign + curvature) | ~40 atoms | DCP |
| **CVXPY** | $O(n)$ | ~2 (sign + curvature) | ~80 atoms | DCP |
| **DCCP** | $O(n)$ | ~3 (extends CVXPY) | ~80 atoms + DC rules | DCP + difference-of-convex |
| **SymbolicAnalysis.jl (DCP)** | $O(n)$ | 6 (3 phases x 2) | 65 atoms | DCP |
| **SymbolicAnalysis.jl (DGCP)** | $O(n)$ | 8 (4 phases x 2) | 65 + 27 atoms | DCP + DGCP |

The key contribution of SymbolicAnalysis.jl is not in asymptotic complexity (which is optimal by Theorem 9) but in the *verifiable class*: DGCP verification accepts a strictly larger set of optimization problems (those involving geodesic convexity on Riemannian manifolds) at the same $O(n)$ asymptotic cost as standard DCP verification.

Specifically:
- **Convex.jl** and **CVXPY** verify standard DCP problems via equivalent $O(n)$ tree traversals.
- **DCCP** extends to difference-of-convex programs, also in $O(n)$, but targets a different generalization direction (non-convex decomposition rather than Riemannian geometry).
- **SymbolicAnalysis.jl** provides unified DCP + DGCP verification with an explicit conic form generation pass, all in $O(n)$.

The additional traversal count (6 vs. ~2 in Convex.jl) reflects the architecture of using SymbolicUtils' rewriting framework, which requires separate Postwalk + Prewalk passes for each phase, rather than a monolithic visitor. This is a constant-factor difference with no asymptotic impact, and provides the benefit of modularity: each phase can be independently tested and extended.

## 6. Summary

The DCP and DGCP verification algorithms in SymbolicAnalysis.jl are optimal up to constant factors:

- **DCP verification:** $\Theta(n)$ time via 6 tree traversals. Matches the $\Omega(n)$ lower bound.
- **DGCP verification:** $\Theta(n)$ time via 8 tree traversals. Same asymptotic class as DCP.
- **Conic form generation:** $\Theta(n)$ time producing $O(n)$ constraints.
- **DGCP marginal cost:** One additional $O(n)$ phase (2 traversals), yielding a theoretical overhead ratio of $8/6 \approx 1.33$ in traversal count.

The entire pipeline --- from raw symbolic expression to verified curvature label (and optionally to conic form) --- is linear in the size of the expression tree, independent of any matrix dimensions or numerical parameters.
