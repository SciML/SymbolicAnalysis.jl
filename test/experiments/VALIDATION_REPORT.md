# DGCP Experiments Validation Report

This report validates whether the experiments in `test/experiments/` properly address the reviewer comments from the DGCP paper revisions.

---

## Experiment 1: non_gconvex_examples.jl

### Reviewer Comment Addressed
**Technical Review #2 (Reviewer 400):** "It would be valuable if the authors could provide explicit examples, using key atoms, to illustrate how the framework recognizes functions that are NOT geodesically convex."

### What the Experiment Does
- Tests 6 expressions that should NOT be verified as g-convex:
  1. `sqrt(X * Y)` - Product of two SPD variables
  2. `X - A` - Matrix subtraction
  3. `tr(X^2)` - Quadratic trace without log transform
  4. `X + Y` - Sum of two matrix variables
  5. `logdet(X)^2` - Square of logdet (not same as 2*logdet)
  6. `logdet(X) * logdet(Y)` - Product of g-linear terms

- Verifies each returns `GUnknownCurvature`
- Includes explanations for WHY each cannot be verified
- Bonus: Demonstrates symbolic non-uniqueness (2*logdet vs logdet^2)

### Validation: PASS
The experiment properly addresses the reviewer's request by:
- Providing explicit examples of non-g-convex/non-DGCP-verifiable functions
- Showing that DGCP correctly returns `GUnknownCurvature` for these cases
- Explaining the mathematical reasoning behind each rejection
- Also addresses Reviewer 385's concern about symbolic non-uniqueness

---

## Experiment 2: dcp_dgcp_comparison.jl

### Reviewer Comment Addressed
**Technical Review #2 (Reviewer 399):** "Under the premise of fair comparison, is there an existing DCP software package that can be directly compared with DGCP? [...] compare their capabilities in performing symbolic analysis and convexity verification"

**Technical Review #2 (Reviewer 400):** "explicitly demonstrate correspondence between DGCP and classical DCP under the assumption of a Euclidean manifold"

### What the Experiment Does
- Compares DGCP results with Euclidean convexity for 7 functions:
  1. `logdet(X)` - Both DCP and DGCP verify
  2. `tr(inv(X))` - G-convex on SPD
  3. `distance(M, A, X)^2` - DGCP only (Riemannian distance)
  4. `S-divergence(X, A)` - DGCP only
  5. `logdet(A' X^{-1} A)` - DGCP only (conjugation)
  6. Tyler's M-Estimator - DGCP only
  7. Karcher Mean - DGCP only

- Reports both DGCP curvature and Euclidean convexity status
- Optionally integrates with Convex.jl for DCP comparison

### Validation: PARTIAL PASS
**Strengths:**
- Shows verification scope difference (what DGCP can verify that DCP cannot)
- Compares Euclidean vs geodesic convexity for each function
- Includes both "both verify" and "DGCP only" examples

**Gaps:**
- Performance comparison is NOT included (reviewer asked about timing comparison for functions both can verify)
- Could strengthen by adding explicit timing benchmarks for `logdet(X)` in both DCP and DGCP

---

## Experiment 3: extended_benchmark.jl

### Reviewer Comment Addressed
**Technical Review #2 (Reviewer 399):** "The experiments in this paper concerning symbolic complexity and verification time remain insufficient. Could the authors design one or more experiments to explore symbolic complexity and verification time in greater depth?"

### What the Experiment Does
- Defines AST complexity metrics:
  - `count_ast_nodes(ex)` - Total nodes in expression tree
  - `ast_depth(ex)` - Maximum depth of expression tree
  - `count_unique_operations(ex)` - Number of unique operations

- Benchmarks 4 problem types across multiple matrix sizes:
  - Tyler's M-Estimator (5-30)
  - Karcher Mean (25-150)
  - Log-Determinant (50-400)
  - Brascamp-Lieb (5-30)

- Records: median time, AST nodes, AST depth, memory allocation
- Creates correlation plots (time vs complexity, time vs size)
- Computes approximate scaling exponents via log-log regression

### Validation: PASS
The experiment fully addresses the reviewer's request by:
- Measuring symbolic complexity (AST nodes, depth)
- Correlating complexity with verification time
- Providing scaling analysis
- Generating visualizations of the relationship

---

## Experiment 4: convergence_comparison.jl

### Reviewer Comment Addressed
**Technical Editor Comment #2:** "it would strengthen the paper to demonstrate the benefits of DGCP by solving the problems as nonconvex using state-of-the-art local nonlinear optimization solvers and also with a Riemannian solver, allowing a comparison that highlights the advantage of certified g-convexity"

### What the Experiment Does
- Compares 3 optimization approaches on Karcher mean problem:
  1. **Euclidean BFGS** (via Optim.jl) - Treats as unconstrained optimization
  2. **Riemannian Gradient Descent** (via Manopt.jl) - Manifold-aware
  3. **Riemannian Conjugate Gradient** (via Manopt.jl) - Faster manifold-aware

- Tests on multiple problem sizes (n=5,10,15; m=10,20,30 data points)
- Tracks:
  - Final objective value
  - Whether result stays on SPD manifold (is_spd check)
  - Computation time
  - Success/failure status

### Validation: PASS
The experiment properly addresses the reviewer's request by:
- Using state-of-the-art Euclidean solver (BFGS)
- Using Riemannian solvers (GD, CG) via Manopt.jl
- Comparing convergence and manifold-feasibility
- Demonstrating that Euclidean solvers may leave the SPD manifold while Riemannian solvers stay on it

---

## Experiment 5: expert_examples.jl

### Reviewer Comment Addressed
**Technical Review #2 (Reviewer 400):** "Can the proposed DGCP framework correctly identify complex cases that challenge even human experts?"

### What the Experiment Does
- Documents 6 complex verification cases with:
  - Mathematical formula
  - Literature reference
  - Estimated difficulty for human experts (Easy/Medium/Hard)
  - DGCP verification result
  - Verification time

- Cases included:
  1. **Tyler's M-Estimator** - Tyler (1987) - Hard
  2. **Brascamp-Lieb Bound** - Sra & Hosseini (2015) - Hard
  3. **Matrix Square Root via S-Divergence** - Sra (2016) - Medium
  4. **Karcher Mean** - Karcher (1977) - Hard
  5. **Diagonal Loading Regularization** - Ledoit & Wolf (2004) - Medium
  6. **Sum of Largest Log-Eigenvalues** - Lewis (1996) - Hard

- For each case, explains what expert verification would require

### Validation: PASS
The experiment properly addresses the reviewer's request by:
- Including genuinely complex cases from the literature
- Providing proper references for each case
- Explaining WHY each case is challenging for human experts
- Demonstrating that DGCP verifies these in milliseconds
- Including 4 "Hard" cases and 2 "Medium" cases

---

## Summary Table

| Experiment | Reviewer Comment | Status | Notes |
|------------|------------------|--------|-------|
| non_gconvex_examples.jl | Non-g-convex identification | **PASS** | Fully addresses with 6 examples + explanations |
| dcp_dgcp_comparison.jl | Fair DCP vs DGCP comparison | **PARTIAL** | Good scope comparison, missing performance comparison |
| extended_benchmark.jl | Symbolic complexity + timing | **PASS** | Full AST metrics + correlation analysis |
| convergence_comparison.jl | Euclidean vs Riemannian solvers | **PASS** | Compares BFGS vs Manopt solvers |
| expert_examples.jl | Complex expert-level cases | **PASS** | 6 cases with proper references |

---

## Recommendations

### For dcp_dgcp_comparison.jl
Add a performance comparison section that:
1. Times DGCP verification of `logdet(X)`
2. Times DCP (Convex.jl) verification of the same expression
3. Reports timing comparison to show DGCP doesn't add significant overhead

### General
- All experiments include proper test sets for automated validation
- References are provided where applicable
- Output formatting is clear and informative

---

## Conclusion

**4 out of 5 experiments fully address their corresponding reviewer comments.** The `dcp_dgcp_comparison.jl` experiment partially addresses the reviewer's request but could be strengthened with explicit performance timing comparisons for functions that both DCP and DGCP can verify.

Overall, the experiments provide strong evidence addressing the major reviewer concerns about:
1. Demonstrating non-g-convex identification
2. Comparing verification scope between DCP and DGCP
3. Analyzing symbolic complexity and its relationship to verification time
4. Showing practical optimization benefits of DGCP-verified problems
5. Handling complex cases that challenge human experts
