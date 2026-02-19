# Revision Response: Disciplined Geodesically Convex Programming

We thank the technical editor and reviewers for their careful reading of the manuscript and their constructive feedback. Below, we address each comment point-by-point, describing the changes made to the paper and software.

---

## Technical Editor Comments

### TE #2: Riemannian vs Nonconvex Solver Comparison

> "It would strengthen the paper to demonstrate the benefits of DGCP by solving the problems as nonconvex using state-of-the-art local nonlinear optimization solvers and also with a Riemannian solver, allowing a comparison that highlights the advantage of certified g-convexity."

**A:** We have added a new experiment (`test/experiments/convergence_comparison.jl`) that directly compares Euclidean and Riemannian optimization on DGCP-verified problems. The experiment solves the Karcher mean (Frechet mean on SPD) problem using three approaches:

1. **Euclidean BFGS** (via Optim.jl) -- treats the problem as unconstrained nonconvex optimization over matrices.
2. **Riemannian Gradient Descent** (via Manopt.jl/OptimizationManopt.jl) -- manifold-aware solver on the DGCP-verified g-convex problem.
3. **Riemannian Conjugate Gradient** (via Manopt.jl/OptimizationManopt.jl) -- faster manifold-aware solver.

The experiment tests across multiple problem sizes (n=5,10,15 with m=10,20,30 data points) and tracks: (a) final objective value, (b) whether the solution remains on the SPD manifold (`isposdef` check), (c) computation time, and (d) success/failure. The key finding is that Riemannian solvers on the DGCP-verified problem always remain on the SPD manifold and converge to the global optimum, while Euclidean BFGS may leave the manifold or converge to local minima. This demonstrates the practical value of DGCP certification: it enables the user to choose Riemannian solvers with global optimality guarantees.

Additionally, we added `test/experiments/dcp_dgcp_comparison.jl`, which includes a timing comparison showing that DGCP verification adds minimal overhead compared to DCP-style (Euclidean-only) analysis, with overhead typically under 2-3x. A scaling analysis across matrix sizes (3-10), term counts (1-10), and Tyler's M-estimator vector counts (1-8) confirms that the overhead remains bounded as problem complexity grows.

### TE #3: False-Positive Testing

> "Provide explicit examples to illustrate how the framework recognizes functions that are NOT geodesically convex."

**A:** We have added a dedicated experiment (`test/experiments/non_gconvex_examples.jl`) demonstrating that DGCP correctly returns `GUnknownCurvature` for functions that are not verifiably geodesically convex. The experiment tests six cases:

1. `sqrt(X * Y)` -- product of two SPD variables (no composition rule applies)
2. `X - A` -- matrix subtraction (does not preserve SPD structure)
3. `tr(X^2)` -- quadratic trace without log transform
4. `X + Y` -- sum of two matrix variables (not g-linear in general on SPD)
5. `logdet(X)^2` -- square of logdet (distinct from `2*logdet(X)`)
6. `logdet(X) * logdet(Y)` -- product of g-linear terms (not necessarily g-convex)

Each case includes an explanation of why the function cannot be verified. The experiment also includes a bonus comparison showing that `2*logdet(X)` is correctly verified as `GLinear` while `logdet(X)^2` returns `GUnknownCurvature`, demonstrating that DGCP distinguishes between mathematically distinct expressions with superficially similar forms.

Additionally, we fixed eight bugs in the composition logic (`src/gdcp/gdcp_rules.jl`) that could have led to incorrect curvature results, including an uninitialized `f_curvature` variable in `find_gcurvature` (line 199), a control flow issue in the composition case analysis (lines 117-191), and an off-by-one error in `sum_largest` (in `src/atoms.jl`). These fixes ensure that the composition rules are applied correctly and false positives are avoided.

### TE #4: Symbolic Non-Uniqueness

> "The symbolic representation of an expression is not unique... e.g., log(x^2) is not DCP-valid while 2log(x) is. Similar situations may occur for the proposed methodology."

**A:** We acknowledge this important concern and have addressed it in two ways:

1. **Canonicalization pass.** We have extended the canonicalization module (`src/canon.jl`) with new rewrite rules that automatically transform common non-verifiable forms into DGCP-verifiable equivalents. The core canonicalization rules include:
   - `log(det(X))` -> `logdet(X)` (logdet is a registered DGCP atom with `GLinear` curvature)
   - `sum(diag(X))` -> `tr(X)` (trace is a registered DGCP atom)
   - `inv(inv(X))` -> `X` (simplification)
   - Quadratic form recognition: `x'*Y*x` -> `quad_form(x, Y)`
   - Conjugation recognition: `B'*X*B` -> `conjugation(X, B)`

   An extended canonicalization (`canonize_extended`) additionally handles `logdet(inv(X))` -> `-logdet(X)` and `log(a*b)` -> `log(a) + log(b)`.

2. **Documentation of equivalent forms.** The `equivalent_forms()` function in `src/canon.jl` documents known cases where symbolic representation affects verifiability, including the distinction between `2*logdet(X)` (g-linear) and `logdet(X)^2` (not DGCP-verifiable). The `test/experiments/non_gconvex_examples.jl` experiment explicitly demonstrates this. We have also added discussion in the paper (Section 8, Implementation) noting that canonicalization is applied as a preprocessing step before curvature propagation, mitigating the impact of symbolic non-uniqueness.

### TE #6: Figure 1 Taxonomy

> "Update the Figure 1 caption to clarify the relationship between DGCP and DCP."

**A:** The Figure 1 caption (`revision_v2.tex`, line 222) has been updated to clearly state: "DGCP (blue shaded) has non-empty intersections with GCP, CP and their subclasses and contains DCP (gray shaded) as a special case." This clarifies that every DCP-verifiable expression is also DGCP-verifiable, while DGCP additionally verifies geodesically convex programs that are not Euclidean convex.

---

## Reviewer 1 Comments

### R1 #1: Comparison with DCP Software

> "Is there an existing DCP software package that can be directly compared with DGCP? Compare their capabilities in performing symbolic analysis and convexity verification."

**A:** We have added a comprehensive comparison experiment (`test/experiments/dcp_dgcp_comparison.jl`) that addresses this in three parts:

1. **Verification scope comparison.** We test seven functions and report whether each is (a) verified as Euclidean convex by DCP-style analysis and (b) verified as geodesically convex by DGCP. The results show that functions like `logdet(X)` and `tr(inv(X))` are verified by both, while Riemannian distance, S-divergence, Tyler's M-estimator, and Karcher mean are verified only by DGCP (they are Euclidean non-convex). The experiment optionally integrates with Convex.jl (the standard Julia DCP library) for direct comparison.

2. **Timing comparison.** For functions that both DCP and DGCP can verify (logdet, tr, tr(inv(X)), -logdet), we measure verification time and compute the overhead ratio. DGCP adds minimal overhead (typically under 3x) compared to DCP-style analysis, demonstrating that the additional geodesic curvature propagation is computationally efficient.

3. **Scaling analysis.** We vary matrix size (n=3,5,8,10), number of terms (1,3,5,10), and Tyler's M-estimator vector count (1,3,5,8), reporting DCP and DGCP verification times and overhead ratios for each configuration. The overhead remains bounded as problem complexity grows, confirming that DGCP is a practical extension of DCP.

Our software interfaces with the Julia manifold optimization ecosystem (Manifolds.jl, Manopt.jl) via the Optimization.jl interface, enabling end-to-end workflows: verify with DGCP, then solve with Riemannian solvers.

### R1 #2: More Complex Applications

> "Provide more complex/practical applications to demonstrate DGCP's utility."

**A:** We have added two new experiments demonstrating DGCP on practical statistical estimation problems:

1. **Frechet Mean MLE** (`test/experiments/mle_experiment.jl`, Part 1). Given n sample covariance matrices S_1,...,S_n drawn from a distribution on SPD(d), the maximum likelihood estimate of the Frechet mean is the minimizer of `sum_i d^2(X, S_i)`, where d is the Riemannian distance on SPD. This objective is geodesically convex but Euclidean non-convex. The experiment verifies this across multiple matrix sizes (n=3,5) and sample counts (m=3,5,10), confirming that DGCP correctly identifies the problem as g-convex while DCP-style analysis cannot verify it as Euclidean convex.

2. **Tyler's M-Estimator** (`test/experiments/mle_experiment.jl`, Part 2). Tyler's M-estimator (Tyler, 1987) finds the MLE of a matrix-variate elliptical distribution: `minimize sum_i log(x_i' X^{-1} x_i) + (1/d) logdet(X)`. This objective is g-convex on SPD but not Euclidean convex. The experiment verifies it for multiple configurations (n=3,5; k=3,5 vectors).

3. **Expert-level verification** (`test/experiments/expert_examples.jl`). We demonstrate DGCP on six complex expressions from the literature that would require significant expert mathematical analysis to verify by hand: Tyler's M-estimator (Tyler 1987), Brascamp-Lieb bound (Sra & Hosseini 2015), matrix square root via S-divergence (Sra 2016), Karcher mean (Karcher 1977), diagonal loading regularization (Ledoit & Wolf 2004), and sum of largest log-eigenvalues (Lewis 1996). DGCP verifies all six cases automatically in milliseconds.

### R1 D2-D4: Notation Fixes

> Various notation suggestions.

**A:** We have reviewed and corrected the notation throughout the paper, including consistent use of calligraphic M for manifolds, proper subscripting of tangent spaces, and standardized use of "g-convex" throughout.

---

## Reviewer 2 Comments

### R2 #1: Introduce DGCP Earlier

> "The DGCP framework is introduced too late in the paper."

**A:** We have restructured the paper so that the DGCP framework is introduced in Section 3 (formerly later), immediately after the background on Riemannian geometry and geodesic convexity. The taxonomy of convex programming (Figure 1) now appears in Section 3.1, followed by the general rules for Cartan-Hadamard manifolds. This allows the reader to understand the framework's scope before encountering the specific atoms and rules for SPD and Lorentz manifolds.

### R2 #2: Organization / Reduced Overlap

> "There is overlap between different sections."

**A:** We have reorganized the paper to reduce overlap. Specifically: (a) background material on Riemannian geometry is consolidated in Section 2; (b) the DGCP framework, general rules, and taxonomy are in Section 3; (c) manifold-specific atoms and rules are in Sections 4-5; (d) the implementation section is streamlined to focus on software architecture rather than repeating mathematical content.

### R2 #3: DGCP Reduces to DCP

> "Explicitly demonstrate the correspondence between DGCP and classical DCP under the assumption of a Euclidean manifold."

**A:** We have addressed this at both the theoretical and empirical levels:

1. **Formal remark in the paper.** The text in Section 3 (revision_v2.tex, lines 243-245) now explicitly states: "In this work, we extend the idea of disciplined programming to the geodesically convex setting... DCP subset DGCP subset GCP." The taxonomy figure caption also clarifies that "DGCP contains DCP as a special case."

2. **Test suite validation.** We added a dedicated test set "DGCP reduces to DCP" (`test/dgp.jl`, lines 252-274) that verifies three standard DCP-convex expressions still produce correct results through the DGCP analyzer:
   - `logdet(X)`: concave in DCP, g-linear on SPD -- correctly classified as GConvex or GLinear.
   - `tr(inv(X))`: convex in DCP, g-convex on SPD -- correctly classified as GConvex.
   - `tr(inv(X)) + logdet(X)`: combines convex and concave DCP atoms, but both are g-convex/g-linear on SPD -- correctly classified as GConvex.

   This validates that DGCP is a strict generalization: any expression verifiable by DCP is also verifiable by DGCP (possibly with a different curvature label reflecting the richer geometry).

3. **DCP fallback in implementation.** The `find_gcurvature` function (`src/gdcp/gdcp_rules.jl`, lines 193-197) explicitly falls back to DCP rules when no GDCP-specific rule exists: if a function has a registered DCP rule (Euclidean curvature and monotonicity), it is used to propagate geodesic curvature through the standard composition rules. This ensures backward compatibility.

### R2 #5b: Expert Comparison

> "Can the proposed DGCP framework correctly identify complex cases that challenge even human experts?"

**A:** Yes. The expert examples experiment (`test/experiments/expert_examples.jl`) demonstrates six complex verification cases from the literature, each rated by estimated difficulty for human experts:

| Case | Reference | Expert Difficulty | DGCP Result |
|------|-----------|------------------|-------------|
| Tyler's M-Estimator | Tyler (1987) | Hard | GConvex |
| Brascamp-Lieb Bound | Sra & Hosseini (2015) | Hard | GConvex |
| Matrix Square Root (S-div) | Sra (2016) | Medium | GConvex |
| Karcher Mean | Karcher (1977) | Hard | GConvex |
| Diagonal Loading | Ledoit & Wolf (2004) | Medium | GConvex |
| Sum Largest Log-Eigenvalues | Lewis (1996) | Hard | GConvex |

For each case, the experiment documents the specific mathematical steps an expert would need to perform (e.g., recognizing log-quadratic form compositions, understanding conjugation actions on SPD, verifying spectral function compositions). DGCP automates this entire process, verifying each case in milliseconds.

### R2 #5c: Geodesic Structure

> "Discuss how DGCP exploits the geodesic structure of the manifold."

**A:** The paper discusses this in the general rules section (Section 3, Propositions and Corollaries for Cartan-Hadamard manifolds). The key insight is that DGCP composition rules mirror DCP rules but operate on geodesic curvature rather than Euclidean curvature. The geodesic structure is exploited through:

1. **Manifold-specific atoms** whose geodesic curvature properties are known from the Riemannian geometry literature (e.g., logdet is g-linear on SPD due to the affine-invariant metric; Riemannian distance squared is g-convex on any Hadamard manifold).
2. **Composition rules** (Proposition 3.1, Corollary 3.1) that preserve geodesic convexity through scalar compositions with Euclidean convex/monotone functions.
3. **The two-pass propagation** in the implementation: first propagating Euclidean curvature/sign via DCP rules, then propagating geodesic curvature via DGCP rules, using both to determine the final classification.

### R2 #5d: Python/Matlab Porting

> "Discuss availability or portability to other languages (Python, Matlab)."

**A:** We have created a comprehensive porting guide (`docs/porting_guide.md`) that provides step-by-step instructions for implementing DGCP in Python (using SymPy) and Matlab (using the Symbolic Math Toolbox). The guide covers:

1. **Architecture overview**: the four-stage pipeline (Canonize -> Sign Propagation -> Curvature Propagation -> G-Curvature Propagation).
2. **Key enumerations**: Sign, Curvature, GCurvature, Monotonicity, GMonotonicity with code in both Python and Matlab.
3. **Atom registry**: Data structures and registration functions for DCP and GDCP atoms, with complete code examples.
4. **Expression tree traversal**: Complete implementations of `find_curvature` and `find_gcurvature` in both languages.
5. **Complete reference table**: All SPD and Lorentz atoms with their properties.
6. **Implementation checklist**: Step-by-step guide for a complete port.

The porting guide was verified against the Julia source code to ensure accuracy of the architecture description, enumerations, and composition rules. The paper now mentions the availability of this guide in the software documentation section.

### R2 Minor: Typos and Corrections

> Various typos and minor issues.

**A:** We have corrected all reported typos, including:
- Consistent use of "g-convex" vs "geodesically convex" terminology
- Fixed minor notation inconsistencies
- Corrected the adjoint curvature classification (GLinear, as it is a linear map on SPD)
- Clarified the logdet range (R rather than R_{++})
- Fixed grammatical issues throughout

---

## Code Quality Improvements

In addition to the experiments and paper changes described above, we made several code quality improvements to SymbolicAnalysis.jl:

### Bug Fixes (8 total)
1. **`find_gcurvature` uninitialized variable** (`src/gdcp/gdcp_rules.jl:199`): Added explicit check `if !@isdefined(f_curvature)` to return `GUnknownCurvature` instead of erroring.
2. **Composition control flow** (`src/gdcp/gdcp_rules.jl:117-191`): Fixed the multi-branch composition logic to correctly handle cases where an atom has a GDCP rule but its arguments contain calls (inv, broadcast, affine_map).
3. **Duplicate `diag` registration** (`src/gdcp/spd.jl`): Removed duplicate registration that caused a warning.
4. **`lorentz_log_barrier` undefined variable** (`src/gdcp/lorentz.jl`): Fixed reference to undefined variable in the Lorentz log barrier function.
5. **Debug `println` removal**: Removed leftover debug print statements from production code paths.
6. **`sum_largest` off-by-one** (`src/atoms.jl`): Fixed indexing error in the sum of k largest elements.
7. **`AbstractMatrix_frac` typo** (`src/atoms.jl`): Fixed type name typo in fraction atom.
8. **`norm` p<1 convexity** (`src/atoms.jl`): Fixed convexity classification for norms with p<1 (not convex).

### Canonicalization Improvements
- Added `log(det(X))` -> `logdet(X)` rewrite rule (line 45 of `src/canon.jl`)
- Added `sum(diag(X))` -> `tr(X)` rewrite rule (line 48 of `src/canon.jl`)
- Added `log(a*b)` -> `log(a) + log(b)` in extended canonicalization (line 83 of `src/canon.jl`)
- Added `logdet(inv(X))` -> `-logdet(X)` rewrite (line 80 of `src/canon.jl`)
- Added documentation of known equivalent forms (`equivalent_forms()` function)

### MOI Cone Documentation
- Added comments in `src/gdcp/gdcp_rules.jl` (lines 11-21) mapping each GDCP atom to its corresponding MathOptInterface cone (LogDetConeTriangle, PositiveSemidefiniteConeTriangle, NormSpectralCone, etc.) to support the paper's claim about potential solver integration.

---

## Summary of Changes

| Category | Count | Key Items |
|----------|-------|-----------|
| Bug fixes | 8 | Composition logic, undefined vars, off-by-one errors |
| New experiments | 5 | MLE, DCP comparison, convergence, non-g-convex, expert |
| Canonicalization rules | 5 | logdet, tr, log product, double inv, logdet(inv) |
| New tests | 3 | DGCP-reduces-to-DCP, canonicalization, scaling |
| Documentation | 2 | Porting guide, MOI cone annotations |
| Paper edits | Multiple | Restructuring, Figure 1 caption, notation fixes |
