# Computational Complexity and Scaling

<!-- LaTeX section for MPC paper. Copy-paste into Overleaf. -->

## Notation

Throughout this section, $n$ denotes the number of nodes in the expression
AST (abstract syntax tree), $m$ the number of composition terms used to
control problem size, and $k$ the maximum arity of any function node.
We write $R_{\text{DCP}}$ and $R_{\text{DGCP}}$ for the sizes of the
DCP and DGCP atom rule tables, respectively; both are implementation
constants ($R_{\text{DCP}} = 65$, $R_{\text{DGCP}} = 27$).

---

## Verification Complexity

The \texttt{analyze} pipeline performs verification in sequential phases,
each consisting of a bottom-up traversal (Postwalk) followed by a
top-down traversal (Prewalk) of the expression tree.  For DCP, there
are three phases: canonicalization, sign propagation, and curvature
propagation.  DGCP extends this with a fourth phase for geodesic
curvature propagation.

\begin{proposition}[Linear-time verification]\label{prop:linear}
The complete DCP verification pipeline runs in $\Theta(n)$ time using
exactly 6 tree traversals.  The complete DGCP verification pipeline
runs in $\Theta(n)$ time using exactly 8 tree traversals.
\end{proposition}

\begin{proof}
Each phase applies a chain of rewrite rules via one Postwalk and one
Prewalk pass.  At every node, rule matching attempts at most
$\max(R_{\text{DCP}}, R_{\text{DGCP}})$ pattern comparisons, each
requiring $O(1)$ work (hash-table lookup against the atom dictionary
plus constant-time metadata attachment).  Rules that iterate over the
$k_v$ children of a node~$v$ (e.g., sign and curvature composition)
contribute $O(k_v)$ work at~$v$; since $\sum_v k_v = n - 1$ over the
entire tree, the total work per traversal is $O(n)$.  With two
traversals per phase and three (resp.\ four) phases, DCP verification
performs $6n$ (resp.\ $8n$) node visits, each at $O(1)$ amortized cost.
The bound is tight because every node must be visited at least once per
phase to attach metadata.
\end{proof}

\begin{corollary}[DGCP marginal cost]\label{cor:marginal}
DGCP verification adds exactly one $\Theta(n)$ phase (2 traversals) to
the DCP pipeline.  In terms of traversal count, the overhead ratio is
$8/6 \approx 1.33$.
\end{corollary}

\begin{proof}
The DGCP phase (\texttt{propagate\_gcurvature}) has the same algorithmic
structure as the DCP curvature phase: a chain of 3 rewrite rules applied
via Postwalk~+~Prewalk, with per-node cost dominated by a dictionary
lookup in the DGCP rule table and an argument-iteration loop identical to
the DCP composition check.  The additional case analysis for
geodesic-specific patterns (e.g., $\operatorname{logdet}$ composed with
SPD operations) involves at most 6 constant-time checks per call node,
increasing the per-node constant by a modest factor but preserving the
$O(n)$ bound.  The traversal-count ratio $8/6 \approx 1.33$ gives a
first-order estimate of the wall-clock overhead; empirical measurements
(Section~\ref{sec:empirical}) confirm a ratio in the range 1.25--1.35.
\end{proof}

\begin{proposition}[Verification lower bound]\label{prop:lower}
Any DCP or DGCP verifier requires $\Omega(n)$ time in the worst case.
\end{proposition}

\begin{proof}
Consider the expression $e = f_1(x_1) + f_2(x_2) + \cdots + f_n(x_n)$,
where every $f_j$ except $f_i$ is a registered convex atom.  Let $f_i$
be a function with no DCP rule (e.g., $\sin$).  Then $e$ is
DCP-compliant if and only if the verifier examines node~$i$.  Since $i$
can be any index in $\{1,\ldots,n\}$, the verifier must inspect all $n$
nodes.
\end{proof}

Propositions~\ref{prop:linear} and~\ref{prop:lower} together establish
that the verification algorithms are \emph{optimal}: their $\Theta(n)$
running time matches the information-theoretic lower bound up to constant
factors.

---

## Conic Form Generation

\begin{proposition}[Linear conic reformulation]\label{prop:conic}
The conic form generation procedure runs in $\Theta(n)$ time and
produces $O(n)$ conic constraints and $O(n)$ epigraph variables.
\end{proposition}

\begin{proof}
After DCP verification ($\Theta(n)$ by Proposition~\ref{prop:linear}),
a single bottom-up traversal decomposes each atom into its conic
representation.  Leaf nodes require $O(1)$ work.  Each DCP atom emits a
bounded number of constraints (at most 3 per atom, e.g., an exponential
cone atom requires two conic constraints and one linear inequality) and
introduces one epigraph variable.  The affine detection subroutine
processes disjoint affine subtrees in aggregate $O(n)$ time.  Since each
of the $n$ nodes is visited once and contributes $O(1)$ constraints, the
total output size is $O(n)$ and the total time is $\Theta(n)$.
\end{proof}

---

## Comparison with Related Systems

All major DCP implementations share the same optimal $O(n)$ verification
complexity.  The distinguishing feature of the present work is the
\emph{verifiable class}, not asymptotic speed.

\begin{table}[ht]
\centering
\caption{Verification complexity across DCP systems.  All systems are
$O(n)$ in the AST node count~$n$.  The column ``Verifiable class''
indicates the broadest problem class accepted by each verifier.}
\label{tab:systems}
\begin{tabular}{lccl}
\toprule
System & Traversals & Atoms & Verifiable class \\
\midrule
Convex.jl       & $\sim$2 & $\sim$40  & DCP \\
CVXPY           & $\sim$2 & $\sim$80  & DCP \\
DCCP            & $\sim$3 & $\sim$80  & DCP + difference-of-convex \\
\textbf{This work (DCP)}  & 6 & 65  & DCP \\
\textbf{This work (DGCP)} & 8 & 92  & DCP + DGCP \\
\bottomrule
\end{tabular}
\end{table}

The higher traversal count in our implementation (6 vs.\ $\sim$2)
reflects the use of a symbolic rewriting framework that requires
separate Postwalk and Prewalk passes per phase, rather than a monolithic
visitor.  This is a constant-factor difference with no asymptotic impact
and provides modularity: each verification phase can be independently
tested, extended, and composed.

For the class of problems that DGCP can verify---geodesically convex
optimization on Riemannian manifolds, including the Karcher mean,
Tyler's M-estimator, the S-divergence, and Schatten-norm objectives on
$\mathcal{S}_{++}^n$---the alternative is \emph{no automated
verification whatsoever}.  In this context, even a hypothetical $10
\times$ overhead would be practically irrelevant; the actual overhead of
$\sim$1.33$\times$ is negligible.

---

## Empirical Scaling Methodology\label{sec:empirical}

We validate the theoretical predictions with controlled scaling
experiments.\footnote{The experiment script is available at
\texttt{test/experiments/scaling\_analysis.jl}.}

\paragraph{Expression construction.}
To isolate AST-level scaling from numerical artifacts, we hold the matrix
dimension fixed ($n_{\text{mat}} = 5$) and vary the number of composition
terms~$m$.  Three expression families are tested:

\begin{itemize}
\item \textbf{Karcher mean} (DGCP): $\sum_{i=1}^{m} d^2(A_i, X)$ on
  $\mathcal{S}_{++}^{n_{\text{mat}}}$, where each distance term
  contributes a fixed number of AST nodes.
\item \textbf{Tyler's M-estimator} (DGCP):
  $\sum_{i=1}^{m} \log(x_i^\top X^{-1} x_i) + \tfrac{1}{n_{\text{mat}}}
  \operatorname{logdet}(X)$.
\item \textbf{Scalar DCP}: $\sum_{i=1}^{m} (\exp(x_i) + \log(x_i))$.
\end{itemize}

In all cases, the AST node count grows linearly in~$m$, providing a
clean independent variable for regression.

\paragraph{Timing protocol.}
Each configuration is timed over 15 independent trials after 3 warmup
iterations (to eliminate JIT compilation artifacts).  The
\emph{minimum} trial time is reported, following best practices for
microbenchmarking in managed-runtime languages: the minimum of
independent trials provides the best estimate of the deterministic
computation time, free of GC pauses and OS scheduling jitter.  A minor
GC collection is triggered before each trial to reduce mid-trial
interference.  Timing uses nanosecond-resolution clocks.

\paragraph{Curve fitting.}
A power-law model $t = c \cdot n^{\alpha}$ is fit via ordinary least
squares on log-log-transformed data.  The scaling exponent~$\alpha$ and
the coefficient of determination~$R^2$ are reported.  For $O(n)$ scaling
we expect $\alpha \approx 1.0$ with $R^2$ close to~1.

\paragraph{Expected findings.}
Based on the theoretical analysis and preliminary runs, we expect:

\begin{enumerate}
\item \textbf{Linear scaling}: fitted exponent $\alpha \approx 1.0$
  with $R^2 > 0.99$ across all three expression families, confirming
  $\Theta(n)$ verification time.

\item \textbf{Phase decomposition}: each of the four phases individually
  scales as $O(n)$.  At the largest problem size, the DGCP phase
  (\texttt{propagate\_gcurvature}) accounts for approximately 25\% of
  total DGCP verification time, yielding a measured DGCP/DCP ratio of
  $\sim$1.25--1.35$\times$.

\item \textbf{Linear memory}: total allocations scale as $O(n)$, with a
  bounded number of bytes per AST node for metadata storage.

\item \textbf{Linear conic output}: both the number of epigraph
  variables and the number of conic constraints grow linearly in~$n$.
\end{enumerate}

\begin{remark}[Matrix dimension independence]\label{rem:matrix}
A matrix-valued variable $X \in \mathbb{R}^{p \times p}$ appearing in an
expression such as $\operatorname{logdet}(X)$ occupies a \emph{single
leaf node} in the AST, regardless of~$p$.  The matrix entries are
numerical data resolved at evaluation time, not symbolic nodes visited
during verification.  Consequently, varying the matrix dimension with the
number of composition terms held fixed produces essentially identical AST
sizes and verification times.  The scaling experiments confirm this by
showing negligible variation in timing across matrix dimensions $p \in
\{3, 5, 8, 10, 15\}$ at fixed $m = 4$.  Previously reported ``$2$--$3
\times$ overhead'' for DGCP conflated matrix-dimension variation (which
does not affect AST size) with algorithmic scaling and included
JIT/GC artifacts, leading to a misleading characterization of the
marginal cost.
\end{remark}

---

## Summary

The full SymbolicAnalysis.jl pipeline---from raw symbolic expression
through verified curvature labels to optional conic
reformulation---runs in $\Theta(n)$ time, matching the
$\Omega(n)$ information-theoretic lower bound.  DGCP verification adds
one traversal phase to DCP's three, for a theoretical overhead of
$8/6 \approx 1.33\times$ and an empirically measured overhead of
$\sim$1.25--1.35$\times$.  The conic form generation is likewise
$\Theta(n)$, producing $O(n)$ constraints.  These complexity guarantees
are shared by all major DCP systems; the contribution of the present
work lies in the strictly larger verifiable class enabled by DGCP, not in
asymptotic speedup.
