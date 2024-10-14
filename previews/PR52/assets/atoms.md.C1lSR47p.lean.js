import{_ as e,c as l,a5 as n,o as d}from"./chunks/framework.DdMmfCjO.js";const o=JSON.parse('{"title":"Atoms for DCP and DGCP","description":"","frontmatter":{},"headers":[],"relativePath":"atoms.md","filePath":"atoms.md","lastUpdated":null}'),i={name:"atoms.md"};function a(s,t,g,f,y,x){return d(),l("div",null,t[0]||(t[0]=[n('<h1 id="Atoms-for-DCP-and-DGCP" tabindex="-1">Atoms for DCP and DGCP <a class="header-anchor" href="#Atoms-for-DCP-and-DGCP" aria-label="Permalink to &quot;Atoms for DCP and DGCP {#Atoms-for-DCP-and-DGCP}&quot;">​</a></h1><p>This page is intended to be a reference for the atoms that are currently implemented in this package with their respective properties. As much as possible atoms are created with functions from base, standard libraries and popular packages, but we also inherit a few functions from the CVX family of packages such as <code>quad_form</code>, <code>quad_over_lin</code> etc. and also introduce some new functions in this package. Description of all such special functions implemented in this package is available in the <a href="/SymbolicAnalysis.jl/previews/PR52/functions#Special-functions">Special functions</a> section of the documentation.</p><h2 id="DCP-Atoms" tabindex="-1">DCP Atoms <a class="header-anchor" href="#DCP-Atoms" aria-label="Permalink to &quot;DCP Atoms {#DCP-Atoms}&quot;">​</a></h2><table tabindex="0"><thead><tr><th style="text-align:left;">Atom</th><th style="text-align:left;">Domain</th><th style="text-align:left;">Sign</th><th style="text-align:left;">Curvature</th><th style="text-align:left;">Monotonicity</th></tr></thead><tbody><tr><td style="text-align:left;">dot</td><td style="text-align:left;">(array_domain(ℝ), array_domain(ℝ))</td><td style="text-align:left;">AnySign</td><td style="text-align:left;">Affine</td><td style="text-align:left;">Increasing</td></tr><tr><td style="text-align:left;">dotsort</td><td style="text-align:left;">(array_domain(ℝ, 1), array_domain(ℝ, 1))</td><td style="text-align:left;">AnySign</td><td style="text-align:left;">Convex</td><td style="text-align:left;">(AnyMono, increasing_if_positive ∘ minimum)</td></tr><tr><td style="text-align:left;">StatsBase.geomean</td><td style="text-align:left;">array_domain(HalfLine{Real,:open}(), 1)</td><td style="text-align:left;">Positive</td><td style="text-align:left;">Concave</td><td style="text-align:left;">Increasing</td></tr><tr><td style="text-align:left;">StatsBase.harmmean</td><td style="text-align:left;">array_domain(HalfLine{Real,:open}(), 1)</td><td style="text-align:left;">Positive</td><td style="text-align:left;">Concave</td><td style="text-align:left;">Increasing</td></tr><tr><td style="text-align:left;">invprod</td><td style="text-align:left;">array_domain(HalfLine{Real,:open}())</td><td style="text-align:left;">Positive</td><td style="text-align:left;">Convex</td><td style="text-align:left;">Decreasing</td></tr><tr><td style="text-align:left;">eigmax</td><td style="text-align:left;">symmetric_domain()</td><td style="text-align:left;">AnySign</td><td style="text-align:left;">Convex</td><td style="text-align:left;">AnyMono</td></tr><tr><td style="text-align:left;">eigmin</td><td style="text-align:left;">symmetric_domain()</td><td style="text-align:left;">AnySign</td><td style="text-align:left;">Concave</td><td style="text-align:left;">AnyMono</td></tr><tr><td style="text-align:left;">eigsummax</td><td style="text-align:left;">(array_domain(ℝ, 2), ℝ)</td><td style="text-align:left;">AnySign</td><td style="text-align:left;">Convex</td><td style="text-align:left;">AnyMono</td></tr><tr><td style="text-align:left;">eigsummin</td><td style="text-align:left;">(array_domain(ℝ, 2), ℝ)</td><td style="text-align:left;">AnySign</td><td style="text-align:left;">Concave</td><td style="text-align:left;">AnyMono</td></tr><tr><td style="text-align:left;">logdet</td><td style="text-align:left;">semidefinite_domain()</td><td style="text-align:left;">AnySign</td><td style="text-align:left;">Concave</td><td style="text-align:left;">AnyMono</td></tr><tr><td style="text-align:left;">LogExpFunctions.logsumexp</td><td style="text-align:left;">array_domain(ℝ, 2)</td><td style="text-align:left;">AnySign</td><td style="text-align:left;">Convex</td><td style="text-align:left;">Increasing</td></tr><tr><td style="text-align:left;">matrix_frac</td><td style="text-align:left;">(array_domain(ℝ, 1), definite_domain())</td><td style="text-align:left;">AnySign</td><td style="text-align:left;">Convex</td><td style="text-align:left;">AnyMono</td></tr><tr><td style="text-align:left;">maximum</td><td style="text-align:left;">array_domain(ℝ)</td><td style="text-align:left;">AnySign</td><td style="text-align:left;">Convex</td><td style="text-align:left;">Increasing</td></tr><tr><td style="text-align:left;">minimum</td><td style="text-align:left;">array_domain(ℝ)</td><td style="text-align:left;">AnySign</td><td style="text-align:left;">Concave</td><td style="text-align:left;">Increasing</td></tr><tr><td style="text-align:left;">norm</td><td style="text-align:left;">(array_domain(ℝ), Interval{:closed, :open}(1, Inf))</td><td style="text-align:left;">Positive</td><td style="text-align:left;">Convex</td><td style="text-align:left;">increasing_if_positive</td></tr><tr><td style="text-align:left;">norm</td><td style="text-align:left;">(array_domain(ℝ), Interval{:closed, :open}(0, 1))</td><td style="text-align:left;">Positive</td><td style="text-align:left;">Convex</td><td style="text-align:left;">increasing_if_positive</td></tr><tr><td style="text-align:left;">perspective(f, x, s)</td><td style="text-align:left;">(function_domain(), ℝ, Positive)</td><td style="text-align:left;">Same as f</td><td style="text-align:left;">Same as f</td><td style="text-align:left;">AnyMono</td></tr><tr><td style="text-align:left;">quad_form</td><td style="text-align:left;">(array_domain(ℝ, 1), semidefinite_domain())</td><td style="text-align:left;">Positive</td><td style="text-align:left;">Convex</td><td style="text-align:left;">(increasing_if_positive, Increasing)</td></tr><tr><td style="text-align:left;">quad_over_lin</td><td style="text-align:left;">(array_domain(ℝ), HalfLine{Real,:open}())</td><td style="text-align:left;">Positive</td><td style="text-align:left;">Convex</td><td style="text-align:left;">(increasing_if_positive, Decreasing)</td></tr><tr><td style="text-align:left;">quad_over_lin</td><td style="text-align:left;">(ℝ, HalfLine{Real,:open}())</td><td style="text-align:left;">Positive</td><td style="text-align:left;">Convex</td><td style="text-align:left;">(increasing_if_positive, Decreasing)</td></tr><tr><td style="text-align:left;">sum</td><td style="text-align:left;">array_domain(ℝ, 2)</td><td style="text-align:left;">AnySign</td><td style="text-align:left;">Affine</td><td style="text-align:left;">Increasing</td></tr><tr><td style="text-align:left;">sum_largest</td><td style="text-align:left;">(array_domain(ℝ, 2), ℤ)</td><td style="text-align:left;">AnySign</td><td style="text-align:left;">Convex</td><td style="text-align:left;">Increasing</td></tr><tr><td style="text-align:left;">sum_smallest</td><td style="text-align:left;">(array_domain(ℝ, 2), ℤ)</td><td style="text-align:left;">AnySign</td><td style="text-align:left;">Concave</td><td style="text-align:left;">Increasing</td></tr><tr><td style="text-align:left;">tr</td><td style="text-align:left;">array_domain(ℝ, 2)</td><td style="text-align:left;">AnySign</td><td style="text-align:left;">Affine</td><td style="text-align:left;">Increasing</td></tr><tr><td style="text-align:left;">trinv</td><td style="text-align:left;">definite_domain()</td><td style="text-align:left;">Positive</td><td style="text-align:left;">Convex</td><td style="text-align:left;">AnyMono</td></tr><tr><td style="text-align:left;">tv</td><td style="text-align:left;">array_domain(ℝ, 1)</td><td style="text-align:left;">Positive</td><td style="text-align:left;">Convex</td><td style="text-align:left;">AnyMono</td></tr><tr><td style="text-align:left;">tv</td><td style="text-align:left;">array_domain(array_domain(ℝ, 2), 1)</td><td style="text-align:left;">Positive</td><td style="text-align:left;">Convex</td><td style="text-align:left;">AnyMono</td></tr><tr><td style="text-align:left;">abs</td><td style="text-align:left;">ℂ</td><td style="text-align:left;">Positive</td><td style="text-align:left;">Convex</td><td style="text-align:left;">increasing_if_positive</td></tr><tr><td style="text-align:left;">conj</td><td style="text-align:left;">ℂ</td><td style="text-align:left;">AnySign</td><td style="text-align:left;">Affine</td><td style="text-align:left;">AnyMono</td></tr><tr><td style="text-align:left;">exp</td><td style="text-align:left;">ℝ</td><td style="text-align:left;">Positive</td><td style="text-align:left;">Convex</td><td style="text-align:left;">Increasing</td></tr><tr><td style="text-align:left;">xlogx</td><td style="text-align:left;">ℝ</td><td style="text-align:left;">AnySign</td><td style="text-align:left;">Convex</td><td style="text-align:left;">AnyMono</td></tr><tr><td style="text-align:left;">huber</td><td style="text-align:left;">(ℝ, HalfLine())</td><td style="text-align:left;">Positive</td><td style="text-align:left;">Convex</td><td style="text-align:left;">increasing_if_positive</td></tr><tr><td style="text-align:left;">imag</td><td style="text-align:left;">ℂ</td><td style="text-align:left;">AnySign</td><td style="text-align:left;">Affine</td><td style="text-align:left;">AnyMono</td></tr><tr><td style="text-align:left;">inv</td><td style="text-align:left;">HalfLine{Real,:open}()</td><td style="text-align:left;">Positive</td><td style="text-align:left;">Convex</td><td style="text-align:left;">Decreasing</td></tr><tr><td style="text-align:left;">log</td><td style="text-align:left;">HalfLine{Real,:open}()</td><td style="text-align:left;">AnySign</td><td style="text-align:left;">Concave</td><td style="text-align:left;">Increasing</td></tr><tr><td style="text-align:left;">log</td><td style="text-align:left;">array_domain(ℝ, 2)</td><td style="text-align:left;">Positive</td><td style="text-align:left;">Concave</td><td style="text-align:left;">Increasing</td></tr><tr><td style="text-align:left;">inv</td><td style="text-align:left;">semidefinite_domain()</td><td style="text-align:left;">AnySign</td><td style="text-align:left;">Convex</td><td style="text-align:left;">Decreasing</td></tr><tr><td style="text-align:left;">sqrt</td><td style="text-align:left;">semidefinite_domain()</td><td style="text-align:left;">Positive</td><td style="text-align:left;">Concave</td><td style="text-align:left;">Increasing</td></tr><tr><td style="text-align:left;">kldivergence</td><td style="text-align:left;">(array_domain(HalfLine{Real,:open}, 1), array_domain(HalfLine{Real,:open}, 1))</td><td style="text-align:left;">Positive</td><td style="text-align:left;">Convex</td><td style="text-align:left;">AnyMono</td></tr><tr><td style="text-align:left;">lognormcdf</td><td style="text-align:left;">ℝ</td><td style="text-align:left;">Negative</td><td style="text-align:left;">Concave</td><td style="text-align:left;">Increasing</td></tr><tr><td style="text-align:left;">log1p</td><td style="text-align:left;">Interval{:open,:open}(-1, Inf)</td><td style="text-align:left;">Negative</td><td style="text-align:left;">Concave</td><td style="text-align:left;">Increasing</td></tr><tr><td style="text-align:left;">logistic</td><td style="text-align:left;">ℝ</td><td style="text-align:left;">Positive</td><td style="text-align:left;">Convex</td><td style="text-align:left;">Increasing</td></tr><tr><td style="text-align:left;">max</td><td style="text-align:left;">(ℝ, ℝ)</td><td style="text-align:left;">AnySign</td><td style="text-align:left;">Convex</td><td style="text-align:left;">Increasing</td></tr><tr><td style="text-align:left;">min</td><td style="text-align:left;">(ℝ, ℝ)</td><td style="text-align:left;">AnySign</td><td style="text-align:left;">Concave</td><td style="text-align:left;">Increasing</td></tr><tr><td style="text-align:left;">^(x, i)</td><td style="text-align:left;">See below</td><td style="text-align:left;">See below</td><td style="text-align:left;">See below</td><td style="text-align:left;">See below</td></tr><tr><td style="text-align:left;">real</td><td style="text-align:left;">ℂ</td><td style="text-align:left;">AnySign</td><td style="text-align:left;">Affine</td><td style="text-align:left;">Increasing</td></tr><tr><td style="text-align:left;">rel_entr</td><td style="text-align:left;">(HalfLine{Real,:open}(), HalfLine{Real,:open}())</td><td style="text-align:left;">AnySign</td><td style="text-align:left;">Convex</td><td style="text-align:left;">(AnyMono, Decreasing)</td></tr><tr><td style="text-align:left;">sqrt</td><td style="text-align:left;">HalfLine()</td><td style="text-align:left;">Positive</td><td style="text-align:left;">Concave</td><td style="text-align:left;">Increasing</td></tr><tr><td style="text-align:left;">xexpx</td><td style="text-align:left;">HalfLine</td><td style="text-align:left;">Positive</td><td style="text-align:left;">Convex</td><td style="text-align:left;">Increasing</td></tr><tr><td style="text-align:left;">conv</td><td style="text-align:left;">(array_domain(ℝ, 1), array_domain(ℝ, 1))</td><td style="text-align:left;">AnySign</td><td style="text-align:left;">Affine</td><td style="text-align:left;">AnyMono</td></tr><tr><td style="text-align:left;">cumsum</td><td style="text-align:left;">array_domain(ℝ)</td><td style="text-align:left;">AnySign</td><td style="text-align:left;">Affine</td><td style="text-align:left;">Increasing</td></tr><tr><td style="text-align:left;">diagm</td><td style="text-align:left;">array_domain(ℝ, 1)</td><td style="text-align:left;">AnySign</td><td style="text-align:left;">Affine</td><td style="text-align:left;">Increasing</td></tr><tr><td style="text-align:left;">diag</td><td style="text-align:left;">array_domain(ℝ, 2)</td><td style="text-align:left;">AnySign</td><td style="text-align:left;">Affine</td><td style="text-align:left;">Increasing</td></tr><tr><td style="text-align:left;">diff</td><td style="text-align:left;">array_domain(ℝ)</td><td style="text-align:left;">AnySign</td><td style="text-align:left;">Affine</td><td style="text-align:left;">Increasing</td></tr><tr><td style="text-align:left;">kron</td><td style="text-align:left;">(array_domain(ℝ, 2), array_domain(ℝ, 2))</td><td style="text-align:left;">AnySign</td><td style="text-align:left;">Affine</td><td style="text-align:left;">Increasing</td></tr></tbody></table><h3 id="Special-Cases-for-(x,-i)" tabindex="-1">Special Cases for ^(x, i) <a class="header-anchor" href="#Special-Cases-for-(x,-i)" aria-label="Permalink to &quot;Special Cases for ^(x, i) {#Special-Cases-for-(x,-i)}&quot;">​</a></h3><table tabindex="0"><thead><tr><th style="text-align:left;">Condition on i</th><th style="text-align:left;">Domain</th><th style="text-align:left;">Sign</th><th style="text-align:left;">Curvature</th><th style="text-align:left;">Monotonicity</th></tr></thead><tbody><tr><td style="text-align:left;">i = 1</td><td style="text-align:left;">ℝ</td><td style="text-align:left;">AnySign</td><td style="text-align:left;">Affine</td><td style="text-align:left;">Increasing</td></tr><tr><td style="text-align:left;">i is even integer</td><td style="text-align:left;">ℝ</td><td style="text-align:left;">Positive</td><td style="text-align:left;">Convex</td><td style="text-align:left;">increasing_if_positive</td></tr><tr><td style="text-align:left;">i is odd integer</td><td style="text-align:left;">HalfLine()</td><td style="text-align:left;">Positive</td><td style="text-align:left;">Convex</td><td style="text-align:left;">Increasing</td></tr><tr><td style="text-align:left;">i ≥ 1</td><td style="text-align:left;">HalfLine()</td><td style="text-align:left;">Positive</td><td style="text-align:left;">Convex</td><td style="text-align:left;">Increasing</td></tr><tr><td style="text-align:left;">0 &lt; i &lt; 1</td><td style="text-align:left;">HalfLine()</td><td style="text-align:left;">Positive</td><td style="text-align:left;">Concave</td><td style="text-align:left;">Increasing</td></tr><tr><td style="text-align:left;">i &lt; 0</td><td style="text-align:left;">HalfLine{Float64,:closed}()</td><td style="text-align:left;">Positive</td><td style="text-align:left;">Convex</td><td style="text-align:left;">Increasing</td></tr></tbody></table><h2 id="DGCP-Atoms-(Symmetric-Positive-Definite)" tabindex="-1">DGCP Atoms (Symmetric Positive Definite) <a class="header-anchor" href="#DGCP-Atoms-(Symmetric-Positive-Definite)" aria-label="Permalink to &quot;DGCP Atoms (Symmetric Positive Definite) {#DGCP-Atoms-(Symmetric-Positive-Definite)}&quot;">​</a></h2><table tabindex="0"><thead><tr><th style="text-align:left;">Atom</th><th style="text-align:left;">Sign</th><th style="text-align:left;">Geodesic Curvature</th><th style="text-align:left;">Monotonicity</th></tr></thead><tbody><tr><td style="text-align:left;">LinearAlgebra.logdet</td><td style="text-align:left;">Positive</td><td style="text-align:left;">GLinear</td><td style="text-align:left;">GIncreasing</td></tr><tr><td style="text-align:left;">conjugation</td><td style="text-align:left;">Positive</td><td style="text-align:left;">GConvex</td><td style="text-align:left;">GIncreasing</td></tr><tr><td style="text-align:left;">LinearAlgebra.tr</td><td style="text-align:left;">Positive</td><td style="text-align:left;">GConvex</td><td style="text-align:left;">GIncreasing</td></tr><tr><td style="text-align:left;">sum</td><td style="text-align:left;">Positive</td><td style="text-align:left;">GConvex</td><td style="text-align:left;">GIncreasing</td></tr><tr><td style="text-align:left;">adjoint</td><td style="text-align:left;">Positive</td><td style="text-align:left;">GLinear</td><td style="text-align:left;">GIncreasing</td></tr><tr><td style="text-align:left;">scalar_mat</td><td style="text-align:left;">Positive</td><td style="text-align:left;">GConvex</td><td style="text-align:left;">GIncreasing</td></tr><tr><td style="text-align:left;">LinearAlgebra.diag</td><td style="text-align:left;">Positive</td><td style="text-align:left;">GConvex</td><td style="text-align:left;">GIncreasing</td></tr><tr><td style="text-align:left;">sdivergence</td><td style="text-align:left;">Positive</td><td style="text-align:left;">GConvex</td><td style="text-align:left;">GIncreasing</td></tr><tr><td style="text-align:left;">Manifolds.distance</td><td style="text-align:left;">Positive</td><td style="text-align:left;">GConvex</td><td style="text-align:left;">GAnyMono</td></tr><tr><td style="text-align:left;">SymbolicAnalysis.quad_form</td><td style="text-align:left;">Positive</td><td style="text-align:left;">GConvex</td><td style="text-align:left;">GIncreasing</td></tr><tr><td style="text-align:left;">LinearAlgebra.eigmax</td><td style="text-align:left;">Positive</td><td style="text-align:left;">GConvex</td><td style="text-align:left;">GIncreasing</td></tr><tr><td style="text-align:left;">log_quad_form</td><td style="text-align:left;">Positive</td><td style="text-align:left;">GConvex</td><td style="text-align:left;">GIncreasing</td></tr><tr><td style="text-align:left;">inv</td><td style="text-align:left;">Positive</td><td style="text-align:left;">GConvex</td><td style="text-align:left;">GDecreasing</td></tr><tr><td style="text-align:left;">diag</td><td style="text-align:left;">Positive</td><td style="text-align:left;">GConvex</td><td style="text-align:left;">GIncreasing</td></tr><tr><td style="text-align:left;">eigsummax</td><td style="text-align:left;">Positive</td><td style="text-align:left;">GConvex</td><td style="text-align:left;">GIncreasing</td></tr><tr><td style="text-align:left;">schatten_norm</td><td style="text-align:left;">Positive</td><td style="text-align:left;">GConvex</td><td style="text-align:left;">GIncreasing</td></tr><tr><td style="text-align:left;">sum_log_eigmax</td><td style="text-align:left;">Positive</td><td style="text-align:left;">GConvex</td><td style="text-align:left;">GIncreasing</td></tr><tr><td style="text-align:left;">affine_map</td><td style="text-align:left;">Positive</td><td style="text-align:left;">GConvex</td><td style="text-align:left;">GIncreasing</td></tr><tr><td style="text-align:left;">hadamard_product</td><td style="text-align:left;">Positive</td><td style="text-align:left;">GConvex</td><td style="text-align:left;">GIncreasing</td></tr></tbody></table>',8)]))}const c=e(i,[["render",a]]);export{o as __pageData,c as default};
