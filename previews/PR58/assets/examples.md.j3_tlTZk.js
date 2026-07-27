import{_ as s,c as n,o as e,ai as p}from"./chunks/framework.CTH-wruX.js";const m=JSON.parse('{"title":"Examples","description":"","frontmatter":{},"headers":[],"relativePath":"examples.md","filePath":"examples.md","lastUpdated":null}'),i={name:"examples.md"};function l(o,a,t,c,r,d){return e(),n("div",null,a[0]||(a[0]=[p(`<h1 id="examples" tabindex="-1">Examples <a class="header-anchor" href="#examples" aria-label="Permalink to &quot;Examples&quot;">​</a></h1><p>Here are some examples demonstrating the use of the <code>analyze</code> function from the <code>SymbolicAnalysis</code> package.</p><h2 id="Basic-Expression-Analysis" tabindex="-1">Basic Expression Analysis <a class="header-anchor" href="#Basic-Expression-Analysis" aria-label="Permalink to &quot;Basic Expression Analysis {#Basic-Expression-Analysis}&quot;">​</a></h2><div class="language-@example vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">@example</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line highlighted"><span>using SymbolicAnalysis, Symbolics, Symbolics.DomainSets</span></span>
<span class="line"><span></span></span>
<span class="line"><span>@variables x</span></span>
<span class="line"><span>ex1 = exp(x) - log(x)</span></span>
<span class="line"><span>result = analyze(ex1)</span></span>
<span class="line"><span>@show result.curvature</span></span></code></pre></div><p>This example analyzes a simple expression <code>exp(x) - log(x)</code>, determining that it&#39;s convex and can have any sign.</p><h2 id="Analysis-on-Manifolds" tabindex="-1">Analysis on Manifolds <a class="header-anchor" href="#Analysis-on-Manifolds" aria-label="Permalink to &quot;Analysis on Manifolds {#Analysis-on-Manifolds}&quot;">​</a></h2><p>We can perform DGCP analysis on the Symmetric Positive Definite (SPD) manifold by passing a manifold from <a href="https://juliamanifolds.github.io/Manifolds.jl/stable/" target="_blank" rel="noreferrer">Manifolds.jl</a> to the <code>analyze</code> function. We consider the Karcher mean problem which involves finding the geometric mean of SPD matrices:</p><div class="language-@example vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">@example</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line highlighted"><span>using SymbolicAnalysis, Symbolics, Manifolds, LinearAlgebra</span></span>
<span class="line"><span></span></span>
<span class="line"><span>@variables X[1:5, 1:5]</span></span>
<span class="line"><span></span></span>
<span class="line"><span>M = SymmetricPositiveDefinite(5)</span></span>
<span class="line"><span></span></span>
<span class="line"><span>As = [rand(5, 5) for i in 1:5]</span></span>
<span class="line"><span>As = [As[i] * As[i]&#39; for i in 1:5]  # Make them SPD</span></span>
<span class="line"><span></span></span>
<span class="line"><span>ex2 = sum(Manifolds.distance(M, As[i], X)^2 for i in 1:5)</span></span>
<span class="line"><span>result = analyze(ex2, M)</span></span>
<span class="line"><span>@show result.curvature</span></span>
<span class="line"><span>@show result.gcurvature</span></span></code></pre></div><p>This analysis shows that the Karcher mean objective function is geodesically convex on the SPD manifold.</p><h3 id="Domain-aware-analysis" tabindex="-1">Domain aware analysis <a class="header-anchor" href="#Domain-aware-analysis" aria-label="Permalink to &quot;Domain aware analysis {#Domain-aware-analysis}&quot;">​</a></h3><p>We can also assert the domain of the variable by assigning <code>VarDomain</code> metadata that takes a <code>Domain</code> from the <a href="https://juliaapproximation.github.io/DomainSets.jl/dev/" target="_blank" rel="noreferrer">DomainSets.jl</a> package.</p><div class="language-@example vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">@example</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line highlighted"><span>@variables x y</span></span>
<span class="line"><span></span></span>
<span class="line"><span>x = setmetadata(</span></span>
<span class="line"><span>    x,</span></span>
<span class="line"><span>    SymbolicAnalysis.VarDomain,</span></span>
<span class="line"><span>    OpenInterval(0,1),</span></span>
<span class="line"><span>)</span></span>
<span class="line"><span></span></span>
<span class="line"><span>y = setmetadata(</span></span>
<span class="line"><span>    y,</span></span>
<span class="line"><span>    SymbolicAnalysis.VarDomain,</span></span>
<span class="line"><span>    OpenInterval(0,1),</span></span>
<span class="line"><span>)</span></span>
<span class="line"><span></span></span>
<span class="line"><span>ex = SymbolicAnalysis.quad_over_lin(x - y, 1 - max(x, y))</span></span>
<span class="line"><span>result = analyze(ex)</span></span>
<span class="line"><span>@show result.curvature</span></span></code></pre></div><p>This example analyzes a quadratic expression over a linear expression, showing that it&#39;s convex.</p><h2 id="Analysis-on-the-Lorentz-Manifold" tabindex="-1">Analysis on the Lorentz Manifold <a class="header-anchor" href="#Analysis-on-the-Lorentz-Manifold" aria-label="Permalink to &quot;Analysis on the Lorentz Manifold {#Analysis-on-the-Lorentz-Manifold}&quot;">​</a></h2><p>We can also perform DGCP analysis on the Lorentz manifold, which is a model of hyperbolic space:</p><div class="language-@example vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">@example</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line highlighted"><span>using SymbolicAnalysis, Symbolics, Manifolds, LinearAlgebra</span></span>
<span class="line"><span></span></span>
<span class="line"><span># Create a Lorentz manifold of dimension 2 (3D ambient space)</span></span>
<span class="line"><span>M = Lorentz(2)</span></span>
<span class="line"><span></span></span>
<span class="line"><span># Define symbolic variables and fixed points</span></span>
<span class="line"><span>@variables p[1:3]</span></span>
<span class="line"><span>q = [0.0, 0.0, 1.0]  # A point on the Lorentz model</span></span>
<span class="line"><span></span></span>
<span class="line"><span># Create a composite function from Lorentz atoms</span></span>
<span class="line"><span>ex = 2.0 * Manifolds.distance(M, q, p) + </span></span>
<span class="line"><span>     SymbolicAnalysis.lorentz_log_barrier(p)</span></span>
<span class="line"><span></span></span>
<span class="line"><span># Analyze the expression</span></span>
<span class="line"><span>result = analyze(ex, M)</span></span>
<span class="line"><span>@show result.gcurvature</span></span></code></pre></div><p>This example shows that the sum of the Lorentz distance function and the log-barrier function is geodesically convex on the Lorentz manifold.</p>`,17)]))}const y=s(i,[["render",l]]);export{m as __pageData,y as default};
