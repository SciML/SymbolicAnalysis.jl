import{_ as s,c as i,o as a,a6 as e}from"./chunks/framework.D_dyxcd-.js";const g=JSON.parse('{"title":"Special functions","description":"","frontmatter":{},"headers":[],"relativePath":"functions.md","filePath":"functions.md","lastUpdated":null}'),t={name:"functions.md"},n=e('<h1 id="Special-functions" tabindex="-1">Special functions <a class="header-anchor" href="#Special-functions" aria-label="Permalink to &quot;Special functions {#Special-functions}&quot;">​</a></h1><p>Since some atoms are not available in the base langugage or other packages we have implemented them here.</p><div style="border-width:1px;border-style:solid;border-color:black;padding:1em;border-radius:25px;"><a id="SymbolicAnalysis.conjugation-Tuple{Any, Any}" href="#SymbolicAnalysis.conjugation-Tuple{Any, Any}">#</a> <b><u>SymbolicAnalysis.conjugation</u></b> — <i>Method</i>. <div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">conjugation</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(X, B)</span></span></code></pre></div><p>Conjugation of a matrix <code>X</code> by a matrix <code>B</code> is defined as <code>B&#39;X*B</code>.</p><p><strong>Arguments</strong></p><div class="language- vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang"></span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span>- `X::Matrix`: A symmetric positive definite matrix.</span></span>\n<span class="line"><span>- `B::Matrix`: A matrix.</span></span></code></pre></div><p><a href="https://github.com/Vaibhavdixit02/SymbolicAnalysis.jl" target="_blank" rel="noreferrer">source</a></p></div><br><div style="border-width:1px;border-style:solid;border-color:black;padding:1em;border-radius:25px;"><a id="SymbolicAnalysis.dotsort-Tuple{AbstractVector, AbstractVector}" href="#SymbolicAnalysis.dotsort-Tuple{AbstractVector, AbstractVector}">#</a> <b><u>SymbolicAnalysis.dotsort</u></b> — <i>Method</i>. <div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">dotsort</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(x, y)</span></span></code></pre></div><p>Sorts <code>x</code> and <code>y</code> and returns the dot product of the sorted vectors.</p><p><strong>Arguments</strong></p><div class="language- vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang"></span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span>- `x::AbstractVector`: A vector.</span></span>\n<span class="line"><span>- `y::AbstractVector`: A vector.</span></span></code></pre></div><p><a href="https://github.com/Vaibhavdixit02/SymbolicAnalysis.jl" target="_blank" rel="noreferrer">source</a></p></div><br><div style="border-width:1px;border-style:solid;border-color:black;padding:1em;border-radius:25px;"><a id="SymbolicAnalysis.eigsummax-Tuple{LinearAlgebra.Symmetric, Int64}" href="#SymbolicAnalysis.eigsummax-Tuple{LinearAlgebra.Symmetric, Int64}">#</a> <b><u>SymbolicAnalysis.eigsummax</u></b> — <i>Method</i>. <div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">eigsummax</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(m</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">::</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">Symmetric</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, k)</span></span></code></pre></div><p>Returns the sum of the <code>k</code> largest eigenvalues of <code>m</code>.</p><p><strong>Arguments</strong></p><div class="language- vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang"></span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span>- `m::Symmetric`: A symmetric matrix.</span></span>\n<span class="line"><span>- `k::Int`: The number of largest eigenvalues to sum.</span></span></code></pre></div><p><a href="https://github.com/Vaibhavdixit02/SymbolicAnalysis.jl" target="_blank" rel="noreferrer">source</a></p></div><br><div style="border-width:1px;border-style:solid;border-color:black;padding:1em;border-radius:25px;"><a id="SymbolicAnalysis.eigsummin-Tuple{LinearAlgebra.Symmetric, Int64}" href="#SymbolicAnalysis.eigsummin-Tuple{LinearAlgebra.Symmetric, Int64}">#</a> <b><u>SymbolicAnalysis.eigsummin</u></b> — <i>Method</i>. <div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">eigsummin</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(m</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">::</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">Symmetric</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, k)</span></span></code></pre></div><p>Returns the sum of the <code>k</code> smallest eigenvalues of <code>m</code>.</p><p><strong>Arguments</strong></p><div class="language- vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang"></span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span>- `m::Symmetric`: A symmetric matrix.</span></span>\n<span class="line"><span>- `k::Int`: The number of smallest eigenvalues to sum.</span></span></code></pre></div><p><a href="https://github.com/Vaibhavdixit02/SymbolicAnalysis.jl" target="_blank" rel="noreferrer">source</a></p></div><br><div style="border-width:1px;border-style:solid;border-color:black;padding:1em;border-radius:25px;"><a id="SymbolicAnalysis.huber" href="#SymbolicAnalysis.huber">#</a> <b><u>SymbolicAnalysis.huber</u></b> — <i>Function</i>. <div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">huber</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(x, M</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">1</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">)</span></span></code></pre></div><p>Returns the Huber loss function of <code>x</code> with threshold <code>M</code>.</p><p><strong>Arguments</strong></p><div class="language- vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang"></span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span>- `x::Number`: A number.</span></span>\n<span class="line"><span>- `M::Number`: The threshold.</span></span></code></pre></div><p><a href="https://github.com/Vaibhavdixit02/SymbolicAnalysis.jl" target="_blank" rel="noreferrer">source</a></p></div><br><div style="border-width:1px;border-style:solid;border-color:black;padding:1em;border-radius:25px;"><a id="SymbolicAnalysis.invprod-Tuple{AbstractVector}" href="#SymbolicAnalysis.invprod-Tuple{AbstractVector}">#</a> <b><u>SymbolicAnalysis.invprod</u></b> — <i>Method</i>. <div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">invprod</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(x</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">::</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">AbstractVector</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">)</span></span></code></pre></div><p>Returns the inverse of the product of the elements of <code>x</code>.</p><p><strong>Arguments</strong></p><div class="language- vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang"></span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span>- `x::AbstractVector`: A vector.</span></span></code></pre></div><p><a href="https://github.com/Vaibhavdixit02/SymbolicAnalysis.jl" target="_blank" rel="noreferrer">source</a></p></div><br><div style="border-width:1px;border-style:solid;border-color:black;padding:1em;border-radius:25px;"><a id="SymbolicAnalysis.lognormcdf-Tuple{Number}" href="#SymbolicAnalysis.lognormcdf-Tuple{Number}">#</a> <b><u>SymbolicAnalysis.lognormcdf</u></b> — <i>Method</i>. <div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">lognormcdf</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(x</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">::</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">Number</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">)</span></span></code></pre></div><p>Returns the log of the normal cumulative distribution function of <code>x</code>.</p><p><strong>Arguments</strong></p><div class="language- vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang"></span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span>- `x::Number`: A number.</span></span></code></pre></div><p><a href="https://github.com/Vaibhavdixit02/SymbolicAnalysis.jl" target="_blank" rel="noreferrer">source</a></p></div><br><div style="border-width:1px;border-style:solid;border-color:black;padding:1em;border-radius:25px;"><a id="SymbolicAnalysis.matrix_frac-Tuple{AbstractVector, AbstractMatrix}" href="#SymbolicAnalysis.matrix_frac-Tuple{AbstractVector, AbstractMatrix}">#</a> <b><u>SymbolicAnalysis.matrix_frac</u></b> — <i>Method</i>. <div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">matrix_frac</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(x</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">::</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">AbstractVector</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, P</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">::</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">AbstractMatrix</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">)</span></span></code></pre></div><p>Returns the quadratic form <code>x&#39; * P^{-1} * x</code>.</p><p><strong>Arguments</strong></p><div class="language- vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang"></span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span>- `x::AbstractVector`: A vector.</span></span>\n<span class="line"><span>- `P::AbstractMatrix`: A matrix.</span></span></code></pre></div><p><a href="https://github.com/Vaibhavdixit02/SymbolicAnalysis.jl" target="_blank" rel="noreferrer">source</a></p></div><br><div style="border-width:1px;border-style:solid;border-color:black;padding:1em;border-radius:25px;"><a id="SymbolicAnalysis.perspective-Tuple{Function, Any, Number}" href="#SymbolicAnalysis.perspective-Tuple{Function, Any, Number}">#</a> <b><u>SymbolicAnalysis.perspective</u></b> — <i>Method</i>. <div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">perspective</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(f</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">::</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">Function</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, x, s</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">::</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">Number</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">)</span></span></code></pre></div><p>Returns the perspective function <code>s * f(x / s)</code>.</p><p><strong>Arguments</strong></p><div class="language- vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang"></span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span>- `f::Function`: A function.</span></span>\n<span class="line"><span>- `x`: A number.</span></span>\n<span class="line"><span>- `s::Number`: A positive number.</span></span></code></pre></div><p><a href="https://github.com/Vaibhavdixit02/SymbolicAnalysis.jl" target="_blank" rel="noreferrer">source</a></p></div><br><div style="border-width:1px;border-style:solid;border-color:black;padding:1em;border-radius:25px;"><a id="SymbolicAnalysis.quad_form-Tuple{AbstractVector, AbstractMatrix}" href="#SymbolicAnalysis.quad_form-Tuple{AbstractVector, AbstractMatrix}">#</a> <b><u>SymbolicAnalysis.quad_form</u></b> — <i>Method</i>. <div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">quad_form</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(x</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">::</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">AbstractVector</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, P</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">::</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">AbstractMatrix</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">)</span></span></code></pre></div><p>Returns the quadratic form <code>x&#39; * P * x</code>.</p><p><strong>Arguments</strong></p><div class="language- vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang"></span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span>- `x::AbstractVector`: A vector.</span></span>\n<span class="line"><span>- `P::AbstractMatrix`: A matrix.</span></span></code></pre></div><p><a href="https://github.com/Vaibhavdixit02/SymbolicAnalysis.jl" target="_blank" rel="noreferrer">source</a></p></div><br><div style="border-width:1px;border-style:solid;border-color:black;padding:1em;border-radius:25px;"><a id="SymbolicAnalysis.quad_over_lin-Tuple{Real, Real}" href="#SymbolicAnalysis.quad_over_lin-Tuple{Real, Real}">#</a> <b><u>SymbolicAnalysis.quad_over_lin</u></b> — <i>Method</i>. <div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">quad_over_lin</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(x, y</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">::</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">Number</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">)</span></span></code></pre></div><p>Returns the quadratic over linear form <code>x^2 / y</code>.</p><p><strong>Arguments</strong></p><div class="language- vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang"></span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span>- `x`: A number or a vector.</span></span>\n<span class="line"><span>- `y::Number`: A positive number.</span></span></code></pre></div><p><a href="https://github.com/Vaibhavdixit02/SymbolicAnalysis.jl" target="_blank" rel="noreferrer">source</a></p></div><br><div style="border-width:1px;border-style:solid;border-color:black;padding:1em;border-radius:25px;"><a id="SymbolicAnalysis.scalar_mat" href="#SymbolicAnalysis.scalar_mat">#</a> <b><u>SymbolicAnalysis.scalar_mat</u></b> — <i>Function</i>. <div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">scalar_mat</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(X, k</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">size</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(X, </span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">1</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">))</span></span></code></pre></div><p>Scalar matrix of a symmetric positive definite matrix <code>X</code> is defined as <code>tr(X)*I(k)</code>.</p><p><strong>Arguments</strong></p><div class="language- vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang"></span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span>- `X::Matrix`: A symmetric positive definite matrix.</span></span>\n<span class="line"><span>- `k::Int`: The size of the identity matrix.</span></span></code></pre></div><p><a href="https://github.com/Vaibhavdixit02/SymbolicAnalysis.jl" target="_blank" rel="noreferrer">source</a></p></div><br><div style="border-width:1px;border-style:solid;border-color:black;padding:1em;border-radius:25px;"><a id="SymbolicAnalysis.sdivergence-Tuple{Any, Any}" href="#SymbolicAnalysis.sdivergence-Tuple{Any, Any}">#</a> <b><u>SymbolicAnalysis.sdivergence</u></b> — <i>Method</i>. <div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">sdivergence</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(X, Y)</span></span></code></pre></div><p>Symmetric divergence of two symmetric positive definite matrices <code>X</code> and <code>Y</code> is defined as <code>logdet((X+Y)/2) - 1/2*logdet(X*Y)</code>.</p><p><strong>Arguments</strong></p><div class="language- vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang"></span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span>- `X::Matrix`: A symmetric positive definite matrix.</span></span>\n<span class="line"><span>- `Y::Matrix`: A symmetric positive definite matrix.</span></span></code></pre></div><p><a href="https://github.com/Vaibhavdixit02/SymbolicAnalysis.jl" target="_blank" rel="noreferrer">source</a></p></div><br><div style="border-width:1px;border-style:solid;border-color:black;padding:1em;border-radius:25px;"><a id="SymbolicAnalysis.sum_largest-Tuple{AbstractMatrix, Integer}" href="#SymbolicAnalysis.sum_largest-Tuple{AbstractMatrix, Integer}">#</a> <b><u>SymbolicAnalysis.sum_largest</u></b> — <i>Method</i>. <div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">sum_largest</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(x</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">::</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">AbstractMatrix</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, k)</span></span></code></pre></div><p>Returns the sum of the <code>k</code> largest elements of <code>x</code>.</p><p><strong>Arguments</strong></p><div class="language- vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang"></span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span>- `x::AbstractMatrix`: A matrix.</span></span>\n<span class="line"><span>- `k::Int`: The number of largest elements to sum.</span></span></code></pre></div><p><a href="https://github.com/Vaibhavdixit02/SymbolicAnalysis.jl" target="_blank" rel="noreferrer">source</a></p></div><br><div style="border-width:1px;border-style:solid;border-color:black;padding:1em;border-radius:25px;"><a id="SymbolicAnalysis.sum_smallest-Tuple{AbstractMatrix, Integer}" href="#SymbolicAnalysis.sum_smallest-Tuple{AbstractMatrix, Integer}">#</a> <b><u>SymbolicAnalysis.sum_smallest</u></b> — <i>Method</i>. <div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">sum_smallest</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(x</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">::</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">AbstractMatrix</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, k)</span></span></code></pre></div><p>Returns the sum of the <code>k</code> smallest elements of <code>x</code>.</p><p><strong>Arguments</strong></p><div class="language- vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang"></span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span>- `x::AbstractMatrix`: A matrix.</span></span>\n<span class="line"><span>- `k::Int`: The number of smallest elements to sum.</span></span></code></pre></div><p><a href="https://github.com/Vaibhavdixit02/SymbolicAnalysis.jl" target="_blank" rel="noreferrer">source</a></p></div><br><div style="border-width:1px;border-style:solid;border-color:black;padding:1em;border-radius:25px;"><a id="SymbolicAnalysis.trinv-Tuple{AbstractMatrix}" href="#SymbolicAnalysis.trinv-Tuple{AbstractMatrix}">#</a> <b><u>SymbolicAnalysis.trinv</u></b> — <i>Method</i>. <div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">trinv</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(x</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">::</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">AbstractMatrix</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">)</span></span></code></pre></div><p>Returns the trace of the inverse of <code>x</code>.</p><p><strong>Arguments</strong></p><div class="language- vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang"></span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span>- `x::AbstractMatrix`: A matrix.</span></span></code></pre></div><p><a href="https://github.com/Vaibhavdixit02/SymbolicAnalysis.jl" target="_blank" rel="noreferrer">source</a></p></div><br><div style="border-width:1px;border-style:solid;border-color:black;padding:1em;border-radius:25px;"><a id="SymbolicAnalysis.tv-Tuple{AbstractVector{&lt;:AbstractMatrix}}" href="#SymbolicAnalysis.tv-Tuple{AbstractVector{&lt;:AbstractMatrix}}">#</a> <b><u>SymbolicAnalysis.tv</u></b> — <i>Method</i>. <div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">tv</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(x</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">::</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">AbstractVector{&lt;:AbstractMatrix}</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">)</span></span></code></pre></div><p>Returns the total variation of <code>x</code>, defined as <code>sum_{i,j} |x_{k+1}[i,j] - x_k[i,j]|</code>.</p><p><strong>Arguments</strong></p><div class="language- vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang"></span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span>- `x::AbstractVector`: A vector of matrices.</span></span></code></pre></div><p><a href="https://github.com/Vaibhavdixit02/SymbolicAnalysis.jl" target="_blank" rel="noreferrer">source</a></p></div><br><div style="border-width:1px;border-style:solid;border-color:black;padding:1em;border-radius:25px;"><a id="SymbolicAnalysis.tv-Tuple{AbstractVector{&lt;:Number}}" href="#SymbolicAnalysis.tv-Tuple{AbstractVector{&lt;:Number}}">#</a> <b><u>SymbolicAnalysis.tv</u></b> — <i>Method</i>. <div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">tv</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(x</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">::</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">AbstractVector{&lt;:Number}</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">)</span></span></code></pre></div><p>Returns the total variation of <code>x</code>, defined as <code>sum_i |x_{i+1} - x_i|</code>.</p><p><strong>Arguments</strong></p><div class="language- vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang"></span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span>- `x::AbstractVector`: A vector.</span></span></code></pre></div><p><a href="https://github.com/Vaibhavdixit02/SymbolicAnalysis.jl" target="_blank" rel="noreferrer">source</a></p></div><br>',38),l=[n];function p(r,o,d,c,h,b){return a(),i("div",null,l)}const k=s(t,[["render",p]]);export{g as __pageData,k as default};
