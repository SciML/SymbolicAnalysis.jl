import{_ as t,c as e,o as a,a6 as o}from"./chunks/framework.D_dyxcd-.js";const u=JSON.parse('{"title":"SymbolicAnalysis.jl","description":"","frontmatter":{},"headers":[],"relativePath":"index.md","filePath":"index.md","lastUpdated":null}'),i={name:"index.md"},r=o('<h1 id="symbolicanalysis-jl" tabindex="-1">SymbolicAnalysis.jl <a class="header-anchor" href="#symbolicanalysis-jl" aria-label="Permalink to &quot;SymbolicAnalysis.jl&quot;">​</a></h1><p>Symbolics-based function property propagation for optimization</p><p>SymbolicAnalysis is a package for implementating the Disciplined Programming approach to optimization, As demonstrated by the <a href="https://dcp.stanford.edu/" target="_blank" rel="noreferrer">DCP framework</a>, and further followups to it for further classes of functions <a href="https://www.cvxpy.org/tutorial/index.html" target="_blank" rel="noreferrer">https://www.cvxpy.org/tutorial/index.html</a> such as DGP, DQP etc, symbolic representation of problems can be leveraged to identify and facilitate building convex (or similar function properties) expressions.</p><p>This package aims to utlize expression graph rewriting and metadata propagation supported by Symbolics.jl, to support propagation of several of these properties - limited right now to Euclidean Convexity and Geodesic Convexity on the Symmetric Positive Definite manifold. This package provides an easy to expand implementation of &quot;atoms&quot; that are functions that have kow properties. This allows users to add atoms to the library more easily than the previous implementations <a href="https://www.cvxpy.org/index.html" target="_blank" rel="noreferrer">CVXPY</a> and <a href="https://github.com/jump-dev/Convex.jl" target="_blank" rel="noreferrer">Convex.jl</a> as well as a more performant implementation of the function property propagation.</p><p>Manifold optimization is supported through the Manopt.jl wrapper of <a href="https://github.com/SciML/Optimization.jl" target="_blank" rel="noreferrer">Optimization.jl</a> that allows utilizing the Geodesic Convexity propagation to be utilized as a certificate of global optimization with riemannian optimization.</p>',5),n=[r];function s(l,p,c,m,d,h){return a(),e("div",null,n)}const _=t(i,[["render",s]]);export{u as __pageData,_ as default};
