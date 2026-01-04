using Documenter, DocumenterVitepress

using SymbolicAnalysis

makedocs(;
    modules = [SymbolicAnalysis],
    authors = "Vaibhav Dixit, Shashi Gowda",
    repo = "https://github.com/SciML/SymbolicAnalysis.jl",
    sitename = "SymbolicAnalysis.jl",
    format = DocumenterVitepress.MarkdownVitepress(
        repo = "https://github.com/SciML/SymbolicAnalysis.jl",
        devurl = "dev"
    ),
    pages = [
        "Home" => "index.md",
        "Examples" => "examples.md",
        "Atoms" => "atoms.md",
        "Special Functions" => "functions.md",
    ],
    warnonly = true
)

deploydocs(; repo = "github.com/SciML/SymbolicAnalysis.jl", push_preview = true)
