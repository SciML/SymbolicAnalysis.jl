#=
Generate REPL-style listing images for Section 4.5 non-g-convex examples.
Produces listing/11.png, listing/12.png, listing/13.png
=#

using CairoMakie

const LISTING_DIR = get(
    ENV,
    "SYMBOLICANALYSIS_LISTING_DIR",
    joinpath(@__DIR__, "..", "..", "_MPC_v2__DGCP", "listing"),
)
mkpath(LISTING_DIR)

function make_listing_image(
        code_lines::Vector{String},
        output_lines::Vector{String},
        filename::String,
    )
    all_lines = vcat(code_lines, output_lines)
    n = length(all_lines)

    fig =
        Figure(size = (800, 30 + 22 * n), fontsize = 13, figure_padding = (15, 15, 10, 10))

    ax = Axis(fig[1, 1], limits = (0, 100, -n, 0.5), yreversed = false)
    hidedecorations!(ax)
    hidespines!(ax)

    for (i, line) in enumerate(all_lines)
        color = i <= length(code_lines) ? :black : RGBf(0.0, 0.5, 0.0)
        text!(
            ax,
            1,
            -(i - 1),
            text = line,
            fontsize = 13,
            font = "JuliaMono",
            color = color,
            align = (:left, :top),
        )
    end

    save(filename, fig, px_per_unit = 3)
    return println("Saved $filename")
end

# Listing 11: Square of logdet
make_listing_image(
    [
        "julia> @variables X[1:3, 1:3]",
        "       M = SymmetricPositiveDefinite(3)",
        "       result = analyze(logdet(X)^2, M)",
        "       println(result.gcurvature)",
    ],
    ["GUnknownCurvature"],
    joinpath(LISTING_DIR, "11.png"),
)

# Listing 12: sin of logdet (non-DCP atom)
make_listing_image(
    ["julia> result = analyze(sin(logdet(X)), M)", "       println(result.gcurvature)"],
    ["GUnknownCurvature"],
    joinpath(LISTING_DIR, "12.png"),
)

# Listing 13: sqrt of trace (concave of convex)
make_listing_image(
    ["julia> result = analyze(sqrt(tr(X)), M)", "       println(result.gcurvature)"],
    ["GUnknownCurvature"],
    joinpath(LISTING_DIR, "13.png"),
)
