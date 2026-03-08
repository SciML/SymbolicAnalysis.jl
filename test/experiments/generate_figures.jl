"""
Generate publication-quality figures from experiment CSV results.

Usage:
    julia --project=test test/experiments/generate_figures.jl

Reads CSVs from test/experiments/results/ and produces PDF + PNG figures.
"""

using CairoMakie
using CSV
using DataFrames

const RESULTS_DIR = joinpath(@__DIR__, "results")

# --------------------------------------------------------------------------- #
# Theme setup
# --------------------------------------------------------------------------- #

# Okabe-Ito colorblind-safe palette
const OI_PALETTE = [
    colorant"#E69F00",
    colorant"#56B4E9",
    colorant"#009E73",
    colorant"#F0E442",
    colorant"#0072B2",
    colorant"#D55E00",
    colorant"#CC79A7",
]

function publication_theme()
    t = Theme(
        fontsize = 10,
        figure_padding = 8,
        Axis = (
            xgridvisible = false,
            ygridvisible = false,
            topspinevisible = false,
            rightspinevisible = false,
            xlabelsize = 11,
            ylabelsize = 11,
            titlesize = 12,
        ),
        Legend = (framevisible = false, labelsize = 9, patchsize = (15, 10)),
    )
    # Try to use a serif font; fall back silently if unavailable
    try
        t = merge(t, Theme(fonts = (; regular = "Times New Roman")))
    catch
    end
    return t
end

set_theme!(publication_theme())

# Helper: save both PDF and PNG
function save_figure(fig, name)
    save(joinpath(RESULTS_DIR, name * ".pdf"), fig)
    save(joinpath(RESULTS_DIR, name * ".png"), fig, px_per_unit = 300 / 72)
    println("  Saved $(name).pdf and $(name).png")
end

# --------------------------------------------------------------------------- #
# Figure 1: DCP vs DGCP Overhead (grouped bar)
# --------------------------------------------------------------------------- #

function figure_timing_overhead()
    df = CSV.read(joinpath(RESULTS_DIR, "timing_comparison.csv"), DataFrame)
    n = nrow(df)
    xs = 1:n

    fig = Figure(size = (504, 288))  # ~7x4 inches at 72 dpi
    ax = Axis(
        fig[1, 1],
        xlabel = "Function",
        ylabel = "Time (us)",
        title = "DCP vs DGCP Verification Time",
        xticks = (collect(xs), df.Function),
        xticklabelrotation = pi / 6,
    )

    w = 0.35
    barplot!(
        ax,
        collect(xs) .- w / 2,
        df.DCP_us;
        width = w,
        color = OI_PALETTE[1],
        label = "DCP",
    )
    barplot!(
        ax,
        collect(xs) .+ w / 2,
        df.DGCP_us;
        width = w,
        color = OI_PALETTE[2],
        label = "DGCP",
    )

    axislegend(ax; position = :lt)

    save_figure(fig, "fig1_timing_overhead")
    return fig
end

# --------------------------------------------------------------------------- #
# Figure 2: Scaling Analysis (2-panel)
# --------------------------------------------------------------------------- #

function figure_scaling()
    df = CSV.read(joinpath(RESULTS_DIR, "scaling_analysis.csv"), DataFrame)

    fig = Figure(size = (720, 288))  # ~10x4 inches

    # Panel (a): time vs Terms for Karcher, MatrixSize==5
    sub_terms = filter(r -> r.Problem == "Karcher" && r.MatrixSize == 5, df)
    ax1 = Axis(
        fig[1, 1],
        xlabel = "Number of terms",
        ylabel = "Time (us)",
        title = "(a) Karcher mean, n = 5",
    )
    scatterlines!(
        ax1,
        sub_terms.Terms,
        sub_terms.DCP_us;
        color = OI_PALETTE[1],
        marker = :circle,
        linewidth = 2,
        label = "DCP",
    )
    scatterlines!(
        ax1,
        sub_terms.Terms,
        sub_terms.DGCP_us;
        color = OI_PALETTE[2],
        marker = :rect,
        linewidth = 2,
        label = "DGCP",
    )
    axislegend(ax1; position = :lt)

    # Panel (b): time vs MatrixSize for Karcher, Terms==3
    sub_size = filter(r -> r.Problem == "Karcher" && r.Terms == 3, df)
    ax2 = Axis(
        fig[1, 2],
        xlabel = "Matrix size n",
        ylabel = "Time (us)",
        title = "(b) Karcher mean, 3 terms",
    )
    scatterlines!(
        ax2,
        sub_size.MatrixSize,
        sub_size.DCP_us;
        color = OI_PALETTE[1],
        marker = :circle,
        linewidth = 2,
        label = "DCP",
    )
    scatterlines!(
        ax2,
        sub_size.MatrixSize,
        sub_size.DGCP_us;
        color = OI_PALETTE[2],
        marker = :rect,
        linewidth = 2,
        label = "DGCP",
    )
    axislegend(ax2; position = :lt)

    save_figure(fig, "fig2_scaling")
    return fig
end

# --------------------------------------------------------------------------- #
# Figure 3: Benchmark Complexity (time vs size)
# --------------------------------------------------------------------------- #

function figure_benchmark()
    df = CSV.read(joinpath(RESULTS_DIR, "extended_benchmark.csv"), DataFrame)

    fig = Figure(size = (504, 288))
    ax = Axis(
        fig[1, 1],
        xlabel = "Matrix size n",
        ylabel = "Time (ms)",
        title = "Verification Time vs Problem Size",
    )

    problems = unique(df.Problem)
    for (i, ptype) in enumerate(problems)
        sub = filter(r -> r.Problem == ptype, df)
        ci = mod1(i, length(OI_PALETTE))
        scatterlines!(
            ax,
            sub.Size,
            sub.Time_ms;
            color = OI_PALETTE[ci],
            marker = :circle,
            linewidth = 2,
            label = ptype,
        )
    end
    axislegend(ax; position = :lt)

    save_figure(fig, "fig3_benchmark")
    return fig
end

# --------------------------------------------------------------------------- #
# Figure 4: Expert Verification (horizontal bars)
# --------------------------------------------------------------------------- #

function figure_expert()
    df = CSV.read(joinpath(RESULTS_DIR, "expert_examples.csv"), DataFrame)
    n = nrow(df)
    ys = 1:n

    colors = [d == "Hard" ? OI_PALETTE[5] : OI_PALETTE[1] for d in df.Difficulty]

    fig = Figure(size = (504, 288))
    ax = Axis(
        fig[1, 1],
        ylabel = "",
        xlabel = "Time (ms)",
        title = "Expert-Level DGCP Verification Time",
        yticks = (collect(ys), df.Case),
    )

    barplot!(ax, collect(ys), df.Time_ms; direction = :x, color = colors)

    # Manual legend entries for difficulty
    elem_hard = PolyElement(color = OI_PALETTE[5])
    elem_med = PolyElement(color = OI_PALETTE[1])
    Legend(
        fig[1, 2],
        [elem_hard, elem_med],
        ["Hard", "Medium"];
        framevisible = false,
        labelsize = 9,
    )

    save_figure(fig, "fig4_expert")
    return fig
end

# --------------------------------------------------------------------------- #
# Figure 5: MLE Comparison (grouped bars)
# --------------------------------------------------------------------------- #

function figure_mle()
    df = CSV.read(joinpath(RESULTS_DIR, "mle_experiment.csv"), DataFrame)
    n = nrow(df)
    xs = 1:n

    labels = df.Problem .* " n=" .* string.(df.n) .* " k=" .* string.(df.Samples)
    dgcp_ms = df.DGCP_s .* 1000
    dcp_ms = df.DCP_s .* 1000

    fig = Figure(size = (576, 288))
    ax = Axis(
        fig[1, 1],
        xlabel = "",
        ylabel = "Time (ms)",
        title = "MLE Verification Time",
        xticks = (collect(xs), labels),
        xticklabelrotation = pi / 4,
    )

    w = 0.35
    barplot!(
        ax,
        collect(xs) .- w / 2,
        dgcp_ms;
        width = w,
        color = OI_PALETTE[1],
        label = "DGCP",
    )
    barplot!(
        ax,
        collect(xs) .+ w / 2,
        dcp_ms;
        width = w,
        color = OI_PALETTE[2],
        label = "DCP",
    )

    axislegend(ax; position = :lt)

    save_figure(fig, "fig5_mle")
    return fig
end

# --------------------------------------------------------------------------- #
# Main
# --------------------------------------------------------------------------- #

println("Generating publication figures from CSVs in $RESULTS_DIR ...")
println()

figure_timing_overhead()
figure_scaling()
figure_benchmark()
figure_expert()
figure_mle()

println()
println("All figures generated.")
