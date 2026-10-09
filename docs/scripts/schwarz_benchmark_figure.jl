# Figure for the Schwarz benchmarks: setup + solve time against problem size on
# the unit square, one line per method. Reads the CSV written by
# schwarz_benchmarks.jl, so the figure can be redrawn without rerunning them.
#
# Run with:  julia --project=docs docs/scripts/schwarz_benchmark_figure.jl
# Writes docs/figures/schwarz/timing.svg and timing.png.

get!(ENV, "GKSwstype", "100")

using Plots

const DIR = joinpath(@__DIR__, "..", "figures", "schwarz")

# One benchmark result, as read back from benchmarks.csv.
struct TimingPoint
    case::String
    dofs::Int
    method::String
    total_s::Float64
    converged::Bool
end

function read_points(path)
    map(readlines(path)[2:end]) do line
        c = match(r"^\"([^\"]*)\",(\d+),(\d+),\"([^\"]*)\",([^,]+),([^,]+),(-?\d+),", line)
        TimingPoint(c[1], parse(Int, c[2]), c[4], parse(Float64, c[5]) + parse(Float64, c[6]),
            parse(Int, c[7]) >= 0)
    end
end

# Categorical slots of the reference palette, in fixed order, one per method,
# with a distinct marker so identity never rests on colour alone. CHOLMOD is
# left out (it tracks the CliqueTrees direct solve within a factor of two and
# stays in the CSV) to keep eight series.
const SERIES = [
    ("direct (ChordalLDLt)", "#2a78d6", :circle),
    ("Galerkin tower alone (one coarse solve)", "#eb6834", :diamond),
    ("Schwarz, one-level (multicolor)", "#1baf7a", :utriangle),
    ("Schwarz, two-level: tower + sweeps", "#eda100", :dtriangle),
    ("Robin Schwarz p* (multicolor)", "#e87ba4", :rect),
    ("Robin Schwarz p* (parallel)", "#008300", :star5),
    ("sheaf ADMM, ρ = p*", "#4a3aa7", :hexagon),
    ("two-level Schwarz CG", "#e34948", :pentagon),
]

points = filter(p -> startswith(p.case, "square") && p.converged, read_points(joinpath(DIR, "benchmarks.csv")))

plt = plot(; xscale=:log10, yscale=:log10,
    xlabel="unknowns", ylabel="setup + solve time (s)",
    title="Poisson on the unit square: 32×32-point subdomains, overlap 2, 8 threads",
    titlefontsize=11, guidefontsize=10, tickfontsize=9, legendfontsize=9,
    legend=:outerright, size=(1150, 540), dpi=150,
    background_color="#fcfcfb", foreground_color_text="#0b0b0b",
    foreground_color_axis="#52514e", foreground_color_border="#52514e",
    gridcolor="#52514e", gridalpha=0.15, gridlinewidth=0.5,
    left_margin=6Plots.mm, bottom_margin=6Plots.mm)
for (method, color, shape) in SERIES
    selected = sort(filter(p -> p.method == method, points); by=p -> p.dofs)
    isempty(selected) && continue
    plot!(plt, [p.dofs for p in selected], [p.total_s for p in selected];
        label=startswith(method, "Galerkin tower alone") ? "Galerkin tower alone (not a solver: 92–94% error)" : method,
        color=color, lw=2, marker=shape, markersize=6,
        markercolor=color, markerstrokecolor="#fcfcfb", markerstrokewidth=1.5)
end
savefig(plt, joinpath(DIR, "timing.svg"))
savefig(plt, joinpath(DIR, "timing.png"))
println("wrote ", joinpath(DIR, "timing.svg"), " and timing.png")
