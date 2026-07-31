# Figures for the Tikhonov Harmonic Tracking feature guide.
# Every curve is computed from the shared harmonic-extension scenario.
#
#   julia --project=docs docs/scripts/tikhonov_figures.jl

# Render GR figures directly to files without opening a GKS window.
get!(ENV, "GKSwstype", "100")

using CellularSheaves
using CellularSheaves.ControlSheaves
using LinearAlgebra
using Plots
using Printf
using Statistics
using TOML

include(joinpath(@__DIR__, "tikhonov_scenario.jl"))
using .TikhonovGuideScenario

const FIG = normpath(joinpath(@__DIR__, "..", "figures", "tikhonov"))
mkpath(FIG)

# Match the plotting language used by the distributed harmonic tracking guide.
default(
    thickness_scaling = 1.35,
    framestyle = :box,
    grid = true,
    gridalpha = 0.18,
    gridstyle = :dot,
    titlefontsize = 10,
    guidefontsize = 9,
    legendfontsize = 8,
    tickfontsize = 8,
    markerstrokewidth = 0,
    size = (720, 340),
)

const PAL = [:steelblue, :darkorange, :green, :crimson]
const INK = :gray20
const TIMES = collect(0.0:0.02:4TARGET_PERIOD)
const EPSILON = 0.2
const BASELINE = planner_rollout(TIMES; epsilon = EPSILON)
const CORRECTED = planner_rollout(TIMES; epsilon = EPSILON, feedforward = true)

saveguide(plt, name) = savefig(plt, joinpath(FIG, name * ".svg"))

function manifold_figure()
    times = collect(0.0:0.02:TARGET_PERIOD)
    ideal = hcat(harmonic_reference.(times)...)
    coordinate = 1
    offsets = ([1.1, 0.8], [-0.9, 0.7], [0.9, -0.8], [-1.0, -0.75])

    plt = plot(
        times, ideal[coordinate, :];
        color = INK,
        linewidth = 3.5,
        label = "moving harmonic reference",
        xlabel = "time (s)",
        ylabel = "one planner coordinate",
        legend = :topright,
        size = (760, 370),
    )

    for (j, offset) in enumerate(offsets)
        filter = TikhonovFilter(ideal[:, 1] + repeat(offset, NA); epsilon = 0.24)
        state = similar(ideal)
        state[:, 1] = filter.x
        for k in 1:(length(times) - 1)
            tikhonov_step!(filter, ideal[:, k], ideal[:, k + 1], times[k + 1] - times[k])
            state[:, k + 1] = filter.x
        end
        plot!(
            plt,
            times, state[coordinate, :];
            color = PAL[j],
            linewidth = 2,
            label = j == 1 ? "filtered planner state" : "",
        )
        scatter!(plt, times[1:1], state[coordinate, 1:1];
            color = PAL[j], marker = :circle, markersize = 4, label = "")
    end

    saveguide(plt, "fig_manifold")
end

function boundary_layer_figure()
    epsilons = (0.4, 0.2, 0.1)
    times = collect(0.0:0.01:2.2)
    reference = harmonic_reference(0.0)
    initial = reference + repeat([0.65, -0.35], NA)
    initial_error = norm(initial - reference)

    plt = plot(
        xlabel = "time (s)",
        ylabel = "normalized planner error",
        yscale = :log10,
        ylims = (1e-9, 1.2),
        legend = :outerright,
        size = (850, 370),
    )

    for (j, epsilon) in enumerate(epsilons)
        filter = TikhonovFilter(initial; epsilon)
        measured = zeros(length(times))
        for k in eachindex(times)
            measured[k] = max(norm(filter.x - reference) / initial_error, 1e-12)
            k < length(times) && tikhonov_step!(
                filter, reference, reference, times[k + 1] - times[k]
            )
        end
        plot!(plt, times, measured; color = PAL[j], linewidth = 2.5,
            label = "epsilon = $epsilon")
        plot!(plt, times, exp.(-times ./ epsilon); color = PAL[j], linewidth = 1.3,
            linestyle = :dash, alpha = 0.75,
            label = j == 1 ? "exp(-t/epsilon) prediction" : "")
    end

    saveguide(plt, "fig_boundary_layer")
end

function lag_geometry_figure()
    ids = findall(t -> t <= TARGET_PERIOD, TIMES)
    times = TIMES[ids]
    coordinate = 1

    plt = plot(
        times, BASELINE.ideal[coordinate, ids];
        color = INK,
        linewidth = 3.5,
        label = "harmonic reference",
        xlabel = "time (s)",
        ylabel = "one planner coordinate",
        legend = :topright,
        size = (760, 370),
    )
    plot!(plt, times, BASELINE.state[coordinate, ids];
        color = :crimson, linewidth = 2.4, linestyle = :dash,
        label = "uncompensated")
    plot!(plt, times, CORRECTED.state[coordinate, ids];
        color = :green, linewidth = 2.2, label = "analytic feedforward")

    saveguide(plt, "fig_lag_geometry")
end

function epsilon_scaling_figure()
    epsilons = [0.025, 0.04, 0.065, 0.1, 0.16, 0.25, 0.4, 0.65]
    times = collect(0.0:0.02:3TARGET_PERIOD)
    tail = times .>= TARGET_PERIOD
    lag = Float64[]
    feedforward = Float64[]
    for epsilon in epsilons
        uncompensated = planner_rollout(times; epsilon)
        corrected = planner_rollout(times; epsilon, feedforward = true)
        push!(lag, sqrt(mean(abs2, uncompensated.error[tail])))
        push!(feedforward, sqrt(mean(abs2, corrected.error[tail])))
    end

    logepsilon = log10.(epsilons)
    loglag = log10.(lag)
    slope = sum((logepsilon .- mean(logepsilon)) .* (loglag .- mean(loglag))) /
        sum(abs2, logepsilon .- mean(logepsilon))
    first_order = lag[1] .* (epsilons ./ epsilons[1])

    plt = plot(
        epsilons, lag;
        color = :crimson,
        marker = :circle,
        linewidth = 2.5,
        xscale = :log10,
        yscale = :log10,
        xlabel = "epsilon",
        ylabel = "tail RMS planner error",
        label = "uncompensated",
        legend = :topleft,
        size = (760, 370),
    )
    plot!(plt, epsilons, feedforward; color = :green, marker = :square,
        linewidth = 2.5, label = "analytic feedforward")
    plot!(plt, epsilons, first_order; color = INK, linewidth = 1.4,
        linestyle = :dash, label = "first-order guide")
    saveguide(plt, "fig_epsilon_scaling")
    return epsilons, lag, feedforward, slope
end

manifold_figure()
boundary_layer_figure()
lag_geometry_figure()
epsilons, lag, feedforward, slope = epsilon_scaling_figure()

tail = TIMES .>= 2TARGET_PERIOD
metrics = Dict(
    "epsilon" => EPSILON,
    "epsilon_sweep" => epsilons,
    "lag_rms" => lag,
    "feedforward_rms" => feedforward,
    "lag_loglog_slope" => slope,
    "baseline_tail_rms" => sqrt(mean(abs2, BASELINE.error[tail])),
    "feedforward_tail_rms" => sqrt(mean(abs2, CORRECTED.error[tail])),
)
open(joinpath(FIG, "metrics.toml"), "w") do io
    TOML.print(io, metrics)
end
@printf(
    "lag slope %.3f; planner RMS %.4e -> %.4e\n",
    slope,
    metrics["baseline_tail_rms"],
    metrics["feedforward_tail_rms"],
)
