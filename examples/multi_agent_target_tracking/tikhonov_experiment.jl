# End-to-end Tikhonov filtering experiment on the existing 13-agent,
# four-target tracking scenario.
#
#   julia --project=examples/multi_agent_target_tracking \
#       examples/multi_agent_target_tracking/tikhonov_experiment.jl

get!(ENV, "GKSwstype", "100")

using CellularSheaves
using CellularSheaves.ControlSheaves.Tikhonov
using JSON3
using LinearAlgebra
using Plots
using Printf
using Statistics

import CellularSheaves.NetworkSheaves.EuclideanSheaves:
    _harmonic_extension_restricted_laplacian

include(joinpath(@__DIR__, "src", "core", "cellular_sheaf.jl"))
using .SheafConsensus

const CONFIG_PATH = joinpath(@__DIR__, "configurations", "config_common.json")
const FIGURE_DIR = normpath(joinpath(@__DIR__, "..", "..", "docs", "figures", "tikhonov"))
const DT = 0.05
const HORIZON = 300
const EPSILON = 0.2
const LOCAL_GAIN = 4.0

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
)

const COLORS = (:steelblue, :darkorange, :green, :crimson)
const INK = :gray20

"Velocity used by the existing figure-eight target scenario."
function figure_eight_velocity(t)
    return [-20sin(t), 20cos(2t), 15cos(t)]
end

function target_offsets(initial)
    d = 15.0
    h = sqrt(6) * d / 3
    desired = hcat(
        initial[:, 1],
        initial[:, 1] + [d, 0, 0],
        initial[:, 1] + [d / 2, sqrt(3) * d / 2, 0],
        initial[:, 1] + [d / 2, sqrt(3) * d / 6, h],
    )
    offsets = zeros(size(initial))
    for i in axes(initial, 2), j in (i + 1):size(initial, 2)
        delta = desired[:, j] - desired[:, i]
        offsets[:, i] .+= delta ./ size(initial, 2)
        offsets[:, j] .-= delta ./ size(initial, 2)
    end
    return offsets
end

function load_scenario()
    cfg = Dict{String, Any}(JSON3.read(read(CONFIG_PATH, String), Dict{String, Any}))
    agents = hcat([Float64.(collect(a["initial_position"])) for a in cfg["agents"]]...)
    targets = hcat([Float64.(collect(t["initial_position"])) for t in cfg["targets"]]...)
    agent_edges = [(Int(e[1]), Int(e[2])) for e in cfg["agent_edge_set"]]
    target_edges = [(Int(e[1]), Int(e[2])) for e in cfg["target_edge_set"]]
    pinning = Float64.(hcat([Float64.(collect(row)) for row in cfg["pinning_matrix"]]...)')
    return (; cfg, agents, targets, agent_edges, target_edges, pinning)
end

function target_rollout(scenario)
    ntargets = size(scenario.targets, 2)
    target_sheaf = SheafConsensus.build_target_sheaf(ntargets, 3, scenario.target_edges)
    coboundary = coboundary_map(target_sheaf)
    offsets = target_offsets(scenario.targets)
    gain = Float64(scenario.cfg["targets_proportional_gain"])
    trajectory = zeros(3, HORIZON, ntargets)
    trajectory[:, 1, :] = scenario.targets

    for k in 1:(HORIZON - 1)
        state = reshape(trajectory[:, k, :] .+ offsets, :)
        disagreement = coboundary' * (coboundary * state)
        velocity = figure_eight_velocity((k - 1) * DT)
        for target in 1:ntargets
            block = (3target - 2):(3target)
            control = -gain .* disagreement[block] ./ 100 .+ velocity
            trajectory[:, k + 1, target] = trajectory[:, k, target] + DT .* control
        end
    end
    return trajectory
end

function harmonic_references(scenario, targets)
    nagents = size(scenario.agents, 2)
    ntargets = size(targets, 3)
    sheaf, _, target_vertices = SheafConsensus.build_agent_sheaf(
        nagents, ntargets, 3, scenario.agent_edges, scenario.pinning
    )
    boundary = Dict(target_vertices[j] => targets[:, 1, j] for j in 1:ntargets)
    _, _, H, B = _harmonic_extension_restricted_laplacian(sheaf, boundary)
    factor = factorize(H)
    references = zeros(3, HORIZON, nagents)
    for k in 1:HORIZON
        references[:, k, :] = reshape(factor \ (-B * vec(targets[:, k, :])), 3, nagents)
    end
    return references
end

function filtered_references(ideal; feedforward=false)
    nagents = size(ideal, 3)
    flat = reshape(permutedims(ideal, (1, 3, 2)), 3nagents, HORIZON)
    filter = TikhonovFilter(flat[:, 1]; epsilon=EPSILON)
    output = similar(flat)
    output[:, 1] = filter.x
    for k in 1:(HORIZON - 1)
        q0, q1 = flat[:, k], flat[:, k + 1]
        if feedforward
            qdot = (q1 - q0) / DT
            u0 = tikhonov_feedforward_reference(q0, qdot, EPSILON)
            u1 = tikhonov_feedforward_reference(q1, qdot, EPSILON)
            tikhonov_step!(filter, u0, u1, DT)
        else
            tikhonov_step!(filter, q0, q1, DT)
        end
        output[:, k + 1] = filter.x
    end
    return permutedims(reshape(output, 3, nagents, HORIZON), (1, 3, 2))
end

function physical_rollout(initial, reference, ideal)
    positions = zeros(size(reference))
    positions[:, 1, :] = initial
    for k in 1:(HORIZON - 1)
        positions[:, k + 1, :] = positions[:, k, :] +
            DT * LOCAL_GAIN .* (reference[:, k, :] - positions[:, k, :])
    end
    planner = vec(sqrt.(mean(abs2, reference - ideal; dims=(1, 3))))
    local_error = vec(sqrt.(mean(abs2, positions - reference; dims=(1, 3))))
    formation = vec(sqrt.(mean(abs2, positions - ideal; dims=(1, 3))))
    return (; positions, references=reference, planner, local_error, formation)
end

scenario = load_scenario()
targets = target_rollout(scenario)
ideal = harmonic_references(scenario, targets)
references = Dict(
    :direct => ideal,
    :lagged => filtered_references(ideal),
    :feedforward => filtered_references(ideal; feedforward=true),
)
rollouts = Dict(name => physical_rollout(scenario.agents, reference, ideal)
                for (name, reference) in references)

tail = (HORIZON ÷ 2 + 1):HORIZON
println("Tikhonov cascade: existing 13-agent/four-target scenario, epsilon=$EPSILON")
@printf("%-14s %13s %13s %16s\n", "planner", "planner RMS", "local RMS", "formation RMS")
for name in (:direct, :lagged, :feedforward)
    result = rollouts[name]
    @printf("%-14s %13.5e %13.5e %16.5e\n", String(name),
        sqrt(mean(abs2, result.planner[tail])),
        sqrt(mean(abs2, result.local_error[tail])),
        sqrt(mean(abs2, result.formation[tail])))
end

mkpath(FIGURE_DIR)
saveguide(plt, name) = savefig(plt, joinpath(FIGURE_DIR, name))

assignments = [begin
    pinned = findfirst(!iszero, scenario.pinning[agent, :])
    pinned === nothing ? argmin([
        norm(ideal[:, 1, agent] - targets[:, 1, target])
        for target in axes(targets, 3)
    ]) : pinned
end for agent in axes(ideal, 3)]

all_x = vcat(vec(targets[1, :, :]), vec(ideal[1, :, :]),
    vec(rollouts[:lagged].positions[1, :, :]), vec(rollouts[:feedforward].positions[1, :, :]))
all_y = vcat(vec(targets[2, :, :]), vec(ideal[2, :, :]),
    vec(rollouts[:lagged].positions[2, :, :]), vec(rollouts[:feedforward].positions[2, :, :]))
x_min, x_max = extrema(all_x)
y_min, y_max = extrema(all_y)
padding = 0.06 * max(x_max - x_min, y_max - y_min)

function tracking_frame(name, title, k)
    result = rollouts[name]
    first_trail = max(1, k - 25)
    plt = plot(
        aspect_ratio=:equal,
        xlims=(x_min - padding, x_max + padding),
        ylims=(y_min - padding, y_max + padding),
        xlabel="horizontal position",
        ylabel="vertical position",
        title=@sprintf("%s   t = %.1f s", title, (k - 1) * DT),
        legend=false,
    )
    for target in axes(targets, 3)
        color = COLORS[target]
        plot!(plt, targets[1, :, target], targets[2, :, target];
            color, linewidth=1.5, linestyle=:dash, alpha=0.25, label="")
        scatter!(plt, targets[1, k:k, target], targets[2, k:k, target];
            color, marker=:star5, markersize=7, label="")
    end
    for agent in axes(ideal, 3)
        color = COLORS[assignments[agent]]
        plot!(plt,
            [result.positions[1, k, agent], ideal[1, k, agent]],
            [result.positions[2, k, agent], ideal[2, k, agent]];
            color, linewidth=0.9, alpha=0.35, label="")
        plot!(plt,
            result.positions[1, first_trail:k, agent],
            result.positions[2, first_trail:k, agent];
            color, linewidth=1.8, alpha=0.75, label="")
        scatter!(plt, ideal[1, k:k, agent], ideal[2, k:k, agent];
            color, marker=:circle, markersize=4.5, markercolor=:white,
            markerstrokecolor=color, markerstrokewidth=1.2, label="")
        scatter!(plt, result.positions[1, k:k, agent], result.positions[2, k:k, agent];
            color, marker=:circle, markersize=3.2, label="")
    end
    return plt
end

animation = @animate for k in 1:5:HORIZON
    plot(
        tracking_frame(:lagged, "Uncompensated", k),
        tracking_frame(:feedforward, "Feedforward", k);
        layout=(1, 2), size=(1000, 430),
        left_margin=4Plots.mm, bottom_margin=4Plots.mm,
    )
end
gif(animation, joinpath(FIGURE_DIR, "fig_deployed_tracking.gif"); fps=12, show_msg=false)

times = (0:(HORIZON - 1)) .* DT
p_planner = plot(times, rollouts[:lagged].planner;
    color=:crimson, linewidth=2.5, label="uncompensated",
    xlabel="time (s)", ylabel="planner RMS", legend=:topright)
plot!(p_planner, times, rollouts[:feedforward].planner;
    color=:green, linewidth=2.5, label="analytic feedforward")

p_formation = plot(times, rollouts[:direct].formation;
    color=INK, linewidth=2.3, label="direct reference",
    xlabel="time (s)", ylabel="formation RMS", legend=:topright)
plot!(p_formation, times, rollouts[:lagged].formation;
    color=:crimson, linewidth=2.3, label="uncompensated")
plot!(p_formation, times, rollouts[:feedforward].formation;
    color=:green, linewidth=2.3, label="analytic feedforward")

errors = plot(p_planner, p_formation; layout=(1, 2), size=(1000, 380),
    left_margin=4Plots.mm, bottom_margin=4Plots.mm)
saveguide(errors, "fig_deployed_errors.svg")
