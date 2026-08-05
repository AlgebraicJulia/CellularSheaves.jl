module CoordinationBenchmarkPlots

using Plots
using CellularSheaves
using Graphs
using Printf
using Random
using CellularSheaves.ControlSheaves.CoordinationBenchmarks:
    CoordinationBenchmark, SettlingResult, CoordinationScenario,
    CoordinationRollout, settle_to_formation, spectral_summary,
    target_trajectory, agent_groups, direct_plan, tree_statistics

const DIFFUSION_COLOR = :steelblue
const DIRECT_COLOR = :darkorange
const GUIDE_COLOR = :gray45

_family_marker(name) = get(Dict(
    :chain => :circle, :ring => :diamond, :grid => :square,
    :rgg => :utriangle, :star => :star5, :twoclique => :hexagon,
    :expander => :dtriangle), name, :circle)

"""
Plot a [`CoordinationBenchmark`](@ref) under one of several metrics.

    plot(result, :settling)       # control ticks to reach the formation vs kappa
    plot(result, :compute)        # total compute, both laws, vs kappa
    plot(result, :communication)  # total radio slots, both laws, vs kappa
    plot(result, :speedup)        # Direct advantage vs kappa, both currencies
    plot(result, :spectrum)       # lambda_min / lambda_max decomposition
"""
@recipe function f(result::CoordinationBenchmark, metric::Symbol = :settling)
    rows = result.rows
    kappa = [r.condition for r in rows]

    xscale --> :log10
    xlabel --> "κ(ℋ)"
    legend --> :topleft
    markerstrokewidth --> 0

    if metric === :settling
        yscale --> :log10
        ylabel --> "control ticks to reach the formation"
        title --> "Diffusion needs O(κ) ticks; Direct needs a constant few"

        @series begin
            seriestype := :scatter
            color := DIFFUSION_COLOR
            markersize := 5
            label := "Diffusion"
            kappa, [r.diffusion.ticks for r in rows]
        end
        @series begin
            seriestype := :scatter
            color := DIRECT_COLOR
            marker := :diamond
            markersize := 5
            label := "Direct"
            kappa, [r.direct.ticks for r in rows]
        end
        @series begin
            seriestype := :path
            color := GUIDE_COLOR
            linestyle := :dash
            label := "∝ κ"
            ks = sort(kappa)
            ks, 0.6 .* ks
        end

    elseif metric === :compute
        yscale --> :log10
        ylabel --> "seconds to reach the formation"
        title --> "Total compute to coordinate"
        @series begin
            seriestype := :scatter
            color := DIFFUSION_COLOR
            markersize := 5
            label := "Diffusion"
            kappa, [r.diffusion.seconds for r in rows]
        end
        @series begin
            seriestype := :scatter
            color := DIRECT_COLOR
            marker := :diamond
            markersize := 5
            label := "Direct"
            kappa, [r.direct.seconds for r in rows]
        end

    elseif metric === :communication
        yscale --> :log10
        ylabel --> "half-duplex slots to reach the formation"
        title --> "Diffusion diffuses every tick; Direct solves once"
        @series begin
            seriestype := :scatter
            color := DIFFUSION_COLOR
            markersize := 5
            label := "Diffusion  (ticks × Δ)"
            kappa, [max(r.diffusion.slots, 1) for r in rows]
        end
        @series begin
            seriestype := :scatter
            color := DIRECT_COLOR
            marker := :diamond
            markersize := 5
            label := "Direct  (2·cp(T), once)"
            kappa, [max(r.direct.slots, 1) for r in rows]
        end

    elseif metric === :speedup
        yscale --> :log10
        ylabel --> "Diffusion cost / Direct cost"
        title --> "Direct's advantage grows with conditioning"
        @series begin
            seriestype := :scatter
            color := DIRECT_COLOR
            markersize := 5
            label := "compute"
            kappa, [r.speedup for r in rows]
        end
        @series begin
            seriestype := :scatter
            color := :seagreen
            marker := :diamond
            markersize := 5
            label := "radio slots"
            kappa, [r.slot_ratio for r in rows]
        end
        @series begin
            seriestype := :hline
            color := GUIDE_COLOR
            linestyle := :dash
            label := "parity"
            [1.0]
        end

    elseif metric === :landscape
        ## the two-currency map: where each formation sits, and who wins there
        yscale := :log10
        xscale := :identity
        xlabel := "treewidth  (Direct's currency)"
        ylabel := "κ(ℋ)  (Diffusion's currency)"
        title --> "Where each formation sits, and who wins"
        legend --> :topright
        tw = [max(r.treewidth, 0) for r in rows]
        kp = [r.condition for r in rows]
        adv = [r.speedup for r in rows]
        won = adv .>= 1
        @series begin
            seriestype := :scatter
            color := DIRECT_COLOR
            markersize := [3 + 3log10(max(a, 1)) for a in adv[won]]
            alpha := 0.75
            label := "Direct faster (marker ∝ margin)"
            tw[won], kp[won]
        end
        if any(.!won)
            @series begin
                seriestype := :scatter
                color := DIFFUSION_COLOR
                marker := :diamond
                markersize := 6
                label := "Diffusion faster"
                tw[.!won], kp[.!won]
            end
        end

    elseif metric === :spectrum
        yscale := :log10
        ylabel := "eigenvalue of ℋ"
        xlabel := "agents"
        xscale := :log10
        title --> "λ_max tracks degree; λ_min tracks target coverage"
        agents = [r.agents for r in rows]
        @series begin
            seriestype := :scatter
            color := DIFFUSION_COLOR
            markersize := 5
            label := "λ_max"
            agents, [r.lambda_max for r in rows]
        end
        @series begin
            seriestype := :scatter
            color := DIRECT_COLOR
            marker := :diamond
            markersize := 5
            label := "λ_min"
            agents, [r.lambda_min for r in rows]
        end
    else
        error("unknown metric $metric; expected one of :settling, :compute, " *
              ":communication, :speedup, :landscape, :spectrum")
    end
end

"""
Plot the settling transient of one or more [`SettlingResult`](@ref)s.

    plot(ecc_result)
    plot!(acc_result)
"""
@recipe function f(res::SettlingResult)
    yscale --> :log10
    xlabel --> "control tick"
    ylabel --> "‖q − q*‖ / √N"
    label --> uppercase(String(res.method))
    linewidth --> 2
    color --> (res.method === :diffusion ? DIFFUSION_COLOR : DIRECT_COLOR)
    markerstrokewidth --> 0
    idx = 1:length(res.error_history)
    idx, max.(res.error_history, 1e-16)
end

"""
    plot_coordination(result; metrics, layout, size)

Multi-panel summary of a [`CoordinationBenchmark`](@ref).
"""
function CellularSheaves.ControlSheaves.CoordinationBenchmarks.plot_coordination(
        result::CoordinationBenchmark;
        metrics = (:settling, :communication),
        layout = (1, length(metrics)),
        size = (1050, 400),
        kwargs...)
    panels = [plot(result, m) for m in metrics]
    return plot(panels...; layout, size, left_margin = 5Plots.mm,
        bottom_margin = 6Plots.mm, kwargs...)
end


# ==========================
# scenario rollouts
# ==========================

const RING_PALETTE = [:steelblue, :darkorange, :green, :crimson]
const TARGET_A, TARGET_B = :gray35, :black
ring_color(g) = RING_PALETTE[mod1(g, length(RING_PALETTE))]
target_color(k) = k == 1 ? TARGET_A : TARGET_B

_ax(r::CoordinationRollout, i, k) = r.positions[(i - 1) * r.scenario.dim + 1, k]
_ay(r::CoordinationRollout, i, k) = r.positions[(i - 1) * r.scenario.dim + 2, k]
_rx(r::CoordinationRollout, i, k) = r.references[(i - 1) * r.scenario.dim + 1, k]
_ry(r::CoordinationRollout, i, k) = r.references[(i - 1) * r.scenario.dim + 2, k]

"""Symmetric square limits covering agents, references and targets."""
function _rollout_limits(rollouts...)
    m = 0.0
    for r in rollouts
        m = max(m, maximum(abs, r.positions), maximum(abs, r.references),
                maximum(abs, r.targets))
    end
    lim = 1.1m
    return (-lim, lim)
end

"""
One frame of a [`CoordinationRollout`](@ref), in the house style of the other
control examples: agents coloured by escort group with a trailing marker path,
faint dotted target orbits, and the targets as large stars.

    plot(rollout, k)
"""
@recipe function f(r::CoordinationRollout, k::Int = length(r.times);
        trail = 30, axis_lims = _rollout_limits(r),
        ring_of = agent_groups(r.scenario), show_reference = true)
    s = r.scenario
    aspect_ratio --> 1
    xlims --> axis_lims
    ylims --> axis_lims
    legend --> :outerright
    size --> (640, 480)
    title --> @sprintf("%s   t = %.2f", uppercasefirst(String(r.method)), r.times[k])
    markerstrokewidth --> 0

    for j in 1:s.ntargets
        @series begin
            seriestype := :path
            color := :gray80
            linewidth := 1
            linestyle := :dot
            label := ""
            r.targets[1, :, j], r.targets[2, :, j]
        end
    end

    lo = max(1, k - trail)
    seen = falses(s.ntargets)
    for i in 1:s.nagents
        g = ring_of[i]
        lab = seen[g] ? "" : "group " * ('A' + g - 1)
        seen[g] = true
        @series begin
            seriestype := :path
            marker := :circle
            markersize := 3
            alpha := 0.6
            linewidth := 1.4
            color := ring_color(g)
            linestyle := g == 1 ? :solid : :dash
            label := lab
            [_ax(r, i, kk) for kk in lo:k], [_ay(r, i, kk) for kk in lo:k]
        end
    end

    if show_reference
        @series begin
            seriestype := :scatter
            marker := :circle
            markersize := 4
            markercolor := :white
            markerstrokecolor := GUIDE_COLOR
            markerstrokewidth := 1.1
            label := "q*"
            [_rx(r, i, k) for i in 1:s.nagents], [_ry(r, i, k) for i in 1:s.nagents]
        end
    end

    for j in 1:s.ntargets
        @series begin
            seriestype := :scatter
            marker := :star5
            markersize := 10
            color := target_color(j)
            label := "T$j"
            [r.targets[1, k, j]], [r.targets[2, k, j]]
        end
    end
end

"""
    animate_coordination(rollouts...; filename, fps, frame_step, trail, frames)

Animate one or more [`CoordinationRollout`](@ref)s side by side, sharing axes so
the two laws are directly comparable frame by frame. `frames` restricts the
animation to the window where something happens.
"""
function CellularSheaves.ControlSheaves.CoordinationBenchmarks.animate_coordination(
        rollouts::CoordinationRollout...;
        filename = "coordination_tracking.gif", fps = 15, frame_step = 3,
        trail = 30, frames = nothing, size = (640 * length(rollouts), 480))
    axis_lims = _rollout_limits(rollouts...)
    rings = [agent_groups(r.scenario) for r in rollouts]
    nsteps = minimum(length(r.times) for r in rollouts)
    span = frames === nothing ? (1:frame_step:nsteps) :
        (first(frames):frame_step:min(last(frames), nsteps))
    anim = @animate for k in span
        panels = [plot(r, k; trail, axis_lims, ring_of = rings[i])
                  for (i, r) in enumerate(rollouts)]
        plot(panels...; layout = (1, length(panels)), size = size)
    end
    return gif(anim, filename; fps = fps, show_msg = false)
end

# ==========================
# topology atlas
# ==========================

const AGENT_NODE = :steelblue
const TARGET_NODE = :crimson

"""Force-directed fallback layout for graphs with no natural drawing."""
function _spring_layout(g; iters = 240, seed = 11)
    n = nv(g)
    rng = MersenneTwister(seed)
    x = [cos(2pi * i / n) + 0.02randn(rng) for i in 1:n]
    y = [sin(2pi * i / n) + 0.02randn(rng) for i in 1:n]
    k = sqrt(1.0 / max(n, 1))
    for it in 1:iters
        temp = 0.1 * (1 - it / iters) + 1e-3
        dx = zeros(n); dy = zeros(n)
        for i in 1:n, j in 1:n
            i == j && continue
            ex, ey = x[i] - x[j], y[i] - y[j]
            d2 = max(ex^2 + ey^2, 1e-6)
            dx[i] += k^2 * ex / d2; dy[i] += k^2 * ey / d2
        end
        for e in edges(g)
            i, j = src(e), dst(e)
            ex, ey = x[i] - x[j], y[i] - y[j]
            d = sqrt(max(ex^2 + ey^2, 1e-9))
            fx, fy = d * ex / k, d * ey / k
            dx[i] -= fx; dy[i] -= fy; dx[j] += fx; dy[j] += fy
        end
        for i in 1:n
            m = sqrt(dx[i]^2 + dy[i]^2)
            if m > 1e-9
                x[i] += dx[i] / m * min(m, temp); y[i] += dy[i] / m * min(m, temp)
            end
        end
    end
    return x, y
end

"""Agent coordinates for drawing a scenario: analytic where the family has an
obvious picture, force-directed otherwise."""
function _scenario_layout(s::CoordinationScenario)
    g, n, name = s.agent_graph, s.nagents, s.name
    if name === :chain
        return collect(range(-1, 1; length = n)), zeros(n)
    elseif name === :ring
        th = range(0, 2pi; length = n + 1)[1:(end - 1)]
        return cos.(th), sin.(th)
    elseif name === :grid
        side = round(Int, sqrt(n)); sc = side > 1 ? 2 / (side - 1) : 1.0
        return [(mod(i - 1, side)) * sc - 1 for i in 1:n],
               [(fld(i - 1, side)) * sc - 1 for i in 1:n]
    elseif name === :star || name === :wheel
        th = range(0, 2pi; length = n)[1:(end - 1)]
        return vcat(0.0, cos.(th)), vcat(0.0, sin.(th))
    elseif name === :complete
        th = range(0, 2pi; length = n + 1)[1:(end - 1)]
        return cos.(th), sin.(th)
    elseif name === :ladder
        half = n ÷ 2; xs = collect(range(-1, 1; length = half))
        return vcat(xs, xs), vcat(fill(0.35, half), fill(-0.35, half))
    elseif name === :prism
        half = n ÷ 2; th = range(0, 2pi; length = half + 1)[1:half]
        return vcat(cos.(th), 0.55cos.(th)), vcat(sin.(th), 0.55sin.(th))
    elseif name === :bipartite
        half = n ÷ 2; ys = collect(range(-1, 1; length = half))
        return vcat(fill(-0.7, half), fill(0.7, half)), vcat(ys, ys)
    elseif name === :tree
        xs = zeros(n); ys = zeros(n)
        for i in 1:n
            d = floor(Int, log2(i)); width = 2.0^d; pos = i - 2^d
            xs[i] = width == 1 ? 0.0 : 2 * (pos + 0.5) / width - 1
            ys[i] = 1 - 2d / max(floor(Int, log2(n)), 1)
        end
        return xs, ys
    elseif name === :caterpillar
        spine = n ÷ 2; xs = collect(range(-1, 1; length = spine))
        return vcat(xs, xs), vcat(zeros(spine), fill(-0.6, spine))
    end
    return _spring_layout(g)
end

_marker_size(n) = n <= 20 ? 7.0 : n <= 60 ? 5.0 : 3.5

function _treewidth_hint(s::CoordinationScenario)
    return tree_statistics(direct_plan(s)).treewidth
end

"""
Draw the structure of a [`CoordinationScenario`](@ref): the agent graph, the
pinned targets, and which agents sense them.

    plot(scenario)

Blue circles are free agents, red diamonds are pinned targets, dashed links are
agent-to-target sensing edges.
"""
@recipe function f(s::CoordinationScenario; show_title = true)
    x, y = _scenario_layout(s)
    aspect_ratio --> 1
    legend --> false
    grid --> false
    axis --> false
    ticks --> false
    markerstrokewidth --> 0
    show_title && (title --> @sprintf("%s   n=%d, tw=%d", s.label, s.nagents,
                                      _treewidth_hint(s)))

    tx = Float64[]; ty = Float64[]
    for seeds in s.sensing
        cx = sum(x[i] for i in seeds) / length(seeds)
        cy = sum(y[i] for i in seeds) / length(seeds)
        r = sqrt(cx^2 + cy^2)
        push!(tx, r > 1e-6 ? cx * (1 + 0.55 / r) : 0.0)
        push!(ty, r > 1e-6 ? cy * (1 + 0.55 / r) : -1.45)
    end

    dense = ne(s.agent_graph) > 4nv(s.agent_graph)
    for e in edges(s.agent_graph)
        @series begin
            seriestype := :path
            color := :gray65
            linewidth := dense ? 0.35 : 0.9
            alpha := dense ? 0.25 : 0.6
            label := ""
            [x[src(e)], x[dst(e)]], [y[src(e)], y[dst(e)]]
        end
    end
    for (k, seeds) in enumerate(s.sensing), i in seeds
        @series begin
            seriestype := :path
            color := TARGET_NODE
            linewidth := 0.7
            linestyle := :dash
            alpha := 0.4
            label := ""
            [x[i], tx[k]], [y[i], ty[k]]
        end
    end
    @series begin
        seriestype := :scatter
        color := AGENT_NODE
        markersize := _marker_size(s.nagents)
        label := ""
        x, y
    end
    @series begin
        seriestype := :scatter
        color := TARGET_NODE
        marker := :diamond
        markersize := 7
        markerstrokecolor := :white
        markerstrokewidth := 0.8
        label := ""
        tx, ty
    end
end

end # module CoordinationBenchmarkPlots

