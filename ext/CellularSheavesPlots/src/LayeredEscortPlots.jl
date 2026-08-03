module LayeredEscortPlots

using Plots

using LinearAlgebra
using Printf
using CellularSheaves.ControlSheaves.DistributedLayeredControl: LayeredSimulationResult

const RING_COLOR = :steelblue
const TARGET_COLOR = :black

# ==========================
# Metric primitives
# ==========================

# By default, plot the whole trajectory if `k` is not provided.
@recipe function f(res::LayeredSimulationResult, metric::Symbol; k=res.problem.steps, state_idx=1, state_name="state")
    prob = res.problem
    sim_d = res.sim_data
    q_d = res.qstar_history
    ts = (1:prob.steps) .* prob.dt
    t_curr = k * prob.dt
    NA = length(prob.sheaf.vertex_stalks) - length(prob.target_nodes)
    TV1 = prob.target_nodes[1]
    
    if metric == :tracking_error
        tracking_dev = zeros(prob.steps)
        for step in 1:prob.steps
            c = [sum(sim_d[step, :, 1])/NA, sum(sim_d[step, :, 2])/NA]
            t = prob.target_trajectory_func(TV1, ts[step])[1:2]
            tracking_dev[step] = norm(c - t)
        end
        xlabel --> "time (s)"
        ylabel --> "error (m)"
        title --> "Target Tracking Error"
        xlims --> (0, ts[end])
        ylims --> (0, maximum(tracking_dev)*1.2 + 0.01)
        linewidth --> 1.5
        color --> :orange
        label --> false
        return ts[1:k], tracking_dev[1:k]

    elseif metric == :formation_error
        formation_dev = zeros(prob.steps)
        r_ring = hasfield(typeof(prob), :r_ring) ? prob.r_ring : 1.0 # fallback for scenarios without r_ring
        for step in 1:prob.steps
            c = [sum(sim_d[step, :, 1])/NA, sum(sim_d[step, :, 2])/NA]
            dev = 0.0
            for i in 1:NA
                dev += abs(norm([sim_d[step, i, 1], sim_d[step, i, 2]] - c) - r_ring)
                next_i = (i % NA) + 1
                dev += abs(norm([sim_d[step, i, 1], sim_d[step, i, 2]] - [sim_d[step, next_i, 1], sim_d[step, next_i, 2]]) - r_ring)
            end
            formation_dev[step] = dev / (2*NA)
        end
        xlabel --> "time (s)"
        ylabel --> "error (m)"
        title --> "Formation Shape Error"
        xlims --> (0, ts[end])
        ylims --> (0, maximum(formation_dev)*1.2 + 0.01)
        linewidth --> 1.5
        color --> :purple
        label --> false
        return ts[1:k], formation_dev[1:k]

    elseif metric == :top_down
        target_orbit_t = range(0, prob.steps*prob.dt; length=100)
        
        all_xs = vcat([sim_d[:, i, 1] for i in 1:NA]...)
        all_ys = vcat([sim_d[:, i, 2] for i in 1:NA]...)
        target_xs = [prob.target_trajectory_func(TV1, t)[1] for t in ts]
        target_ys = [prob.target_trajectory_func(TV1, t)[2] for t in ts]
        min_x, max_x = minimum(vcat(all_xs, target_xs)), maximum(vcat(all_xs, target_xs))
        min_y, max_y = minimum(vcat(all_ys, target_ys)), maximum(vcat(all_ys, target_ys))

        cx = (min_x + max_x) / 2
        cy = (min_y + max_y) / 2
        span = max(max_x - min_x, max_y - min_y) / 2 + max(max(max_x - min_x, max_y - min_y) * 0.1, 0.5)

        aspect_ratio --> 1
        xlims --> (cx - span, cx + span)
        ylims --> (cy - span, cy + span)
        xlabel --> "x position (m)"
        ylabel --> "y position (m)"
        title --> "World Top-Down View (x-y Plane)"
        
        @series begin
            color := :gray80
            linestyle := :dot
            linewidth := 1
            label := false
            [prob.target_trajectory_func(TV1, t)[1] for t in target_orbit_t], [prob.target_trajectory_func(TV1, t)[2] for t in target_orbit_t]
        end
        @series begin
            seriestype := :scatter
            marker := :star5
            markersize := 10
            color := TARGET_COLOR
            label := false
            [prob.target_trajectory_func(TV1, t_curr)[1]], [prob.target_trajectory_func(TV1, t_curr)[2]]
        end
        
        c_curr_x = sum(sim_d[k, :, 1])/NA
        c_curr_y = sum(sim_d[k, :, 2])/NA
        @series begin
            seriestype := :scatter
            marker := :square
            markersize := 5
            color := :red
            label := false
            [c_curr_x], [c_curr_y]
        end
        
        for i in 1:NA
            @series begin
                seriestype := :path
                marker := :circle
                markersize := 3
                alpha := 0.6
                linewidth := 1.4
                color := RING_COLOR
                label := false
                sim_d[1:k, i, 1], sim_d[1:k, i, 2]
            end
            @series begin
                seriestype := :scatter
                marker := :square
                markersize := 4
                color := :purple
                alpha := 0.3
                label := false
                [q_d[k, i, 1]], [q_d[k, i, 2]]
            end
        end
        
        ring_x = [sim_d[k, i, 1] for i in 1:NA]
        push!(ring_x, sim_d[k, 1, 1])
        ring_y = [sim_d[k, i, 2] for i in 1:NA]
        push!(ring_y, sim_d[k, 1, 2])
        @series begin
            color := :gray80
            linestyle := :dot
            linewidth := 1
            label := false
            ring_x, ring_y
        end

    elseif metric == :target_centered
        r_ring = hasfield(typeof(prob), :r_ring) ? prob.r_ring : 1.0
        all_rel_xs = vcat([[sim_d[step, i, 1] - prob.target_trajectory_func(TV1, ts[step])[1] for step in 1:prob.steps] for i in 1:NA]...)
        all_rel_ys = vcat([[sim_d[step, i, 2] - prob.target_trajectory_func(TV1, ts[step])[2] for step in 1:prob.steps] for i in 1:NA]...)
        all_rel_q_xs = vcat([[q_d[step, i, 1] - prob.target_trajectory_func(TV1, ts[step])[1] for step in 1:prob.steps] for i in 1:NA]...)
        all_rel_q_ys = vcat([[q_d[step, i, 2] - prob.target_trajectory_func(TV1, ts[step])[2] for step in 1:prob.steps] for i in 1:NA]...)
        
        max_rel = maximum(vcat(abs.(all_rel_xs), abs.(all_rel_ys), abs.(all_rel_q_xs), abs.(all_rel_q_ys), [r_ring * 1.2]))
        span_rel = max_rel + max(max_rel * 0.1, 0.2)
        lims_rel = (-span_rel, span_rel)
        circ_ang = range(0, 2π; length = 100)
        
        aspect_ratio --> 1
        xlims --> lims_rel
        ylims --> lims_rel
        xlabel --> "rel x to target (m)"
        ylabel --> "rel y to target (m)"
        title --> "Target-Centered View ($(r_ring)m)"
        
        @series begin
            color := :gray80
            linestyle := :dash
            linewidth := 1
            label := false
            r_ring .* cos.(circ_ang), r_ring .* sin.(circ_ang)
        end
        @series begin
            seriestype := :scatter
            marker := :star5
            markersize := 10
            color := TARGET_COLOR
            label := false
            [0.0], [0.0]
        end
        
        for i in 1:NA
            rel_x_hist = [sim_d[step, i, 1] - prob.target_trajectory_func(TV1, ts[step])[1] for step in 1:k]
            rel_y_hist = [sim_d[step, i, 2] - prob.target_trajectory_func(TV1, ts[step])[2] for step in 1:k]
            @series begin
                seriestype := :path
                marker := :circle
                markersize := 3
                alpha := 0.6
                linewidth := 1.4
                color := RING_COLOR
                label := false
                rel_x_hist, rel_y_hist
            end
            ref_rel_x = q_d[k, i, 1] - prob.target_trajectory_func(TV1, t_curr)[1]
            ref_rel_y = q_d[k, i, 2] - prob.target_trajectory_func(TV1, t_curr)[2]
            @series begin
                seriestype := :scatter
                marker := :square
                markersize := 4
                color := :purple
                alpha := 0.3
                label := false
                [ref_rel_x], [ref_rel_y]
            end
        end
        
        ring_rel_x = [sim_d[k, i, 1] - prob.target_trajectory_func(TV1, t_curr)[1] for i in 1:NA]
        push!(ring_rel_x, sim_d[k, 1, 1] - prob.target_trajectory_func(TV1, t_curr)[1])
        ring_rel_y = [sim_d[k, i, 2] - prob.target_trajectory_func(TV1, t_curr)[2] for i in 1:NA]
        push!(ring_rel_y, sim_d[k, 1, 2] - prob.target_trajectory_func(TV1, t_curr)[2])
        @series begin
            color := :gray80
            linestyle := :dot
            linewidth := 1
            label := false
            ring_rel_x, ring_rel_y
        end

    elseif metric == :state_timeseries
        xlabel --> "time (s)"
        ylabel --> state_name
        title --> "$(titlecase(state_name)) Dynamics"
        xlims --> (0, ts[end])
        
        # calculate bounds across all agents
        max_val = maximum(sim_d[:, :, state_idx])
        min_val = minimum(sim_d[:, :, state_idx])
        pad = (max_val - min_val) * 0.2
        if pad == 0.0
            pad = 1.0
        end
        ylims --> (min_val - pad, max_val + pad)
        
        for i in 1:NA
            @series begin
                linewidth := 1.2
                color := (i == 1 ? :steelblue : (i == 2 ? :darkorange : RING_COLOR))
                label := (i <= 2 ? "A$i Actual" : false) # label the first two for legends if needed
                ts[1:k], sim_d[1:k, i, state_idx]
            end
            if size(q_d, 3) >= state_idx
                @series begin
                    alpha := 0.5
                    ls := :dash
                    marker := :circle
                    ms := 2
                    color := (i == 1 ? :steelblue : (i == 2 ? :darkorange : RING_COLOR))
                    label := (i <= 2 ? "A$i q*" : false)
                    ts[1:k], q_d[1:k, i, state_idx]
                end
            end
        end
        
        # Also plot targets if target_trajectory_func provides it up to state_idx length
        # We assume targets have positions/states in the same index order
        NT = length(prob.target_nodes)
        for j in 1:NT
            TNode = prob.target_nodes[j]
            t_traj = [prob.target_trajectory_func(TNode, t) for t in ts[1:k]]
            if length(t_traj[1]) >= state_idx
                @series begin
                    ls := :dot
                    color := (j == 1 ? :gray : :black)
                    label := "Target $j"
                    ts[1:k], [traj[state_idx] for traj in t_traj]
                end
            end
        end

    else
        error("Unknown metric $metric")
    end
end


# ==========================
# Scenario 5 frame recipe  (2-agent planar tracking)
# ==========================

"""
Helper: precompute full-trajectory bounding box and static limits for scenario5 panels.
Returns NamedTuple with all precomputed extents.
"""
function _scenario5_precompute(res::LayeredSimulationResult)
    prob = res.problem
    sim  = res.sim_data
    qh   = res.qstar_history
    ts   = (1:prob.steps) .* prob.dt
    NA   = size(sim, 2)
    T1, T2 = prob.target_nodes[1], prob.target_nodes[2]

    t1_y = [prob.target_trajectory_func(T1, t)[1] for t in ts]
    t1_z = [prob.target_trajectory_func(T1, t)[2] for t in ts]
    t2_y = [prob.target_trajectory_func(T2, t)[1] for t in ts]
    t2_z = [prob.target_trajectory_func(T2, t)[2] for t in ts]

    all_y  = vcat(sim[:, :, 1]..., t1_y, t2_y)
    all_z  = vcat(sim[:, :, 2]..., t1_z, t2_z)
    min_y, max_y = minimum(all_y), maximum(all_y)
    min_z, max_z = minimum(all_z), maximum(all_z)
    cy  = (min_y + max_y) / 2;  cz = (min_z + max_z) / 2
    span_yz = max(max_y - min_y, max_z - min_z) / 2 + 0.4

    all_theta = sim[:, :, 3]
    max_th    = maximum(abs.(all_theta))
    pad_th    = max(max_th * 0.2, 0.05)

    terr_all  = vcat([[norm(sim[step, i, 1:2] - qh[step, i, :]) for step in 1:prob.steps] for i in 1:NA]...)
    max_terr  = maximum(terr_all)

    return (
        ts=ts, NA=NA, T1=T1, T2=T2,
        t1_y=t1_y, t1_z=t1_z, t2_y=t2_y, t2_z=t2_z,
        cy=cy, cz=cz, span_yz=span_yz,
        max_th=max_th, pad_th=pad_th,
        max_terr=max_terr,
    )
end

"""
Build the y-z position-plane panel for frame k of a Scenario 5 simulation.
"""
function _scenario5_yz_panel(res::LayeredSimulationResult, pre::NamedTuple, k::Int)
    prob = res.problem
    sim  = res.sim_data
    qh   = res.qstar_history
    ts   = pre.ts
    NA   = pre.NA
    T1, T2 = pre.T1, pre.T2
    t_curr = k * prob.dt

    p = plot(;
        title        = "Agent Positions [y-z Plane]",
        xlabel       = "y position [m]",
        ylabel       = "z position [m]",
        aspect_ratio = 1,
        xlims        = (pre.cy - pre.span_yz, pre.cy + pre.span_yz),
        ylims        = (pre.cz - pre.span_yz, pre.cz + pre.span_yz),
        legend       = :topleft,
        bottom_margin = 5Plots.mm,
        left_margin   = 5Plots.mm,
    )
    # Target orbit trails
    plot!(p, pre.t1_y, pre.t1_z; ls=:dot, lw=1, color=:gray60, label=false)
    plot!(p, pre.t2_y, pre.t2_z; ls=:dot, lw=1, color=:gray30, label=false)
    # Current target positions
    scatter!(p, [prob.target_trajectory_func(T1, t_curr)[1]],
                [prob.target_trajectory_func(T1, t_curr)[2]];
             marker=:star5, ms=10, color=:gray60, label="Target 1")
    scatter!(p, [prob.target_trajectory_func(T2, t_curr)[1]],
                [prob.target_trajectory_func(T2, t_curr)[2]];
             marker=:star5, ms=10, color=:gray30, label="Target 2")
    # Agent trails and harmonic extension squares
    for i in 1:NA
        c = i == 1 ? RING_COLOR : :darkorange
        plot!(p, sim[1:k, i, 1], sim[1:k, i, 2];
              seriestype=:path, marker=:circle, ms=3, lw=1.4, alpha=0.7,
              color=c, label="Agent $i")
        scatter!(p, [qh[k, i, 1]], [qh[k, i, 2]];
                 marker=:square, ms=5, color=:purple, alpha=0.3,
                 label=(i == 1 ? "Harmonic ext." : false))
    end
    return p
end

"""
Build the tilt-angle panel for frame k of a Scenario 5 simulation.
"""
function _scenario5_tilt_panel(res::LayeredSimulationResult, pre::NamedTuple, k::Int)
    sim = res.sim_data
    ts  = pre.ts

    p = plot(;
        title  = "Quadrotor Tilt [θ vs Time]",
        xlabel = "time [s]",
        ylabel = "tilt angle θ [rad]",
        xlims  = (0, ts[end]),
        ylims  = (-pre.max_th - pre.pad_th, pre.max_th + pre.pad_th),
        legend = :topleft,
        bottom_margin = 5Plots.mm,
        left_margin   = 5Plots.mm,
    )
    plot!(p, ts[1:k], sim[1:k, 1, 3]; lw=1.2, color=RING_COLOR,    label="Agent 1")
    plot!(p, ts[1:k], sim[1:k, 2, 3]; lw=1.2, color=:darkorange, label="Agent 2")
    return p
end

"""
Build the tracking-error panel for frame k of a Scenario 5 simulation.
"""
function _scenario5_error_panel(res::LayeredSimulationResult, pre::NamedTuple, k::Int)
    sim = res.sim_data
    qh  = res.qstar_history
    ts  = pre.ts

    p = plot(;
        title  = "Tracking Error [‖pos − q*‖]",
        xlabel = "time [s]",
        ylabel = "position error [m]",
        yscale = :log10,
        xlims  = (0, ts[end]),
        ylims  = (1e-4, max(pre.max_terr * 2.0, 1.0)),
        legend = :topleft,
        bottom_margin = 5Plots.mm,
        left_margin   = 5Plots.mm,
    )
    err1 = [norm(sim[step, 1, 1:2] - qh[step, 1, :]) for step in 1:k]
    err2 = [norm(sim[step, 2, 1:2] - qh[step, 2, :]) for step in 1:k]
    plot!(p, ts[1:k], max.(err1, 1e-6); lw=1.5, color=RING_COLOR,    label="Agent 1")
    plot!(p, ts[1:k], max.(err2, 1e-6); lw=1.5, color=:darkorange, label="Agent 2")
    return p
end

# ==========================
# Main layout recipe (for escort mission)
# ==========================
@recipe function f(res::LayeredSimulationResult, k::Int)
    prob = res.problem
    t_curr = k * prob.dt

    layout := (2, 3)
    size --> (1200, 700)
    plot_title := @sprintf("6-Agent SE(3) Moving Escort Ring (t = %.2f s)", t_curr)
    legend := false

    @series begin
        subplot := 1
        k := k
        res, :top_down
    end

    @series begin
        subplot := 2
        k := k
        res, :target_centered
    end

    @series begin
        subplot := 3
        k := k
        state_idx := 4
        # Note: in original code roll/pitch were deg, but here they are rad if we don't scale.
        # Let's override y data by scaling? The easiest way is to let the user pass `yguide` 
        # But wait, original code scaled it with rad2deg... Let's just plot radians for now 
        # since `state_timeseries` is generic.
        state_name := "roll angle ϕ (rad)"
        res, :state_timeseries
    end

    @series begin
        subplot := 4
        k := k
        state_idx := 5
        state_name := "pitch angle θ (rad)"
        res, :state_timeseries
    end

    @series begin
        subplot := 5
        k := k
        res, :formation_error
    end

    @series begin
        subplot := 6
        k := k
        res, :tracking_error
    end
end

end # module
