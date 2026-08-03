module LayeredEscortPlots

using Plots
using LinearAlgebra
using Printf
using CellularSheaves.ControlSheaves.DistributedLayeredControl: LayeredSimulationResult

const RING_COLOR = :steelblue
const TARGET_COLOR = :black

@recipe function f(res::LayeredSimulationResult, k::Int)
    prob = res.problem
    sim_d = res.sim_data
    q_d = res.qstar_history
    ts = (1:prob.steps) .* prob.dt
    t_curr = k * prob.dt
    NA = length(prob.sheaf.vertex_stalks) - length(prob.target_nodes)
    TV1 = prob.target_nodes[1]
    
    # Compute axis limits and metrics based on full trajectory (cached implicitly by plotting once, though ideally could be pre-computed)
    all_xs = vcat([sim_d[:, i, 1] for i in 1:NA]...)
    all_ys = vcat([sim_d[:, i, 2] for i in 1:NA]...)
    target_xs = [prob.target_trajectory_func(TV1, t)[1] for t in ts]
    target_ys = [prob.target_trajectory_func(TV1, t)[2] for t in ts]
    min_xy = minimum(vcat(all_xs, all_ys, target_xs, target_ys))
    max_xy = maximum(vcat(all_xs, all_ys, target_xs, target_ys))
    pad_xy = (max_xy - min_xy) * 0.15 + 0.2
    lims_xy = (min_xy - pad_xy, max_xy + pad_xy)
    lims_rel = (-prob.r_ring*1.8, prob.r_ring*1.8)
    
    max_roll = maximum(abs.(rad2deg.(sim_d[:, :, 4])))
    max_pitch = maximum(abs.(rad2deg.(sim_d[:, :, 5])))
    lims_roll = (-max_roll*1.2 - 1.0, max_roll*1.2 + 1.0)
    lims_pitch = (-max_pitch*1.2 - 1.0, max_pitch*1.2 + 1.0)

    formation_dev = zeros(prob.steps)
    tracking_dev = zeros(prob.steps)
    for step in 1:prob.steps
        c = [sum(sim_d[step, :, 1])/NA, sum(sim_d[step, :, 2])/NA]
        t = prob.target_trajectory_func(TV1, ts[step])[1:2]
        tracking_dev[step] = norm(c - t)
        
        dev = 0.0
        for i in 1:NA
            dev += abs(norm([sim_d[step, i, 1], sim_d[step, i, 2]] - c) - prob.r_ring)
            next_i = (i % NA) + 1
            dev += abs(norm([sim_d[step, i, 1], sim_d[step, i, 2]] - [sim_d[step, next_i, 1], sim_d[step, next_i, 2]]) - prob.r_ring)
        end
        formation_dev[step] = dev / (2*NA)
    end
    
    layout := (2, 3)
    size --> (1200, 700)
    plot_title := @sprintf("6-Agent SE(3) Moving Escort Ring (t = %.2f s)", t_curr)
    legend := false

    # ==========================
    # Panel 1: World Top-Down View
    # ==========================
    target_orbit_t = range(0, prob.steps*prob.dt; length=100)
    @series begin
        subplot := 1
        aspect_ratio := 1
        xlims := lims_xy
        ylims := lims_xy
        xlabel := "x position (m)"
        ylabel := "y position (m)"
        title := "World Top-Down View (x-y Plane)"
        color := :gray80
        linestyle := :dot
        linewidth := 1
        [prob.target_trajectory_func(TV1, t)[1] for t in target_orbit_t], [prob.target_trajectory_func(TV1, t)[2] for t in target_orbit_t]
    end
    @series begin
        subplot := 1
        seriestype := :scatter
        marker := :star5
        markersize := 10
        color := TARGET_COLOR
        [prob.target_trajectory_func(TV1, t_curr)[1]], [prob.target_trajectory_func(TV1, t_curr)[2]]
    end
    
    c_curr_x = sum(sim_d[k, :, 1])/NA
    c_curr_y = sum(sim_d[k, :, 2])/NA
    @series begin
        subplot := 1
        seriestype := :scatter
        marker := :square
        markersize := 5
        color := :red
        [c_curr_x], [c_curr_y]
    end
    
    for i in 1:NA
        @series begin
            subplot := 1
            seriestype := :path
            marker := :circle
            markersize := 3
            alpha := 0.6
            linewidth := 1.4
            color := RING_COLOR
            sim_d[1:k, i, 1], sim_d[1:k, i, 2]
        end
        @series begin
            subplot := 1
            seriestype := :scatter
            marker := :square
            markersize := 4
            color := :purple
            alpha := 0.3
            [q_d[k, i, 1]], [q_d[k, i, 2]]
        end
    end
    
    ring_x = [sim_d[k, i, 1] for i in 1:NA]
    push!(ring_x, sim_d[k, 1, 1])
    ring_y = [sim_d[k, i, 2] for i in 1:NA]
    push!(ring_y, sim_d[k, 1, 2])
    @series begin
        subplot := 1
        color := :gray80
        linestyle := :dot
        linewidth := 1
        ring_x, ring_y
    end

    # ==========================
    # Panel 2: Target-Centered Escort Ring
    # ==========================
    circ_ang = range(0, 2π; length = 100)
    @series begin
        subplot := 2
        aspect_ratio := 1
        xlims := lims_rel
        ylims := lims_rel
        xlabel := "rel x to target (m)"
        ylabel := "rel y to target (m)"
        title := "Target-Centered Escort Ring ($(prob.r_ring)m)"
        color := :gray80
        linestyle := :dash
        linewidth := 1
        prob.r_ring .* cos.(circ_ang), prob.r_ring .* sin.(circ_ang)
    end
    @series begin
        subplot := 2
        seriestype := :scatter
        marker := :star5
        markersize := 10
        color := TARGET_COLOR
        [0.0], [0.0]
    end
    
    for i in 1:NA
        rel_x_hist = [sim_d[step, i, 1] - prob.target_trajectory_func(TV1, ts[step])[1] for step in 1:k]
        rel_y_hist = [sim_d[step, i, 2] - prob.target_trajectory_func(TV1, ts[step])[2] for step in 1:k]
        @series begin
            subplot := 2
            seriestype := :path
            marker := :circle
            markersize := 3
            alpha := 0.6
            linewidth := 1.4
            color := RING_COLOR
            rel_x_hist, rel_y_hist
        end
        ref_rel_x = q_d[k, i, 1] - prob.target_trajectory_func(TV1, t_curr)[1]
        ref_rel_y = q_d[k, i, 2] - prob.target_trajectory_func(TV1, t_curr)[2]
        @series begin
            subplot := 2
            seriestype := :scatter
            marker := :square
            markersize := 4
            color := :purple
            alpha := 0.3
            [ref_rel_x], [ref_rel_y]
        end
    end
    
    ring_rel_x = [sim_d[k, i, 1] - prob.target_trajectory_func(TV1, t_curr)[1] for i in 1:NA]
    push!(ring_rel_x, sim_d[k, 1, 1] - prob.target_trajectory_func(TV1, t_curr)[1])
    ring_rel_y = [sim_d[k, i, 2] - prob.target_trajectory_func(TV1, t_curr)[2] for i in 1:NA]
    push!(ring_rel_y, sim_d[k, 1, 2] - prob.target_trajectory_func(TV1, t_curr)[2])
    @series begin
        subplot := 2
        color := :gray80
        linestyle := :dot
        linewidth := 1
        ring_rel_x, ring_rel_y
    end

    # ==========================
    # Panel 3: Roll Tilt Dynamics
    # ==========================
    @series begin
        subplot := 3
        xlabel := "time (s)"
        ylabel := "roll angle ϕ (deg)"
        title := "Roll Tilt Dynamics ϕ(t)"
        xlims := (0, ts[end])
        ylims := lims_roll
        # Just to set the axis limits correctly
        [0.0], [0.0]
    end
    for i in 1:NA
        @series begin
            subplot := 3
            linewidth := 1.2
            color := RING_COLOR
            ts[1:k], rad2deg.(sim_d[1:k, i, 4])
        end
    end

    # ==========================
    # Panel 4: Pitch Tilt Dynamics
    # ==========================
    @series begin
        subplot := 4
        xlabel := "time (s)"
        ylabel := "pitch angle θ (deg)"
        title := "Pitch Tilt Dynamics θ(t)"
        xlims := (0, ts[end])
        ylims := lims_pitch
        [0.0], [0.0]
    end
    for i in 1:NA
        @series begin
            subplot := 4
            linewidth := 1.2
            color := RING_COLOR
            ts[1:k], rad2deg.(sim_d[1:k, i, 5])
        end
    end

    # ==========================
    # Panel 5: Formation Shape Error
    # ==========================
    @series begin
        subplot := 5
        xlabel := "time (s)"
        ylabel := "error (m)"
        title := "Formation Shape Error"
        xlims := (0, ts[end])
        ylims := (0, maximum(formation_dev)*1.2 + 0.01)
        linewidth := 1.5
        color := :purple
        ts[1:k], formation_dev[1:k]
    end

    # ==========================
    # Panel 6: Target Tracking Error
    # ==========================
    @series begin
        subplot := 6
        xlabel := "time (s)"
        ylabel := "error (m)"
        title := "Target Tracking Error"
        xlims := (0, ts[end])
        ylims := (0, maximum(tracking_dev)*1.2 + 0.01)
        linewidth := 1.5
        color := :orange
        ts[1:k], tracking_dev[1:k]
    end
end

end # module
