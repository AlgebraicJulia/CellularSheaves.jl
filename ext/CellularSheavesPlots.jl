module CellularSheavesPlots

using CellularSheaves
using CellularSheaves.AsynchSheaves
using CellularSheaves.ControlSheaves.MultiAgentTracking: ScenarioResult
using Plots

function CellularSheaves.AsynchSheaves.empty_experiment_plot(; kwargs...)
    return plot(title="", thickness_scaling=2.0, yformatter=:plain, xformatter=:plain; kwargs...)
end

function CellularSheaves.AsynchSheaves.plot_loss_curve!(plt, losses, label; kwargs...)
    return plot!(plt, losses, label=label, linewidth=2; kwargs...)
end

function CellularSheaves.AsynchSheaves.plot_log_loss_curve!(plt, losses, label; kwargs...)
    return plot!(plt, yscale=:log10, losses, label=label, linewidth=2; kwargs...)
end

# ---------------------------------------------------------------------------
# Plots.jl recipe for ScenarioResult
# Renders a ScenarioResult as a two-panel y(t) / z(t) figure.
# Agent trajectories use solid/dashed lines; target reference trajectories
# are shown as dotted lines for comparison.
# ---------------------------------------------------------------------------

@recipe function f(sr::ScenarioResult)
    layout := (1, 2)
    size := (800, 380)
    plot_title := "$(sr.label): null dim = $(sr.null_dim),  ||dz|| = $(round(sr.residual; sigdigits=3))"
    agent_colors  = [:steelblue, :darkorange, :green, :crimson]
    agent_styles  = [:solid, :dash, :dashdot, :dot]
    target_colors = [:gray, :black, :darkgreen, :purple]
    for (i, traj) in enumerate(sr.agent_trajs)
        @series begin
            subplot   := 1
            title     := "y(t)"
            xlabel    := "t (s)"
            ylabel    := "y (m)"
            label     := "A$i"
            lw        := 2
            linecolor := agent_colors[i]
            linestyle := agent_styles[i]
            sr.times, traj[:, sr.y_col]
        end
        @series begin
            subplot   := 2
            title     := "z(t)"
            xlabel    := "t (s)"
            ylabel    := "z (m)"
            label     := "A$i"
            lw        := 2
            linecolor := agent_colors[i]
            linestyle := agent_styles[i]
            sr.times, traj[:, sr.z_col]
        end
    end
    for (j, traj_j) in enumerate(sr.target_trajs)
        @series begin
            subplot   := 1
            label     := "T$j"
            lw        := 1
            linecolor := target_colors[j]
            linestyle := :dot
            sr.times, getindex.(traj_j, sr.y_col)
        end
        @series begin
            subplot   := 2
            label     := "T$j"
            lw        := 1
            linecolor := target_colors[j]
            linestyle := :dot
            sr.times, getindex.(traj_j, sr.z_col)
        end
    end
end

end # module CellularSheavesPlots
