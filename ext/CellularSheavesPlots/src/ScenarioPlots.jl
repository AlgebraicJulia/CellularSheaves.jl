module ScenarioPlots

using Plots
using CellularSheaves.ControlSheaves.MultiAgentTracking: ScenarioResult
using ..Utils

# ----------------------------------------------------------------------
# Plots.jl recipe for `ScenarioResult`
# ----------------------------------------------------------------------
@recipe function f(sr::ScenarioResult;
                   layout=(1, 2),
                   size=(800, 380),
                   legend=:topright,
                   legendfontsize=7)

    title := "$(sr.label): null dim = $(sr.null_dim), ||dz|| = $(round(sr.residual; sigdigits=3))"
    xlabel := "t (s)"
    xlims := (0.0, min(10.0, Float64(sr.times[end])))

    for (i, traj) in enumerate(sr.agent_trajs)
        @series begin
            subplot := 1
            ylabel := "y (m)"
            title := "y(t)"
            label := "A$(i)"
            linecolor := pick_agent_color(i)
            linestyle := pick_agent_style(i)
            lw := 2
            sr.times, traj[:, sr.y_col]
        end
    end

    for (i, traj) in enumerate(sr.agent_trajs)
        @series begin
            subplot := 2
            ylabel := "z (m)"
            title := "z(t)"
            label := "A$(i)"
            linecolor := pick_agent_color(i)
            linestyle := pick_agent_style(i)
            lw := 2
            sr.times, traj[:, sr.z_col]
        end
    end

    for (j, traj) in enumerate(sr.target_trajs)
        ty = [safe_coord(v, sr.y_col) for v in traj]
        tz = [safe_coord(v, sr.z_col) for v in traj]

        @series begin
            subplot := 1
            label := "T$(j)"
            lw := 1
            linecolor := pick_target_color(j)
            linestyle := :dot
            sr.times, ty
        end

        @series begin
            subplot := 2
            label := "T$(j)"
            lw := 1
            linecolor := pick_target_color(j)
            linestyle := :dot
            sr.times, tz
        end
    end
end

end # module
