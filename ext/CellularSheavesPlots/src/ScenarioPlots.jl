module ScenarioPlots

using Plots, Statistics, ..Utils
export ScenarioResultRecipe

# Colours and styles that match the existing examples
const AGENT_COLORS  = [:steelblue, :darkorange, :green, :crimson]
const AGENT_STYLES  = [:solid, :dash, :dashdot, :dot]
const TARGET_COLORS = [:gray, :black, :darkgreen, :purple]

# ----------------------------------------------------------------------
# Plots.jl recipe for `ScenarioResult`
# ----------------------------------------------------------------------
@recipe function f(sr::ScenarioResult;
                   layout   = (1,2),
                   size     = (800,380),
                   legend   = :topright,
                   legendfontsize = 7)

    # Global attributes ---------------------------------------------------
    title      := "$(sr.label): null dim = $(sr.null_dim),  ||dz|| = $(round(sr.residual; sigdigits=3))"
    xlabel     := "t (s)"
    xlims      := (0.0, min(10.0, Float64(sr.times[end])))

    # ---- y‑panel --------------------------------------------------------
    @series begin
        subplot   := 1
        ylabel    := "y (m)"
        title     := "y(t)"
        for (i, traj) in enumerate(sr.agent_trajs)
            label     := "A$(i)"
            linecolor := AGENT_COLORS[i]
            linestyle := AGENT_STYLES[i]
            lw        := 2
            sr.times, traj[:, sr.y_col]
        end
    end

    # ---- z‑panel --------------------------------------------------------
    @series begin
        subplot   := 2
        ylabel    := "z (m)"
        title     := "z(t)"
        for (i, traj) in enumerate(sr.agent_trajs)
            label     := "A$(i)"
            linecolor := AGENT_COLORS[i]
            linestyle := AGENT_STYLES[i]
            lw        := 2
            sr.times, traj[:, sr.z_col]
        end
    end

    # ---- target trajectories (appear on both panels) ----------------------
    for (j, traj) in enumerate(sr.target_trajs)
        # y‑panel
        @series begin
            subplot   := 1
            label     := "T$(j)"
            lw        := 1
            linecolor := TARGET_COLORS[j]
            linestyle := :dot
            sr.times, getindex.(traj, IDX_Y)
        end
        # z‑panel
        @series begin
            subplot   := 2
            label     := "T$(j)"
            lw        := 1
            linecolor := TARGET_COLORS[j]
            linestyle := :dot
            sr.times, getindex.(traj, IDX_Z)
        end
    end
end

end # module
