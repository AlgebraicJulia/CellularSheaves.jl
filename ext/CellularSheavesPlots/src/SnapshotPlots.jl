module SnapshotPlots

using Plots
using CellularSheaves.ControlSheaves.MultiAgentTracking: ScenarioResult
using ..Utils

export animate_tracking_xy

function animate_tracking_xy(
    result::ScenarioResult;
    fps::Int=15,
    filename::AbstractString="./tracking_animation.gif",
    x_col::Int=1,
    y_col::Int=2,
    xlims::Tuple{<:Real,<:Real}=(-2.0, 4.0),
    ylims::Tuple{<:Real,<:Real}=(-2.0, 4.0),
)
    fps >= 1 || throw(ArgumentError("fps must be >= 1"))

    anim = @animate for (step, t) in enumerate(result.times)
        p = plot()

        for (i, traj) in enumerate(result.agent_trajs)
            tx, ty = target_path(result.target_trajs[i], x_col, y_col, step)
            plot!(
                p,
                tx,
                ty;
                marker=:star5,
                color=pick_target_color(i),
                label="T$(i)",
                ms=4,
                alpha=0.6,
            )

            plot!(
                p,
                traj[1:step, x_col],
                traj[1:step, y_col];
                alpha=0.6,
                color=pick_agent_color(i),
                linestyle=pick_agent_style(i),
                label="A$(i)",
                marker=:circle,
                ms=4,
                xlim=xlims,
                ylim=ylims,
            )
        end

        title!(p, "t = $(round(t; digits=2))")
        p
    end

    gif(anim, filename, fps=fps)
    return anim
end

end # module SnapshotPlots
