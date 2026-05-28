module SnapshotPlots

using Plots
using CellularSheaves.ControlSheaves.MultiAgentTracking: ScenarioResult
using ..Utils

export animate_tracking_xy

struct TrackingXYFrame{R<:ScenarioResult,TX,TY,L<:Tuple{<:Real,<:Real},M<:Tuple{<:Real,<:Real}}
    result::R
    step::Int
    x_col::TX
    y_col::TY
    xlims::L
    ylims::M
end

@recipe function f(frame::TrackingXYFrame)
    sr = frame.result
    step = frame.step
    x_col = frame.x_col
    y_col = frame.y_col

    legend := :outerright
    xlims := frame.xlims
    ylims := frame.ylims
    title := "t = $(round(sr.times[step]; digits=2))"

    for (i, traj) in enumerate(sr.target_trajs)
        tx, ty = target_path(traj, x_col, y_col, step)

        @series begin
            seriestype := :scatter
            marker := :star5
            markersize := 4
            alpha := 0.6
            label := "T$(i)"
            color := pick_target_color(i)
            tx, ty
        end
    end

    for (i, traj) in enumerate(sr.agent_trajs)
        @series begin
            seriestype := :path
            marker := :circle
            markersize := 4
            alpha := 0.6
            label := "A$(i)"
            color := pick_agent_color(i)
            linestyle := pick_agent_style(i)
            traj[1:step, x_col], traj[1:step, y_col]
        end
    end
end

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

    anim = @animate for (step, _) in enumerate(result.times)
        plot(TrackingXYFrame(result, step, x_col, y_col, xlims, ylims))
    end

    gif(anim, filename, fps=fps)
    return anim
end

end # module SnapshotPlots
