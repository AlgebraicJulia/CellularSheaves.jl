module Trajectory3D

using Plots
using CellularSheaves.ControlSheaves.MultiAgentTracking: ScenarioResult
using ..Utils

export plot_trajectories_3d

struct Trajectory3DData{R<:ScenarioResult}
    result::R
    x_col::Int
    y_col::Int
    z_col::Int
end

@recipe function f(data::Trajectory3DData;
                   camera=(45, 25))
    sr = data.result
    x_col = data.x_col
    y_col = data.y_col
    z_col = data.z_col

    xlabel := "x (m)"
    ylabel := "y (m)"
    zlabel := "z (m)"
    title := "Trajectories"
    legend := :outerleft
    size := (1100, 800)
    camera := camera

    for (i, traj) in enumerate(sr.agent_trajs)
        @series begin
            seriestype := :path3d
            label := "A$(i)"
            color := pick_agent_color(i)
            linestyle := pick_agent_style(i)
            linewidth := 2
            traj[:, x_col], traj[:, y_col], traj[:, z_col]
        end
    end

    for (i, traj) in enumerate(sr.target_trajs)
        tx, ty = target_path(traj, x_col, y_col, length(traj))
        tz = [safe_coord(v, z_col) for v in traj]
        @series begin
            seriestype := :path3d
            label := "T$(i)"
            color := pick_target_color(i)
            linestyle := :dot
            linewidth := 2
            tx, ty, tz
        end
    end
end

function plot_trajectories_3d(
    result::ScenarioResult;
    x_col::Int=1,
    y_col::Int=result.y_col,
    z_col::Int=result.z_col,
    camera::Tuple{<:Real,<:Real}=(45, 25),
)
    return plot(Trajectory3DData(result, x_col, y_col, z_col); camera=camera)
end

end # module Trajectory3D
