module Trajectory3D

using Plots
using CellularSheaves.ControlSheaves.MultiAgentTracking: ScenarioResult
using ..Utils

export plot_trajectories_3d

function plot_trajectories_3d(
    result::ScenarioResult;
    x_col::Int=1,
    y_col::Int=result.y_col,
    z_col::Int=result.z_col,
    camera::Tuple{<:Real,<:Real}=(45, 25),
)
    p = plot(
        xlabel="x (m)",
        ylabel="y (m)",
        zlabel="z (m)",
        title="Trajectories",
        camera=camera,
        legend=:outerleft,
        size=(1100, 800),
    )

    for (i, traj) in enumerate(result.agent_trajs)
        plot!(
            p,
            traj[:, x_col],
            traj[:, y_col],
            traj[:, z_col];
            label="A$(i)",
            color=pick_agent_color(i),
            linestyle=pick_agent_style(i),
            linewidth=2,
        )
    end

    for (i, traj) in enumerate(result.target_trajs)
        tx, ty = target_path(traj, x_col, y_col, length(traj))
        tz = [safe_coord(v, z_col) for v in traj]
        plot!(
            p,
            tx,
            ty,
            tz;
            label="T$(i)",
            color=pick_target_color(i),
            linestyle=:dot,
            linewidth=2,
        )
    end

    return p
end

end # module Trajectory3D
