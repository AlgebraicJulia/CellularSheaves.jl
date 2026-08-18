module CellularSheavesPlots

using CellularSheaves
using CellularSheaves.AsynchSheaves
using CellularSheaves.ControlSheaves.MultiAgentTracking: ScenarioResult
using Plots
using Printf

include("CellularSheavesPlots/src/Utils.jl")
include("CellularSheavesPlots/src/Style.jl")
include("CellularSheavesPlots/src/ScenarioPlots.jl")
include("CellularSheavesPlots/src/ErrorNormPlots.jl")
include("CellularSheavesPlots/src/Trajectory3D.jl")
include("CellularSheavesPlots/src/SnapshotPlots.jl")
include("CellularSheavesPlots/src/SpyPlots.jl")
include("CellularSheavesPlots/src/CoordinationBenchmarkPlots.jl")

function CellularSheaves.AsynchSheaves.empty_experiment_plot(; kwargs...)
    return plot(title="", thickness_scaling=2.0, yformatter=:plain, xformatter=:plain; kwargs...)
end

function CellularSheaves.AsynchSheaves.plot_loss_curve!(plt, losses, label; kwargs...)
    return plot!(plt, losses, label=label, linewidth=2; kwargs...)
end

function CellularSheaves.AsynchSheaves.plot_log_loss_curve!(plt, losses, label; kwargs...)
    return plot!(plt, yscale=:log10, losses, label=label, linewidth=2; kwargs...)
end

function CellularSheaves.ControlSheaves.MultiAgentTracking.animate_tracking_xy(
    result::ScenarioResult;
    fps::Int=15,
    filename::AbstractString="./tracking_animation.gif",
    x_col::Int=1,
    y_col::Int=2,
    xlims::Tuple{<:Real,<:Real}=(-2.0, 4.0),
    ylims::Tuple{<:Real,<:Real}=(-2.0, 4.0),
)
    return SnapshotPlots.animate_tracking_xy(
        result;
        fps=fps,
        filename=filename,
        x_col=x_col,
        y_col=y_col,
        xlims=xlims,
        ylims=ylims,
    )
end
include("CellularSheavesPlots/src/LayeredEscortPlots.jl")

function CellularSheaves.ControlSheaves.DistributedLayeredControl.animate_layered_escort(res::CellularSheaves.ControlSheaves.DistributedLayeredControl.LayeredSimulationResult; fps::Int=10, frame_step::Int=4, filename::AbstractString="layered_escort.gif", kwargs...)
    steps = res.problem.steps
    anim = @animate for k in 1:frame_step:steps
        plot(res, k; kwargs...)
    end
    gif(anim, filename; fps=fps)
    return anim
end

function CellularSheaves.ControlSheaves.DistributedLayeredControl.animate_scenario5(
    res::CellularSheaves.ControlSheaves.DistributedLayeredControl.LayeredSimulationResult;
    fps::Int=10,
    frame_step::Int=4,
    filename::AbstractString="scenario5.gif",
    label_suffix::String="",
    kwargs...
)
    pre   = LayeredEscortPlots._scenario5_precompute(res)
    steps = res.problem.steps
    anim  = @animate for k in 1:frame_step:steps
        t_curr = k * res.problem.dt
        p1 = LayeredEscortPlots._scenario5_yz_panel(res, pre, k)
        p2 = LayeredEscortPlots._scenario5_tilt_panel(res, pre, k)
        p3 = LayeredEscortPlots._scenario5_error_panel(res, pre, k)
        plot(p1, p2, p3;
            layout     = (1, 3),
            size       = (1200, 420),
            plot_title = @sprintf("Scenario 5 [%s]  t = %.2f s", label_suffix, t_curr),
            kwargs...)
    end
    gif(anim, filename; fps=fps)
    return anim
end

function CellularSheaves.ControlSheaves.Layered.animate_comprehensive_escort(
    sim_data, qstar_history, target_history, time_grid, r1, r2, r3, filename;
    frame_step::Int=4, fps::Int=10, title_prefix::String=""
)
    return LayeredEscortPlots.animate_comprehensive_escort(
        sim_data, qstar_history, target_history, time_grid, r1, r2, r3, filename;
        frame_step=frame_step, fps=fps, title_prefix=title_prefix
    )
end

end # module CellularSheavesPlots
