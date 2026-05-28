module CellularSheavesPlots

using CellularSheaves
using CellularSheaves.AsynchSheaves
using CellularSheaves.ControlSheaves.MultiAgentTracking: ScenarioResult
using Plots

include("CellularSheavesPlots/src/Utils.jl")
include("CellularSheavesPlots/src/Style.jl")
include("CellularSheavesPlots/src/ScenarioPlots.jl")
include("CellularSheavesPlots/src/ErrorNormPlots.jl")
include("CellularSheavesPlots/src/Trajectory3D.jl")
include("CellularSheavesPlots/src/SnapshotPlots.jl")
include("CellularSheavesPlots/src/SpyPlots.jl")

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

end # module CellularSheavesPlots
