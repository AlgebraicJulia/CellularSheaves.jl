module CellularSheavesPlots

using CellularSheaves
using CellularSheaves.AsynchSheaves
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

end # module CellularSheavesPlots
