module ErrorNormPlots

using Plots

export rms, plot_error_norm

function rms(x::AbstractVector{<:Real})
    isempty(x) && return 0.0
    return sqrt(sum(abs2, x) / length(x))
end

function plot_error_norm(times::AbstractVector{<:Real}, errors::AbstractVector{<:Real};
    label::AbstractString="error",
    title::AbstractString="Tracking error norm",
    xlabel::AbstractString="Time (s)",
    ylabel::AbstractString="Error",
    yscale=:log10,
)
    p = plot(
        times,
        errors;
        label=label,
        title=title,
        xlabel=xlabel,
        ylabel=ylabel,
        yscale=yscale,
        linewidth=1.8,
    )
    return p
end

end # module ErrorNormPlots
