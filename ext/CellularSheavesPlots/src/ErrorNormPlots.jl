module ErrorNormPlots

using Plots

export rms, plot_error_norm

struct ErrorNormData{T<:AbstractVector{<:Real},E<:AbstractVector{<:Real}}
    times::T
    errors::E
end

@recipe function f(data::ErrorNormData;
                   label="error",
                   title="Tracking error norm",
                   xlabel="Time (s)",
                   ylabel="Error",
                   yscale=:log10)
    seriestype := :path
    linewidth := 1.8
    label := label
    title := title
    xlabel := xlabel
    ylabel := ylabel
    yscale := yscale
    data.times, data.errors
end

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
    return plot(
        ErrorNormData(times, errors);
        label=label,
        title=title,
        xlabel=xlabel,
        ylabel=ylabel,
        yscale=yscale,
    )
end

end # module ErrorNormPlots
