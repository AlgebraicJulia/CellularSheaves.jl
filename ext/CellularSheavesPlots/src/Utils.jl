module Utils

export AGENT_COLORS, AGENT_STYLES, TARGET_COLORS
export pick_agent_color, pick_agent_style, pick_target_color
export safe_coord, target_path

const AGENT_COLORS = [:steelblue, :darkorange, :green, :crimson]
const AGENT_STYLES = [:solid, :dash, :dashdot, :dot]
const TARGET_COLORS = [:gray, :black, :darkgreen, :purple]

pick_agent_color(i::Int) = AGENT_COLORS[mod1(i, length(AGENT_COLORS))]
pick_agent_style(i::Int) = AGENT_STYLES[mod1(i, length(AGENT_STYLES))]
pick_target_color(i::Int) = TARGET_COLORS[mod1(i, length(TARGET_COLORS))]

function safe_coord(v::AbstractVector{<:Real}, idx::Int)
    return 1 <= idx <= length(v) ? Float64(v[idx]) : NaN
end

function target_path(traj::Vector{<:AbstractVector{<:Real}}, x_col::Int, y_col::Int, upto::Int)
    n = min(length(traj), upto)
    xs = Vector{Float64}(undef, n)
    ys = Vector{Float64}(undef, n)
    for i in 1:n
        xs[i] = safe_coord(traj[i], x_col)
        ys[i] = safe_coord(traj[i], y_col)
    end
    return xs, ys
end

end # module Utils
