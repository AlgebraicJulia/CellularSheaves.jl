# A recorded rollout of a filtered system, and the plotting entry points that consume it.
# The implementations live in the CellularSheavesPlots extension, following the same
# arrangement as `animate_layered_escort`.

"""
    FilterRollout

One recorded run of a filtered system, holding the quantities the layered architecture is
validated against: how close the agents came to one another and to any obstacle, how large a
command was demanded, how the Lyapunov function evolved, and when the barrier was binding.

`positions` is `steps × agents × dims`. `attitude`, when supplied, is `steps × agents` of
roll angles, and makes the animations draw each agent as a tilting airframe rather than a
dot — which is what is needed when the interesting motion is in attitude rather than in
position. `clearance` is `nothing` when the scenario has no
obstacle. Plot it with the recipes keyed by metric — `:top_down`, `:separation`,
`:clearance`, `:command`, `:lyapunov` — or animate it with [`animate_filter_rollout`](@ref).
"""
struct FilterRollout
    times::Vector{Float64}
    positions::Array{Float64, 3}
    separation::Vector{Float64}
    clearance::Union{Nothing, Vector{Float64}}
    command::Vector{Float64}
    lyapunov::Vector{Float64}
    slack::Vector{Float64}
    barrier_active::Vector{Bool}
    saturated::Vector{Bool}
    separation_floor::Float64
    keep_out::Float64
    command_cap::Float64
    attitude::Union{Nothing, Matrix{Float64}}
    label::String
end

function FilterRollout(times, positions; separation, clearance = nothing, command,
                       lyapunov, slack = zeros(length(times)),
                       barrier_active = falses(length(times)),
                       separation_floor = 0.0, keep_out = 0.0, command_cap = Inf,
                       attitude = nothing, label = "")
    @argcheck size(positions, 1) == length(times)
    @argcheck attitude === nothing || size(attitude, 1) == length(times)
    cmd = collect(command)
    # Where the cone is binding: that, and not the barrier, is the constraint that competes
    # with stability at an agent whose own program carries no barrier row.
    saturated = isfinite(command_cap) ? cmd .>= 0.999 * command_cap : falses(length(cmd))
    return FilterRollout(collect(times), positions, collect(separation),
                         clearance === nothing ? nothing : collect(clearance),
                         cmd, collect(lyapunov), collect(slack),
                         collect(barrier_active), saturated, separation_floor, keep_out,
                         command_cap, attitude === nothing ? nothing : Matrix{Float64}(attitude),
                         label)
end

"""
    barrier_window(rollout)
    saturation_window(rollout)

First and last time at which the barrier, respectively the actuator cone, was binding, or
`nothing` if it never was. Used to shade the interval over which that family, rather than
stability, was setting the command — and so the interval over which the Lyapunov function is
permitted to rise.
"""
function barrier_window(rollout::FilterRollout)
    idx = findall(rollout.barrier_active)
    isempty(idx) && return nothing
    return (rollout.times[first(idx)], rollout.times[last(idx)])
end

@doc (@doc barrier_window)
function saturation_window(rollout::FilterRollout)
    idx = findall(rollout.saturated)
    isempty(idx) && return nothing
    return (rollout.times[first(idx)], rollout.times[last(idx)])
end

"""
    relaxation_window(rollout; threshold = 1e-2)

First and last time at which the CLF relaxation was meaningfully open, or `nothing` if it
never was. This is the honest indicator of where stability was yielded: the certified
inequality is ``\\dot{V} \\le -c_3 V + \\delta``, so ``V`` may rise exactly where
``\\delta > 0`` and nowhere else. A quadratic penalty leaves ``\\delta`` at an ``O(1/p)``
floor rather than at zero, so a threshold distinguishes a genuine relaxation from that floor.
"""
function relaxation_window(rollout::FilterRollout; threshold::Real = 1e-2)
    idx = findall(>(threshold), rollout.slack)
    isempty(idx) && return nothing
    return (rollout.times[first(idx)], rollout.times[last(idx)])
end

"""
    animate_filter_rollout(unfiltered, filtered; filename, fps, frame_step, kwargs...)

Animate two rollouts side by side in the world frame. `mode` selects the semantics drawn,
so that a figure never asserts something its scenario does not enforce: `:safety` draws
keep-out regions and flags a breach, `:actuator` draws the commanded vector against its cap,
`:stability` draws the airframe attitude, and `:contended` draws both hard bounds together
with a raster of which constraint family is acting. Requires `Plots` to be loaded.
"""
function animate_filter_rollout end
