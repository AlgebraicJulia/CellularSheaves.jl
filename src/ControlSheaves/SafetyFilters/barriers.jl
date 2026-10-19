# Barrier functions. Each type answers one question: given a displacement to another object,
# what is the barrier value, and what linear row in the command enforces
# ``\dot{h} \ge -\gamma h`` for this agent's share of it?

const _RELATIVE_DEGREE_TOL = 1e-10

# Near contact the derivative of sqrt(4a(r-d)) diverges. The denominator is floored so the
# row stays finite; the effect is to cap how hard the barrier can demand deceleration in the
# last millimetre before the keep-out radius.
const _BRAKING_FLOOR = 1e-3

"""
    AbstractBarrier

Supertype of the collision barriers. A barrier supplies its value ``h`` and the row
``a^\\top u \\ge b`` that enforces ``\\dot{h} \\ge -\\gamma h`` for one agent, through
[`barrier_row`](@ref).
"""
abstract type AbstractBarrier end

"""
    DistanceBarrier(d_safe)

``h = \\|P x_i - P x_j\\|^2 - d_{\\text{safe}}^2``.

Relative degree one exactly when the configuration is directly actuated, i.e. when
``P g \\ne 0``. For anything else — a double integrator, a quadrotor — use
[`BrakingBarrier`](@ref); see [`RelativeDegreeError`](@ref).
"""
struct DistanceBarrier <: AbstractBarrier
    d_safe::Float64
end

"""
    BrakingBarrier(d_safe, a_max)

``h = \\sqrt{4 a_{\\max} (\\|\\Delta p\\| - d_{\\text{safe}})} + \\Delta p^\\top \\Delta v / \\|\\Delta p\\|``,
the safety barrier certificate of Wang, Ames, and Egerstedt (T-RO 2017): the pair is safe
while the closing speed stays within what the two agents can brake away before reaching the
keep-out radius.

Because ``h`` involves relative velocity its derivative reaches an acceleration input
directly, so the barrier has relative degree one wherever velocity is actuated — one better
than [`DistanceBarrier`](@ref). The price is that each agent must sense a neighbour's
velocity as well as its position, and that an obstacle is assumed to hold its velocity over
the control period.
"""
struct BrakingBarrier <: AbstractBarrier
    d_safe::Float64
    a_max::Float64
end

separation(b::AbstractBarrier) = b.d_safe

"""
    barrier_row(barrier, ctx, dp, dv, target_rate, share, extra_radius, gamma, disturbance)

Return `(h, coefficients, bound)`, the barrier value and the row
``\\text{coefficients}^\\top u \\ge \\text{bound}`` enforcing this agent's share of
``\\dot{h} \\ge -\\gamma h``.

`share` is `2` for a neighbouring agent, which enforces the other half of the pair itself,
and `1` against an obstacle, which does not. `extra_radius` inflates the keep-out beyond
`d_safe`, and `disturbance` is the drift budget this agent reserves for itself.
"""
function barrier_row(barrier::DistanceBarrier, ctx, dp::Vector{Float64},
                     dv::Vector{Float64}, target_rate::Vector{Float64},
                     share::Float64, extra_radius::Float64, gamma::Float64,
                     disturbance::Float64)
    inflation = barrier.d_safe + extra_radius
    h = dot(dp, dp) - inflation^2
    coefficients = 2.0 .* (ctx.Pg' * dp)
    bound = -(gamma / share) * h - 2.0 * dot(dp, ctx.Pf) + 2.0 * dot(dp, target_rate) +
            2.0 * norm(dp) * disturbance
    return h, coefficients, bound
end

function barrier_row(barrier::BrakingBarrier, ctx, dp::Vector{Float64},
                     dv::Vector{Float64}, target_rate::Vector{Float64},
                     share::Float64, extra_radius::Float64, gamma::Float64,
                     disturbance::Float64)
    inflation = barrier.d_safe + extra_radius
    r = norm(dp)
    @argcheck r > 0
    n = dp ./ r
    rdot = dot(n, dv)
    s = sqrt(4.0 * barrier.a_max * max(r - inflation, 0.0))
    h = s + rdot
    # d/dt of the two terms of h, holding the input aside.
    d_sqrt = 2.0 * barrier.a_max * rdot / max(s, _BRAKING_FLOOR)
    curvature = (dot(dv, dv) - rdot^2) / r
    coefficients = ctx.Pvg' * n
    bound = -(gamma / share) * h - (d_sqrt + curvature) / share - dot(n, ctx.Pvf) +
            disturbance
    return h, coefficients, bound
end

"""
    barrier_value(barrier, ctx, dp, dv, extra_radius)

The barrier value alone, without building a row. Used to report margins for targets whose
row was dropped, and by post-hoc verification of a rollout.
"""
function barrier_value(barrier::DistanceBarrier, ctx, dp::Vector{Float64},
                       dv::Vector{Float64}, extra_radius::Float64)
    return dot(dp, dp) - (barrier.d_safe + extra_radius)^2
end

function barrier_value(barrier::BrakingBarrier, ctx, dp::Vector{Float64},
                       dv::Vector{Float64}, extra_radius::Float64)
    r = norm(dp)
    @argcheck r > 0
    n = dp ./ r
    return sqrt(4.0 * barrier.a_max * max(r - barrier.d_safe - extra_radius, 0.0)) +
           dot(n, dv)
end
