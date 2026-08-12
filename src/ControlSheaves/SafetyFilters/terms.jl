# The three constraint families of the filter, each independent of the others. A term reads
# the context and emits rows; it never inspects another term. Anything that genuinely
# depends on two families — dropping rows the actuator cone makes unreachable, breaking a
# head-on deadlock — belongs to the composer in `filter.jl`, not to a term.

"""
    FilterContext(model, x, ref; ref_velocity, others, obstacles, obstacle_velocities)

Everything a filter term may read about the current instant: the agent's model and state,
the reference it is tracking, and the neighbours and obstacles around it. The drift, input
map, and their configuration and velocity projections are evaluated once here and shared by
every term.
"""
struct FilterContext
    model::ControlAffineModel
    x::Vector{Float64}
    ref::Vector{Float64}
    ref_velocity::Vector{Float64}
    others::Vector{Vector{Float64}}
    obstacles::Vector{Vector{Float64}}
    obstacle_velocities::Union{Nothing, Vector{Vector{Float64}}}
    nu::Int
    fx::Vector{Float64}
    gx::Matrix{Float64}
    px::Vector{Float64}
    Pf::Vector{Float64}
    Pg::Matrix{Float64}
    vx::Vector{Float64}
    Pvf::Vector{Float64}
    Pvg::Matrix{Float64}
end

function FilterContext(model::ControlAffineModel, x::AbstractVector{<:Real},
                       ref::AbstractVector{<:Real}, nu::Integer;
                       ref_velocity::AbstractVector{<:Real} = zeros(length(x)),
                       others = Vector{Vector{Float64}}(),
                       obstacles = Vector{Vector{Float64}}(),
                       obstacle_velocities = nothing)
    nx = state_dimension(model)
    @argcheck length(x) == nx
    @argcheck length(ref) == nx
    @argcheck length(ref_velocity) == nx
    @argcheck obstacle_velocities === nothing ||
              length(obstacle_velocities) == length(obstacles)

    state = Vector{Float64}(x)
    fx = Vector{Float64}(model.drift(state))
    gx = Matrix{Float64}(model.input(state))
    @argcheck length(fx) == nx
    @argcheck size(gx) == (nx, nu)

    P = model.position
    has_velocity = model.velocity !== nothing
    Pv = has_velocity ? model.velocity : zeros(0, nx)
    return FilterContext(model, state, Vector{Float64}(ref), Vector{Float64}(ref_velocity),
                         [Vector{Float64}(o) for o in others],
                         [Vector{Float64}(o) for o in obstacles],
                         obstacle_velocities === nothing ? nothing :
                             [Vector{Float64}(v) for v in obstacle_velocities],
                         Int(nu), fx, gx, P * state, P * fx, P * gx,
                         Pv * state, Pv * fx, Pv * gx)
end

"""
    AbstractFilterTerm

Supertype of the constraint families. A term implements [`constraints`](@ref), returning the
rows it contributes, an optional cone, and whatever diagnostics it can report.
"""
abstract type AbstractFilterTerm end

"""
    constraints(term, ctx) -> (rows, cone, diagnostics)

Rows are [`LinearRow`](@ref)s in the ``a^\\top u \\ge b`` convention. `cone`, when not
`nothing`, is a named tuple `(metric, bound)` describing ``\\|W u\\|_2 \\le u_{\\max}``.
"""
function constraints end

"""
    StabilityTerm(; metric = I, rate = 2.0, penalty = 1.0, disturbance_bound = 0.0)

The ISS-CLF family. With ``V = \\tfrac{1}{2} e^\\top M e`` for ``e = x - r``, it contributes
the single relaxable row

```math
(g^\\top M e)^\\top u \\le -c_3 V - e^\\top M f + e^\\top M \\dot{r} - \\|M e\\| \\Delta + \\delta .
```

The relaxation ``\\delta \\ge 0`` is what keeps the composed program feasible when the
barrier or the actuator cone forbids the commanded decay rate. Because the penalty is
quadratic it is inexact: ``\\delta`` settles at ``O(1/p)`` rather than zero, so the
certified inequality is ``\\dot{V} \\le -c_3 V + O(1/p)``.
"""
Base.@kwdef struct StabilityTerm <: AbstractFilterTerm
    metric::Union{UniformScaling{Bool}, Matrix{Float64}} = I
    rate::Float64 = 2.0
    penalty::Float64 = 1.0
    disturbance_bound::Float64 = 0.0
end

function constraints(term::StabilityTerm, ctx::FilterContext)
    @argcheck term.rate > 0
    @argcheck term.penalty > 0
    @argcheck term.disturbance_bound >= 0
    nx = length(ctx.x)
    M = term.metric isa UniformScaling ? Matrix{Float64}(term.metric, nx, nx) : term.metric
    @argcheck size(M) == (nx, nx)

    e = ctx.x - ctx.ref
    Me = M * e
    V = 0.5 * dot(e, Me)
    coefficients = ctx.gx' * Me
    bound = -term.rate * V - dot(Me, ctx.fx) + dot(Me, ctx.ref_velocity) -
            norm(Me) * term.disturbance_bound
    # a'u >= b - delta form: the inequality above is an upper bound on (g'Me)'u.
    rows = LinearRow[LinearRow(-coefficients, -bound; penalty = term.penalty)]
    return (rows = rows, cone = nothing, cones = ConeRow[], diagnostics = (lyapunov = V,))
end

"""
    SafetyTerm(barrier; gamma = 6.0, sense = 3.0, obstacle_radius = 0.0,
               disturbance_bound = 0.0)

The control-barrier family. It contributes one hard row per neighbour and per obstacle
within `sense`, each enforcing this agent's share of ``\\dot{h} \\ge -\\gamma h`` for the
given [`AbstractBarrier`](@ref).

Responsibility for a pair of agents is split evenly (Wang, Ames, and Egerstedt, T-RO 2017):
each enforces half, and the two halves sum to the pairwise condition without any
communication in the execution loop. An obstacle enforces nothing, so the agent carries the
full ``\\gamma`` and the obstacle's own velocity is carried into the row.

These rows are never relaxed. The composed program can therefore be infeasible, which
[`safety_filter`](@ref) reports rather than papering over.
"""
Base.@kwdef struct SafetyTerm{B <: AbstractBarrier} <: AbstractFilterTerm
    barrier::B
    gamma::Float64 = 6.0
    sense::Float64 = 3.0
    obstacle_radius::Float64 = 0.0
    disturbance_bound::Float64 = 0.0
end

SafetyTerm(barrier::AbstractBarrier; kwargs...) = SafetyTerm(; barrier = barrier, kwargs...)

function constraints(term::SafetyTerm, ctx::FilterContext)
    @argcheck term.gamma > 0
    @argcheck term.sense > 0
    @argcheck term.obstacle_radius >= 0
    @argcheck term.disturbance_bound >= 0
    if term.barrier isa BrakingBarrier && ctx.model.velocity === nothing
        throw(ArgumentError("the braking barrier needs a velocity selector: build the " *
                            "model from an AbstractAgentDynamics or with " *
                            "double_integrator_model"))
    end

    P = ctx.model.position
    Pv = ctx.model.velocity
    rows = LinearRow[]
    values = Float64[]

    for (targets, share, extra, drifts) in
            ((ctx.others, 2.0, 0.0, nothing),
             (ctx.obstacles, 1.0, term.obstacle_radius, ctx.obstacle_velocities))
        for (index, target) in enumerate(targets)
            @argcheck length(target) == length(ctx.x)
            dp = ctx.px - P * target
            dot(dp, dp) > term.sense^2 && continue

            target_rate = drifts === nothing ? zeros(length(dp)) : P * drifts[index]
            # Relative velocity: against a neighbour it comes from that agent's own state,
            # against an obstacle from the supplied velocity.
            dv = if Pv === nothing
                Float64[]
            elseif drifts === nothing
                ctx.vx - Pv * target
            else
                ctx.vx - target_rate
            end

            h, coefficients, bound = barrier_row(term.barrier, ctx, dp, dv, target_rate,
                                                 share, extra, term.gamma,
                                                 term.disturbance_bound)
            # A vanishing row means the barrier has relative degree greater than one with
            # respect to the filtered input: no choice of u changes hdot instantaneously,
            # so the row is vacuous and the program is degenerate rather than infeasible.
            if norm(coefficients) <= _RELATIVE_DEGREE_TOL * max(1.0, norm(dp))
                throw(RelativeDegreeError(term.barrier isa BrakingBarrier ?
                                          norm(ctx.Pvg) : norm(ctx.Pg)))
            end
            push!(values, h)
            push!(rows, LinearRow(coefficients, bound))
        end
    end
    return (rows = rows, cone = nothing, cones = ConeRow[], diagnostics = (barriers = values,))
end

"""
    ActuatorTerm(; metric = I, bound = Inf)

The actuator family: the second-order cone ``\\|W u\\|_2 \\le u_{\\max}``. It contributes no
linear rows, only the cone, and is inert when `bound` is infinite.
"""
Base.@kwdef struct ActuatorTerm <: AbstractFilterTerm
    metric::Union{UniformScaling{Bool}, Matrix{Float64}} = I
    bound::Float64 = Inf
end

function constraints(term::ActuatorTerm, ctx::FilterContext)
    @argcheck term.bound > 0
    W = term.metric isa UniformScaling ? Matrix{Float64}(term.metric, ctx.nu, ctx.nu) :
        term.metric
    @argcheck size(W, 2) == ctx.nu
    cone = isfinite(term.bound) ? (metric = W, bound = term.bound) : nothing
    return (rows = LinearRow[], cone = cone, cones = ConeRow[], diagnostics = NamedTuple())
end

"""
    DiscreteStabilityTerm(; Ad, Bd, metric, dt, rate)

The stability family enforced on the *step* rather than on the derivative. In terms of the
tracking error ``e = x - r``, held over one period,

```math
V(A_d e + B_d u) \\le (1 - c_3 \\Delta t) V(e), \\qquad
V(z) = \\tfrac{1}{2} z^\\top M z .
```

Writing ``M = L^\\top L``, this is exactly the second-order cone

```math
\\|L(A_d e + B_d u)\\|_2 \\le \\sqrt{1 - c_3 \\Delta t}\\, \\|L e\\|_2 ,
```

which [`ConeRow`](@ref) expresses with ``W = L B_d`` and offset ``L A_d e``.

[`StabilityTerm`](@ref) instead imposes the continuous inequality ``\\dot{V} \\le -c_3 V``
at the sample instants. The two agree as ``\\Delta t \\to 0``, but they are not the same
constraint: satisfying the derivative condition at a sample says nothing about where an
explicit step of length ``\\Delta t`` actually lands, and the discrepancy grows with the
size of the state change over a period. Against a benign nominal command that difference is
invisible; against a nominal that is itself destabilising it dominates, and the continuous
form can certify every step of a trajectory that diverges. Prefer this term whenever the
nominal command is not already known to be stabilising: a learned policy, a transferred
gain, a delayed or saturated loop.

The condition is hard and carries no relaxation, so it can be infeasible — no command may
contract ``V`` by the demanded factor within one period. Pair it with a modest `rate`, and
treat an uncertified result as the report that it is. An actuator limit is a separate
family: compose this term with an [`ActuatorTerm`](@ref) rather than bounding the command
here.
"""
Base.@kwdef struct DiscreteStabilityTerm <: AbstractFilterTerm
    Ad::Matrix{Float64}
    Bd::Matrix{Float64}
    metric::Matrix{Float64}
    dt::Float64
    rate::Float64
end

function constraints(term::DiscreteStabilityTerm, ctx::FilterContext)
    @argcheck term.dt > 0
    @argcheck term.rate > 0
    @argcheck term.rate * term.dt < 1
    nx = length(ctx.x)
    @argcheck size(term.Ad) == (nx, nx)
    @argcheck size(term.Bd) == (nx, ctx.nu)
    @argcheck size(term.metric) == (nx, nx)

    L = cholesky(Symmetric(term.metric)).U
    e = ctx.x - ctx.ref
    V = 0.5 * dot(e, term.metric * e)
    contraction = sqrt(1 - term.rate * term.dt)
    # The successor error is A_d e + B_d u when the reference is held over the period.
    bound = contraction * norm(L * e)
    cone = ConeRow(L * term.Bd, max(bound, 1e-12); offset = L * (term.Ad * e))
    return (rows = LinearRow[], cone = nothing, cones = ConeRow[cone],
            diagnostics = (lyapunov = V,))
end
