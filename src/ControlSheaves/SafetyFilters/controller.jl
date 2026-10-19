# The filter as an AbstractAgentController, so a filtered agent slots into the same
# layered-control tooling as an unconstrained one.

"""
    CBFCLFController(K, params, model) <: AbstractAgentController
    CBFCLFController(inner::LQRController, params, model)

An [`AbstractAgentController`](@ref) that wraps a nominal feedback gain with the local
CBF--CLF--SOCP filter, so that a safety-filtered agent slots into the same layered-control
tooling as [`LQRController`](@ref).

`model` is the [`ControlAffineModel`](@ref) seen by the filter. Filtering a physical
actuator command requires the configuration barrier to have relative degree one with respect
to that command; otherwise [`safety_filter`](@ref) raises a [`RelativeDegreeError`](@ref).
"""
struct CBFCLFController{S} <: AbstractAgentController
    K::Matrix{Float64}
    params::CBFCLFParams{S}
    model::ControlAffineModel
end

CBFCLFController(inner::LQRController, params::CBFCLFParams, model::ControlAffineModel) =
    CBFCLFController(inner.K, params, model)

"""
    step_agent!(w::AgentState, ctrl::CBFCLFController, qstar_target, others, dt;
                qstar_dot_target=nothing, obstacles=[])

Advance one agent by a single control period under the safety filter, and return
`(state, reference, result)` where `result` is the [`SafetyFilterResult`](@ref) for the step.

The reference is filtered by the agent's existing Tikhonov filter exactly as in the
unconstrained [`step_agent!`](@ref) methods; only the command is modified. `others` holds the
states of the neighbors visible to this agent.

When the program is infeasible or numerically unresolved the nominal command is applied and
`result.certified` is `false`, so a caller can count uncertified steps rather than silently
treating them as safe. Applications that must never apply an uncertified command should
inspect `result` and take their own fallback action.
"""
function step_agent!(w::AgentState, ctrl::CBFCLFController, qstar_target::AbstractVector,
                     others::AbstractVector{<:AbstractVector}, dt::Real;
                     qstar_dot_target = nothing,
                     obstacles::AbstractVector{<:AbstractVector{<:Real}} = Vector{Vector{Float64}}())
    nx = length(w.x)
    ref_dim = min(nx, length(qstar_target))
    qstar_local = zeros(nx)
    qstar_local[1:ref_dim] = qstar_target[1:ref_dim]

    # The velocity target enters `tikhonov_step!` in the velocity slots of the state, which
    # is the convention the unconstrained `step_agent!` methods use.
    velocity_target = zeros(nx)
    v_idxs = velocity_indices(w.dyn)
    if qstar_dot_target === nothing
        velocity_target .= qstar_local
    else
        v_dim = min(length(v_idxs), length(qstar_dot_target))
        for k in 1:v_dim
            velocity_target[v_idxs[k]] = qstar_dot_target[k]
        end
    end
    previous_ref = copy(w.filter.x)
    previous_vel = w.filter isa JointTikhonovFilter ? copy(w.filter.v) : previous_ref
    tikhonov_step!(w.filter, qstar_local, velocity_target, dt)
    x_ref = copy(w.filter.x)
    if qstar_dot_target !== nothing && w.filter isa JointTikhonovFilter
        x_ref[v_idxs] .= w.filter.v[v_idxs]
    end

    # The stability row needs the rate of change of the reference the agent actually tracks.
    # Two candidates differ, and the distinction matters: the *instantaneous* derivative
    # (target - x_ref)/eps describes the reference at the sample instant, while the *mean*
    # rate over the period describes how far it really moves before the next sample. They
    # agree only when the reference filter is slow compared with the control period. Here it
    # is not — the escort runs eps = 0.02 against dt = 0.05, so the filter very nearly
    # converges within one period and the two differ by a factor of four. Certifying against
    # the instantaneous rate then bounds a reference motion that is not the one the plant
    # experiences, and the dissipation inequality fails to hold over the step. The mean rate
    # is used, and it is exact rather than approximated because the filter update is a closed
    # form of the state before and after the step.
    reference_rate = (w.filter.x .- previous_ref) ./ dt
    if qstar_dot_target !== nothing && w.filter isa JointTikhonovFilter
        reference_rate[v_idxs] .= (w.filter.v[v_idxs] .- previous_vel[v_idxs]) ./ dt
    end

    u_nom = -ctrl.K * (w.x - x_ref)
    result = safety_filter(ctrl.params, ctrl.model, u_nom, w.x, x_ref;
                           ref_velocity = reference_rate, others = others,
                           obstacles = obstacles)
    u = result.certified ? result.command : u_nom
    w.x .= w.Ad * w.x .+ w.Bd * u
    return (copy(w.x), x_ref, result)
end
