# A single bundle of knobs covering the common case: one stability term, one safety term,
# one actuator term. Convenient, and the shape most call sites want, but the three families
# are independent and can be built and exercised separately; see `terms.jl`.

"""
    CBFCLFParams(; kwargs...)

Tuning for the usual composition of one [`StabilityTerm`](@ref), one [`SafetyTerm`](@ref),
and one [`ActuatorTerm`](@ref). [`filter_terms`](@ref) shows exactly which term each keyword
reaches.

  - `d_safe`, `obstacle_radius`, `gamma`, `sense`, `cbf_disturbance_bound` — the safety term.
    `gamma` is the class-``\\mathcal{K}`` gain in ``\\dot{h} \\ge -\\gamma h``; larger is more
    permissive and intervenes later.
  - `braking_acceleration` — selects [`BrakingBarrier`](@ref) over [`DistanceBarrier`](@ref).
    Set it whenever the configuration is not directly actuated.
  - `clf_metric`, `clf_rate`, `clf_penalty`, `disturbance_bound` — the stability term.
  - `actuator_metric`, `actuator_cap` — the actuator cone ``\\|W u\\|_2 \\le u_{\\max}``.
  - `deadlock_bias` — composition-level symmetry breaking, disabled at `0.0`. Two agents
    approaching exactly head-on have the barrier gradient antiparallel to the CLF descent
    direction, so every admissible command makes no progress and the pair stalls at contact
    distance. A small positive value perturbs the nominal command tangentially before the
    solve; the result is then minimally invasive with respect to the perturbed command
    rather than the original, so it is off by default.
"""
Base.@kwdef struct CBFCLFParams{S}
    d_safe::Float64 = 0.6
    obstacle_radius::Float64 = 0.0
    gamma::Float64 = 6.0
    sense::Float64 = 3.0
    clf_rate::Float64 = 2.0
    clf_penalty::Float64 = 1.0
    disturbance_bound::Float64 = 0.0
    cbf_disturbance_bound::Float64 = 0.0
    actuator_cap::Float64 = Inf
    actuator_metric::Union{UniformScaling{Bool}, Matrix{Float64}} = I
    clf_metric::Union{UniformScaling{Bool}, Matrix{Float64}} = I
    deadlock_bias::Float64 = 0.0
    braking_acceleration::Union{Nothing, Float64} = nothing
    settings::S = safety_filter_settings()
end

"""
    filter_terms(params) -> (stability, safety, actuator)

The three independent terms a [`CBFCLFParams`](@ref) stands for. Use this to exercise one
family at a time: pass a single term to [`safety_filter`](@ref) and the others simply are
not there.
"""
function filter_terms(params::CBFCLFParams)
    @argcheck params.d_safe > 0
    barrier = params.braking_acceleration === nothing ?
        DistanceBarrier(params.d_safe) :
        BrakingBarrier(params.d_safe, params.braking_acceleration)
    stability = StabilityTerm(metric = params.clf_metric, rate = params.clf_rate,
                              penalty = params.clf_penalty,
                              disturbance_bound = params.disturbance_bound)
    safety = SafetyTerm(barrier; gamma = params.gamma, sense = params.sense,
                        obstacle_radius = params.obstacle_radius,
                        disturbance_bound = params.cbf_disturbance_bound)
    actuator = ActuatorTerm(metric = params.actuator_metric, bound = params.actuator_cap)
    return stability, safety, actuator
end

"""
    safety_filter(params, model, u_nom, x, ref; ref_velocity, others, obstacles,
                  obstacle_velocities)

Solve the local CBF--CLF--SOCP for one agent and return a [`SafetyFilterResult`](@ref).

This is the bundled form of the term-wise interface: it builds the three terms of
[`filter_terms`](@ref) and a [`FilterContext`](@ref), then composes them.
"""
function safety_filter(params::CBFCLFParams, model::ControlAffineModel,
                       u_nom::AbstractVector{<:Real}, x::AbstractVector{<:Real},
                       ref::AbstractVector{<:Real};
                       ref_velocity::AbstractVector{<:Real} = zeros(length(x)),
                       others::AbstractVector{<:AbstractVector{<:Real}} = Vector{Vector{Float64}}(),
                       obstacles::AbstractVector{<:AbstractVector{<:Real}} = Vector{Vector{Float64}}(),
                       obstacle_velocities = nothing)
    ctx = FilterContext(model, x, ref, length(u_nom); ref_velocity = ref_velocity,
                        others = others, obstacles = obstacles,
                        obstacle_velocities = obstacle_velocities)
    return safety_filter(filter_terms(params), ctx, u_nom;
                         settings = params.settings, deadlock_bias = params.deadlock_bias)
end

"""
    safety_filter_all(params, model, u_noms, xs, refs; ref_velocities, obstacles,
                      obstacle_velocities)

Apply [`safety_filter`](@ref) to every agent of an ensemble, passing each agent the states of
all the others. The solves are independent by construction: coordination has already been
discharged by the sheaf Laplacian solve that produced `refs`, so no communication is required
inside the execution loop.
"""
function safety_filter_all(params::CBFCLFParams, model::ControlAffineModel,
                           u_noms::AbstractVector, xs::AbstractVector, refs::AbstractVector;
                           ref_velocities = nothing,
                           obstacles::AbstractVector{<:AbstractVector{<:Real}} = Vector{Vector{Float64}}(),
                           obstacle_velocities = nothing)
    @argcheck length(u_noms) == length(xs) == length(refs)
    velocities = ref_velocities === nothing ? [zeros(length(x)) for x in xs] : ref_velocities
    @argcheck length(velocities) == length(xs)
    return [safety_filter(params, model, u_noms[i], xs[i], refs[i];
                          ref_velocity = velocities[i],
                          others = [xs[j] for j in eachindex(xs) if j != i],
                          obstacles = obstacles,
                          obstacle_velocities = obstacle_velocities)
            for i in eachindex(xs)]
end

"""
    barrier_values(model, x, others; d_safe, sense = Inf)

Pairwise [`DistanceBarrier`](@ref) values against every entry of `others` within `sense`.
Shared by the filter and by post-hoc verification, so that a rollout is checked against
exactly the quantity the filter constrained.
"""
function barrier_values(model::ControlAffineModel, x::AbstractVector{<:Real},
                        others::AbstractVector{<:AbstractVector{<:Real}};
                        d_safe::Real, sense::Real = Inf)
    @argcheck d_safe > 0
    @argcheck sense > 0
    px = model.position * Vector{Float64}(x)
    values = Float64[]
    for xj in others
        displacement = px - model.position * Vector{Float64}(xj)
        radius2 = dot(displacement, displacement)
        radius2 > sense^2 && continue
        push!(values, radius2 - d_safe^2)
    end
    return values
end
