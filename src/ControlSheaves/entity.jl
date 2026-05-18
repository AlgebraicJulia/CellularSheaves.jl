# src/core/entity.jl

abstract type AbstractEntity{D <: AbstractDynamics} end

mutable struct Entity{D <: AbstractDynamics}
    id::String
    num_states::Int
    time_step_delta::Float64
    positions::Matrix{Float64}
    velocities::Matrix{Float64}
    control_output::Vector{Float64}
    synchronization_error::Vector{Float64}
    offsets::Vector{Float64}
    neighbors::Vector{AbstractEntity}
end

mutable struct Target{D <: AbstractDynamics} <: AbstractEntity{D}
    entity::Entity{D}
    k1::Float64
    is_centroid::Bool
    tracking_type::String
end

mutable struct Agent{D <: AbstractDynamics} <: AbstractEntity{D}
    entity::Entity{D}
    k1::Float64
    pin_row::Vector{Float64}
    targets::Vector{<:AbstractEntity}
end

function Entity{D}(initial_position::Vector{Float64}, time_steps::Int, id::String,
                   num_states::Int, time_step_delta::Float64) where {D <: AbstractDynamics}
    positions = zeros(Float64, num_states, time_steps)
    velocities = zeros(Float64, num_states, time_steps)
    positions[:, 1] .= initial_position
    control_output = zeros(Float64, num_states)
    synchronization_error = zeros(Float64, num_states)
    offsets = zeros(Float64, num_states)
    neighbors = AbstractEntity[]
    return Entity{D}(id, num_states, time_step_delta, positions, velocities,
                     control_output, synchronization_error, offsets, neighbors)
end

function Target{D}(initial_position::Vector{Float64}, time_steps::Int, id::String,
                   num_states::Int, time_step_delta::Float64,
                   targets_proportional_gain::Float64;
                   is_centroid::Bool=false,
                   tracking_type::AbstractString="") where {D <: AbstractDynamics}
    entity = Entity{D}(initial_position, time_steps, id, num_states, time_step_delta)
    return Target{D}(entity, targets_proportional_gain, is_centroid, String(tracking_type))
end

function Agent{D}(initial_position::Vector{Float64}, time_steps::Int, id::String,
                  num_states::Int, time_step_delta::Float64,
                  agents_proportional_gain::Float64,
                  targets::Vector{<:AbstractEntity}) where {D <: AbstractDynamics}
    entity = Entity{D}(initial_position, time_steps, id, num_states, time_step_delta)
    pin_row = zeros(Float64, length(targets))
    return Agent{D}(entity, agents_proportional_gain, pin_row, targets)
end

# Delegate selected field accesses to the inner Entity.
function Base.getproperty(entity::AbstractEntity, name::Symbol)
    if name in (:id, :num_states, :time_step_delta, :positions, :velocities,
                :control_output, :synchronization_error, :offsets, :neighbors)
        return getfield(getfield(entity, :entity), name)
    else
        return getfield(entity, name)
    end
end

function run_dynamics(entity::AbstractEntity{D}, input) where {D <: AbstractDynamics}
    return run_dynamics(D, input)
end

function run_dynamics!(entity::AbstractEntity{D}, output, input) where {D <: AbstractDynamics}
    return run_dynamics!(D, output, input)
end

# edit function so save computations for basic dynamics NoDynamics, CurrentAgentDynamics, CurrentTargetDynamics
# function update_dynamics!(entity::AbstractEntity, step::Int)
#     vel      = view(entity.velocities, :, step)
#     pos_prev = view(entity.positions,  :, step - 1)

#     run_dynamics!(entity, vel, pos_prev)
#     vel .+= entity.control_output

#     result = integrate_step(pos_prev, step - 2, entity.time_step_delta) do _t, pos
#         return run_dynamics(entity, pos) + entity.control_output
#     end

#     entity.positions[:, step] .= result
#     return entity
# end

function _direct_update_dynamics!(entity::AbstractEntity, step::Int)
    vel      = view(entity.velocities, :, step)
    pos_prev = view(entity.positions,  :, step - 1)

    run_dynamics!(entity, vel, pos_prev)
    vel .+= entity.control_output

    entity.positions[:, step] .= pos_prev .+ entity.time_step_delta .* vel
    return entity
end

function update_dynamics!(entity::AbstractEntity{NoDynamics}, step::Int)
    return _direct_update_dynamics!(entity, step)
end

function update_dynamics!(entity::AbstractEntity{CurrentAgentDynamics}, step::Int)
    return _direct_update_dynamics!(entity, step)
end

function update_dynamics!(entity::AbstractEntity{CurrentTargetDynamics}, step::Int)
    return _direct_update_dynamics!(entity, step)
end

function update_dynamics!(entity::AbstractEntity{D}, step::Int) where {D <: AbstractDynamics}
    vel      = view(entity.velocities, :, step)
    pos_prev = view(entity.positions,  :, step - 1)

    run_dynamics!(entity, vel, pos_prev)
    vel .+= entity.control_output

    result = integrate_step(pos_prev, step - 2, entity.time_step_delta) do _t, pos
        return run_dynamics(entity, pos) + entity.control_output
    end

    entity.positions[:, step] .= result
    return entity
end

# ─────────────────────────────────────────────────────────────────────────────
# BASELINE control laws (explicit neighbor-sum; unchanged from original)
# ─────────────────────────────────────────────────────────────────────────────

function compute_control_output!(agent::Agent, step::Int)
    position = view(agent.positions, :, step - 1)

    neighborhood_consensus_term = zeros(Float64, agent.num_states)
    self_total = position .+ agent.offsets

    for neighbor in agent.neighbors
        if neighbor isa Agent
            neighbor_pos   = view(neighbor.positions, :, step - 1)
            neighbor_total = neighbor_pos .+ neighbor.offsets
            neighborhood_consensus_term .+= neighbor_total .- self_total
        end
    end

    target_consensus_term = zeros(Float64, agent.num_states)
    if !isempty(agent.targets) && !isempty(agent.pin_row)
        for (index, weight) in enumerate(agent.pin_row)
            if weight != 0.0
                tgt     = agent.targets[index]
                tgt_pos = view(tgt.positions, :, step - 1)
                target_consensus_term .+= weight .* (tgt_pos .- self_total)
            end
        end
    end

    agent.synchronization_error .= target_consensus_term .+ neighborhood_consensus_term
    agent.control_output         .= agent.k1 .* agent.synchronization_error

    if agent.num_states >= 3
        agent.control_output[3] = 0.0
    end
    return agent
end

function compute_control_output!(target::Target, step::Int)
    time = (step - 1) * target.time_step_delta
    desired_velocity = zeros(Float64, target.num_states)

    if target.is_centroid && target.tracking_type == "f8_dynamics"
        desired_velocity .= f8_dynamics(time)
    end

    neighborhood_consensus_term = zeros(Float64, target.num_states)
    self_total = view(target.positions, :, step - 1) .+ target.offsets

    for neighbor in target.neighbors
        if neighbor isa Target
            neighbor_total = view(neighbor.positions, :, step - 1) .+ neighbor.offsets
            neighborhood_consensus_term .+= neighbor_total .- self_total
        end
    end

    target.synchronization_error .= neighborhood_consensus_term
    target.control_output         .= target.k1 .* target.synchronization_error ./ 100.0 .+ desired_velocity
    return target
end

# ─────────────────────────────────────────────────────────────────────────────
# SHEAF-BASED control laws
# ─────────────────────────────────────────────────────────────────────────────

function compute_control_output_sheaf!(agent::Agent, i::Int,
                                        control_block::Matrix{Float64}, step::Int)
    agent.control_output .= control_block[:, i]
    if agent.num_states >= 3
        agent.control_output[3] = 0.0
    end
    return agent
end

"""
    compute_control_output_sheaf!(target, j, control_block, step)

Write the j-th column of `control_block` (nd×n_targets) into
`target.control_output` and add the f8 desired velocity for centroid targets.
"""
function compute_control_output_sheaf!(target::Target, j::Int,
                                        control_block::Matrix{Float64}, step::Int)
    time = (step - 1) * target.time_step_delta
    desired_velocity = zeros(Float64, target.num_states)
    if target.is_centroid && target.tracking_type == "f8_dynamics"
        desired_velocity .= f8_dynamics(time)
    end
    target.control_output .= control_block[:, j] ./ 100.0 .+ desired_velocity
    return target
end
