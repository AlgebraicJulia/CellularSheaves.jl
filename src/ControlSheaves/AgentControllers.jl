module AgentControllers

using LinearAlgebra
using ArgCheck
using ...NetworkSheaves.TrajectorySheaves: continuous_to_discrete_zoh
using ..Tikhonov: AbstractTikhonovFilter, TikhonovFilter, JointTikhonovFilter, tikhonov_step!

export AbstractAgentDynamics, QuadrotorDynamics, PlanarQuadrotorDynamics,
       AbstractAgentController, LQRController,
       AbstractAgentState, AgentState,
       solve_dare, step_agent!,
       position_indices, velocity_indices, state_dim, initial_state

abstract type AbstractAgentDynamics end

"""
    QuadrotorDynamics <: AbstractAgentDynamics

10D Quadrotor dynamics linearized around hover. 
State x = [x, y, z, phi, theta, x_dot, y_dot, z_dot, phi_dot, theta_dot]
"""
Base.@kwdef struct QuadrotorDynamics <: AbstractAgentDynamics
    g::Float64 = 9.81
    m::Float64 = 0.5
    Ixx::Float64 = 0.01
    Iyy::Float64 = 0.01
end

function continuous_matrices(dyn::QuadrotorDynamics)
    Ac = zeros(10, 10)
    Ac[1, 6] = 1.0
    Ac[2, 7] = 1.0
    Ac[3, 8] = 1.0
    Ac[4, 9] = 1.0
    Ac[5, 10] = 1.0
    Ac[6, 5] = dyn.g
    Ac[7, 4] = -dyn.g

    Bc = zeros(10, 3)
    Bc[8, 1] = 1.0 / dyn.m       # net thrust deviation
    Bc[9, 2] = 1.0 / dyn.Ixx     # roll moment
    Bc[10, 3] = 1.0 / dyn.Iyy    # pitch moment
    
    return Ac, Bc
end

"""
    PlanarQuadrotorDynamics <: AbstractAgentDynamics

6D Planar Quadrotor dynamics.
State x = [y, z, theta, y_dot, z_dot, theta_dot]
"""
Base.@kwdef struct PlanarQuadrotorDynamics <: AbstractAgentDynamics
    g::Float64 = 9.81
    m::Float64 = 0.5
    I_quad::Float64 = 0.01
    ell::Float64 = 0.25
end

function continuous_matrices(dyn::PlanarQuadrotorDynamics)
    Ac = zeros(6, 6)
    Ac[1, 4] = 1.0
    Ac[2, 5] = 1.0
    Ac[3, 6] = 1.0
    Ac[4, 3] = -dyn.g
    
    Bc = zeros(6, 2)
    Bc[5, 1] = 1.0 / dyn.m
    Bc[5, 2] = 1.0 / dyn.m
    Bc[6, 1] = dyn.ell / (2*dyn.I_quad)
    Bc[6, 2] = -dyn.ell / (2*dyn.I_quad)
    return Ac, Bc
end

function discrete_matrices(dyn::AbstractAgentDynamics, dt::Float64)
    Ac, Bc = continuous_matrices(dyn)
    continuous_to_discrete_zoh(Ac, Bc, dt)
end

"""
    velocity_indices(dyn::QuadrotorDynamics) -> 6:8
    velocity_indices(dyn::PlanarQuadrotorDynamics) -> 4:5

Returns the velocity state indices in the full state vector for the given agent dynamics.
"""
velocity_indices(dyn::QuadrotorDynamics) = 6:8
velocity_indices(dyn::PlanarQuadrotorDynamics) = 4:5

"""
    position_indices(dyn::QuadrotorDynamics) -> 1:3
    position_indices(dyn::PlanarQuadrotorDynamics) -> 1:2

Returns the position state indices in the full state vector for the given agent dynamics.
"""
position_indices(dyn::QuadrotorDynamics) = 1:3
position_indices(dyn::PlanarQuadrotorDynamics) = 1:2

"""
    state_dim(dyn::QuadrotorDynamics) -> 10
    state_dim(dyn::PlanarQuadrotorDynamics) -> 6

Returns the full state vector dimension for the given agent dynamics.
"""
state_dim(dyn::QuadrotorDynamics) = 10
state_dim(dyn::PlanarQuadrotorDynamics) = 6

"""
    initial_state(dyn::AbstractAgentDynamics, position::AbstractVector)
    initial_state(dyn::AbstractAgentDynamics, position::AbstractVector, velocity::AbstractVector)
    initial_state(dyn::AbstractAgentDynamics, position::AbstractVector, velocity::AbstractVector, acceleration::AbstractVector)

Constructs a full state vector for dynamics model `dyn` starting at world `position`
(length `position_indices(dyn)`), with `velocity` (length `velocity_indices(dyn)`,
defaults to rest) and `acceleration` (defaults to steady/zero). When `acceleration` is
supplied, the initial attitude is trimmed using the same flat-output feedforward formula
used by `step_agent!`, so agents seeded with the trajectory's initial acceleration start
already banked into the turn instead of level.
"""
function initial_state(dyn::AbstractAgentDynamics, position::AbstractVector)
    initial_state(dyn, position, zeros(length(velocity_indices(dyn))))
end

function initial_state(dyn::AbstractAgentDynamics, position::AbstractVector, velocity::AbstractVector)
    initial_state(dyn, position, velocity, zeros(length(position)))
end

function initial_state(dyn::QuadrotorDynamics, position::AbstractVector, velocity::AbstractVector, acceleration::AbstractVector)
    @argcheck length(position) == length(position_indices(dyn))
    @argcheck length(velocity) == length(velocity_indices(dyn))
    x = zeros(state_dim(dyn))
    x[position_indices(dyn)] .= position
    x[velocity_indices(dyn)] .= velocity
    if length(acceleration) >= 2
        g = dyn.g
        ax, ay = acceleration[1], acceleration[2]
        x[4] = -ay / g
        x[5] =  ax / g
    end
    return x
end

function initial_state(dyn::PlanarQuadrotorDynamics, position::AbstractVector, velocity::AbstractVector, acceleration::AbstractVector)
    @argcheck length(position) == length(position_indices(dyn))
    @argcheck length(velocity) == length(velocity_indices(dyn))
    x = zeros(state_dim(dyn))
    x[position_indices(dyn)] .= position
    x[velocity_indices(dyn)] .= velocity
    if length(acceleration) >= 1
        g = dyn.g
        ay = acceleration[1]
        x[3] = -ay / g
    end
    return x
end

"""
    solve_dare(A, B, Q, R)

Solves the Discrete Algebraic Riccati Equation for LQR gain.
"""
function solve_dare(A::AbstractMatrix, B::AbstractMatrix, Q::AbstractMatrix, R::AbstractMatrix)
    P = copy(Q)
    for i in 1:200
        P_next = A' * P * A - (A' * P * B) * ((R + B' * P * B) \ (B' * P * A)) + Q
        if norm(P_next - P) < 1e-6
            break
        end
        P = P_next
    end
    return (R + B' * P * B) \ (B' * P * A)
end

abstract type AbstractAgentController end

struct LQRController <: AbstractAgentController
    K::Matrix{Float64}
end

function LQRController(dyn::AbstractAgentDynamics, dt::Float64, Q::AbstractMatrix, R::AbstractMatrix)
    Ad, Bd = discrete_matrices(dyn, dt)
    K = solve_dare(Ad, Bd, Q, R)
    LQRController(K)
end

abstract type AbstractAgentState end

"""
    AgentState <: AbstractAgentState

Holds the current state vector `x`, an `AbstractTikhonovFilter` (either position or joint position/velocity),
LQR gain `K_lqr`, discrete system matrices `(Ad, Bd)`, and dynamics model `dyn`.
"""
mutable struct AgentState <: AbstractAgentState
    x::Vector{Float64}
    filter::AbstractTikhonovFilter{Float64, Vector{Float64}}
    K_lqr::Matrix{Float64}
    Ad::Matrix{Float64}
    Bd::Matrix{Float64}
    dyn::AbstractAgentDynamics
end

function AgentState(x0::Vector{Float64}, dyn::AbstractAgentDynamics, dt::Float64, K_lqr::Matrix{Float64}, eps::Float64; use_velocity::Bool=false)
    flt = use_velocity ? JointTikhonovFilter(zeros(length(x0)); epsilon = eps) : TikhonovFilter(zeros(length(x0)); epsilon = eps)
    Ad, Bd = discrete_matrices(dyn, dt)
    AgentState(copy(x0), flt, copy(K_lqr), Ad, Bd, dyn)
end

"""
    step_agent!(w::AgentState, qstar_target::Vector{Float64}, dt::Float64)

Steps standard feedback agent dynamics using position-only reference filtering.
"""
function step_agent!(w::AgentState, qstar_target::Vector{Float64}, dt::Float64)
    nx = length(w.x)
    ref_dim = min(nx, length(qstar_target))
    
    qstar_local = zeros(nx)
    qstar_local[1:ref_dim] = qstar_target[1:ref_dim]
    
    tikhonov_step!(w.filter, qstar_local, qstar_local, dt)
    x_ref = w.filter.x
    
    u = -w.K_lqr * (w.x - x_ref)
    w.x .= w.Ad * w.x .+ w.Bd * u
    
    return (copy(w.x), copy(x_ref))
end

"""
    step_agent!(w::AgentState, qstar_target::Vector{Float64}, qstar_dot_target::Vector{Float64}, dt::Float64)

Steps feedback agent dynamics using joint reference position AND reference velocity filtering.
Reference velocity is mapped into the state velocity components, eliminating tracking lag via standard LQR feedback.
"""
function step_agent!(w::AgentState, qstar_target::Vector{Float64}, qstar_dot_target::Vector{Float64}, dt::Float64)
    nx = length(w.x)
    ref_dim = min(nx, length(qstar_target))
    
    qstar_local = zeros(nx)
    qstar_local[1:ref_dim] = qstar_target[1:ref_dim]
    
    qstar_dot_local = zeros(nx)
    v_idxs = velocity_indices(w.dyn)
    v_dim = min(length(v_idxs), length(qstar_dot_target))
    for k in 1:v_dim
        qstar_dot_local[v_idxs[k]] = qstar_dot_target[k]
    end
    
    tikhonov_step!(w.filter, qstar_local, qstar_dot_local, dt)
    x_ref = copy(w.filter.x)
    if w.filter isa JointTikhonovFilter
        x_ref[v_idxs] .= w.filter.v[v_idxs]
    end

    u = -w.K_lqr * (w.x - x_ref)
    w.x .= w.Ad * w.x .+ w.Bd * u
    
    return (copy(w.x), copy(x_ref))
end

"""
    step_agent!(w::AgentState, qstar_target::Vector{Float64}, qstar_dot_target::Vector{Float64}, qstar_ddot_target::Vector{Float64}, dt::Float64)

Steps agent dynamics using position, velocity, and acceleration references.
Calculates acceleration-based attitude references (pitch/roll angles) and feedforward thrust,
enabling exact differential flatness trajectory tracking without lag.
"""
function step_agent!(w::AgentState, qstar_target::Vector{Float64}, qstar_dot_target::Vector{Float64}, qstar_ddot_target::Vector{Float64}, dt::Float64)
    nx = length(w.x)
    ref_dim = min(nx, length(qstar_target))
    
    qstar_local = zeros(nx)
    qstar_local[1:ref_dim] = qstar_target[1:ref_dim]
    
    qstar_dot_local = zeros(nx)
    v_idxs = velocity_indices(w.dyn)
    v_dim = min(length(v_idxs), length(qstar_dot_target))
    for k in 1:v_dim
        qstar_dot_local[v_idxs[k]] = qstar_dot_target[k]
    end
    
    if w.dyn isa QuadrotorDynamics && length(qstar_ddot_target) >= 3
        g = w.dyn.g
        ax, ay = qstar_ddot_target[1], qstar_ddot_target[2]
        qstar_local[4] = -ay / g
        qstar_local[5] =  ax / g
    elseif w.dyn isa PlanarQuadrotorDynamics && length(qstar_ddot_target) >= 2
        g = w.dyn.g
        ay = qstar_ddot_target[1]
        qstar_local[3] = -ay / g
    end

    tikhonov_step!(w.filter, qstar_local, qstar_dot_local, dt)
    x_ref = copy(w.filter.x)
    if w.filter isa JointTikhonovFilter
        x_ref[v_idxs] .= w.filter.v[v_idxs]
    end
    
    u_fb = -w.K_lqr * (w.x - x_ref)
    
    u_ff = zeros(size(w.Bd, 2))
    if w.dyn isa QuadrotorDynamics && length(qstar_ddot_target) >= 3
        u_ff[1] = w.dyn.m * qstar_ddot_target[3]
    elseif w.dyn isa PlanarQuadrotorDynamics && length(qstar_ddot_target) >= 2
        u_ff[1] = w.dyn.m * qstar_ddot_target[2] / 2.0
        u_ff[2] = w.dyn.m * qstar_ddot_target[2] / 2.0
    end

    u = u_fb + u_ff
    w.x .= w.Ad * w.x .+ w.Bd * u
    
    return (copy(w.x), copy(x_ref))
end

end # module
