module AgentControllers

using LinearAlgebra
using ...NetworkSheaves.TrajectorySheaves: continuous_to_discrete_zoh
using ..Tikhonov: TikhonovFilter, tikhonov_step!

export AbstractAgentDynamics, QuadrotorDynamics, LQRController, AgentState, solve_dare, step_agent!

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

function discrete_matrices(dyn::QuadrotorDynamics, dt::Float64)
    Ac, Bc = continuous_matrices(dyn)
    continuous_to_discrete_zoh(Ac, Bc, dt)
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

struct LQRController
    K::Matrix{Float64}
end

"""
    LQRController(dyn::QuadrotorDynamics, dt::Float64, Q::AbstractMatrix, R::AbstractMatrix)

Constructs an LQR controller for the given dynamics and LQR cost matrices.
"""
function LQRController(dyn::QuadrotorDynamics, dt::Float64, Q::AbstractMatrix, R::AbstractMatrix)
    Ad, Bd = discrete_matrices(dyn, dt)
    K = solve_dare(Ad, Bd, Q, R)
    LQRController(K)
end

"""
    AgentState

Holds the current state vector `x`, a `TikhonovFilter`, and references to the dynamics and controller. 
This struct is instantiated on the worker processes.
"""
mutable struct AgentState
    x::Vector{Float64}
    filter::TikhonovFilter{Float64, Vector{Float64}}
    K_lqr::Matrix{Float64}
    Ad::Matrix{Float64}
    Bd::Matrix{Float64}
end

"""
    AgentState(x0::Vector{Float64}, dyn::QuadrotorDynamics, dt::Float64, K_lqr::Matrix{Float64}, eps::Float64)

Initializes the agent's flight computer state.
"""
function AgentState(x0::Vector{Float64}, dyn::QuadrotorDynamics, dt::Float64, K_lqr::Matrix{Float64}, eps::Float64)
    flt = TikhonovFilter(zeros(length(x0)); epsilon = eps)
    Ad, Bd = discrete_matrices(dyn, dt)
    AgentState(copy(x0), flt, copy(K_lqr), Ad, Bd)
end

"""
    step_agent!(w::AgentState, qstar_target::Vector{Float64}, dt::Float64)

Steps the agent dynamics. Extracts the necessary top-level reference signal from `qstar_target` 
(based on the agent state dimension) and applies the LQR control.
"""
function step_agent!(w::AgentState, qstar_target::Vector{Float64}, dt::Float64)
    # The reference is dynamically padded/extracted to match the state dimension of this agent
    nx = length(w.x)
    ref_dim = min(nx, length(qstar_target))
    
    qstar_local = zeros(nx)
    qstar_local[1:ref_dim] = qstar_target[1:ref_dim]
    
    # Run the Tikhonov filter for the reference signal
    tikhonov_step!(w.filter, qstar_local, qstar_local, dt)
    x_ref = w.filter.x
    
    # Calculate control effort
    u = -w.K_lqr * (w.x - x_ref)
    
    # Step dynamics
    w.x .= w.Ad * w.x .+ w.Bd * u
    
    return (copy(w.x), copy(x_ref))
end

end # module
