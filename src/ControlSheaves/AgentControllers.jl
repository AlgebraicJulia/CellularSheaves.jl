module AgentControllers

using LinearAlgebra
using ...NetworkSheaves.TrajectorySheaves: continuous_to_discrete_zoh
using ..Tikhonov: AbstractTikhonovFilter, TikhonovFilter, JointTikhonovFilter, tikhonov_step!

export AbstractAgentDynamics, QuadrotorDynamics, PlanarQuadrotorDynamics,
       AbstractAgentController, LQRController, FeedforwardLQRController,
       AbstractAgentState, AgentState, FeedforwardAgentState,
       solve_dare, step_agent!

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

"""
    FeedforwardLQRController <: AbstractAgentController

Holds state feedback gain `K`, continuous system matrix `Ac`, and pseudoinverse `Bc_pinv`
for calculating feedforward control effort:

    u_ff = B^dagger * (xdot_ref - A_c * x)
"""
struct FeedforwardLQRController <: AbstractAgentController
    K::Matrix{Float64}
    Ac::Matrix{Float64}
    Bc_pinv::Matrix{Float64}
end

function FeedforwardLQRController(dyn::AbstractAgentDynamics, dt::Float64, Q::AbstractMatrix, R::AbstractMatrix)
    Ac, Bc = continuous_matrices(dyn)
    Ad, Bd = discrete_matrices(dyn, dt)
    K = solve_dare(Ad, Bd, Q, R)
    Bc_pinv = pinv(Bc)
    FeedforwardLQRController(K, Ac, Bc_pinv)
end

abstract type AbstractAgentState end

"""
    AgentState <: AbstractAgentState

Holds the current state vector `x`, a `TikhonovFilter`, and references to dynamics and controller.
"""
mutable struct AgentState <: AbstractAgentState
    x::Vector{Float64}
    filter::TikhonovFilter{Float64, Vector{Float64}}
    K_lqr::Matrix{Float64}
    Ad::Matrix{Float64}
    Bd::Matrix{Float64}
end

function AgentState(x0::Vector{Float64}, dyn::AbstractAgentDynamics, dt::Float64, K_lqr::Matrix{Float64}, eps::Float64)
    flt = TikhonovFilter(zeros(length(x0)); epsilon = eps)
    Ad, Bd = discrete_matrices(dyn, dt)
    AgentState(copy(x0), flt, copy(K_lqr), Ad, Bd)
end

"""
    velocity_indices(dyn::QuadrotorDynamics) -> 6:8
    velocity_indices(dyn::PlanarQuadrotorDynamics) -> 4:5

Returns the velocity state indices in the full state vector for the given agent dynamics.
"""
velocity_indices(dyn::QuadrotorDynamics) = 6:8
velocity_indices(dyn::PlanarQuadrotorDynamics) = 4:5

"""
    FeedforwardAgentState <: AbstractAgentState

Holds state vector `x`, a `JointTikhonovFilter` for joint reference and velocity filtering,
and matrices required for feedforward control.
"""
mutable struct FeedforwardAgentState <: AbstractAgentState
    x::Vector{Float64}
    filter::JointTikhonovFilter{Float64, Vector{Float64}}
    K_lqr::Matrix{Float64}
    Ac::Matrix{Float64}
    Bc_pinv::Matrix{Float64}
    Ad::Matrix{Float64}
    Bd::Matrix{Float64}
    dyn::AbstractAgentDynamics
end

function FeedforwardAgentState(x0::Vector{Float64}, dyn::AbstractAgentDynamics, dt::Float64, ctrl::FeedforwardLQRController, eps::Float64)
    flt = JointTikhonovFilter(zeros(length(x0)); epsilon = eps)
    Ac, _ = continuous_matrices(dyn)
    Ad, Bd = discrete_matrices(dyn, dt)
    FeedforwardAgentState(copy(x0), flt, copy(ctrl.K), copy(Ac), copy(ctrl.Bc_pinv), Ad, Bd, dyn)
end

function FeedforwardAgentState(x0::Vector{Float64}, dyn::AbstractAgentDynamics, dt::Float64, K_lqr::Matrix{Float64}, eps::Float64)
    Ac, Bc = continuous_matrices(dyn)
    Ad, Bd = discrete_matrices(dyn, dt)
    flt = JointTikhonovFilter(zeros(length(x0)); epsilon = eps)
    FeedforwardAgentState(copy(x0), flt, copy(K_lqr), copy(Ac), pinv(Bc), Ad, Bd, dyn)
end

"""
    step_agent!(w::AgentState, qstar_target::Vector{Float64}, dt::Float64)

Steps standard feedback agent dynamics.
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
    step_agent!(w::FeedforwardAgentState, qstar_target::Vector{Float64}, qstar_dot_target::Vector{Float64}, dt::Float64)

Steps feedforward-enhanced agent dynamics using feedforward control signal:
    u = -K*(x - x_ref) + B^dagger * (v_ref - A_c * x)
"""
function step_agent!(w::FeedforwardAgentState, qstar_target::Vector{Float64}, qstar_dot_target::Vector{Float64}, dt::Float64)
    nx = length(w.x)
    ref_dim = min(nx, length(qstar_target))
    
    qstar_local = zeros(nx)
    qstar_local[1:ref_dim] = qstar_target[1:ref_dim]
    
    # Map spatial velocity reference to the agent's velocity state indices
    qstar_dot_local = zeros(nx)
    v_idxs = velocity_indices(w.dyn)
    v_dim = min(length(v_idxs), length(qstar_dot_target))
    for k in 1:v_dim
        qstar_dot_local[v_idxs[k]] = qstar_dot_target[k]
    end
    
    tikhonov_step!(w.filter, qstar_local, qstar_dot_local, dt)
    x_ref = w.filter.x
    v_ref = w.filter.v
    
    u_fb = -w.K_lqr * (w.x - x_ref)
    u_ff = w.Bc_pinv * (v_ref - w.Ac * w.x)
    u = u_fb + u_ff
    
    w.x .= w.Ad * w.x .+ w.Bd * u
    
    return (copy(w.x), copy(x_ref))
end

end # module
