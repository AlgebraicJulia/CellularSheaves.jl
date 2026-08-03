# # 12-Agent Escort Worker Component
#
# This file contains the onboard vehicle state and control implementation for each agent's worker process.

using CellularSheaves
using CellularSheaves.ControlSheaves.Tikhonov
using LinearAlgebra

"""
    AgentWorkerState

State container held on each agent's dedicated Julia worker process.
Encapsulates physical state `x` (10D: x, y, z, phi, theta, dx, dy, dz, dphi, dtheta),
10D Tikhonov reference filter, discrete state-space dynamics (Ad, Bd),
and the 10D Discrete LQR feedback matrix `K_lqr`.
"""
mutable struct AgentWorkerState
    x::Vector{Float64}
    filter::TikhonovFilter{Float64, Vector{Float64}}
    K_lqr::Matrix{Float64}
    Ad::Matrix{Float64}
    Bd::Matrix{Float64}
end

const LOCAL_AGENT = Ref{Union{Nothing, AgentWorkerState}}(nothing)

"""
    init_worker_agent!(x0, K, A, B, eps)

Initialize the agent's flight computer state on its local worker process.
"""
function init_worker_agent!(x0::Vector{Float64}, K::Matrix{Float64}, A::Matrix{Float64}, B::Matrix{Float64}, eps::Float64)
    flt = TikhonovFilter(zeros(length(x0)); epsilon = eps)
    LOCAL_AGENT[] = AgentWorkerState(copy(x0), flt, copy(K), copy(A), copy(B))
    return nothing
end

"""
    step_worker_agent!(qstar_3d_i, dt)

Step the agent's onboard controller and vehicle dynamics for one time step `dt`.
1. Embeds 3D spatial target `qstar_3d_i = [x*, y*, z*]` into 10D reference `[x*, y*, z*, 0, 0, 0, 0, 0, 0, 0]`.
2. Filters the 10D reference through local Tikhonov filter.
3. Evaluates 10D Discrete LQR control law `u = -K_lqr * (x - x_ref)`.
4. Integrates vehicle state `x_{t+1} = Ad * x_t + Bd * u_t`.
Returns `(x_actual, x_ref)`.
"""
function step_worker_agent!(qstar_3d_i::Vector{Float64}, dt::Float64)
    w = LOCAL_AGENT[]
    qstar_10d = zeros(10)
    qstar_10d[1:3] = qstar_3d_i
    tikhonov_step!(w.filter, qstar_10d, qstar_10d, dt)
    x_ref = w.filter.x
    u = -w.K_lqr * (w.x - x_ref)
    w.x .= w.Ad * w.x .+ w.Bd * u
    return (copy(w.x), copy(x_ref))
end
