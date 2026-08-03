# # Scenario 5 Worker Agent Component
#
# This file contains the onboard vehicle state and control implementation for each agent's worker process.

using CellularSheaves
using CellularSheaves.ControlSheaves.Tikhonov
using LinearAlgebra

"""
    AgentWorkerState

State container held on each agent's dedicated Julia worker process.
Encapsulates physical state `x`, Tikhonov reference filter, discrete state-space dynamics (Ad, Bd),
and the discrete LQR feedback matrix `K_lqr`.
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
    step_worker_agent!(qstar_i, dt)

Step the agent's onboard controller and vehicle dynamics for one time step `dt`.
1. Filters high-level harmonic cochain `qstar_i` through Tikhonov reference filter.
2. Evaluates Discrete LQR control law `u = -K_lqr * (x - x_ref)`.
3. Integrates vehicle state `x_{t+1} = Ad * x_t + Bd * u_t`.
Returns `(x_actual, x_ref)`.
"""
function step_worker_agent!(qstar_i::Vector{Float64}, dt::Float64)
    w = LOCAL_AGENT[]
    tikhonov_step!(w.filter, qstar_i, qstar_i, dt)
    x_ref = w.filter.x
    u = -w.K_lqr * (w.x - x_ref)
    w.x .= w.Ad * w.x .+ w.Bd * u
    return (copy(w.x), copy(x_ref))
end
