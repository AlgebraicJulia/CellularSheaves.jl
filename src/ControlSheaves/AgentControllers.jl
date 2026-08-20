module AgentControllers

using LinearAlgebra
using ArgCheck
using ...NetworkSheaves.TrajectorySheaves: continuous_to_discrete_zoh
using ..Tikhonov: AbstractTikhonovFilter, TikhonovFilter, JointTikhonovFilter, tikhonov_step!

export AbstractAgentDynamics, QuadrotorDynamics, PlanarQuadrotorDynamics,
       AbstractControlAffine, SingleIntegrator, DoubleIntegrator,
       AbstractAgentController, LQRController,
       AbstractAgentState, AgentState,
       solve_dare, step_agent!,
       position_indices, velocity_indices, state_dim, initial_state,
       continuous_matrices, discrete_matrices,
       SE3QuadrotorDynamics, SE3Reference, GeometricSE3Controller,
       SE3Errors, SE3Certificate, SE3AgentState,
       geometric_control, se3_lyapunov, se3_certificate, se3_tune_lyapunov, se3_autotune,
       se3_error_vector, se3_measured_rate, se3_unpack, se3_pack, se3_derivative, se3_step,
       attitude_indices, angular_velocity_indices,
       hat_so3, vee_so3, expm_so3, project_SO3, deriv_unit_vector

abstract type AbstractAgentDynamics end

"""
    AbstractControlAffine <: AbstractAgentDynamics

Agents whose configuration is *directly actuated*: the input reaches it in one derivative
for a [`SingleIntegrator`](@ref) and two for a [`DoubleIntegrator`](@ref), rather than
through an attitude loop as in [`QuadrotorDynamics`](@ref). That is what decides whether a
configuration barrier is usable, so it is what the hierarchy records.
"""
abstract type AbstractControlAffine <: AbstractAgentDynamics end

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

"""
    continuous_matrices(dyn::AbstractAgentDynamics)

Continuous-time state and input matrices `(Ac, Bc)` of the linearized agent model,
``\\dot{x} = A_c x + B_c u``.
"""
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

"""
    SingleIntegrator(n; input = I)

Single integrator on ``\\mathbb{R}^n``: state is the configuration alone, ``\\dot{q} = B u``
with ``B = I`` by default. The configuration is the whole state, so a configuration barrier
has relative degree one and the plain `DistanceBarrier` applies.

`input` is the ``n \\times m`` map from command to configuration rate; a non-identity `input`
describes an agent that cannot move equally freely in every direction.
"""
struct SingleIntegrator <: AbstractControlAffine
    n::Int
    input::Matrix{Float64}
end

function SingleIntegrator(n::Integer; input = I)
    @argcheck n > 0
    B = input isa UniformScaling ? Matrix{Float64}(input, n, n) : Matrix{Float64}(input)
    @argcheck size(B, 1) == n
    return SingleIntegrator(Int(n), B)
end

"""
    DoubleIntegrator(n; damping = nothing, input = I)

Double integrator on ``n`` position variables, in the Darboux coordinates ``(q, p)`` of
``\\mathbb{R}^{2n}``: the state is a configuration stacked on its rate, and

```math
A = \\begin{pmatrix} 0 & I \\\\ 0 & D \\end{pmatrix}, \\qquad
B = \\begin{pmatrix} 0 \\\\ B_p \\end{pmatrix} .
```

Only the blocks acting on the rate are specified; the kinematic ``\\dot{q} = p`` coupling and
the zero rows are padded on. `damping` is ``D``, defaulting to the *free* double integrator;
`input` is ``B_p``, defaulting to the identity so the command is an acceleration.

The configuration is not directly actuated, so a configuration barrier has relative degree
two and needs the braking form; see `BrakingBarrier` and `RelativeDegreeError`.
"""
struct DoubleIntegrator <: AbstractControlAffine
    n::Int
    damping::Matrix{Float64}
    input::Matrix{Float64}
end

function DoubleIntegrator(n::Integer; damping = nothing, input = I)
    @argcheck n > 0
    D = damping === nothing ? zeros(n, n) :
        damping isa UniformScaling ? Matrix{Float64}(damping, n, n) : Matrix{Float64}(damping)
    B = input isa UniformScaling ? Matrix{Float64}(input, n, n) : Matrix{Float64}(input)
    @argcheck size(D) == (n, n)
    @argcheck size(B, 1) == n
    return DoubleIntegrator(Int(n), D, B)
end

function continuous_matrices(dyn::SingleIntegrator)
    return zeros(dyn.n, dyn.n), copy(dyn.input)
end

# The Darboux padding: the caller gave the blocks acting on the rate, and the kinematic
# coupling and zero rows are supplied here rather than at every call site.
function continuous_matrices(dyn::DoubleIntegrator)
    n = dyn.n
    Ac = zeros(2n, 2n)
    Ac[1:n, (n + 1):(2n)] = Matrix{Float64}(I, n, n)
    Ac[(n + 1):(2n), (n + 1):(2n)] = dyn.damping
    Bc = vcat(zeros(n, size(dyn.input, 2)), dyn.input)
    return Ac, Bc
end

position_indices(dyn::SingleIntegrator) = 1:(dyn.n)
# A single integrator has no rate state at all, so the velocity selector is empty. This is
# what makes `ControlAffineModel` hand the filter a `nothing` velocity, and so what makes the
# braking barrier report that it cannot be used here.
velocity_indices(dyn::SingleIntegrator) = 1:0
state_dim(dyn::SingleIntegrator) = dyn.n

position_indices(dyn::DoubleIntegrator) = 1:(dyn.n)
velocity_indices(dyn::DoubleIntegrator) = (dyn.n + 1):(2 * dyn.n)
state_dim(dyn::DoubleIntegrator) = 2 * dyn.n

# There is no attitude to trim, so the acceleration argument carries no information here.
function initial_state(dyn::AbstractControlAffine, position::AbstractVector,
                       velocity::AbstractVector, acceleration::AbstractVector)
    @argcheck length(position) == length(position_indices(dyn))
    @argcheck length(velocity) == length(velocity_indices(dyn))
    x = zeros(state_dim(dyn))
    x[position_indices(dyn)] .= position
    x[velocity_indices(dyn)] .= velocity
    return x
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
    _default_initial_position(dyn::AbstractAgentDynamics, i::Int) -> Vector{Float64}

Default "airstrip" starting position for agent `i` when no explicit initial position is
supplied: agents are lined up along the first position coordinate at a fixed spacing, at a
common hover altitude (the last position coordinate). This is deliberately distinct from the
target formation so that convergence into formation is visible in a simulation.
"""
function _default_initial_position(dyn::AbstractAgentDynamics, i::Int)
    pos = zeros(length(position_indices(dyn)))
    pos[1] = (i - 1) * 0.5   # spacing between agents along the line
    pos[end] = 1.5           # common hover altitude
    return pos
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

"""
    AbstractAgentController

Supertype for the local control laws an agent applies to track its reference. Implementations
supply a [`step_agent!`](@ref) method advancing an [`AbstractAgentState`](@ref) by one
control period.
"""
abstract type AbstractAgentController end

"""
    LQRController(K)
    LQRController(dyn, dt, Q, R)

Infinite-horizon LQR feedback ``u = -K (x - x_{\\text{ref}})``. The second form solves the
discrete algebraic Riccati equation for `dyn` discretized at `dt`; see [`solve_dare`](@ref).
"""
struct LQRController <: AbstractAgentController
    K::Matrix{Float64}
end

function LQRController(dyn::AbstractAgentDynamics, dt::Float64, Q::AbstractMatrix, R::AbstractMatrix)
    Ad, Bd = discrete_matrices(dyn, dt)
    K = solve_dare(Ad, Bd, Q, R)
    LQRController(K)
end

"""
    AbstractAgentState

Supertype for the per-agent state an [`AbstractAgentController`](@ref) advances. An
implementation carries whatever a [`step_agent!`](@ref) method needs to take one control
period: at minimum the state vector itself and the reference filter feeding it. See
[`AgentState`](@ref) for the concrete type used throughout.
"""
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

include("geometric_se3.jl")

end # module
