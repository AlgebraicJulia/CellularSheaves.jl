"""
    MechanicalHJB

Optimal feedback for a torque-limited two-link arm (a double pendulum with a
motor at each joint) from the Hamilton–Jacobi–Bellman equation in canonical
coordinates, solved on a grid by policy iteration with the matrix-free
domain-decomposition solver of [`DoubleIntegratorHJB`](@ref CellularSheaves.ControlSheaves.DoubleIntegratorHJB)
(`GridPolicyIteration`: CPU threads, GPUs, MPI boxes).

**The mechanical system.** Two uniform links in a vertical plane, joined by
revolute joints. The configuration is ``q = (θ_1, θ_2)``: ``θ_1`` the angle of
the first link from the upward vertical, ``θ_2`` the angle of the second link
relative to the first, so ``q = 0`` is the arm balanced upright and
``q = (π, 0)`` hangs down. The kinetic energy is ``\\tfrac12 \\dot q^\\mathsf{T} M(q) \\dot q``
with the mass matrix

```math
M(q) = \\begin{pmatrix} a_1 + 2a_2\\cos θ_2 & a_3 + a_2 \\cos θ_2 \\\\ a_3 + a_2\\cos θ_2 & a_3 \\end{pmatrix},
\\quad a_1 = I_1 + I_2 + m_1 c_1^2 + m_2(\\ell_1^2 + c_2^2),\\; a_2 = m_2 \\ell_1 c_2,\\; a_3 = I_2 + m_2 c_2^2,
```

(masses ``m_i``, lengths ``\\ell_i``, centres of mass ``c_i`` from the joint,
moments of inertia ``I_i`` about the centres; Spong, Hutchinson and
Vidyasagar, *Robot Modeling and Control*, 2006, §7.4) and the potential energy
is ``U(q) = b_1 \\cos θ_1 + b_2 \\cos(θ_1 + θ_2)`` with ``b_1 = (m_1 c_1 + m_2 \\ell_1) g``,
``b_2 = m_2 c_2 g``. In the canonical coordinates ``x = (q, p)``, ``p = M(q)\\dot q``
the generalized momenta, the arm follows Hamilton's equations with the joint
torques ``τ`` and viscous joint friction ``β``:

```math
H(q, p) = \\tfrac12 p^\\mathsf{T} M(q)^{-1} p + U(q), \\qquad
\\dot q = \\frac{∂H}{∂p}, \\qquad \\dot p = -\\frac{∂H}{∂q} - β \\dot q + τ,
\\qquad |τ_k| \\le \\bar τ_k .
```

The torques enter the momentum equations directly (the arm is fully
actuated), but they are bounded: with ``\\bar τ`` below the gravity torque the
arm cannot simply lift itself upright, it has to pump energy by swinging
(Tedrake, *Underactuated Robotics*, ch. 3).

**The control problem.** Minimize the discounted cost

```math
V(x) = \\min_{τ(\\cdot)} \\int_0^\\infty e^{-ρt}\\Big(w \\sum_k \\big(1 - \\cos(θ_k - θ_k^\\star)\\big)
    + c\\, T(q, p) + \\tfrac r2 \\|τ\\|^2\\Big)\\,dt
```

of being away from the target configuration ``q^\\star`` (by default upright),
of kinetic energy ``T = \\tfrac12 p^\\mathsf{T} M^{-1} p`` and of torque. ``V``
solves ``ρV = \\min_τ [\\ell(x, τ) + ∇V \\cdot f(x, τ)]``, minimized pointwise
over the box of torques.

**Discretization.** The state space is the torus of joint angles times a box
of momenta ``|p_k| \\le P_k``. The grid is periodic in the angles (an even
number of points, so that red–black Gauss–Seidel sweeps alternate colours
around the circle) and the momentum box is a state constraint: the scheme
drops flux leaving it (a reflecting boundary), so pick ``P`` large enough that
optimal motions stay inside. The upwind scheme of Kushner and Dupuis gives a
monotone M-matrix for every policy, for any drift, so policy iteration
(Howard; Bokanowski, Maroso and Zidani, 2009) converges as for the double
integrator.
"""
module MechanicalHJB

using ArgCheck
using CommonSolve
using CommonSolve: solve
using LinearAlgebra
using ..DoubleIntegratorHJB: StateGrid, HJBSolution, GridPolicyIteration
import ..DoubleIntegratorHJB: value_at, control_at, _state_grid, _periodic, _control_count, _discount, _kernel_model,
    _drifts, _cost, _dirichlet, _boundary_value, _argmin_hamiltonian, _initial_value, _initial_control, _box_argmin,
    _grid_solve

export TwoLinkArm, mass_matrix, potential_energy, hamiltonian, hamiltonian_vector_field, MechanicalHJBProblem,
    closed_loop, value_at, control_at, GridPolicyIteration

# ===========================================================================
# the arm
# ===========================================================================

"""
    TwoLinkArm(; masses = (1, 1), lengths = (1, 1), centers = lengths ./ 2,
               inertias = masses .* lengths .^ 2 ./ 12, gravity = 9.81, damping = (0, 0))

A planar two-link arm with revolute joints, in SI units: link masses (kg),
lengths (m), distances of the centres of mass from the joints (m), moments of
inertia about the centres of mass (kg m², default: uniform rods), gravitational
acceleration (m/s²) and viscous joint friction (N m s). See the module
documentation for the equations of motion.
"""
struct TwoLinkArm
    masses::NTuple{2,Float64}
    lengths::NTuple{2,Float64}
    centers::NTuple{2,Float64}
    inertias::NTuple{2,Float64}
    gravity::Float64
    damping::NTuple{2,Float64}
end

function TwoLinkArm(; masses = (1.0, 1.0), lengths = (1.0, 1.0), centers = lengths ./ 2,
        inertias = masses .* lengths .^ 2 ./ 12, gravity = 9.81, damping = (0.0, 0.0))
    @argcheck all(>(0), masses) && all(>(0), lengths) "masses and lengths must be positive"
    @argcheck all(0 .<= centers .<= lengths) "the centres of mass must lie on the links"
    @argcheck all(>(0), inertias) "the moments of inertia must be positive"
    @argcheck gravity >= 0 && all(>=(0), damping)
    f(t) = Float64.(Tuple(t))
    return TwoLinkArm(f(masses), f(lengths), f(centers), f(inertias), Float64(gravity), f(damping))
end

# The constants of M(q) and U(q) (module docstring).
function _inertia(arm::TwoLinkArm)
    (m1, m2), (l1, _), (c1, c2), (I1, I2) = arm.masses, arm.lengths, arm.centers, arm.inertias
    return (I1 + I2 + m1 * c1^2 + m2 * (l1^2 + c2^2), m2 * l1 * c2, I2 + m2 * c2^2)
end

function _gravity(arm::TwoLinkArm)
    (m1, m2), (l1, _), (c1, c2) = arm.masses, arm.lengths, arm.centers
    return ((m1 * c1 + m2 * l1) * arm.gravity, m2 * c2 * arm.gravity)
end

"""
    mass_matrix(arm, q) -> Matrix

The mass (inertia) matrix ``M(q)`` of the arm at the joint angles `q`.
"""
function mass_matrix(arm::TwoLinkArm, q::AbstractVector)
    a1, a2, a3 = _inertia(arm)
    c = cos(q[2])
    return [a1+2a2*c a3+a2*c; a3+a2*c a3]
end

"""
    potential_energy(arm, q)

The gravitational potential energy ``U(q)`` (zero with both links horizontal;
largest upright, at ``q = 0``).
"""
function potential_energy(arm::TwoLinkArm, q::AbstractVector)
    b1, b2 = _gravity(arm)
    return b1 * cos(q[1]) + b2 * cos(q[1] + q[2])
end

"""
    hamiltonian(arm, x)

The total energy ``H(q, p) = \\tfrac12 p^\\mathsf{T} M(q)^{-1} p + U(q)`` at the
state `x = [q; p]`.
"""
function hamiltonian(arm::TwoLinkArm, x::AbstractVector)
    q, p = x[1:2], x[3:4]
    return dot(p, mass_matrix(arm, q) \ p) / 2 + potential_energy(arm, q)
end

"""
    hamiltonian_vector_field(arm, x, τ) -> Vector

The time derivative ``(\\dot q, \\dot p)`` of the state `x = [q; p]` under the
joint torques `τ`: Hamilton's equations with friction (module docstring).
"""
hamiltonian_vector_field(arm::TwoLinkArm, x::AbstractVector, τ::AbstractVector) =
    collect(_flow(_ArmModel(arm), Tuple(x), (Float64(τ[1]), Float64(τ[2]))))

# ===========================================================================
# the kernel model
# ===========================================================================

# Everything the grid kernels need, as plain bits (runs on GPUs).
struct _ArmModel
    a1::Float64
    a2::Float64
    a3::Float64
    b1::Float64
    b2::Float64
    β1::Float64
    β2::Float64
    target::NTuple{2,Float64}
    configuration_weight::Float64
    kinetic_weight::Float64
    control_weight::Float64
    torque_bound::NTuple{2,Float64}
end

_ArmModel(arm::TwoLinkArm; target = (0.0, 0.0), configuration_weight = 0.0, kinetic_weight = 0.0,
          control_weight = 0.0, torque_bound = (0.0, 0.0)) =
    _ArmModel(_inertia(arm)..., _gravity(arm)..., arm.damping..., target, configuration_weight, kinetic_weight,
        control_weight, torque_bound)

# The joint velocities q̇ = M⁻¹p and the torque-free momentum equations
# ṗ + τ = -∂H/∂q - βq̇ (with ∂T/∂θ₂ = -½ q̇ᵀ (∂M/∂θ₂) q̇ at fixed p).
@inline function _free(m::_ArmModel, x)
    θ1, θ2, p1, p2 = x[1], x[2], x[3], x[4]
    s2, c2 = sincos(θ2)
    M11 = m.a1 + 2m.a2 * c2
    M12 = m.a3 + m.a2 * c2
    M22 = m.a3
    det = M11 * M22 - M12^2
    v1 = (M22 * p1 - M12 * p2) / det
    v2 = (M11 * p2 - M12 * p1) / det
    s12 = sin(θ1 + θ2)
    g1 = m.b1 * sin(θ1) + m.b2 * s12 - m.β1 * v1
    g2 = -m.a2 * s2 * (v1^2 + v1 * v2) + m.b2 * s12 - m.β2 * v2
    return v1, v2, g1, g2
end

@inline function _flow(m::_ArmModel, x, τ)
    v1, v2, g1, g2 = _free(m, x)
    return (v1, v2, g1 + τ[1], g2 + τ[2])
end

@inline _drifts(m::_ArmModel, x, u) = _flow(m, x, u)

@inline function _cost(m::_ArmModel, x, u)
    v1, v2, _, _ = _free(m, x)
    configuration = 2 - cos(x[1] - m.target[1]) - cos(x[2] - m.target[2])
    kinetic = (x[3] * v1 + x[4] * v2) / 2
    return m.configuration_weight * configuration + m.kinetic_weight * kinetic + m.control_weight * (u[1]^2 + u[2]^2) / 2
end

_dirichlet(::_ArmModel) = false
@inline _boundary_value(::_ArmModel, x) = 0.0             # not used: reflecting momentum boundary
@inline _initial_value(::_ArmModel, x) = 0.0
@inline _initial_control(::_ArmModel, x) = (0.0, 0.0)     # start from the unactuated arm

@inline function _argmin_hamiltonian(m::_ArmModel, x, Dp::NTuple{2}, Dm::NTuple{2})
    _, _, g1, g2 = _free(m, x)
    r = m.control_weight
    return (_box_argmin(r, m.torque_bound[1], g1, Dp[1], Dm[1]), _box_argmin(r, m.torque_bound[2], g2, Dp[2], Dm[2]))
end

# ===========================================================================
# the problem
# ===========================================================================

"""
    MechanicalHJBProblem(arm::TwoLinkArm; angle_points = 32, momentum_points = 33,
                         momentum_bound = nothing, target = (0, 0), torque_bound = (5, 2.5),
                         configuration_weight = 1, kinetic_weight = 0.01, control_weight = 0.02,
                         discount = 0.2)

The discounted HJB problem of the module docstring for `arm`, on a grid of the
state ``x = (θ_1, θ_2, p_1, p_2)``: `angle_points` (even) points around each
joint circle ``[-π, π)``, `momentum_points` points on ``[-P_k, P_k]``. The
default `momentum_bound` ``P`` is 1.2 times the largest momenta the
unactuated, frictionless arm reaches falling from upright to hanging down
(energy ``2(b_1 + b_2)``; ``|p_k| \\le \\sqrt{2 T M_{kk}}``). `target` is the
configuration ``q^\\star`` held at no cost (default upright, the unstable
equilibrium), `torque_bound` the motor limits ``\\bar τ`` (N m). The default
discount ``ρ = 0.2`` (a horizon of about 5 s) makes a swing-up worth its
cost: hanging down costs ``2w`` per second for ever, ``2w/ρ`` discounted, and a
swing-up of a few seconds less.
With a much larger ``ρ`` the optimal policy from hanging down is to stay there.

Solve with `solve(problem, GridPolicyIteration(...))`; query the result with
[`value_at`](@ref), [`control_at`](@ref) and [`closed_loop`](@ref).
"""
struct MechanicalHJBProblem
    arm::TwoLinkArm
    grid::StateGrid
    target::NTuple{2,Float64}
    torque_bound::NTuple{2,Float64}
    configuration_weight::Float64
    kinetic_weight::Float64
    control_weight::Float64
    discount::Float64
end

function MechanicalHJBProblem(arm::TwoLinkArm; angle_points::Integer = 32, momentum_points::Integer = 33,
        momentum_bound = nothing, target = (0.0, 0.0), torque_bound = (5.0, 2.5), configuration_weight::Real = 1.0,
        kinetic_weight::Real = 0.01, control_weight::Real = 0.02, discount::Real = 0.2)
    @argcheck angle_points >= 4 && iseven(angle_points) "need an even number (≥ 4) of angle points"
    @argcheck momentum_points >= 3
    @argcheck all(>=(0), torque_bound) "torque bounds must be nonnegative"
    @argcheck configuration_weight >= 0 && kinetic_weight >= 0 && control_weight > 0
    @argcheck discount > 0 "the discount rate must be positive"
    P = momentum_bound === nothing ? _default_momentum_bound(arm) : Float64.(Tuple(momentum_bound))
    @argcheck all(>(0), P) "momentum bounds must be positive"
    h = 2π / angle_points
    grid = StateGrid([-π, -π, -P[1], -P[2]], [π - h, π - h, P[1], P[2]],
        [angle_points, angle_points, momentum_points, momentum_points])
    return MechanicalHJBProblem(arm, grid, Float64.(Tuple(target)), Float64.(Tuple(torque_bound)),
        Float64(configuration_weight), Float64(kinetic_weight), Float64(control_weight), Float64(discount))
end

function _default_momentum_bound(arm::TwoLinkArm)
    a1, a2, a3 = _inertia(arm)
    b1, b2 = _gravity(arm)
    T = 2 * (b1 + b2)
    return (1.2 * sqrt(2T * (a1 + 2a2)), 1.2 * sqrt(2T * a3))
end

_state_grid(prob::MechanicalHJBProblem) = prob.grid
_periodic(::MechanicalHJBProblem) = (true, true, false, false)
_control_count(::MechanicalHJBProblem) = 2
_discount(prob::MechanicalHJBProblem) = prob.discount
_kernel_model(prob::MechanicalHJBProblem) =
    _ArmModel(prob.arm; target = prob.target, configuration_weight = prob.configuration_weight,
        kinetic_weight = prob.kinetic_weight, control_weight = prob.control_weight, torque_bound = prob.torque_bound)

CommonSolve.solve(prob::MechanicalHJBProblem, alg::GridPolicyIteration) = _grid_solve(prob, alg)

# ===========================================================================
# the solution
# ===========================================================================

# Multilinear interpolation of the columns of F (or of a vector) at x: around
# the joint circles (angles wrap), and clamped to the momentum box.
function _interpolate(g::StateGrid, F::AbstractVecOrMat, x::AbstractVector)
    lo = zeros(Int, 4)
    hi = zeros(Int, 4)
    t = zeros(4)
    for k in 1:4
        n = g.points[k]
        s = (x[k] - g.lower[k]) / g.spacing[k]
        if k <= 2
            s = mod(s, n)
            lo[k] = min(floor(Int, s), n - 1)
            hi[k] = mod(lo[k] + 1, n)
        else
            s = clamp(s, 0, n - 1)
            lo[k] = min(floor(Int, s), n - 2)
            hi[k] = lo[k] + 1
        end
        t[k] = s - lo[k]
    end
    acc = F isa AbstractVector ? 0.0 : zeros(size(F, 1))
    for corner in 0:15
        weight = 1.0
        n = 1
        for k in 1:4
            bit = (corner >> (k - 1)) & 1
            weight *= bit == 1 ? t[k] : 1 - t[k]
            n += (bit == 1 ? hi[k] : lo[k]) * g.strides[k]
        end
        iszero(weight) && continue
        if F isa AbstractVector
            acc += weight * F[n]
        else
            acc .+= weight .* view(F, :, n)
        end
    end
    return acc
end

"""
    value_at(solution, x) -> Float64

For a [`MechanicalHJBProblem`](@ref): the value function at the state
`x = [θ₁, θ₂, p₁, p₂]`, interpolated multilinearly (angles taken modulo 2π,
momenta clamped to the grid's box).
"""
value_at(sol::HJBSolution{MechanicalHJBProblem}, x::AbstractVector) = _interpolate(sol.problem.grid, sol.values, x)

"""
    control_at(solution, x) -> Vector

For a [`MechanicalHJBProblem`](@ref): the optimal joint torques at the state
`x`, the grid policy interpolated multilinearly and clamped to the torque
bounds.
"""
function control_at(sol::HJBSolution{MechanicalHJBProblem}, x::AbstractVector)
    τ = _interpolate(sol.problem.grid, sol.controls, x)
    return clamp.(τ, .-collect(sol.problem.torque_bound), collect(sol.problem.torque_bound))
end

"""
    closed_loop(solution, x0; duration = 10, dt = 0.005) -> (; times, states, torques, cost)

Simulate the arm from the state `x0 = [θ₁, θ₂, p₁, p₂]` under the feedback
of `solution`, sampled and held over each step of length `dt` (a digital
controller), the dynamics integrated by the classical Runge–Kutta method.
Returns the times, the states (`4 × N`, angles unwrapped), the torques
(`2 × N`, the last column repeated) and the discounted cost accumulated over
`duration` by the trapezoidal rule.
"""
function closed_loop(sol::HJBSolution{MechanicalHJBProblem}, x0::AbstractVector; duration::Real = 10.0,
        dt::Real = 0.005)
    @argcheck length(x0) == 4 && duration > 0 && dt > 0
    prob = sol.problem
    m = _kernel_model(prob)
    steps = ceil(Int, duration / dt)
    times = collect((0:steps) .* dt)
    states = zeros(4, steps + 1)
    torques = zeros(2, steps + 1)
    states[:, 1] .= x0
    f(x, τ) = _flow(m, x, τ)
    cost = 0.0
    for i in 1:steps
        x = Tuple(states[:, i])
        τ = Tuple(control_at(sol, collect(x)))
        torques[:, i] .= τ
        k1 = f(x, τ)
        k2 = f(x .+ dt / 2 .* k1, τ)
        k3 = f(x .+ dt / 2 .* k2, τ)
        k4 = f(x .+ dt .* k3, τ)
        xn = x .+ dt / 6 .* (k1 .+ 2 .* k2 .+ 2 .* k3 .+ k4)
        states[:, i + 1] .= xn
        cost += dt / 2 * (exp(-prob.discount * times[i]) * _cost(m, x, τ) + exp(-prob.discount * times[i + 1]) * _cost(m, xn, τ))
    end
    torques[:, end] .= torques[:, end - 1]
    return (; times, states, torques, cost)
end

end # module MechanicalHJB
