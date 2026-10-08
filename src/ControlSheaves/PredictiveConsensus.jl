"""
    PredictiveConsensus

Optimal consensus control of a fleet of planar double integrators, solved by
exchanging *predicted trajectories* between neighbours.

Each agent ``i`` has position ``q_i \\in \\mathbb R^2``, velocity
``v_i \\in \\mathbb R^2`` and acceleration input ``u_i \\in \\mathbb R^2``. The
fleet is coupled only through a pinned coordination sheaf (a
[`CoordinationScenario`](@ref CellularSheaves.ControlSheaves.CoordinationBenchmarks.CoordinationScenario)):
its Laplacian block ``\\mathcal H`` and target coupling ``B`` define the
formation ``q^\\star = \\mathcal H^{-1} B p`` (the harmonic extension of the
target positions ``p``) and the formation energy
``\\tfrac12 (q - q^\\star)^\\top \\mathcal H (q - q^\\star)``.
The team solves the finite-horizon linear–quadratic problem

```math
\\min_{x, u} \\sum_{t=1}^{T} \\Big[\\tfrac{w}{2}\\|q(t) - q^\\star\\|_{\\mathcal H}^2
+ \\tfrac{c}{2}\\|v(t)\\|^2\\Big] + \\sum_{t=0}^{T-1} \\tfrac{r}{2}\\|u(t)\\|^2
\\quad\\text{s.t.}\\quad x_i(t+1) = A x_i(t) + B_u u_i(t),\\; \\|u_i(t)\\| \\le \\bar u .
```

Without the bound this is a linear–quadratic regulator on the stacked state,
solved exactly by a Riccati recursion ([`JointRiccati`](@ref)). The
distributed method, [`PredictedTrajectorySweeps`](@ref), is block coordinate
descent on agents: agent ``i`` holds its neighbours' predicted position
trajectories ``\\hat q_j(\\cdot)`` fixed and solves its own trajectory QP,
whose only dependence on the neighbours is the linear term
``w\\,\\mathcal H_{ij}\\hat q_j(t)``. In sheaf language the predicted
trajectories are a 0-cochain of the coordination sheaf tensored with
``\\mathbb R^{T}`` (one copy per time step), and the update is the
zero-overlap multiplicative Schwarz method of
[`SchwarzMethods`](@ref CellularSheaves.NetworkSheaves.SchwarzMethods) on the
agents' trajectory spaces. Each local problem is a conic QP solved by the
Mumblebee interior-point method; neighbour predictions change only its linear
term and the initial state only its right-hand side, so one factorization
structure is set up per agent and reused with `reinit!`.

Because the team objective is convex and the constraints are separable by
agent, the colored (Gauss–Seidel) sweep decreases the team cost monotonically
and converges to the team optimum. The fully parallel (Jacobi) sweep converges
when damped by ``\\theta \\le 1/\\chi``, where ``\\chi`` is the number of
colors, since ``\\mathcal H \\preceq \\chi\\,\\mathrm{blockdiag}(\\mathcal H)``.

Two closed-loop baselines are provided for comparison:
[`RiccatiFeedback`](@ref) (the centralized infinite-horizon LQR, the optimum)
and [`SecondOrderDiffusion`](@ref), the second-order consensus law
``u_i = -k_p \\eta_i - k_v v_i`` with ``\\eta = \\mathcal H (q - q^\\star)``
the sheaf disagreement (Ren & Beard, *Distributed Consensus in Multi-vehicle
Cooperative Control*, 2008; arXiv:2512.24886 for the sheaf form).
"""
module PredictiveConsensus

using ..CoordinationBenchmarks: CoordinationScenario
using CellularSheaves.NetworkSheaves.SchwarzMethods: SchwarzSweep, MulticolorSweep, ParallelSweep
using ArgCheck
using CommonSolve
using CommonSolve: init, solve, solve!
using Graphs
using LinearAlgebra
using Mumblebee.IPM: IPMProblem, IPMSettings, AbstractCone, CofreeCone, SecondOrderCone,
    OPTIMAL, NEAR_OPTIMAL
using Mumblebee.IPM: reinit!
using SparseArrays

export PlanarDoubleIntegrator, ConsensusLQ, ConsensusProblem, ConsensusPlan,
    formation, team_cost, coasting_plan,
    JointRiccati, CentralizedQP, PredictedTrajectorySweeps, SweepWorkspace,
    ConsensusController, RiccatiFeedback, SecondOrderDiffusion, RecedingHorizon,
    ClosedLoop, rollout

# ===========================================================================
# problem data
# ===========================================================================

"""
    PlanarDoubleIntegrator(dt)

The planar double integrator ``\\ddot q = u``, ``q, u \\in \\mathbb R^2``,
discretized exactly under a zero-order hold of length `dt`. The state is
``x = (q_1, q_2, v_1, v_2)`` and

```math
A = \\begin{bmatrix} I & \\Delta t\\, I \\\\ 0 & I \\end{bmatrix}, \\qquad
B_u = \\begin{bmatrix} \\tfrac12 \\Delta t^2 I \\\\ \\Delta t\\, I \\end{bmatrix}.
```
"""
struct PlanarDoubleIntegrator
    dt::Float64
    A::Matrix{Float64}
    B::Matrix{Float64}
end

function PlanarDoubleIntegrator(dt::Real)
    @argcheck dt > 0
    I2 = Matrix{Float64}(I, 2, 2)
    A = [I2 dt*I2; zeros(2, 2) I2]
    B = [dt^2 / 2 * I2; dt * I2]
    return PlanarDoubleIntegrator(dt, A, B)
end

"""
    ConsensusLQ(scenario, plant, horizon; targets, position_weight = 1,
                velocity_weight = 0.1, control_weight = 0.1, control_bound = Inf)

A consensus problem for a fleet of [`PlanarDoubleIntegrator`](@ref)s on the
pinned coordination sheaf of `scenario` (stalk dimension 2), over `horizon`
steps. `targets` is the `2 × ntargets` matrix of (fixed) target positions;
the formation is ``q^\\star = \\mathcal H^{-1} B p``. The weights are ``w``,
``c`` and ``r`` in the module docstring, and `control_bound` is the
Euclidean thrust limit ``\\bar u`` (a second-order cone constraint).
"""
struct ConsensusLQ
    scenario::CoordinationScenario
    plant::PlanarDoubleIntegrator
    horizon::Int
    targets::Matrix{Float64}
    position_weight::Float64
    velocity_weight::Float64
    control_weight::Float64
    control_bound::Float64
    formation::Vector{Float64}                       # q⋆, 2 per agent
    coupling::Vector{Vector{Pair{Int, Matrix{Float64}}}}   # agent i => [j => H_ij]
end

function ConsensusLQ(scenario::CoordinationScenario, plant::PlanarDoubleIntegrator, horizon::Integer;
        targets::AbstractMatrix, position_weight::Real = 1.0, velocity_weight::Real = 0.1,
        control_weight::Real = 0.1, control_bound::Real = Inf)
    @argcheck scenario.dim == 2 "the planar double integrator needs a stalk of dimension 2"
    @argcheck size(targets) == (2, scenario.ntargets)
    @argcheck horizon >= 1
    @argcheck position_weight > 0 && velocity_weight >= 0 && control_weight > 0
    @argcheck control_bound > 0
    H = scenario.H
    qstar = cholesky(Symmetric(Matrix(H))) \ (scenario.Bmat * vec(Matrix{Float64}(targets)))
    coupling = map(1:scenario.nagents) do i
        ri = _pos(i)
        [j => Matrix(H[ri, _pos(j)]) for j in [i; scenario.agent_nbrs[i]]]
    end
    return ConsensusLQ(scenario, plant, horizon, Matrix{Float64}(targets), position_weight,
        velocity_weight, control_weight, control_bound, qstar, coupling)
end

_pos(i) = (2i - 1):(2i)
_bounded(lq::ConsensusLQ) = isfinite(lq.control_bound)
nagents(lq::ConsensusLQ) = lq.scenario.nagents

"""
    formation(lq) -> Matrix

The formation ``q^\\star`` as a `2 × nagents` matrix.
"""
formation(lq::ConsensusLQ) = reshape(copy(lq.formation), 2, nagents(lq))

"""
    ConsensusProblem(lq, x0)

The consensus problem `lq` from initial states `x0`, a `4 × nagents` matrix
whose columns are ``(q_i, v_i)``.
"""
struct ConsensusProblem
    lq::ConsensusLQ
    x0::Matrix{Float64}
    function ConsensusProblem(lq::ConsensusLQ, x0::AbstractMatrix)
        @argcheck size(x0) == (4, nagents(lq))
        return new(lq, Matrix{Float64}(x0))
    end
end

"""
    ConsensusPlan

An open-loop plan for the whole fleet.

# Fields
- `states`: `4 × (T+1) × nagents`, the state ``x_i(t)`` at `[:, t+1, i]`.
- `controls`: `2 × T × nagents`, the input ``u_i(t)`` at `[:, t+1, i]`.
- `cost`: the team cost ([`team_cost`](@ref)).
- `iterations`: sweeps used (zero for direct methods).
- `residuals`: per sweep, the largest change in any agent's plan.
- `costs`: per sweep, the team cost after the sweep.
- `converged`: whether the residual fell below the tolerance.
"""
struct ConsensusPlan
    states::Array{Float64, 3}
    controls::Array{Float64, 3}
    cost::Float64
    iterations::Int
    residuals::Vector{Float64}
    costs::Vector{Float64}
    converged::Bool
end

"""
    team_cost(lq, states, controls) -> Float64

The team objective of the module docstring for a fleet plan with `states`
`4 × (T+1) × N` and `controls` `2 × T × N` (any `T`, so closed-loop runs are
scored the same way).
"""
function team_cost(lq::ConsensusLQ, states::AbstractArray{<:Real, 3}, controls::AbstractArray{<:Real, 3})
    w, c, r = lq.position_weight, lq.velocity_weight, lq.control_weight
    N = nagents(lq)
    J = 0.0
    e = zeros(2N)
    for t in 2:size(states, 2)
        for i in 1:N
            e[_pos(i)] .= view(states, 1:2, t, i) .- view(lq.formation, _pos(i))
            J += c / 2 * sum(abs2, view(states, 3:4, t, i))
        end
        J += w / 2 * dot(e, lq.scenario.H, e)
    end
    return J + r / 2 * sum(abs2, controls)
end

"""
    coasting_plan(problem) -> (states, controls)

The plan with zero input: each agent coasts at its initial velocity. Used as
the initial predicted trajectories.
"""
function coasting_plan(prob::ConsensusProblem)
    lq = prob.lq
    T, N, A = lq.horizon, nagents(lq), lq.plant.A
    states = zeros(4, T + 1, N)
    for i in 1:N
        states[:, 1, i] .= view(prob.x0, :, i)
        for t in 1:T
            states[:, t + 1, i] .= A * view(states, :, t, i)
        end
    end
    return states, zeros(2, T, N)
end

# ===========================================================================
# the trajectory QP over a set of agents
# ===========================================================================

# Variables of the trajectory QP over the agents `agents` (local index a):
# for each agent, T stalks (x(t), u(t)) of size 6 and one stalk x(T) of
# size 4; then, if the input is bounded, one second-order-cone slack
# (ū, u(t)) of size 3 per agent and step. Rows: per agent, x(0) = x0 (4),
# the dynamics (4 per step) and, if bounded, the slack identities (3 per step).
struct TrajectoryLayout
    agents::Vector{Int}
    horizon::Int
    bounded::Bool
end

_block(l::TrajectoryLayout) = 6l.horizon + 4
_x(l::TrajectoryLayout, a, t) = (a - 1) * _block(l) + 6t .+ (1:4)
_q(l::TrajectoryLayout, a, t) = (a - 1) * _block(l) + 6t .+ (1:2)
_v(l::TrajectoryLayout, a, t) = (a - 1) * _block(l) + 6t .+ (3:4)
_u(l::TrajectoryLayout, a, t) = (a - 1) * _block(l) + 6t .+ (5:6)
_s(l::TrajectoryLayout, a, t) = length(l.agents) * _block(l) + ((a - 1) * l.horizon + t) * 3 .+ (1:3)
_rows(l::TrajectoryLayout) = 4 + 4l.horizon + (l.bounded ? 3l.horizon : 0)
_row0(l::TrajectoryLayout, a) = (a - 1) * _rows(l) .+ (1:4)
_rowdyn(l::TrajectoryLayout, a, t) = (a - 1) * _rows(l) + 4 + 4t .+ (1:4)
_rowsoc(l::TrajectoryLayout, a, t) = (a - 1) * _rows(l) + 4 + 4l.horizon + 3t .+ (1:3)

function _stalks(l::TrajectoryLayout)
    k, T = length(l.agents), l.horizon
    sizes = repeat([fill(6, T); 4], k)
    cones = AbstractCone[CofreeCone() for _ in sizes]
    if l.bounded
        append!(sizes, fill(3, k * T))
        append!(cones, [SecondOrderCone() for _ in 1:(k * T)])
    end
    return sizes, cones
end

# The trajectory QP for `agents`, with the predicted positions of every other
# agent entering the linear term. Returns the problem and its layout.
function _trajectory_qp(lq::ConsensusLQ, agents::Vector{Int}, x0::AbstractMatrix, states::AbstractArray)
    T = lq.horizon
    l = TrajectoryLayout(agents, T, _bounded(lq))
    local_index = Dict(i => a for (a, i) in enumerate(agents))
    w, c, r = lq.position_weight, lq.velocity_weight, lq.control_weight
    sizes, cones = _stalks(l)
    n = sum(sizes)
    m = length(agents) * _rows(l)
    Qi, Qj, Qv = Int[], Int[], Float64[]
    entry!(ii, jj, vv, rows, cols, M) = for (b, cj) in enumerate(cols), (a, ri) in enumerate(rows)
        push!(ii, ri); push!(jj, cj); push!(vv, M[a, b])
    end
    for k in 1:n
        push!(Qi, k); push!(Qj, k); push!(Qv, 0.0)          # every stalk keeps its diagonal
    end
    I2 = Matrix{Float64}(I, 2, 2)
    for (a, i) in enumerate(agents)
        for t in 1:T
            for (j, Hij) in lq.coupling[i]
                b = get(local_index, j, 0)
                b == 0 || entry!(Qi, Qj, Qv, _q(l, a, t), _q(l, b, t), w * Hij)
            end
            entry!(Qi, Qj, Qv, _v(l, a, t), _v(l, a, t), c * I2)
            entry!(Qi, Qj, Qv, _u(l, a, t - 1), _u(l, a, t - 1), r * I2)
        end
    end
    Q = sparse(Qi, Qj, Qv, n, n)

    Bi, Bj, Bv = Int[], Int[], Float64[]
    A, Bu = lq.plant.A, lq.plant.B
    I4 = Matrix{Float64}(I, 4, 4)
    for a in eachindex(agents)
        entry!(Bi, Bj, Bv, _row0(l, a), _x(l, a, 0), I4)
        for t in 0:(T - 1)
            rows = _rowdyn(l, a, t)
            entry!(Bi, Bj, Bv, rows, _x(l, a, t + 1), I4)
            entry!(Bi, Bj, Bv, rows, _x(l, a, t), -A)
            entry!(Bi, Bj, Bv, rows, _u(l, a, t), -Bu)
            if l.bounded
                rows = _rowsoc(l, a, t)
                s = _s(l, a, t)
                entry!(Bi, Bj, Bv, rows, s, Matrix{Float64}(I, 3, 3))
                entry!(Bi, Bj, Bv, rows[2:3], _u(l, a, t), -I2)
            end
        end
    end
    B = sparse(Bi, Bj, Bv, m, n)
    f = zeros(n)
    _linear_term!(f, lq, l, states)
    g = zeros(m)
    _rhs!(g, lq, l, x0)
    return IPMProblem(Q, B, f, g, 0.0, cones, sizes), l
end

# f = w (B p − Σ_{j outside} H_ij q̂_j(t)) on the positions, using H q⋆ = B p.
function _linear_term!(f, lq::ConsensusLQ, l::TrajectoryLayout, states)
    fill!(f, 0.0)
    inside = Set(l.agents)
    w = lq.position_weight
    for (a, i) in enumerate(l.agents)
        Hqstar = zeros(2)
        for (j, Hij) in lq.coupling[i]
            Hqstar .+= Hij * view(lq.formation, _pos(j))
        end
        for t in 1:l.horizon
            ft = w .* Hqstar
            for (j, Hij) in lq.coupling[i]
                j in inside && continue
                mul!(ft, Hij, view(states, 1:2, t + 1, j), -w, 1.0)
            end
            f[_q(l, a, t)] .= ft
        end
    end
    return f
end

function _rhs!(g, lq::ConsensusLQ, l::TrajectoryLayout, x0)
    fill!(g, 0.0)
    for (a, i) in enumerate(l.agents)
        g[_row0(l, a)] .= view(x0, :, i)
        if l.bounded
            for t in 0:(l.horizon - 1)
                g[first(_rowsoc(l, a, t))] = lq.control_bound
            end
        end
    end
    return g
end

function _extract!(states, controls, l::TrajectoryLayout, p)
    for (a, i) in enumerate(l.agents)
        for t in 0:l.horizon
            states[:, t + 1, i] .= view(p, _x(l, a, t))
            t < l.horizon && (controls[:, t + 1, i] .= view(p, _u(l, a, t)))
        end
    end
    return states, controls
end

function _check(result)
    result.status in (OPTIMAL, NEAR_OPTIMAL) ||
        error("trajectory QP not solved: status $(result.status)")
    return result
end

# ===========================================================================
# direct solvers
# ===========================================================================

"""
    JointRiccati()

The exact team optimum of an unbounded [`ConsensusLQ`](@ref) by the
finite-horizon Riccati recursion on the stacked state of all agents,
``O(T (4N)^3)`` work. The centralized reference.
"""
struct JointRiccati end

# Stacked LQR data in the shifted coordinates (q − q⋆, v), agent-major.
function _stacked(lq::ConsensusLQ)
    N = nagents(lq)
    A = kron(Matrix{Float64}(I, N, N), lq.plant.A)
    B = kron(Matrix{Float64}(I, N, N), lq.plant.B)
    Q = zeros(4N, 4N)
    for i in 1:N
        for (j, Hij) in lq.coupling[i]
            Q[4(i - 1) .+ (1:2), 4(j - 1) .+ (1:2)] .= lq.position_weight .* Hij
        end
        Q[4(i - 1) .+ (3:4), 4(i - 1) .+ (3:4)] .= lq.velocity_weight .* Matrix{Float64}(I, 2, 2)
    end
    R = lq.control_weight * Matrix{Float64}(I, 2N, 2N)
    return A, B, Q, R
end

function _shift(lq::ConsensusLQ, x::AbstractMatrix)
    X = vec(copy(x))
    for i in 1:nagents(lq)
        X[4(i - 1) .+ (1:2)] .-= view(lq.formation, _pos(i))
    end
    return X
end

function _riccati_step(A, B, Q, R, P)
    K = (R + B' * P * B) \ (B' * P * A)
    Acl = A - B * K
    return K, Symmetric(Q + K' * R * K + Acl' * P * Acl)
end

function CommonSolve.solve(prob::ConsensusProblem, ::JointRiccati)
    lq = prob.lq
    @argcheck !_bounded(lq) "JointRiccati needs an unbounded input; use CentralizedQP"
    T, N = lq.horizon, nagents(lq)
    A, B, Q, R = _stacked(lq)
    gains = Vector{Matrix{Float64}}(undef, T)
    P = Symmetric(Q)
    for t in T:-1:1
        gains[t], P = _riccati_step(A, B, Q, R, Matrix(P))
    end
    states = zeros(4, T + 1, N)
    controls = zeros(2, T, N)
    X = _shift(lq, prob.x0)
    for t in 1:T
        U = -gains[t] * X
        controls[:, t, :] .= reshape(U, 2, N)
        states[:, t, :] .= reshape(X, 4, N)
        X = A * X + B * U
    end
    states[:, T + 1, :] .= reshape(X, 4, N)
    for i in 1:N, t in 1:(T + 1)
        states[1:2, t, i] .+= view(lq.formation, _pos(i))
    end
    return ConsensusPlan(states, controls, team_cost(lq, states, controls), 0, Float64[], Float64[], true)
end

"""
    CentralizedQP(; settings = IPMSettings{Float64}())

The team optimum as one conic QP over all agents' trajectories, solved by the
Mumblebee interior-point method. Handles the control bound; the reference for
bounded problems.
"""
Base.@kwdef struct CentralizedQP{S <: IPMSettings}
    settings::S = IPMSettings{Float64}()
end

function CommonSolve.solve(prob::ConsensusProblem, alg::CentralizedQP)
    lq = prob.lq
    states, controls = coasting_plan(prob)
    qp, l = _trajectory_qp(lq, collect(1:nagents(lq)), prob.x0, states)
    result = _check(solve(qp, alg.settings))
    _extract!(states, controls, l, result.p)
    return ConsensusPlan(states, controls, team_cost(lq, states, controls), 0, Float64[], Float64[], true)
end

# ===========================================================================
# predicted-trajectory sweeps
# ===========================================================================

"""
    PredictedTrajectorySweeps(; sweep = MulticolorSweep(), damping = nothing,
                              tol = 1e-8, maxiter = 500,
                              settings = IPMSettings{Float64}(), warm_start = false)

Distributed solver: every agent repeatedly re-plans its own trajectory against
its neighbours' latest predicted trajectories.

- `sweep = MulticolorSweep()`: agents of one color of the communication graph
  re-plan in parallel, colors in turn (block Gauss–Seidel). One neighbour
  exchange per color.
- `sweep = ParallelSweep()`: all agents re-plan at once and move a fraction
  `damping` (default ``1/\\chi``, ``\\chi`` the number of colors) of the way
  to their new plan (damped block Jacobi). One exchange per sweep.

Each local problem is solved by the Mumblebee IPM, reusing its symbolic
factorization across sweeps (`reinit!` with the new linear term). With
`warm_start` the previous local solution is passed as the starting point;
this is off by default because with a control bound that solution lies on the
cone boundary, where the interior-point iteration fails. Stops when no agent's plan
changes by more than `tol` (relative to the plan's scale) in a sweep.
"""
Base.@kwdef struct PredictedTrajectorySweeps{W <: SchwarzSweep, S <: IPMSettings}
    sweep::W = MulticolorSweep()
    damping::Union{Nothing, Float64} = nothing
    tol::Float64 = 1e-8
    maxiter::Int = 500
    settings::S = IPMSettings{Float64}()
    warm_start::Bool = false
end

"""
    SweepWorkspace

The state of a [`PredictedTrajectorySweeps`](@ref) run: one IPM solver per
agent (set up once), the current fleet plan (each agent's latest prediction),
and the coloring. Created by `init(problem, algorithm)` and advanced by
`solve!(workspace)`; [`RecedingHorizon`](@ref) reuses it across time steps.
"""
mutable struct SweepWorkspace{A <: PredictedTrajectorySweeps, V}
    lq::ConsensusLQ
    algorithm::A
    x0::Matrix{Float64}
    solvers::Vector{V}
    layouts::Vector{TrajectoryLayout}
    f::Vector{Vector{Float64}}
    g::Vector{Vector{Float64}}
    primal::Vector{Vector{Float64}}
    states::Array{Float64, 3}
    controls::Array{Float64, 3}
    next_states::Array{Float64, 3}
    next_controls::Array{Float64, 3}
    colors::Vector{Vector{Int}}
    damping::Float64
end

function CommonSolve.init(prob::ConsensusProblem, alg::PredictedTrajectorySweeps;
        states = nothing, controls = nothing)
    lq = prob.lq
    N = nagents(lq)
    s0, c0 = coasting_plan(prob)
    states = states === nothing ? s0 : copy(states)
    controls = controls === nothing ? c0 : copy(controls)
    built = [_trajectory_qp(lq, [i], prob.x0, states) for i in 1:N]
    solvers = [init(qp, alg.settings) for (qp, _) in built]
    layouts = [l for (_, l) in built]
    f = [zeros(length(qp.f)) for (qp, _) in built]
    g = [zeros(length(qp.g)) for (qp, _) in built]
    coloring = Graphs.greedy_color(lq.scenario.agent_graph)
    colors = [findall(==(k), coloring.colors) for k in 1:coloring.num_colors]
    damping = something(alg.damping, 1 / length(colors))
    @argcheck 0 < damping <= 1
    return SweepWorkspace(lq, alg, copy(prob.x0), solvers, layouts, f, g,
        [Float64[] for _ in 1:N], states, controls, copy(states), copy(controls), colors, damping)
end

CommonSolve.solve(prob::ConsensusProblem, alg::PredictedTrajectorySweeps) = solve!(init(prob, alg))

# Mumblebee's `reinit!` (at 9d308da) keeps the barrier Hessian H of the last
# iterate, and with an empty history the KKT augmentation is scaled by 1/‖H‖.
# After a solve that ended near a cone boundary ‖H‖ is huge and the next solve
# fails. A fresh solver has H = Q; restore that until upstream does.
_reset_hessian!(solver) = copyto!(solver.H, solver.Q)

# Re-plan agent i against the current predictions, writing into (states, controls).
function _replan!(ws::SweepWorkspace, i, states, controls)
    l = ws.layouts[i]
    _linear_term!(ws.f[i], ws.lq, l, ws.states)
    _rhs!(ws.g[i], ws.lq, l, ws.x0)
    p0 = ws.algorithm.warm_start && !isempty(ws.primal[i]) ? ws.primal[i] : nothing
    _reset_hessian!(ws.solvers[i])
    reinit!(ws.solvers[i]; f = ws.f[i], g = ws.g[i], p0)
    result = _check(solve!(ws.solvers[i]))
    ws.primal[i] = result.p
    _extract!(states, controls, l, result.p)
    return nothing
end

_change(a, b) = maximum(abs, a - b; init = 0.0)

function _sweep!(ws::SweepWorkspace, ::MulticolorSweep)
    old_states, old_controls = copy(ws.states), copy(ws.controls)
    for color in ws.colors
        Threads.@threads for i in color
            _replan!(ws, i, ws.states, ws.controls)
        end
    end
    return max(_change(ws.states, old_states), _change(ws.controls, old_controls))
end

function _sweep!(ws::SweepWorkspace, ::ParallelSweep)
    Threads.@threads for i in 1:nagents(ws.lq)
        _replan!(ws, i, ws.next_states, ws.next_controls)
    end
    θ = ws.damping
    change = max(_change(ws.states, ws.next_states), _change(ws.controls, ws.next_controls))
    ws.states .+= θ .* (ws.next_states .- ws.states)
    ws.controls .+= θ .* (ws.next_controls .- ws.controls)
    return θ * change
end

function CommonSolve.solve!(ws::SweepWorkspace; maxiter::Integer = ws.algorithm.maxiter)
    residuals = Float64[]
    costs = Float64[]
    converged = false
    for _ in 1:maxiter
        change = _sweep!(ws, ws.algorithm.sweep)
        scale = max(1.0, maximum(abs, ws.states))
        push!(residuals, change / scale)
        push!(costs, team_cost(ws.lq, ws.states, ws.controls))
        if residuals[end] <= ws.algorithm.tol
            converged = true
            break
        end
    end
    cost = isempty(costs) ? team_cost(ws.lq, ws.states, ws.controls) : costs[end]
    return ConsensusPlan(copy(ws.states), copy(ws.controls), cost, length(residuals),
        residuals, costs, converged)
end

# ===========================================================================
# closed loop
# ===========================================================================

"""
    ConsensusController

A feedback law for the fleet, applied by [`rollout`](@ref). Implementations:
[`RiccatiFeedback`](@ref), [`SecondOrderDiffusion`](@ref) and
[`RecedingHorizon`](@ref).
"""
abstract type ConsensusController end

"""
    RiccatiFeedback(lq; tol = 1e-12, maxiter = 100_000)

The centralized infinite-horizon LQR for the stacked fleet (the discrete
algebraic Riccati equation, solved by fixed-point iteration),
``u = -K\\,(q - q^\\star, v)``. Optimal for the unbounded problem and fully
centralized: every input depends on every agent's state. With a control bound
the input is clipped to the disc.
"""
struct RiccatiFeedback <: ConsensusController
    lq::ConsensusLQ
    gain::Matrix{Float64}
end

function RiccatiFeedback(lq::ConsensusLQ; tol = 1e-12, maxiter = 100_000)
    A, B, Q, R = _stacked(lq)
    P = Matrix(Q)
    K = zeros(size(B, 2), size(A, 1))
    for _ in 1:maxiter
        K, Pnew = _riccati_step(A, B, Q, R, P)
        done = norm(Pnew - P) <= tol * norm(Pnew)
        P = Matrix(Pnew)
        done && break
    end
    return RiccatiFeedback(lq, K)
end

"""
    SecondOrderDiffusion(lq, kp, kv)

The decentralized second-order consensus law
``u_i = -k_p \\eta_i - k_v v_i``, with ``\\eta = \\mathcal H (q - q^\\star)``
the sheaf disagreement, which each agent computes from its neighbours'
current positions only. With a control bound the input is clipped to the disc.
"""
struct SecondOrderDiffusion <: ConsensusController
    lq::ConsensusLQ
    kp::Float64
    kv::Float64
end

"""
    RecedingHorizon(problem, algorithm; sweeps = 1)

Distributed model-predictive control: at every step the fleet re-plans over
the horizon with [`PredictedTrajectorySweeps`](@ref), starting from the
previous plan shifted by one step, running at most `sweeps` sweeps, and
applies the first input. The local IPM solvers are set up once.
"""
mutable struct RecedingHorizon{W <: SweepWorkspace} <: ConsensusController
    workspace::W
    sweeps::Int
end

function RecedingHorizon(prob::ConsensusProblem, alg::PredictedTrajectorySweeps; sweeps::Integer = 1)
    @argcheck sweeps >= 1
    return RecedingHorizon(init(prob, alg), Int(sweeps))
end

function _clip!(u, bound)
    isfinite(bound) || return u
    for i in axes(u, 2)
        nu = norm(view(u, :, i))
        nu > bound && (u[:, i] .*= bound / nu)
    end
    return u
end

function _control!(u, ctrl::RiccatiFeedback, x, step)
    u .= reshape(-ctrl.gain * _shift(ctrl.lq, x), 2, :)
    return _clip!(u, ctrl.lq.control_bound)
end

function _control!(u, ctrl::SecondOrderDiffusion, x, step)
    lq = ctrl.lq
    for i in 1:nagents(lq)
        η = zeros(2)
        for (j, Hij) in lq.coupling[i]
            η .+= Hij * (view(x, 1:2, j) .- view(lq.formation, _pos(j)))
        end
        u[:, i] .= -ctrl.kp .* η .- ctrl.kv .* view(x, 3:4, i)
    end
    return _clip!(u, lq.control_bound)
end

function _control!(u, ctrl::RecedingHorizon, x, step)
    ws = ctrl.workspace
    if step > 1
        A = ws.lq.plant.A
        ws.states[:, 1:(end - 1), :] .= ws.states[:, 2:end, :]
        ws.controls[:, 1:(end - 1), :] .= ws.controls[:, 2:end, :]
        ws.controls[:, end, :] .= 0
        for i in axes(ws.states, 3)
            ws.states[:, end, i] .= A * view(ws.states, :, size(ws.states, 2) - 1, i)
        end
    end
    ws.x0 .= x
    ws.states[:, 1, :] .= x
    solve!(ws; maxiter = ctrl.sweeps)
    u .= view(ws.controls, :, 1, :)
    return u
end

"""
    ClosedLoop

The result of [`rollout`](@ref): `states` `4 × (steps+1) × N`, `controls`
`2 × steps × N`, the team `cost` of the run, and `controller_seconds`, the
wall time spent computing inputs.
"""
struct ClosedLoop
    states::Array{Float64, 3}
    controls::Array{Float64, 3}
    cost::Float64
    controller_seconds::Float64
end

"""
    rollout(lq, controller, x0, steps) -> ClosedLoop

Simulate the fleet under `controller` for `steps` steps from `x0`.
"""
function rollout(lq::ConsensusLQ, ctrl::ConsensusController, x0::AbstractMatrix, steps::Integer)
    N = nagents(lq)
    @argcheck size(x0) == (4, N)
    A, B = lq.plant.A, lq.plant.B
    states = zeros(4, steps + 1, N)
    controls = zeros(2, steps, N)
    states[:, 1, :] .= x0
    u = zeros(2, N)
    seconds = 0.0
    for k in 1:steps
        x = states[:, k, :]
        seconds += @elapsed _control!(u, ctrl, x, k)
        controls[:, k, :] .= u
        for i in 1:N
            states[:, k + 1, i] .= A * view(x, :, i) .+ B * view(u, :, i)
        end
    end
    return ClosedLoop(states, controls, team_cost(lq, states, controls), seconds)
end

end # module PredictiveConsensus
