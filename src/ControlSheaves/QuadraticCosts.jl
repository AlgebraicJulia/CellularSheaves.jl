"""
    QuadraticCosts.jl

Utility functions for building and solving a quadratic control cost on the
expanded tracking sheaf.  The idea is to penalise only the control components of
each agent's stalk (the `nu` last entries), optionally using a user‑provided
function `cost_func(agent, time) -> Q_local` that returns a positive‑semidefinite
matrix of size `nu × nu`.

The main workflow is:
1. `build_control_cost_matrix(prob, cost)` constructs a
   global sparse matrix `Q` of dimension `total_dim × total_dim` where `total_dim`
   is the length of the vector of all stalk variables.  All entries are zero
   except on the control sub‑blocks.
2. `solve_quadratic_on_basis(z_harmonic, null_basis, Q)` reduces the
   quadratic problem onto the nullspace of the Laplacian constraints and returns
   the optimal full‑state vector.

These functions are deliberately pure and have no side‑effects, making them
easy to unit‑test.
"""

module QuadraticCosts

using SparseArrays
using LinearAlgebra
using ..MultiAgentTracking: TrackingProblem, agent_vertex, extract_state_trajectories,
    extract_target_trajectories, build_time_expanded_tracking_sheaf
using ....NetworkSheaves: sheaf_laplacian_matrix_direct

export build_control_cost_matrix, solve_quadratic_on_basis

"""
    build_control_cost_matrix(prob::TrackingProblem, stalks, λ::Number) -> SparseMatrixCSC{Float64,Int}

Construct a global quadratic cost matrix `Q` that penalises the control
components of each agent with a uniform weight `λ`.  Internally this builds a
block‑diagonal matrix where each control block is `λ * I(nu)`.
"""
function build_control_cost_matrix(prob::TrackingProblem, stalks::AbstractVector{<:Integer}, λ::Number=1.0)
    return build_control_cost_matrix(prob, stalks,
        (agent_idx::Int, t_step::Int) -> λ * I(size(prob.Bd[agent_idx], 2)))
end

"""
    build_control_cost_matrix(prob::TrackingProblem, stalks, cost_func::Function) -> SparseMatrixCSC{Float64,Int}

Construct a global quadratic cost matrix `Q` using a user‑provided function
`cost_func(agent_idx, t_step) -> Q_local`.

* `prob` – the `TrackingProblem` describing the agents, dynamics, etc.
* `stalks` – the sheaf's `vertex_stalks`, giving the DOF layout `Q` is built for.
* `cost_func` – a function `(agent_idx, t_step) -> Q_local` where
  `Q_local` is a `nu × nu` positive semidefinite matrix for that agent at that timestep.

The returned matrix has size `total_dim × total_dim` where `total_dim = sum(stalks)`.
The control block is read from `stalks` (the trailing `nu` DOFs of each agent
vertex), so an agent vertex that carries state only — e.g. the `t = 0` vertex of
the MPC-window sheaf — has no control DOFs and is skipped automatically.
"""
function build_control_cost_matrix(prob::TrackingProblem, stalks::AbstractVector{<:Integer}, cost_func::Function)
    offsets = [0; cumsum(stalks)]
    Q = spzeros(Float64, last(offsets), last(offsets))
    for i in 1:prob.n_agents
        nx = size(prob.Ad[i], 1); nu = size(prob.Bd[i], 2)
        for t_step in 0:prob.k
            v = agent_vertex(prob, i, t_step)
            stalks[v] == nx + nu || continue   # state-only vertex: no control to penalise
            rng = (offsets[v] + nx + 1):(offsets[v] + nx + nu)
            Q_local = cost_func(i, t_step)
            @assert size(Q_local) == (nu, nu) "cost_func must return a $(nu)×$(nu) matrix for agent $i at time $t_step"
            Q[rng, rng] = Q_local
        end
    end
    return Q
end

"""
    solve_quadratic_on_basis(point::AbstractVector, basis::AbstractMatrix,
                                 Q::AbstractMatrix) -> Vector{Float64}

Given a solution `point`, a matrix `N` whose columns span
the space of constraints, and a global quadratic cost matrix
`Q`, solve the reduced system

```
(N' * Q * N) * α = -N' * Q * point
```

and return the optimal full‑state vector `z_opt = point + N * α`.
If the nullspace dimension is zero the original solution is returned unchanged.
"""
function solve_quadratic_on_basis(
    point::AbstractVector,
    basis::AbstractMatrix,
    Q::AbstractMatrix,
)
    nd = size(basis, 2)
    if nd == 0
        return copy(point)
    end
    N = basis
    QN = Matrix(N' * Q * N)
    rhs = -(N' * Q * point)
    α_opt = svd(QN) \ rhs
    return point + N * α_opt
end

end # module QuadraticCosts
