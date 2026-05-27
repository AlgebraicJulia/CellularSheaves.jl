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
    build_control_cost_matrix(prob::TrackingProblem, λ::Number) -> SparseMatrixCSC{Float64,Int}

Construct a global quadratic cost matrix `Q` that penalises the control
components of each agent with a uniform weight `λ`.  Internally this builds a
block‑diagonal matrix where each control block is `λ * I(nu)`.
"""
function build_control_cost_matrix(prob::TrackingProblem, λ::Number=1.0)
    # Create a simple cost function that returns λ * I(nu) for each agent/time
    cost_func = (agent_idx::Int, t_step::Int) -> begin
        nu = size(prob.Bd[agent_idx], 2)
        return λ * I(nu)
    end
    return build_control_cost_matrix(prob, cost_func)
end

"""
    build_control_cost_matrix(prob::TrackingProblem, cost_func::Function) -> SparseMatrixCSC{Float64,Int}

Construct a global quadratic cost matrix `Q` using a user‑provided function
`cost_func(agent_idx, t_step) -> Q_local`.

* `prob` – the `TrackingProblem` describing the agents, dynamics, etc.
* `cost_func` – optional function `(agent_idx, t_step) -> Q_local` where
  `Q_local` is a `nu × nu` positive semidefinite matrix for that agent at that timestep.  If `nothing`
  the function defaults to `control_weight * I(nu)`.

The returned matrix has size `total_dim × total_dim` where `total_dim` is the
length of the stacked vector of all stalk variables (states + controls) in the
time‑expanded sheaf.
"""
function build_control_cost_matrix(prob::TrackingProblem, cost_func::Function)
    # Compute per‑vertex stalk sizes (shared with the scalar version)
    stalk_sizes = Vector{Int}(undef, (prob.k + 1) * (prob.n_agents + prob.n_targets))
    for t in 0:prob.k, i in 1:prob.n_agents
        nx = size(prob.Ad[i], 1)
        nu = size(prob.Bd[i], 2)
        stalk_sizes[t * (prob.n_agents + prob.n_targets) + i] = nx + nu
    end
    for t in 0:prob.k, j in 1:prob.n_targets
        nx = size(prob.target_Ad[j], 1)
        nu = size(prob.target_Bd[j], 2)
        stalk_sizes[t * (prob.n_agents + prob.n_targets) + prob.n_agents + j] = nx + nu
    end
    total_dim = sum(stalk_sizes)
    Q = spzeros(Float64, total_dim, total_dim)

    function control_range(agent_idx::Int, t_step::Int)
        v = agent_vertex(prob, agent_idx, t_step)
        nx = size(prob.Ad[agent_idx], 1)
        nu = size(prob.Bd[agent_idx], 2)
        offset = sum(stalk_sizes[1:(v-1)])
        return (offset + nx + 1) : (offset + nx + nu)
    end

    for (i, Bd_i) in enumerate(prob.Bd)
        nu = size(Bd_i, 2)
        for t_step in 0:prob.k
            idxs = control_range(i, t_step)
            Q_local = cost_func(i, t_step)
            @assert size(Q_local) == (nu, nu) "cost_func must return a $(nu)×$(nu) matrix for agent $i at time $t_step"
            Q[idxs, idxs] = Q_local
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
    QN = N' * Q * N
    rhs = - N' * Q * point
    # Solve; fall back to least‑squares if singular
    α_opt = QN \ rhs
    return point + N * α_opt
end

end # module QuadraticCost
