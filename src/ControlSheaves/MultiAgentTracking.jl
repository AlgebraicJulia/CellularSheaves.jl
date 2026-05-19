"""
    MultiAgentTracking

Utilities for multi-agent, multi-target coordination via time-expanded cellular sheaves.

Exports types for problem specification (`TrackingEdge`, `TrackingProblem`,
`BobbingTarget`, `ScenarioResult`) and functions to build the sheaf, generate
reference trajectories, solve, and extract results.
"""
module MultiAgentTracking

using LinearAlgebra
using BlockArrays
using ArgCheck

using ...NetworkSheaves: EuclideanSheaf, add_sheaf_edge!, harmonic_extension,
                          sheaf_laplacian_matrix_direct

export TrackingEdge, TrackingProblem, BobbingTarget, ScenarioResult
export trajectory
export selector_matrix, state_projection_matrix
export agent_vertex, target_vertex
export build_time_expanded_tracking_sheaf
export generate_reference_trajectory, extract_state_trajectories, run_scenario

# ---------------------------------------------------------------------------
# Problem-specification types
# ---------------------------------------------------------------------------

"""
    TrackingEdge

One agent–target tracking relationship with its own restriction maps.

Fields:
- `agent_index`: index of the agent in `1:n_agents`.
- `target_index`: index of the target in `1:n_targets`.
- `agent_restriction`: matrix applied to the agent stalk.
- `target_restriction`: matrix applied to the target stalk.

Using a vector of `TrackingEdge`s (rather than a single assignment vector)
allows many-to-many relationships where each pairing can use different
projection matrices.
"""
struct TrackingEdge
    agent_index::Int
    target_index::Int
    agent_restriction::Matrix{Float64}
    target_restriction::Matrix{Float64}
end

"""
    TrackingProblem

All parameters of a multi-agent, multi-target tracking problem encoded as a
time-expanded cellular sheaf.

Fields:
- `n_agents`, `n_targets`: fleet sizes.
- `k`: horizon length; the trajectory has `k+1` timesteps `t = 0, …, k`.
- `Ad`, `Bd`: ZOH-discretised state-space matrices (size `nx × nx`, `nx × nu`).
- `agent_edges`: undirected agent–agent consensus pairs `(i, j)`.
- `tracking_edges`: many-to-many agent–target edges (see `TrackingEdge`).
- `consensus_restriction`: projection used for all agent–agent consensus edges.
- `consensus_timesteps`: timesteps at which consensus edges are added.
- `tracking_timesteps`: timesteps at which tracking edges are added.
- `include_target_dynamics`: add temporal dynamics edges between target vertices.
- `consensus_weight`, `tracking_weight`: Laplacian weights (the restriction maps
  are pre-scaled by their square roots).
"""
struct TrackingProblem
    n_agents::Int
    n_targets::Int
    k::Int
    Ad::Matrix{Float64}
    Bd::Matrix{Float64}
    agent_edges::Vector{Tuple{Int,Int}}
    tracking_edges::Vector{TrackingEdge}
    consensus_restriction::Matrix{Float64}
    consensus_timesteps::Vector{Int}
    tracking_timesteps::Vector{Int}
    include_target_dynamics::Bool
    consensus_weight::Float64
    tracking_weight::Float64
end

# ---------------------------------------------------------------------------
# Target entity type
# ---------------------------------------------------------------------------

"""
    BobbingTarget

Parameters of a vertically oscillating ("bobbing") target.

Fields:
- `y_fixed`: fixed lateral position (m).
- `z_center`: mean altitude (m).
- `z_amplitude`: oscillation amplitude (m).
- `omega`: angular frequency (rad/s).

Call `trajectory(bt, t_range, h, nx, nu)` to materialise the stalk sequence
over any integer time interval.
"""
struct BobbingTarget
    y_fixed::Float64
    z_center::Float64
    z_amplitude::Float64
    omega::Float64
end

"""
    trajectory(bt::BobbingTarget, t_range, h, nx, nu, y_idx, z_idx, zdot_idx)

Return a `Vector` of length `length(t_range)` stalk vectors for `bt`.

Each stalk has dimension `nx + nu`.  The lateral state component at `y_idx`
is set to `bt.y_fixed`; the altitude component at `z_idx` follows a sinusoid
centred at `bt.z_center` with amplitude `bt.z_amplitude` and angular frequency
`bt.omega`; the altitude-rate component at `zdot_idx` is set to the analytic
derivative.  All other components are zero.

`t_range` may be any integer range, enabling evaluation over an arbitrary
sub-interval without recomputing from scratch.
"""
function trajectory(
    bt::BobbingTarget,
    t_range::AbstractVector{<:Integer},
    h::Float64,
    nx::Int,
    nu::Int,
    y_idx::Int,
    z_idx::Int,
    zdot_idx::Int,
)
    @argcheck 1 <= y_idx    <= nx "y index out of range 1:$nx"
    @argcheck 1 <= z_idx    <= nx "z index out of range 1:$nx"
    @argcheck 1 <= zdot_idx <= nx "zdot index out of range 1:$nx"
    return map(t_range) do t
        t_s    = t * h
        z_t    = bt.z_center + bt.z_amplitude * sin(bt.omega * t_s)
        zdot_t = bt.z_amplitude * bt.omega    * cos(bt.omega * t_s)
        x = zeros(nx)
        x[y_idx]    = bt.y_fixed
        x[z_idx]    = z_t
        x[zdot_idx] = zdot_t
        vcat(x, zeros(nu))
    end
end

# ---------------------------------------------------------------------------
# Result type
# ---------------------------------------------------------------------------

"""
    ScenarioResult

Bundle of trajectory data and statistics for one solved scenario.

Fields:
- `label`: scenario name.
- `times`: discrete time vector (length `k+1`).
- `agent_trajs`: `Vector` of `(k+1) × nx` state matrices, one per agent.
- `target_trajs`: `Vector` of stalk-vector sequences, one per target.
- `null_dim`: nullspace dimension of the restricted Laplacian.
- `residual`: Laplacian energy `sqrt(z' * L * z)`.
- `y_col`, `z_col`: which columns of the state matrix correspond to the y and z
  coordinates used for plotting.
"""
struct ScenarioResult
    label::String
    times::Vector{Float64}
    agent_trajs::Vector{Matrix{Float64}}
    target_trajs::Vector{Vector{Vector{Float64}}}
    null_dim::Int
    residual::Float64
    y_col::Int
    z_col::Int
end

# ---------------------------------------------------------------------------
# Helper functions
# ---------------------------------------------------------------------------

"""
    selector_matrix(indices, n)

Return the `length(indices) × n` row-selection matrix `S` such that
`S * v == v[indices]`.  Throws `ArgumentError` if any index is outside `1:n`.
"""
function selector_matrix(indices::AbstractVector{<:Integer}, n::Int)
    if !(all(i -> i >= 1, indices) && all(i -> i <= n, indices))
        throw(ArgumentError("All indices must be in 1:$n, got $indices"))
    end
    S = zeros(Float64, length(indices), n)
    for (row, col) in enumerate(indices)
        S[row, col] = 1.0
    end
    return S
end

"""
    state_projection_matrix(state_indices, nx, nu)

Return a `length(state_indices) × (nx + nu)` matrix that selects the given
state coordinates from an augmented `(nx + nu)`-stalk.
"""
function state_projection_matrix(state_indices::AbstractVector{<:Integer}, nx::Int, nu::Int)
    sel = selector_matrix(state_indices, nx)
    return hcat(sel, zeros(length(state_indices), nu))
end

"""
    agent_vertex(prob, i, t)

Index of agent `i` at timestep `t` in the sheaf built from `prob`.
"""
agent_vertex(prob::TrackingProblem, i::Int, t::Int) =
    t * (prob.n_agents + prob.n_targets) + i

"""
    target_vertex(prob, j, t)

Index of target `j` at timestep `t` in the sheaf built from `prob`.
"""
target_vertex(prob::TrackingProblem, j::Int, t::Int) =
    t * (prob.n_agents + prob.n_targets) + prob.n_agents + j

# ---------------------------------------------------------------------------
# Sheaf construction
# ---------------------------------------------------------------------------

"""
    build_time_expanded_tracking_sheaf(prob::TrackingProblem)

Construct the time-expanded `EuclideanSheaf` from a `TrackingProblem`.

Four edge families are added in order:

1. **Agent dynamics** — `agent(i,t) ↔ agent(i,t+1)` for every agent `i` and
   step `t = 0, …, k-1`, encoding `x_{t+1} = Ad * x_t + Bd * u_t`.
2. **Target dynamics** (when `include_target_dynamics`) — same structure for
   each target vertex.
3. **Consensus** — `agent(i,t) ↔ agent(j,t)` for each `(i,j)` in `agent_edges`,
   active at each `t` in `consensus_timesteps`.
4. **Tracking** — for each `TrackingEdge(i, j, Ra, Rt)`, adds
   `agent(i,t) ↔ target(j,t)` active at each `t` in `tracking_timesteps`,
   using restriction maps `Ra` (agent side) and `Rt` (target side).

Restriction maps for consensus and tracking edges are pre-scaled by
`sqrt(consensus_weight)` and `sqrt(tracking_weight)` respectively so that the
resulting Laplacian is `weight * R' * R`.
"""
function build_time_expanded_tracking_sheaf(prob::TrackingProblem)
    nx_loc  = size(prob.Ad, 1)
    nu_loc  = size(prob.Bd, 2)
    nt      = nx_loc + nu_loc
    n_per_t = prob.n_agents + prob.n_targets
    sheaf   = EuclideanSheaf{Float64}(fill(nt, (prob.k + 1) * n_per_t))
    now_map  = hcat(Matrix(prob.Ad), Matrix(prob.Bd))
    next_map = hcat(Matrix{Float64}(I, nx_loc, nx_loc), zeros(nx_loc, nu_loc))
    for t in 0:prob.k-1, i in 1:prob.n_agents
        add_sheaf_edge!(sheaf,
            agent_vertex(prob, i, t), agent_vertex(prob, i, t + 1),
            now_map, next_map)
    end
    if prob.include_target_dynamics
        for t in 0:prob.k-1, j in 1:prob.n_targets
            add_sheaf_edge!(sheaf,
                target_vertex(prob, j, t), target_vertex(prob, j, t + 1),
                now_map, next_map)
        end
    end
    cscale = sqrt(prob.consensus_weight)
    tscale = sqrt(prob.tracking_weight)
    for t in prob.consensus_timesteps, (i, j) in prob.agent_edges
        add_sheaf_edge!(sheaf,
            agent_vertex(prob, i, t), agent_vertex(prob, j, t),
            cscale * prob.consensus_restriction, cscale * prob.consensus_restriction)
    end
    for t in prob.tracking_timesteps, te in prob.tracking_edges
        add_sheaf_edge!(sheaf,
            agent_vertex(prob, te.agent_index, t),
            target_vertex(prob, te.target_index, t),
            tscale * te.agent_restriction, tscale * te.target_restriction)
    end
    return sheaf
end

# ---------------------------------------------------------------------------
# Reference trajectory
# ---------------------------------------------------------------------------

"""
    generate_reference_trajectory(x0, xk, k, Ad, Bd, nx, nu)

Compute the minimum-energy trajectory from `x0` to `xk` in `k` steps by
harmonic extension on a single-agent dynamics sheaf.

Each stalk has dimension `nx + nu`.  Interior control components are in
general nonzero—they are the minimum-energy inputs satisfying the linear
dynamics.  This function is designed to produce fully-pinned target boundary
data; the control components of the returned stalks do not affect the
multi-agent solution once the targets are pinned as boundary vertices.
"""
function generate_reference_trajectory(x0, xk, k, Ad, Bd, nx, nu)
    ref_sheaf = EuclideanSheaf{Float64}(fill(nx + nu, k + 1))
    now_map  = hcat(Matrix(Ad), Matrix(Bd))
    next_map = hcat(Matrix{Float64}(I, nx, nx), zeros(nx, nu))
    for t in 1:k
        add_sheaf_edge!(ref_sheaf, t, t + 1, now_map, next_map)
    end
    bnd = Dict{Int,Vector{Float64}}(
        1     => vcat(Vector{Float64}(x0), zeros(nu)),
        k + 1 => vcat(Vector{Float64}(xk), zeros(nu)),
    )
    z_ref, _ = harmonic_extension(ref_sheaf, bnd)
    return [Array(z_ref[Block(t)]) for t in 1:k+1]
end

# ---------------------------------------------------------------------------
# Extraction and solving
# ---------------------------------------------------------------------------

"""
    extract_state_trajectories(z_harmonic, prob)

Return one `(k+1) × nx` matrix per agent, pre-allocated and filled directly
from the harmonic-extension output (no intermediate concatenations).
"""
function extract_state_trajectories(z_harmonic, prob::TrackingProblem)
    nx_loc = size(prob.Ad, 1)
    return [
        begin
            X = Matrix{Float64}(undef, prob.k + 1, nx_loc)
            for t in 0:prob.k
                X[t + 1, :] = Array(z_harmonic[Block(agent_vertex(prob, i, t))])[1:nx_loc]
            end
            X
        end
        for i in 1:prob.n_agents
    ]
end

"""
    run_scenario(label, prob, boundary, times; target_trajs, y_col, z_col)

Build the sheaf from `prob`, run `harmonic_extension`, compute the Laplacian
energy `sqrt(z' * L * z)`, and return a `ScenarioResult`.

`sheaf_laplacian_matrix_direct` is used for the energy (rather than
`coboundary_map`) because the coboundary matrix only spans vertices appearing
in at least one edge; isolated target vertices (when
`include_target_dynamics = false`) would cause a dimension mismatch.

`y_col` and `z_col` are the column indices of the y and z coordinates in
the state vector; they are stored in the result for plotting.  Defaults are 1
and 2, matching the planar-quadrotor model convention.
"""
function run_scenario(
    label::String,
    prob::TrackingProblem,
    boundary::Dict{Int,Vector{Float64}},
    times::AbstractVector{<:Real};
    target_trajs::Vector{Vector{Vector{Float64}}} = Vector{Vector{Vector{Float64}}}(),
    y_col::Int = 1,
    z_col::Int = 2,
)
    sheaf = build_time_expanded_tracking_sheaf(prob)
    z_harmonic, null_basis = harmonic_extension(sheaf, boundary)
    L  = sheaf_laplacian_matrix_direct(sheaf)
    xv = Array(z_harmonic)
    r  = sqrt(max(0.0, dot(xv, L * xv)))
    nd = size(null_basis, 2)
    println("$label:  ||dz|| = $(round(r; sigdigits=4)),   null dim = $nd")
    trajs = extract_state_trajectories(z_harmonic, prob)
    return ScenarioResult(label, collect(times), trajs, target_trajs, nd, r, y_col, z_col)
end

end # module MultiAgentTracking
