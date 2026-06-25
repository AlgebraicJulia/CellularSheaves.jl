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

Specification of a heterogeneous multi-agent / multi-target tracking problem
encoded as a time-expanded cellular sheaf.

Fields:
- `n_agents`, `n_targets`: fleet sizes.
- `k`: horizon length; the trajectory has `k+1` timesteps `t = 0, …, k`.
- `Ad`: ZOH-discretised dynamics matrices, one per agent (length `n_agents`).
- `Bd`: ZOH-discretised control matrices, one per agent (length `n_agents`).
- `target_Ad`: dynamics matrices, one per target (length `n_targets`).
- `target_Bd`: control matrices, one per target (length `n_targets`).
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
    Ad::Vector{Matrix{Float64}}          # length = n_agents
    Bd::Vector{Matrix{Float64}}          # length = n_agents
    target_Ad::Vector{Matrix{Float64}}   # length = n_targets
    target_Bd::Vector{Matrix{Float64}}   # length = n_targets
    agent_edges::Vector{Tuple{Int,Int}}
    tracking_edges::Vector{TrackingEdge}
    consensus_restriction::Matrix{Float64}
    consensus_timesteps::Vector{Int}
    tracking_timesteps::Vector{Int}
    include_target_dynamics::Bool
    consensus_weight::Float64
    tracking_weight::Float64
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

# Restriction maps for the agent dynamics edge agent(i,t) ↔ agent(i,t+1).
# Default:            source [A B], destination [I 0]  ⟹  x_{t+1} = A x_t + B u_t.
# state_only_initial: controls reindexed so u_{t+1} produces x_{t+1} — source
#   [A 0] (or just [A] at the state-only t=0 vertex), destination [I −B].
function _agent_dynamics_maps(Ad, Bd, t::Int, state_only_initial::Bool)
    nx, nu = size(Ad, 1), size(Bd, 2)
    if state_only_initial
        src_map = t == 0 ? Matrix(Ad) : hcat(Matrix(Ad), zeros(nx, nu))
        return src_map, hcat(Matrix{Float64}(I, nx, nx), -Matrix(Bd))
    else
        return hcat(Matrix(Ad), Matrix(Bd)), hcat(Matrix{Float64}(I, nx, nx), zeros(nx, nu))
    end
end

"""
    build_time_expanded_tracking_sheaf(prob::TrackingProblem; state_only_initial=false)

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

When `state_only_initial = true` the initial agent vertices (`t = 0`) carry
**state only** (no control), and the control is reindexed so that `u_t` produces
`x_t`: the agent dynamics edge reads `A·x_t` from the source and `x_{t+1} − B·u_{t+1}`
from the destination, encoding `x_{t+1} = A·x_t + B·u_{t+1}`. Consensus/tracking
restriction maps at `t = 0` are then truncated to the agent's state columns
(exact, since those columns multiplied the now-absent control block by zero).
This makes the initial condition a clean full-vertex pin and is the form used by
the MPC window solver; the resulting optimal trajectory is identical to the
default form, which keeps `(x_0, u_0)` on the first vertex.
"""
function build_time_expanded_tracking_sheaf(prob::TrackingProblem; state_only_initial::Bool=false)
    @argcheck prob.consensus_weight >= 0 "consensus_weight must be non-negative, got $(prob.consensus_weight)"
    @argcheck prob.tracking_weight  >= 0 "tracking_weight must be non-negative, got $(prob.tracking_weight)"
    @argcheck all(t -> 0 <= t <= prob.k, prob.consensus_timesteps) "all consensus_timesteps must be in 0:$(prob.k), got $(prob.consensus_timesteps)"
    @argcheck all(t -> 0 <= t <= prob.k, prob.tracking_timesteps)  "all tracking_timesteps must be in 0:$(prob.k), got $(prob.tracking_timesteps)"
    @argcheck all(((i,j),) -> 1 <= i <= prob.n_agents && 1 <= j <= prob.n_agents, prob.agent_edges) "agent_edges indices must be in 1:$(prob.n_agents), got $(prob.agent_edges)"
    @argcheck all(te -> 1 <= te.agent_index  <= prob.n_agents,  prob.tracking_edges) "TrackingEdge agent_index must be in 1:$(prob.n_agents)"
    @argcheck all(te -> 1 <= te.target_index <= prob.n_targets, prob.tracking_edges) "TrackingEdge target_index must be in 1:$(prob.n_targets)"
    @argcheck length(prob.Ad)        == prob.n_agents  "Ad vector length must equal n_agents ($(prob.n_agents)), got $(length(prob.Ad))"
    @argcheck length(prob.Bd)        == prob.n_agents  "Bd vector length must equal n_agents ($(prob.n_agents)), got $(length(prob.Bd))"
    @argcheck length(prob.target_Ad) == prob.n_targets "target_Ad vector length must equal n_targets ($(prob.n_targets)), got $(length(prob.target_Ad))"
    @argcheck length(prob.target_Bd) == prob.n_targets "target_Bd vector length must equal n_targets ($(prob.n_targets)), got $(length(prob.target_Bd))"
    for i in 1:prob.n_agents
        @argcheck size(prob.Ad[i], 1) == size(prob.Ad[i], 2) "Ad[$i] must be square, got size $(size(prob.Ad[i]))"
        @argcheck size(prob.Bd[i], 1) == size(prob.Ad[i], 1) "Bd[$i] rows must match Ad[$i] rows, got $(size(prob.Bd[i],1)) vs $(size(prob.Ad[i],1))"
    end
    for j in 1:prob.n_targets
        @argcheck size(prob.target_Ad[j], 1) == size(prob.target_Ad[j], 2) "target_Ad[$j] must be square, got size $(size(prob.target_Ad[j]))"
        @argcheck size(prob.target_Bd[j], 1) == size(prob.target_Ad[j], 1) "target_Bd[$j] rows must match target_Ad[$j] rows, got $(size(prob.target_Bd[j],1)) vs $(size(prob.target_Ad[j],1))"
    end

    # Allocate heterogeneous stalk sizes for every vertex (agents + targets).
    # With state_only_initial the t=0 agent vertices hold state only (no control).
    n_per_t = prob.n_agents + prob.n_targets
    stalk_sizes = Vector{Int}(undef, (prob.k + 1) * n_per_t)
    for t in 0:prob.k, i in 1:prob.n_agents
        nx_i = size(prob.Ad[i], 1); nu_i = size(prob.Bd[i], 2)
        stalk_sizes[t * n_per_t + i] = (state_only_initial && t == 0) ? nx_i : nx_i + nu_i
    end
    for t in 0:prob.k, j in 1:prob.n_targets
        stalk_sizes[t * n_per_t + prob.n_agents + j] = size(prob.target_Ad[j], 1) + size(prob.target_Bd[j], 2)
    end
    sheaf = EuclideanSheaf{Float64}(stalk_sizes)

    # Agent dynamics edges (per-agent)
    for i in 1:prob.n_agents
        for t in 0:prob.k-1
            src_map, dst_map = _agent_dynamics_maps(prob.Ad[i], prob.Bd[i], t, state_only_initial)
            add_sheaf_edge!(sheaf,
                agent_vertex(prob, i, t), agent_vertex(prob, i, t + 1),
                src_map, dst_map)
        end
    end

    # Target dynamics edges (per-target, if requested)
    if prob.include_target_dynamics
        for j in 1:prob.n_targets
            At_j = prob.target_Ad[j]; Bt_j = prob.target_Bd[j]
            nx_j = size(At_j, 1); nu_j = size(Bt_j, 2)
            now_map  = hcat(Matrix(At_j), Matrix(Bt_j))
            next_map = hcat(Matrix{Float64}(I, nx_j, nx_j), zeros(nx_j, nu_j))
            for t in 0:prob.k-1
                add_sheaf_edge!(sheaf,
                    target_vertex(prob, j, t), target_vertex(prob, j, t + 1),
                    now_map, next_map)
            end
        end
    end

    # In window mode the t=0 agent stalk is state-only, so its restriction maps
    # drop the (now absent) control columns; otherwise they are used as given.
    agent_restr(R, i, t) = (state_only_initial && t == 0) ? R[:, 1:size(prob.Ad[i], 1)] : R
    cscale = sqrt(prob.consensus_weight)
    tscale = sqrt(prob.tracking_weight)
    for t in prob.consensus_timesteps, (i, j) in prob.agent_edges
        add_sheaf_edge!(sheaf,
            agent_vertex(prob, i, t), agent_vertex(prob, j, t),
            cscale * agent_restr(prob.consensus_restriction, i, t),
            cscale * agent_restr(prob.consensus_restriction, j, t))
    end
    for t in prob.tracking_timesteps, te in prob.tracking_edges
        add_sheaf_edge!(sheaf,
            agent_vertex(prob, te.agent_index, t),
            target_vertex(prob, te.target_index, t),
            tscale * agent_restr(te.agent_restriction, te.agent_index, t),
            tscale * te.target_restriction)
    end
    return sheaf
end


# ---------------------------------------------------------------------------
# Extraction and solving
# ---------------------------------------------------------------------------

"""
    extract_state_trajectories(z_harmonic, prob)

Return one `(k+1) × nx_i` matrix per agent, where `nx_i` is the state dimension
for agent `i`.  Pre-allocated and filled directly from the harmonic-extension
output (no intermediate concatenations).
"""
function extract_state_trajectories(z_harmonic, prob::TrackingProblem)
    return map(1:prob.n_agents) do i
        nx_i = size(prob.Ad[i], 1)
        X = Matrix{Float64}(undef, prob.k + 1, nx_i)
        for t in 0:prob.k
            X[t + 1, :] = Array(z_harmonic[Block(agent_vertex(prob, i, t))])[1:nx_i]
        end
        X
    end
end

"""
    extract_target_trajectories(z_harmonic, prob)

Return a `Vector` of length `n_targets`, where each element is a
`Vector{Vector{Float64}}` of length `k+1` containing the full stalk
(state + control) for that target at each timestep.
"""
function extract_target_trajectories(z_harmonic, prob::TrackingProblem)
    return map(1:prob.n_targets) do j
        traj = Vector{Vector{Float64}}(undef, prob.k + 1)
        for t in 0:prob.k
            traj[t + 1] = Array(z_harmonic[Block(target_vertex(prob, j, t))])
        end
        traj
    end
end