# # Multi-Agent Target Tracking via a Time-Expanded Cellular Sheaf
#
# This example uses a time-expanded sheaf to study how coordination constraints
# determine the space of feasible agent trajectories.  Four scenarios are compared;
# the nullspace dimension of the restricted sheaf Laplacian measures how many
# free degrees of freedom remain after constraints are imposed.
#
# **Vehicle model.** Each agent or target is a planar quadrotor.  The state is
# ```math
# x = [y,\, z,\, \varphi,\, \dot y,\, \dot z,\, \dot\varphi]^\top \in \mathbb{R}^6
# ```
# where ``y`` is the lateral position, ``z`` is the altitude, and ``\varphi`` is
# the roll angle.  The control input is ``u = [T_1, T_2]^\top \in \mathbb{R}^2``.
#
# **Sheaf structure.** Each vertex ``(\text{entity}, t)`` carries a stalk
# ``\mathbb{R}^{n_x + n_u} = \mathbb{R}^8``.  Three edge families encode:
# 1. *Dynamics*: ``x_{t+1} = A_d x_t + B_d u_t``.
# 2. *Consensus*: inter-agent agreement on a chosen coordinate subspace.
# 3. *Tracking*: agent–target alignment in a chosen coordinate subspace.
#
# By varying which coordinates are constrained and at which timesteps, we obtain
# four qualitatively different solution spaces whose nullspace dimensions tell
# the story of the engineering constraints.

using CellularSheaves
using CellularSheaves.TrajectorySheaves: continuous_to_discrete_zoh
using LinearAlgebra
using BlockArrays
using Plots

# ## Planar Quadrotor Dynamics
#
# State: ``x = [y, z, \varphi, \dot y, \dot z, \dot\varphi]``.
# The linearisation is taken around the hover trim condition.

g = 9.81
m_veh = 0.5
I_quad = 0.01
ell = 0.25

Ac = [0.0  0.0   0.0   1.0  0.0  0.0;
      0.0  0.0   0.0   0.0  1.0  0.0;
      0.0  0.0   0.0   0.0  0.0  1.0;
      0.0  0.0  -g     0.0  0.0  0.0;
      0.0  0.0   0.0   0.0  0.0  0.0;
      0.0  0.0   0.0   0.0  0.0  0.0]

Bc = [0.0               0.0;
      0.0               0.0;
      0.0               0.0;
      0.0               0.0;
      1.0 / m_veh       1.0 / m_veh;
      ell / (2I_quad)  -ell / (2I_quad)]

h = 0.05
Ad, Bd = continuous_to_discrete_zoh(Ac, Bc, h)
nx = size(Ad, 1)
nu = size(Bd, 2)

# State-index constants used when constructing projection matrices.
const IDX_Y   = 1  # lateral position y
const IDX_Z   = 2  # altitude z
const IDX_PHI = 3  # roll angle φ

# ## Reusable Helper Functions
#
# All four scenarios share the same building blocks.  The functions below are
# designed to be composed without modification.

"""
    selector_matrix(indices, n)

Return the ``|\\text{indices}| \\times n`` row-selection matrix ``S`` such that
``S v = v[\\text{indices}]``.

Throws `ArgumentError` if any index is outside `1:n`.
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

Return a ``|\\text{state\\_indices}| \\times (n_x + n_u)`` projection matrix that
selects the given state coordinates from an ``(n_x + n_u)``-dimensional stalk.
The control block is mapped to zero.
"""
function state_projection_matrix(state_indices::AbstractVector{<:Integer}, nx::Int, nu::Int)
    sel = selector_matrix(state_indices, nx)
    return hcat(sel, zeros(length(state_indices), nu))
end

"""
    build_time_expanded_tracking_sheaf(
        n_agents, n_targets, k, Ad, Bd;
        agent_edges, assignment,
        consensus_restriction, agent_tracking_restriction, target_tracking_restriction,
        consensus_timesteps, tracking_timesteps,
        include_target_dynamics,
        consensus_weight, tracking_weight)

Construct a time-expanded `EuclideanSheaf` encoding multi-agent dynamics,
inter-agent consensus, and agent–target tracking.

# Vertex indexing

Vertices are indexed as follows (total ``(k+1)(n_a + n_t)`` vertices):
- `agent(i, t) = t * n_per_t + i`  for `i ∈ 1:n_agents`, `t ∈ 0:k`
- `target(j, t) = t * n_per_t + n_agents + j`  for `j ∈ 1:n_targets`, `t ∈ 0:k`

where `n_per_t = n_agents + n_targets`.

# Edge families

1. **Agent dynamics** `agent(i,t) ↔ agent(i,t+1)`, `t = 0,…,k−1`:
   restriction maps ``[A_d \\mid B_d]`` and ``[I \\mid 0]`` encode
   ``x_{t+1} = A_d x_t + B_d u_t``.
2. **Target dynamics** (when `include_target_dynamics = true`):
   same dynamics structure for target vertices.
3. **Consensus** `agent(i,t) ↔ agent(j,t)` for `(i,j) ∈ agent_edges`,
   `t ∈ consensus_timesteps`: both sides use `consensus_restriction`.
4. **Tracking** `agent(i,t) ↔ target(assign(i),t)` for `t ∈ tracking_timesteps`:
   agent side uses `agent_tracking_restriction`,
   target side uses `target_tracking_restriction`.

# Arguments
- `n_agents`, `n_targets`: fleet sizes.
- `k`: horizon length; trajectory has `k+1` time steps ``t = 0, \\ldots, k``.
- `Ad`, `Bd`: ZOH-discretised state-space matrices.
- `agent_edges`: list of `(i,j)` pairs for undirected agent–agent consensus edges.
- `assignment`: `assignment[i]` is the target index tracked by agent `i`
  (must satisfy `1 ≤ assignment[i] ≤ n_targets`).
- `consensus_restriction`: restriction map for agent–agent edges
  (maps ``(n_x+n_u)``-stalk to shared-coordinate space).
- `agent_tracking_restriction`, `target_tracking_restriction`: restriction maps
  for the agent and target sides of agent–target edges.
- `consensus_timesteps`: timesteps at which consensus edges are active (default `0:k`).
- `tracking_timesteps`: timesteps at which tracking edges are active (default `0:k`).
- `include_target_dynamics`: add target temporal dynamics edges (default `true`).
- `consensus_weight`, `tracking_weight`: Laplacian weights for the two edge families.
"""
function build_time_expanded_tracking_sheaf(
    n_agents::Int,
    n_targets::Int,
    k::Int,
    Ad::AbstractMatrix,
    Bd::AbstractMatrix;
    agent_edges::Vector{Tuple{Int,Int}},
    assignment::Vector{Int},
    consensus_restriction::AbstractMatrix{Float64},
    agent_tracking_restriction::AbstractMatrix{Float64},
    target_tracking_restriction::AbstractMatrix{Float64},
    consensus_timesteps::AbstractVector{<:Integer} = 0:k,
    tracking_timesteps::AbstractVector{<:Integer} = 0:k,
    include_target_dynamics::Bool = true,
    consensus_weight::Float64 = 1.0,
    tracking_weight::Float64 = 5.0,
)
    if length(assignment) != n_agents ||
            !(all(j -> j >= 1, assignment) && all(j -> j <= n_targets, assignment))
        throw(ArgumentError(
            "assignment must be a length-$n_agents vector with all entries in " *
            "1:$n_targets, got $assignment"))
    end
    nx_loc = size(Ad, 1)
    nu_loc = size(Bd, 2)
    nt = nx_loc + nu_loc
    n_per_t = n_agents + n_targets
    vertex_dims = fill(nt, (k + 1) * n_per_t)
    sheaf = EuclideanSheaf{Float64}(vertex_dims)
    vid_agent(i, t) = t * n_per_t + i
    vid_target(j, t) = t * n_per_t + n_agents + j
    now_map  = hcat(Matrix(Ad), Matrix(Bd))
    next_map = hcat(Matrix{Float64}(I, nx_loc, nx_loc), zeros(nx_loc, nu_loc))
    for t in 0:k-1, i in 1:n_agents
        add_sheaf_edge!(sheaf, vid_agent(i, t), vid_agent(i, t + 1), now_map, next_map)
    end
    if include_target_dynamics
        for t in 0:k-1, j in 1:n_targets
            add_sheaf_edge!(sheaf, vid_target(j, t), vid_target(j, t + 1), now_map, next_map)
        end
    end
    cscale = sqrt(consensus_weight)
    tscale = sqrt(tracking_weight)
    for t in consensus_timesteps, (i, j) in agent_edges
        add_sheaf_edge!(sheaf, vid_agent(i, t), vid_agent(j, t),
            cscale * consensus_restriction, cscale * consensus_restriction)
    end
    for t in tracking_timesteps, i in 1:n_agents
        j = assignment[i]
        add_sheaf_edge!(sheaf, vid_agent(i, t), vid_target(j, t),
            tscale * agent_tracking_restriction, tscale * target_tracking_restriction)
    end
    index = (
        agent   = vid_agent,
        target  = vid_target,
        n_per_t = n_per_t,
        nx      = nx_loc,
        nu      = nu_loc,
        nt      = nt,
    )
    return sheaf, index
end

"""
    bobbing_target_trajectory(y_fixed, z_center, z_amplitude, k, h, nx, nu; n_periods=2)

Return a `Vector` of `k+1` stalk vectors for a sinusoidally oscillating ("bobbing") target.

Altitude follows ``z(t) = z_c + A\\sin(2\\pi n_p\\, t/(kh))``,
with ``\\dot z`` set to the analytic derivative.  Lateral position is `y_fixed`;
roll angle and all other states are zero.
"""
function bobbing_target_trajectory(
    y_fixed::Float64,
    z_center::Float64,
    z_amplitude::Float64,
    k::Int,
    h::Float64,
    nx::Int,
    nu::Int;
    n_periods::Int = 2,
)
    T_total = k * h
    omega = 2π * n_periods / T_total
    return map(0:k) do t
        t_s    = t * h
        z_t    = z_center + z_amplitude * sin(omega * t_s)
        zdot_t = z_amplitude * omega * cos(omega * t_s)
        x = zeros(nx)
        x[1] = y_fixed
        x[2] = z_t
        x[5] = zdot_t
        vcat(x, zeros(nu))
    end
end

"""
    generate_reference_trajectory(x0, xk, k, Ad, Bd, nx, nu)

Compute the minimum-energy trajectory from `x0` to `xk` in `k` steps by
harmonic extension on a single-agent sheaf.  Returns a `Vector` of `k+1`
stalk vectors (state concatenated with zero controls).
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

"""
    extract_state_trajectories(z_harmonic, idx, k, n_agents, nx)

Return a `Vector` of `n_agents` matrices of size `(k+1) × nx`, each containing
the state trajectory of one agent extracted from the block vector `z_harmonic`.
"""
function extract_state_trajectories(z_harmonic, idx, k::Int, n_agents::Int, nx::Int)
    return [
        reduce(vcat,
            [transpose(Array(z_harmonic[Block(idx.agent(i, t))])[1:nx]) for t in 0:k])
        for i in 1:n_agents
    ]
end

"""
    run_scenario(label, sheaf, boundary, idx, k, n_agents, nx)

Run `harmonic_extension` on `sheaf` with the given `boundary` conditions.
Prints the Laplacian energy ``\\sqrt{z^\\top L z}`` (equal to ``\\|dz\\|`` for all
edges in the sheaf) and the nullspace dimension, and returns
`(state_trajectories, nullspace_dim, residual)`.

`sheaf_laplacian_matrix_direct` is used instead of `coboundary_map` because
`coboundary_map` produces a matrix whose **column count** equals the sum of stalk
dimensions only for vertices that appear in at least one edge.  When targets are
isolated (i.e. `include_target_dynamics = false`), those vertices are absent from
the coboundary matrix, making `d * Array(z_harmonic)` a size mismatch.
`sheaf_laplacian_matrix_direct` always spans all ``n_{\\text{total}}`` stalk dimensions.
"""
function run_scenario(label::String, sheaf, boundary, idx, k::Int, n_agents::Int, nx::Int)
    z_harmonic, null_basis = harmonic_extension(sheaf, boundary)
    L  = sheaf_laplacian_matrix_direct(sheaf)
    xv = Array(z_harmonic)
    r  = sqrt(max(0.0, dot(xv, L * xv)))
    nd = size(null_basis, 2)
    println("$label:  ||dz|| = $(round(r; sigdigits=4)),   null dim = $nd")
    trajs = extract_state_trajectories(z_harmonic, idx, k, n_agents, nx)
    return trajs, nd, r
end

# ## Common Setup
#
# All four scenarios use two agents and two targets over a 40-step horizon (2 s).
# The projection matrices are defined once and reused.
#
# - ``R_{yz}``: projects the ``(n_x + n_u)``-stalk onto ``(y, z)`` (used in Scenario 1).
# - ``R_y``: projects onto ``y`` only (consensus in Scenarios 2–4).
# - ``R_z``: projects onto ``z`` only (tracking in Scenarios 2–4).

n_agents  = 2
n_targets = 2
k         = 40

times = h .* (0:k)

R_yz = state_projection_matrix([IDX_Y, IDX_Z], nx, nu)  # 2×8: projects onto (y, z)
R_y  = state_projection_matrix([IDX_Y],         nx, nu)  # 1×8: projects onto y
R_z  = state_projection_matrix([IDX_Z],         nx, nu)  # 1×8: projects onto z

# Bobbing reference trajectories (used in Scenarios 2–4).
# Target 1 oscillates around z = 1.0 m; Target 2 around z = 2.0 m.
# Both complete two full vertical cycles over the 2-second horizon.
traj_bt1 = bobbing_target_trajectory(0.0, 1.0, 0.3, k, h, nx, nu; n_periods=2)
traj_bt2 = bobbing_target_trajectory(0.0, 2.0, 0.3, k, h, nx, nu; n_periods=2)

# ## Scenario 1: Full (y,z) Coordination at Every Timestep
#
# Both the agent–agent consensus edges and the agent–target tracking edges use
# the ``(y, z)`` projection and are active at **every** timestep.  Two agents
# travel to distinct goals so the problem is non-trivial; their reference
# trajectories are pre-computed by harmonic extension so they satisfy the dynamics.
#
# This is the most constrained scenario.

sheaf1, idx1 = build_time_expanded_tracking_sheaf(
    n_agents, n_targets, k, Ad, Bd;
    agent_edges                 = [(1, 2)],
    assignment                  = [1, 2],
    consensus_restriction       = R_yz,
    agent_tracking_restriction  = R_yz,
    target_tracking_restriction = R_yz,
    consensus_timesteps         = collect(0:k),
    tracking_timesteps          = collect(0:k),
    include_target_dynamics     = true,
    consensus_weight            = 1.0,
    tracking_weight             = 5.0,
)

# Agent 1 travels from ``(y, z) = (0, 1)`` to ``(2, 1.5)``;
# Agent 2 from ``(y, z) = (0, 2)`` to ``(2, 2.5)``.
# Target endpoints match their assigned agent endpoints so tracking is satisfiable.

x0_a1_s1 = [0.0, 1.0, 0.0, 0.0, 0.0, 0.0]
xk_a1_s1 = [2.0, 1.5, 0.0, 0.0, 0.0, 0.0]
x0_a2_s1 = [0.0, 2.0, 0.0, 0.0, 0.0, 0.0]
xk_a2_s1 = [2.0, 2.5, 0.0, 0.0, 0.0, 0.0]

traj_t1_s1 = generate_reference_trajectory(x0_a1_s1, xk_a1_s1, k, Ad, Bd, nx, nu)
traj_t2_s1 = generate_reference_trajectory(x0_a2_s1, xk_a2_s1, k, Ad, Bd, nx, nu)

# The entire target trajectory is pinned as boundary data.

bnd1 = Dict{Int,Vector{Float64}}()
bnd1[idx1.agent(1, 0)] = vcat(x0_a1_s1, zeros(nu))
bnd1[idx1.agent(2, 0)] = vcat(x0_a2_s1, zeros(nu))
bnd1[idx1.agent(1, k)] = vcat(xk_a1_s1, zeros(nu))
bnd1[idx1.agent(2, k)] = vcat(xk_a2_s1, zeros(nu))
for t in 0:k
    bnd1[idx1.target(1, t)] = traj_t1_s1[t + 1]
    bnd1[idx1.target(2, t)] = traj_t2_s1[t + 1]
end

trajs1, ndim1, res1 = run_scenario("Scenario 1", sheaf1, bnd1, idx1, k, n_agents, nx)

p1_y = plot(times, trajs1[1][:, 1]; lw=2, label="A1",
    title="Scenario 1 — y(t)", xlabel="t (s)", ylabel="y (m)", legend=:topleft)
plot!(p1_y, times, trajs1[2][:, 1]; lw=2, ls=:dash, label="A2")
plot!(p1_y, times, getindex.(traj_t1_s1, 1); lw=1, ls=:dot, color=:gray, label="T1")
plot!(p1_y, times, getindex.(traj_t2_s1, 1); lw=1, ls=:dot, color=:black, label="T2")

p1_z = plot(times, trajs1[1][:, 2]; lw=2, label="A1",
    title="Scenario 1 — z(t)", xlabel="t (s)", ylabel="z (m)", legend=:topleft)
plot!(p1_z, times, trajs1[2][:, 2]; lw=2, ls=:dash, label="A2")
plot!(p1_z, times, getindex.(traj_t1_s1, 2); lw=1, ls=:dot, color=:gray, label="T1")
plot!(p1_z, times, getindex.(traj_t2_s1, 2); lw=1, ls=:dot, color=:black, label="T2")

plot(p1_y, p1_z; layout=(1, 2), size=(800, 350),
    plot_title="Scenario 1: null dim = $ndim1,  ||dz|| = $(round(res1; sigdigits=3))")

# ## Scenario 2: y-Consensus and z-Tracking, All Timesteps, Aligned Initial Conditions
#
# Targets oscillate vertically (the "bobbing" motion). Consensus edges constrain
# agents to agree only on *lateral position* ``y``; tracking edges constrain each
# agent to match its assigned target only in *altitude* ``z``.  Both edge families
# remain active at **every** timestep.
#
# The agents start **aligned**: they share the same ``y``-coordinate and each starts
# at the altitude of its assigned target.  Because only one coordinate per edge type
# is constrained (versus two in Scenario 1), a larger nullspace is expected.
#
# Since all target vertices are pinned as boundary data, target dynamics edges are
# omitted; their inclusion would introduce residuals without affecting the agent
# solution.

sheaf2, idx2 = build_time_expanded_tracking_sheaf(
    n_agents, n_targets, k, Ad, Bd;
    agent_edges                 = [(1, 2)],
    assignment                  = [1, 2],
    consensus_restriction       = R_y,
    agent_tracking_restriction  = R_z,
    target_tracking_restriction = R_z,
    consensus_timesteps         = collect(0:k),
    tracking_timesteps          = collect(0:k),
    include_target_dynamics     = false,
    consensus_weight            = 1.0,
    tracking_weight             = 5.0,
)

# Aligned initial conditions: both agents at y = 0, each at its target's initial altitude.

x0_a1_s2 = [0.0, traj_bt1[1][2], 0.0, 0.0, traj_bt1[1][5], 0.0]
x0_a2_s2 = [0.0, traj_bt2[1][2], 0.0, 0.0, traj_bt2[1][5], 0.0]

bnd2 = Dict{Int,Vector{Float64}}()
bnd2[idx2.agent(1, 0)] = vcat(x0_a1_s2, zeros(nu))
bnd2[idx2.agent(2, 0)] = vcat(x0_a2_s2, zeros(nu))
for t in 0:k
    bnd2[idx2.target(1, t)] = traj_bt1[t + 1]
    bnd2[idx2.target(2, t)] = traj_bt2[t + 1]
end

trajs2, ndim2, res2 = run_scenario("Scenario 2", sheaf2, bnd2, idx2, k, n_agents, nx)

p2_y = plot(times, trajs2[1][:, 1]; lw=2, label="A1",
    title="Scenario 2 — y(t)", xlabel="t (s)", ylabel="y (m)", legend=:topright)
plot!(p2_y, times, trajs2[2][:, 1]; lw=2, ls=:dash, label="A2")
plot!(p2_y, times, getindex.(traj_bt1, 1); lw=1, ls=:dot, color=:gray, label="T1 y")
plot!(p2_y, times, getindex.(traj_bt2, 1); lw=1, ls=:dot, color=:black, label="T2 y")

p2_z = plot(times, trajs2[1][:, 2]; lw=2, label="A1",
    title="Scenario 2 — z(t)", xlabel="t (s)", ylabel="z (m)", legend=:topright)
plot!(p2_z, times, trajs2[2][:, 2]; lw=2, ls=:dash, label="A2")
plot!(p2_z, times, getindex.(traj_bt1, 2); lw=1, ls=:dot, color=:gray, label="T1 z")
plot!(p2_z, times, getindex.(traj_bt2, 2); lw=1, ls=:dot, color=:black, label="T2 z")

plot(p2_y, p2_z; layout=(1, 2), size=(800, 350),
    plot_title="Scenario 2: null dim = $ndim2,  ||dz|| = $(round(res2; sigdigits=3))")

# ## Scenario 3: y-Consensus and z-Tracking at the Last Timestep Only
#
# The same projection matrices are used as in Scenario 2, but consensus and
# tracking edges are now present **only at** ``t = k``.  Agents evolve freely
# (subject only to dynamics) for ``t = 0, \ldots, k-1`` and must satisfy the
# coordination conditions at the terminal time.
#
# Both targets are fixed at ``y = 0`` (aligned in ``y``).  Agents start
# misaligned with each other and with the targets in both ``y`` and ``z``.
#
# Removing ``k`` timesteps of coordination edges drastically increases the
# nullspace: there are many dynamics-consistent trajectories that begin at the
# specified initial conditions and converge to the terminal constraint.

sheaf3, idx3 = build_time_expanded_tracking_sheaf(
    n_agents, n_targets, k, Ad, Bd;
    agent_edges                 = [(1, 2)],
    assignment                  = [1, 2],
    consensus_restriction       = R_y,
    agent_tracking_restriction  = R_z,
    target_tracking_restriction = R_z,
    consensus_timesteps         = [k],
    tracking_timesteps          = [k],
    include_target_dynamics     = false,
    consensus_weight            = 1.0,
    tracking_weight             = 5.0,
)

# Agents start misaligned with each other and with the targets (both y and z).
# Targets share y = 0 (aligned in y) but bob at different altitudes.

x0_a1_s3 = [-0.5, 0.5, 0.0, 0.0, 0.0, 0.0]
x0_a2_s3 = [ 1.0, 1.5, 0.0, 0.0, 0.0, 0.0]

bnd3 = Dict{Int,Vector{Float64}}()
bnd3[idx3.agent(1, 0)] = vcat(x0_a1_s3, zeros(nu))
bnd3[idx3.agent(2, 0)] = vcat(x0_a2_s3, zeros(nu))
for t in 0:k
    bnd3[idx3.target(1, t)] = traj_bt1[t + 1]
    bnd3[idx3.target(2, t)] = traj_bt2[t + 1]
end

trajs3, ndim3, res3 = run_scenario("Scenario 3", sheaf3, bnd3, idx3, k, n_agents, nx)

p3_y = plot(times, trajs3[1][:, 1]; lw=2, label="A1",
    title="Scenario 3 — y(t)", xlabel="t (s)", ylabel="y (m)", legend=:topright)
plot!(p3_y, times, trajs3[2][:, 1]; lw=2, ls=:dash, label="A2")
plot!(p3_y, times, getindex.(traj_bt1, 1); lw=1, ls=:dot, color=:gray, label="T1 y")

p3_z = plot(times, trajs3[1][:, 2]; lw=2, label="A1",
    title="Scenario 3 — z(t)", xlabel="t (s)", ylabel="z (m)", legend=:topright)
plot!(p3_z, times, trajs3[2][:, 2]; lw=2, ls=:dash, label="A2")
plot!(p3_z, times, getindex.(traj_bt1, 2); lw=1, ls=:dot, color=:gray, label="T1 z")
plot!(p3_z, times, getindex.(traj_bt2, 2); lw=1, ls=:dot, color=:black, label="T2 z")

plot(p3_y, p3_z; layout=(1, 2), size=(800, 350),
    plot_title="Scenario 3: null dim = $ndim3,  ||dz|| = $(round(res3; sigdigits=3))")

# ## Scenario 4: Last-Timestep Constraints, Targets Not Aligned in y or z
#
# Identical sheaf structure to Scenario 3, but Target 2 is moved to a different
# lateral position (``y = 1.5`` instead of ``y = 0``).  The targets are now
# **not** initially aligned with each other in either ``y`` or ``z``.
#
# Because the sheaf topology is the same as Scenario 3, the **nullspace
# dimension is identical**.  The difference lies in the harmonic extension
# itself: the two agents must still converge to the same ``y`` (consensus)
# and to their respective target's ``z`` (tracking) at ``t = k``, but the
# terminal target configuration is now laterally offset.  This illustrates
# that nullspace dimension is a structural property of the coordination
# architecture, not of the specific target configuration.

traj_bt2_s4 = bobbing_target_trajectory(1.5, 2.0, 0.3, k, h, nx, nu; n_periods=2)

sheaf4, idx4 = build_time_expanded_tracking_sheaf(
    n_agents, n_targets, k, Ad, Bd;
    agent_edges                 = [(1, 2)],
    assignment                  = [1, 2],
    consensus_restriction       = R_y,
    agent_tracking_restriction  = R_z,
    target_tracking_restriction = R_z,
    consensus_timesteps         = [k],
    tracking_timesteps          = [k],
    include_target_dynamics     = false,
    consensus_weight            = 1.0,
    tracking_weight             = 5.0,
)

x0_a1_s4 = [-0.5, 0.5, 0.0, 0.0, 0.0, 0.0]
x0_a2_s4 = [ 1.0, 1.5, 0.0, 0.0, 0.0, 0.0]

bnd4 = Dict{Int,Vector{Float64}}()
bnd4[idx4.agent(1, 0)] = vcat(x0_a1_s4, zeros(nu))
bnd4[idx4.agent(2, 0)] = vcat(x0_a2_s4, zeros(nu))
for t in 0:k
    bnd4[idx4.target(1, t)] = traj_bt1[t + 1]
    bnd4[idx4.target(2, t)] = traj_bt2_s4[t + 1]
end

trajs4, ndim4, res4 = run_scenario("Scenario 4", sheaf4, bnd4, idx4, k, n_agents, nx)

p4_y = plot(times, trajs4[1][:, 1]; lw=2, label="A1",
    title="Scenario 4 — y(t)", xlabel="t (s)", ylabel="y (m)", legend=:topright)
plot!(p4_y, times, trajs4[2][:, 1]; lw=2, ls=:dash, label="A2")
plot!(p4_y, times, getindex.(traj_bt1,    1); lw=1, ls=:dot, color=:gray,  label="T1 y")
plot!(p4_y, times, getindex.(traj_bt2_s4, 1); lw=1, ls=:dot, color=:black, label="T2 y")

p4_z = plot(times, trajs4[1][:, 2]; lw=2, label="A1",
    title="Scenario 4 — z(t)", xlabel="t (s)", ylabel="z (m)", legend=:topright)
plot!(p4_z, times, trajs4[2][:, 2]; lw=2, ls=:dash, label="A2")
plot!(p4_z, times, getindex.(traj_bt1,    2); lw=1, ls=:dot, color=:gray,  label="T1 z")
plot!(p4_z, times, getindex.(traj_bt2_s4, 2); lw=1, ls=:dot, color=:black, label="T2 z")

plot(p4_y, p4_z; layout=(1, 2), size=(800, 350),
    plot_title="Scenario 4: null dim = $ndim4,  ||dz|| = $(round(res4; sigdigits=3))")

# ## Comparison: How Constraints Shape the Solution Space
#
# The four scenarios are summarised below.  The nullspace dimension and
# Laplacian residual ``\|dz\|`` are printed at runtime.
#
# | Scenario | Consensus coords | Tracking coords | Active timesteps | Initial alignment |
# |----------|-----------------|-----------------|-----------------|-------------------|
# | 1        | ``(y, z)``      | ``(y, z)``      | all ``0:k``     | distinct endpoints|
# | 2        | ``y``           | ``z``           | all ``0:k``     | aligned with targets|
# | 3        | ``y``           | ``z``           | terminal ``k``  | unaligned; targets share ``y``|
# | 4        | ``y``           | ``z``           | terminal ``k``  | unaligned; targets offset in ``y``|
#
# Three observations stand out:
#
# **Constraint density governs nullspace size.**
# Scenario 1 (``(y,z)`` constrained at every step, two distinct agent goals)
# produces a unique solution (null dim = 0): the competing constraints over-specify
# the system and the harmonic extension finds a single energy-minimising compromise,
# reflected in the nonzero ``\|dz\|``.
# Scenario 2 (only ``y`` and ``z`` constrained separately, same agents aligned at
# start) allows a family of solutions (null dim > 0) because fewer coordinates are
# constrained per edge.
# Restricting coordination to only the terminal timestep (Scenarios 3–4) removes
# ``k`` timesteps of constraints and dramatically enlarges the feasible space.
#
# **Nullspace dimension is a topological invariant.**
# Scenarios 3 and 4 differ only in the *values* of the target trajectories
# (targets sharing ``y = 0`` vs.\ laterally offset).  Because the sheaf topology —
# which vertices exist and which edges connect them — is identical, the restricted
# Laplacian has the same nullspace dimension in both cases.
# The two harmonic extensions produce different trajectories, but they live in
# equally large affine solution spaces.
#
# **Feasibility is separate from trajectory shape.**
# In Scenario 2, the agents start already satisfying the all-time constraints, so
# the harmonic solution closely follows the bobbing targets.  In Scenarios 3–4,
# agents start misaligned; the harmonic extension finds the minimum-energy path that
# converges to the terminal constraint, while the large nullspace admits many
# equally valid deformations in the unconstrained directions.

println("Summary: null dims = [$ndim1, $ndim2, $ndim3, $ndim4]")

bar(["S1\n(y,z) all-t", "S2\ny/z all-t", "S3\ny/z last-t", "S4\ny/z last-t (mis)"],
    [ndim1, ndim2, ndim3, ndim4];
    ylabel="Nullspace dimension",
    title="Nullspace dimension vs. coordination architecture",
    legend=false, color=[:steelblue, :darkorange, :green, :crimson])
