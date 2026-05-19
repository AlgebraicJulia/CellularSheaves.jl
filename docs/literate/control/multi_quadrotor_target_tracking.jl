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

# State-index constants.
const IDX_Y   = 1   # lateral position y
const IDX_Z   = 2   # altitude z
const IDX_PHI = 3   # roll angle φ
const IDX_YDT = 4   # ẏ
const IDX_ZDT = 5   # ż
const IDX_PHDT = 6  # φ̇

# ## Data Types
#
# All parameters of the multi-agent tracking problem live in two structs.
# `TrackingEdge` represents one directed tracking relationship between a specific
# agent and a specific target, carrying its own pair of restriction maps.
# This generalises the old single `assignment` vector to a **many-to-many**
# relationship: any agent can track any subset of targets (or none), and
# each pairing can project onto different coordinates.

"""
    TrackingEdge

One agent–target tracking relationship with its own restriction maps.

Fields:
- `agent_index`: agent index in `1:n_agents`.
- `target_index`: target index in `1:n_targets`.
- `agent_restriction`: projection applied to the agent stalk.
- `target_restriction`: projection applied to the target stalk.
"""
struct TrackingEdge
    agent_index::Int
    target_index::Int
    agent_restriction::Matrix{Float64}
    target_restriction::Matrix{Float64}
end

"""
    TrackingProblem

Encapsulates all parameters of a multi-agent, multi-target tracking problem
as a time-expanded cellular sheaf.

Fields:
- `n_agents`, `n_targets`: fleet sizes.
- `k`: horizon length (trajectory has `k+1` timesteps, `t = 0, …, k`).
- `Ad`, `Bd`: ZOH-discretised state-space matrices.
- `agent_edges`: undirected agent–agent consensus edge pairs.
- `tracking_edges`: many-to-many agent–target edges (see `TrackingEdge`).
- `consensus_restriction`: projection used for all agent–agent consensus edges.
- `consensus_timesteps`, `tracking_timesteps`: timesteps at which each edge family is active.
- `include_target_dynamics`: whether to add temporal dynamics edges for targets.
- `consensus_weight`, `tracking_weight`: Laplacian weights for the two coordination families.
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

# ## Target Entity Type
#
# Rather than a standalone function, a `BobbingTarget` struct stores all
# physical parameters of a sinusoidally oscillating target.  Call
# `trajectory(bt, t_range, h, nx, nu)` to materialise the stalk sequence
# over any integer time range.

"""
    BobbingTarget

Parameters defining a vertically oscillating ("bobbing") target.

Fields:
- `y_fixed`: fixed lateral position (m).
- `z_center`: mean altitude (m).
- `z_amplitude`: altitude oscillation amplitude (m).
- `omega`: angular frequency (rad/s).
"""
struct BobbingTarget
    y_fixed::Float64
    z_center::Float64
    z_amplitude::Float64
    omega::Float64
end

"""
    trajectory(bt::BobbingTarget, t_range::AbstractVector{<:Integer}, h, nx, nu)

Return a `Vector` of `length(t_range)` stalk vectors for `bt`, one per integer timestep.

Altitude follows ``z(t) = z_c + A\\sin(\\omega t_s)`` where ``t_s = t \\cdot h``;
``\\dot z`` is set to the analytic derivative.  Lateral position is `y_fixed`;
all other state components and all controls are zero.

`t_range` can be any integer range, e.g. `0:k` or `10:20`, enabling trajectory
evaluation over an arbitrary sub-interval without recomputing from scratch.
"""
function trajectory(bt::BobbingTarget, t_range::AbstractVector{<:Integer}, h::Float64, nx::Int, nu::Int)
    return map(t_range) do t
        t_s = t * h
        z_t    = bt.z_center + bt.z_amplitude * sin(bt.omega * t_s)
        zdot_t = bt.z_amplitude * bt.omega    * cos(bt.omega * t_s)
        x = zeros(nx)
        x[IDX_Y]   = bt.y_fixed
        x[IDX_Z]   = z_t
        x[IDX_ZDT] = zdot_t
        vcat(x, zeros(nu))
    end
end

# ## Plotting Recipe
#
# `ScenarioResult` bundles a solved scenario for plotting.
# The `@recipe` below produces a two-panel ``y(t)`` / ``z(t)`` figure
# without repeating the same 10-line block for each scenario.

"""
    ScenarioResult

Bundle of trajectory data and statistics for one solved scenario.
Plot via `plot(result::ScenarioResult)`.

Fields:
- `label`: scenario name used in the plot title.
- `times`: discrete time vector (length `k+1`).
- `agent_trajs`: `Vector` of `(k+1) × nx` state matrices, one per agent.
- `target_trajs`: `Vector` of stalk-vector sequences, one per target.
- `null_dim`: nullspace dimension of the restricted Laplacian.
- `residual`: Laplacian energy ``\\sqrt{z^\\top L z}``.
"""
struct ScenarioResult
    label::String
    times::Vector{Float64}
    agent_trajs::Vector{Matrix{Float64}}
    target_trajs::Vector{Vector{Vector{Float64}}}
    null_dim::Int
    residual::Float64
end

@recipe function f(sr::ScenarioResult)
    layout := (1, 2)
    size := (800, 380)
    plot_title := "$(sr.label): null dim = $(sr.null_dim),  ||dz|| = $(round(sr.residual; sigdigits=3))"
    agent_colors  = [:steelblue, :darkorange, :green, :crimson]
    agent_styles  = [:solid, :dash, :dashdot, :dot]
    target_colors = [:gray, :black, :darkgreen, :purple]
    for (i, traj) in enumerate(sr.agent_trajs)
        @series begin
            subplot   := 1
            title     := "y(t)"
            xlabel    := "t (s)"
            ylabel    := "y (m)"
            label     := "A$i"
            lw        := 2
            linecolor := agent_colors[i]
            linestyle := agent_styles[i]
            sr.times, traj[:, IDX_Y]
        end
        @series begin
            subplot   := 2
            title     := "z(t)"
            xlabel    := "t (s)"
            ylabel    := "z (m)"
            label     := "A$i"
            lw        := 2
            linecolor := agent_colors[i]
            linestyle := agent_styles[i]
            sr.times, traj[:, IDX_Z]
        end
    end
    for (j, traj_j) in enumerate(sr.target_trajs)
        @series begin
            subplot   := 1
            label     := "T$j"
            lw        := 1
            linecolor := target_colors[j]
            linestyle := :dot
            sr.times, getindex.(traj_j, IDX_Y)
        end
        @series begin
            subplot   := 2
            label     := "T$j"
            lw        := 1
            linecolor := target_colors[j]
            linestyle := :dot
            sr.times, getindex.(traj_j, IDX_Z)
        end
    end
end

# ## Helper Functions

"""
    selector_matrix(indices, n)

Return the ``|\\text{indices}| \\times n`` row-selection matrix ``S`` such that
``S v = v[\\text{indices}]``.  Throws `ArgumentError` if any index is outside `1:n`.
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

Return a ``|\\text{state\\_indices}| \\times (n_x + n_u)`` projection that selects
the given state coordinates from an augmented ``(n_x + n_u)``-stalk.
"""
function state_projection_matrix(state_indices::AbstractVector{<:Integer}, nx::Int, nu::Int)
    sel = selector_matrix(state_indices, nx)
    return hcat(sel, zeros(length(state_indices), nu))
end

"""
    agent_vertex(prob, i, t)  /  target_vertex(prob, j, t)

Vertex index of agent `i` (or target `j`) at timestep `t` in the sheaf
built from `prob`.  Vertices are laid out as
`t * n_per_t + 1 … t * n_per_t + n_agents` (agents)
followed by `t * n_per_t + n_agents + 1 … t * n_per_t + n_agents + n_targets` (targets).
"""
agent_vertex(prob::TrackingProblem, i::Int, t::Int) =
    t * (prob.n_agents + prob.n_targets) + i

target_vertex(prob::TrackingProblem, j::Int, t::Int) =
    t * (prob.n_agents + prob.n_targets) + prob.n_agents + j

"""
    build_time_expanded_tracking_sheaf(prob::TrackingProblem) -> EuclideanSheaf

Construct the time-expanded `EuclideanSheaf` from a `TrackingProblem`.

Edge families:
1. Agent temporal dynamics: `agent(i,t) ↔ agent(i,t+1)` encodes `x_{t+1} = A_d x_t + B_d u_t`.
2. Target temporal dynamics (when `include_target_dynamics`): same structure.
3. Consensus: `agent(i,t) ↔ agent(j,t)` for `(i,j) ∈ agent_edges`, active at `consensus_timesteps`.
4. Tracking: each `TrackingEdge(i,j,Rₐ,Rₜ)` adds `agent(i,t) ↔ target(j,t)`,
   active at `tracking_timesteps`, using restriction maps `Rₐ` (agent side) and `Rₜ` (target side).
"""
function build_time_expanded_tracking_sheaf(prob::TrackingProblem)
    nx_loc = size(prob.Ad, 1)
    nu_loc = size(prob.Bd, 2)
    nt     = nx_loc + nu_loc
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

"""
    generate_reference_trajectory(x0, xk, k, Ad, Bd, nx, nu)

Compute the minimum-energy trajectory from `x0` to `xk` in `k` steps by harmonic
extension on a single-agent dynamics sheaf.

The stalk at each step is ``(x_t, u_t) \\in \\mathbb{R}^{n_x + n_u}``.  Interior
controls ``u_t`` for ``t = 1, \\ldots, k-1`` will in general be **nonzero**: they
are the minimum-energy inputs that drive the state from `x0` to `xk` while
satisfying the linear dynamics.  Since this function is used to produce fully-pinned
target boundary data, the control components do not affect the multi-agent solution.
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
    extract_state_trajectories(z_harmonic, prob) -> Vector{Matrix{Float64}}

Return one ``(k+1) \\times n_x`` matrix per agent.  Matrices are pre-allocated;
no intermediate concatenations are performed.
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
    run_scenario(label, prob, boundary, times; target_trajs=[]) -> ScenarioResult

Build the sheaf from `prob`, run `harmonic_extension`, compute the Laplacian energy
``\\sqrt{z^\\top L z}``, and return a `ScenarioResult` ready to plot.

`sheaf_laplacian_matrix_direct` is used for the residual (rather than `coboundary_map`)
because the coboundary matrix only spans vertices that appear in at least one edge;
isolated target vertices (when `include_target_dynamics = false`) cause a column-count
mismatch.  `sheaf_laplacian_matrix_direct` always covers all stalk dimensions.
"""
function run_scenario(
    label::String,
    prob::TrackingProblem,
    boundary::Dict{Int,Vector{Float64}},
    times::AbstractVector{<:Real};
    target_trajs::Vector{Vector{Vector{Float64}}} = Vector{Vector{Vector{Float64}}}(),
)
    sheaf = build_time_expanded_tracking_sheaf(prob)
    z_harmonic, null_basis = harmonic_extension(sheaf, boundary)
    L  = sheaf_laplacian_matrix_direct(sheaf)
    xv = Array(z_harmonic)
    r  = sqrt(max(0.0, dot(xv, L * xv)))
    nd = size(null_basis, 2)
    println("$label:  ||dz|| = $(round(r; sigdigits=4)),   null dim = $nd")
    trajs = extract_state_trajectories(z_harmonic, prob)
    return ScenarioResult(label, collect(times), trajs, target_trajs, nd, r)
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
times     = h .* collect(0:k)

R_yz = state_projection_matrix([IDX_Y, IDX_Z], nx, nu)  # 2×8: (y, z)
R_y  = state_projection_matrix([IDX_Y],         nx, nu)  # 1×8: y only
R_z  = state_projection_matrix([IDX_Z],         nx, nu)  # 1×8: z only

# Bobbing targets used in Scenarios 2–4.
# Both complete two full vertical cycles over the 2-second horizon.
omega_2periods = 2π * 2 / (k * h)
bt1 = BobbingTarget(0.0, 1.0, 0.3, omega_2periods)   # target 1: y = 0, z ≈ 1 m
bt2 = BobbingTarget(0.0, 2.0, 0.3, omega_2periods)   # target 2: y = 0, z ≈ 2 m
traj_bt1 = trajectory(bt1, 0:k, h, nx, nu)
traj_bt2 = trajectory(bt2, 0:k, h, nx, nu)

# Named edge configurations reused across scenarios.
edges_12 = [(1, 2)]  # consensus: agents 1 and 2 must agree
te_yz = [TrackingEdge(1, 1, R_yz, R_yz), TrackingEdge(2, 2, R_yz, R_yz)]  # Scenario 1
te_z  = [TrackingEdge(1, 1, R_z,  R_z),  TrackingEdge(2, 2, R_z,  R_z)]   # Scenarios 2–4

# ## Scenario 1: Full (y,z) Coordination at Every Timestep
#
# This scenario imposes the maximum coordination: both consensus edges and tracking
# edges use the ``(y, z)`` projection, active at **every** timestep.
#
# The two agents travel to *distinct* goals while each tracking its own target in
# full ``(y, z)``.  But the consensus edge demands that agents agree in ``(y, z)``
# at every step.  Since the two targets travel *different* ``(y, z)`` trajectories,
# the combination "A1 tracks T1", "A2 tracks T2", and "A1 agrees with A2" is
# **infeasible** unless the targets coincide.  The sheaf Laplacian has no kernel
# (null dim = 0) and `harmonic_extension` returns the least-squares compromise
# with a nonzero residual ``\|dz\| > 0``.  This is Scenario 1's design intent:
# to show that over-determined coordination constraints collapse the solution space
# to a single point (the energy minimum), not a feasible trajectory.

x0_a1_s1 = [0.0, 1.0, 0.0, 0.0, 0.0, 0.0]
xk_a1_s1 = [2.0, 1.5, 0.0, 0.0, 0.0, 0.0]
x0_a2_s1 = [0.0, 2.0, 0.0, 0.0, 0.0, 0.0]
xk_a2_s1 = [2.0, 2.5, 0.0, 0.0, 0.0, 0.0]

traj_t1_s1 = generate_reference_trajectory(x0_a1_s1, xk_a1_s1, k, Ad, Bd, nx, nu)
traj_t2_s1 = generate_reference_trajectory(x0_a2_s1, xk_a2_s1, k, Ad, Bd, nx, nu)

prob1 = TrackingProblem(
    n_agents, n_targets, k, Matrix(Ad), Matrix(Bd),
    edges_12, te_yz, R_yz,
    collect(0:k), collect(0:k),
    true, 1.0, 5.0,
)

bnd1 = Dict{Int,Vector{Float64}}()
bnd1[agent_vertex(prob1, 1, 0)] = vcat(x0_a1_s1, zeros(nu))
bnd1[agent_vertex(prob1, 2, 0)] = vcat(x0_a2_s1, zeros(nu))
bnd1[agent_vertex(prob1, 1, k)] = vcat(xk_a1_s1, zeros(nu))
bnd1[agent_vertex(prob1, 2, k)] = vcat(xk_a2_s1, zeros(nu))
for t in 0:k
    bnd1[target_vertex(prob1, 1, t)] = traj_t1_s1[t + 1]
    bnd1[target_vertex(prob1, 2, t)] = traj_t2_s1[t + 1]
end

result1 = run_scenario("Scenario 1", prob1, bnd1, times;
    target_trajs = [traj_t1_s1, traj_t2_s1])

plot(result1)

# ## Scenario 2: y-Consensus and z-Tracking, All Timesteps, Aligned Initial Conditions
#
# Targets bob vertically.  Consensus edges constrain agents to agree only on
# **lateral position** ``y``; tracking edges constrain each agent to match its
# assigned target in **altitude** ``z``.  Both edge families are active at every
# timestep.
#
# Each agent–target pair now uses different restriction maps (``R_y`` vs.\ ``R_z``),
# illustrating the many-to-many `TrackingEdge` API.  Because the two targets share
# ``y = 0`` and the consensus only constrains ``y``, the system is consistent:
# agents can agree in ``y`` while independently tracking different altitudes.
# A non-trivial nullspace (null dim > 0) arises because ``z`` (and all other state
# components) remain unconstrained by the consensus edges.
#
# Since all target vertices are pinned as boundary data, target dynamics edges
# are omitted.

prob2 = TrackingProblem(
    n_agents, n_targets, k, Matrix(Ad), Matrix(Bd),
    edges_12, te_z, R_y,
    collect(0:k), collect(0:k),
    false, 1.0, 5.0,
)

x0_a1_s2 = [0.0, traj_bt1[1][IDX_Z], 0.0, 0.0, traj_bt1[1][IDX_ZDT], 0.0]
x0_a2_s2 = [0.0, traj_bt2[1][IDX_Z], 0.0, 0.0, traj_bt2[1][IDX_ZDT], 0.0]

bnd2 = Dict{Int,Vector{Float64}}()
bnd2[agent_vertex(prob2, 1, 0)] = vcat(x0_a1_s2, zeros(nu))
bnd2[agent_vertex(prob2, 2, 0)] = vcat(x0_a2_s2, zeros(nu))
for t in 0:k
    bnd2[target_vertex(prob2, 1, t)] = traj_bt1[t + 1]
    bnd2[target_vertex(prob2, 2, t)] = traj_bt2[t + 1]
end

result2 = run_scenario("Scenario 2", prob2, bnd2, times;
    target_trajs = [traj_bt1, traj_bt2])

plot(result2)

# ## Scenario 3: y-Consensus and z-Tracking at the Last Timestep Only
#
# The same projection matrices as Scenario 2, but coordination edges are
# active **only at** ``t = k``.  Agents evolve freely (dynamics only) for
# ``t = 0, \ldots, k-1`` and must converge to the terminal constraint.
#
# Both targets share ``y = 0`` (aligned laterally), but agents start misaligned
# with each other and with the targets in both ``y`` and ``z``.
# Removing ``k`` timesteps of coordination edges dramatically enlarges the
# feasible space (null dim >> Scenario 2).

prob3 = TrackingProblem(
    n_agents, n_targets, k, Matrix(Ad), Matrix(Bd),
    edges_12, te_z, R_y,
    [k], [k],
    false, 1.0, 5.0,
)

x0_a1_s3 = [-0.5, 0.5, 0.0, 0.0, 0.0, 0.0]
x0_a2_s3 = [ 1.0, 1.5, 0.0, 0.0, 0.0, 0.0]

bnd3 = Dict{Int,Vector{Float64}}()
bnd3[agent_vertex(prob3, 1, 0)] = vcat(x0_a1_s3, zeros(nu))
bnd3[agent_vertex(prob3, 2, 0)] = vcat(x0_a2_s3, zeros(nu))
for t in 0:k
    bnd3[target_vertex(prob3, 1, t)] = traj_bt1[t + 1]
    bnd3[target_vertex(prob3, 2, t)] = traj_bt2[t + 1]
end

result3 = run_scenario("Scenario 3", prob3, bnd3, times;
    target_trajs = [traj_bt1, traj_bt2])

plot(result3)

# ## Scenario 4: Last-Timestep Constraints, Targets Not Aligned in y or z
#
# Identical sheaf topology to Scenario 3, but Target 2 is now at ``y = 1.5`` m
# instead of ``y = 0``.  The targets are **not** initially aligned with each other
# in either ``y`` or ``z``.
#
# Because the sheaf topology (which edges exist at which timesteps) is unchanged,
# the restricted Laplacian's nullspace has the **same dimension** as Scenario 3.
# The two harmonic extensions produce different trajectories, but they live in
# equally large affine solution spaces, confirming that nullspace dimension is a
# structural property of the coordination architecture, not of the boundary data.

bt2_s4 = BobbingTarget(1.5, 2.0, 0.3, omega_2periods)
traj_bt2_s4 = trajectory(bt2_s4, 0:k, h, nx, nu)

prob4 = TrackingProblem(
    n_agents, n_targets, k, Matrix(Ad), Matrix(Bd),
    edges_12, te_z, R_y,
    [k], [k],
    false, 1.0, 5.0,
)

x0_a1_s4 = [-0.5, 0.5, 0.0, 0.0, 0.0, 0.0]
x0_a2_s4 = [ 1.0, 1.5, 0.0, 0.0, 0.0, 0.0]

bnd4 = Dict{Int,Vector{Float64}}()
bnd4[agent_vertex(prob4, 1, 0)] = vcat(x0_a1_s4, zeros(nu))
bnd4[agent_vertex(prob4, 2, 0)] = vcat(x0_a2_s4, zeros(nu))
for t in 0:k
    bnd4[target_vertex(prob4, 1, t)] = traj_bt1[t + 1]
    bnd4[target_vertex(prob4, 2, t)] = traj_bt2_s4[t + 1]
end

result4 = run_scenario("Scenario 4", prob4, bnd4, times;
    target_trajs = [traj_bt1, traj_bt2_s4])

plot(result4)

# ## Comparison: How Constraints Shape the Solution Space
#
# The four scenarios are summarised below.  The nullspace dimension and
# Laplacian residual ``\|dz\|`` are printed at runtime.
#
# | Scenario | Consensus | Tracking | Active timesteps | Initial alignment |
# |----------|-----------|----------|-----------------|-------------------|
# | 1 | ``(y,z)`` | ``(y,z)`` | all ``0:k`` | distinct endpoints |
# | 2 | ``y``     | ``z``     | all ``0:k`` | aligned with targets |
# | 3 | ``y``     | ``z``     | terminal ``k`` | unaligned; targets share ``y`` |
# | 4 | ``y``     | ``z``     | terminal ``k`` | unaligned; targets offset in ``y`` |
#
# Three observations stand out:
#
# **Constraint density governs nullspace size.**
# Scenario 1 (``(y,z)`` constrained at every step, two distinct agent goals)
# is **infeasible**: the consensus edge forces agents to share the same ``(y,z)``
# trajectory, but tracking edges demand they follow different targets.  The result
# is null dim = 0 with a nonzero least-squares residual.
# Scenario 2 (only ``y`` consensus + ``z`` tracking, all timesteps) is consistent
# and admits a family of solutions (null dim > 0) because fewer coordinates per
# edge are constrained.
# Restricting to the terminal timestep only (Scenarios 3–4) removes ``k`` sets of
# constraints and dramatically enlarges the solution space.
#
# **Nullspace dimension is a topological invariant.**
# Scenarios 3 and 4 share the same sheaf topology (identical edge structure) and
# therefore the same nullspace dimension, even though their target trajectories
# differ (targets aligned vs. laterally offset).  The harmonic extensions produce
# different specific trajectories, but they live in equally large affine spaces.
#
# **Feasibility is separate from trajectory shape.**
# In Scenario 2, agents start consistent with all-time constraints, so the harmonic
# solution closely tracks the bobbing targets.  In Scenarios 3–4, agents start
# misaligned; the minimum-energy harmonic path converges to the terminal constraint
# while the large nullspace deforms it freely in the unconstrained directions.

println("Summary: null dims = [$(result1.null_dim), $(result2.null_dim), $(result3.null_dim), $(result4.null_dim)]")

bar(["S1\n(y,z) all-t", "S2\ny/z all-t", "S3\ny/z last-t", "S4\ny/z last-t\n(mis-tgt)"],
    [result1.null_dim, result2.null_dim, result3.null_dim, result4.null_dim];
    ylabel = "Nullspace dimension",
    title  = "Nullspace dimension vs. coordination architecture",
    legend = false,
    color  = [:steelblue, :darkorange, :green, :crimson])
