"""
    MultiAgentTracking

Utilities for multi-agent, multi-target coordination via time-expanded cellular sheaves.

Exports types for problem specification (`TrackingEdge`, `TrackingProblem`,
`BobbingTarget`, `ScenarioResult`) and functions to build the sheaf, generate
reference trajectories, solve, and extract results.
"""
module MultiAgentTracking

using LinearAlgebra
using SparseArrays
using BlockArrays
using ArgCheck
using CliqueTrees.Multifrontal

using ...NetworkSheaves: EuclideanSheaf, add_sheaf_edge!, harmonic_extension,
                          sheaf_laplacian_matrix,
                          ldlt_pseudoinverse_and_null, ldlt_pinv_solve

export TrackingEdge, TrackingProblem, BobbingTarget, ScenarioResult
export trajectory
export selector_matrix, state_projection_matrix
export agent_vertex, target_vertex
export build_time_expanded_tracking_sheaf
export generate_reference_trajectory, extract_state_trajectories, extract_target_trajectories, run_scenario
export animate_tracking_xy
export run_mpc_scenario
export WindowSolverCache, tracking_extension_operator
export window_targets, window_problem


# ---------------------------------------------------------------------------
# Helper functions
# ---------------------------------------------------------------------------

"""
    selector_matrix(indices, n) -> Matrix{Float64}

`length(indices) × n` row-selection matrix `S` such that `S * v == v[indices]`.
"""
function selector_matrix(indices::AbstractVector{<:Integer}, n::Int)
    @argcheck all(i -> 1 <= i <= n, indices) "All indices must be in 1:$n, got $indices"
    S = zeros(Float64, length(indices), n)
    for (row, col) in enumerate(indices); S[row, col] = 1.0; end
    return S
end

"""
    state_projection_matrix(state_indices, nx, nu) -> Matrix{Float64}

`length(state_indices) × (nx + nu)` matrix selecting state coordinates from an augmented stalk.
"""
state_projection_matrix(state_indices, nx::Int, nu::Int) =
    hcat(selector_matrix(state_indices, nx), zeros(length(state_indices), nu))

include("Targets.jl")
include("TrackingProblems.jl")
include("ScenarioResults.jl")
include("QuadraticCosts.jl")
using .QuadraticCosts: build_control_cost_matrix, solve_quadratic_on_basis


"""
    run_scenario(label, prob, boundary, times; target_trajs, y_col, z_col)

Build the sheaf from `prob`, run `harmonic_extension`, compute the Laplacian
energy `sqrt(z' * L * z)`, and return a `ScenarioResult`.

`sheaf_laplacian_matrix` is used for the energy (rather than
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
    target_trajs::Union{Nothing,Vector{Vector{Vector{Float64}}}} = nothing,
    y_col::Int = 1,
    z_col::Int = 2,
    cost::Union{Number,Function}=1.0   # (agent, time) -> Q_control matrix or a number
)
    # -------------------------------------------------------------------
    # Build the sheaf and obtain the harmonic‑extension solution
    # -------------------------------------------------------------------
    sheaf = build_time_expanded_tracking_sheaf(prob)
    z_harmonic, null_basis = harmonic_extension(sheaf, boundary)
    L  = sheaf_laplacian_matrix(sheaf)
    z_harmonic_array = Array(z_harmonic)
    nd = size(null_basis, 2)

    # -------------------------------------------------------------------
    # Optional quadratic control cost
    # -------------------------------------------------------------------
    # Build the global quadratic cost matrix (only on control components) and
    # solve the reduced problem on the nullspace of the harmonic‑extension
    # constraints.  The heavy lifting is delegated to `QuadraticCost` utilities.
    if cost != 0.0
        Q = build_control_cost_matrix(prob, sheaf.vertex_stalks, cost)
        z_opt = solve_quadratic_on_basis(z_harmonic_array, null_basis, Q)
    else
        # No extra cost – just keep the harmonic solution.
        z_opt = z_harmonic_array
    end
    residual = sqrt(max(0.0, dot(z_opt, L * z_opt)))

    # -------------------------------------------------------------------
    # Extract trajectories from the (possibly) optimized solution
    # -------------------------------------------------------------------
    # Re‑wrap the vector as a BlockVector so the existing extraction utilities
    # can operate unchanged.
    z_block = BlockArray(z_opt, sheaf.vertex_stalks)
    trajs = extract_state_trajectories(z_block, prob)
    tt = isnothing(target_trajs) ? extract_target_trajectories(z_block, prob) : target_trajs
    return ScenarioResult(label, collect(times), trajs, tt, nd, residual, y_col, z_col)
end

"""
    window_targets(target_trajs, t0, t1)

Return target trajectories sliced to integer timesteps `t0:t1`.
"""
window_targets(target_trajs, t0::Int, t1::Int) = [traj[t0+1:t1+1] for traj in target_trajs]

"""
    window_problem(prob, t, t_end) -> TrackingProblem

Extract the sub-problem for timesteps `[t, t_end]` from a full `TrackingProblem`.
Activation timesteps are shifted to be relative to the window start.
"""
function window_problem(prob::TrackingProblem, t::Int, t_end::Int)
    @argcheck 0 <= t <= t_end <= prob.k "window [$t, $t_end] must be within [0, $(prob.k)]"
    W = t_end - t
    return TrackingProblem(
        prob.n_agents, prob.n_targets, W,
        prob.Ad, prob.Bd, prob.target_Ad, prob.target_Bd,
        prob.agent_edges, prob.tracking_edges, prob.consensus_restriction,
        [ts - t for ts in prob.consensus_timesteps if t <= ts <= t_end],
        [ts - t for ts in prob.tracking_timesteps  if t <= ts <= t_end],
        prob.include_target_dynamics, prob.consensus_weight, prob.tracking_weight,
    )
end

"""
    _assemble_boundary(prob, x_agents, target_window) -> Dict{Int,Vector{Float64}}

Build the vertex-indexed boundary dictionary for `harmonic_extension` over the
MPC-window sheaf (built with `state_only_initial = true`):
- Agent `i`'s state-only vertex at the window's first timestep is pinned to `x_agents[i]`.
- Target `j`'s full stalk at each timestep is pinned to `target_window[j][t+1]`.

Both are full-vertex pins, so this returns a `Dict{Int,Vector{Float64}}`; the
target stalk dimensions are checked by `harmonic_extension`.
`target_window` must have `prob.k + 1` samples per target.
"""
function _assemble_boundary(
    prob::TrackingProblem,
    x_agents::Vector{Vector{Float64}},
    target_window::Vector{Vector{Vector{Float64}}},
)
    pinned = Dict{Int,Vector{Float64}}()
    for i in 1:prob.n_agents
        nx_i = size(prob.Ad[i], 1)
        @argcheck length(x_agents[i]) == nx_i "x_agents[$i] must have length nx=$nx_i, got $(length(x_agents[i]))"
        pinned[agent_vertex(prob, i, 0)] = x_agents[i]
    end
    for j in 1:prob.n_targets, ts in 0:prob.k
        pinned[target_vertex(prob, j, ts)] = target_window[j][ts + 1]
    end
    return pinned
end

"""
    _apply_first_control(window_prob, z_opt, stalks) -> Vector{Vector{Float64}}

Advance each agent one step using the first decision control of the window.
In the reindexed window sheaf the initial vertex (`t=0`) is state-only and the
control `u_1` on vertex `t=1` drives the first transition, so this computes
`x_{t+1} = Ad_i * x_0 + Bd_i * u_1`.  Both `x_0` (the pinned initial state) and
`u_1` are read from the section.
"""
function _apply_first_control(
    window_prob::TrackingProblem,
    z_opt::AbstractVector,
    stalks::Vector{Int},
)
    z_block = BlockArray(z_opt, stalks)
    return map(1:window_prob.n_agents) do i
        nx_i = size(window_prob.Ad[i], 1)
        nu_i = size(window_prob.Bd[i], 2)
        x_0  = Array(z_block[Block(agent_vertex(window_prob, i, 0))])[1:nx_i]
        u_1  = Array(z_block[Block(agent_vertex(window_prob, i, 1))])[nx_i+1 : nx_i+nu_i]
        window_prob.Ad[i] * x_0 + window_prob.Bd[i] * u_1
    end
end

# Single source of truth for the window-sheaf boundary ordering.  Visits every
# pinned block in canonical order — each agent's initial state, then every target
# stalk by timestep — calling `f(dofs, dest, vals)`: `dofs` are the block's
# columns in the global Laplacian, `dest` its slice of the packed boundary vector
# `x_B`, and `vals` the pinned values (`nothing` at operator-build time, when only
# the DOF layout is needed).  Both `_window_boundary_dofs` and
# `_fill_boundary_vector!` route through this, so their orderings cannot drift.
function _foreach_window_boundary(f, prob::TrackingProblem, stalks::AbstractVector{<:Integer};
                                  x_agents=nothing, targets=nothing)
    offsets = [0; cumsum(stalks)]
    pos = 0
    for i in 1:prob.n_agents
        v   = agent_vertex(prob, i, 0)
        nxi = size(prob.Ad[i], 1)
        vals = isnothing(x_agents) ? nothing : view(x_agents[i], 1:nxi)
        f(offsets[v]+1 : offsets[v]+nxi, pos+1 : pos+nxi, vals)
        pos += nxi
    end
    for ts in 0:prob.k, j in 1:prob.n_targets
        v  = target_vertex(prob, j, ts)
        nd = stalks[v]
        vals = isnothing(targets) ? nothing : targets[j][ts + 1]
        f(offsets[v]+1 : offsets[v]+nd, pos+1 : pos+nd, vals)
        pos += nd
    end
    return nothing
end

# Boundary DOFs of a window sheaf in canonical order (agent initial states, then
# target stalks by timestep), derived from `_foreach_window_boundary` so it stays
# aligned with the `x_B` packing by construction.
function _window_boundary_dofs(prob::TrackingProblem, stalks::AbstractVector{<:Integer})
    dofs = Int[]
    _foreach_window_boundary((d, _dest, _vals) -> append!(dofs, d), prob, stalks)
    return dofs
end

# Pack the boundary values into `x_B` in the same canonical order, sharing
# `_foreach_window_boundary` with `_window_boundary_dofs` so the operator's
# columns and `x_B` stay aligned.
function _fill_boundary_vector!(x_B, prob::TrackingProblem, stalks, x_agents, targets)
    _foreach_window_boundary(prob, stalks; x_agents=x_agents, targets=targets) do _dofs, dest, vals
        copyto!(view(x_B, dest), vals)
    end
    return x_B
end

# Shared tail of both solver paths: window Laplacian energy of the section and the
# first applied control advancing each agent one step.
function _finish_step(local_prob::TrackingProblem, z, L, stalks)
    residual = sqrt(max(0.0, dot(z, L * z)))
    x_now    = _apply_first_control(local_prob, z, stalks)
    return x_now, residual
end


# ---------------------------------------------------------------------------
# Cost-optimal harmonic extension operator (cached MPC core)
# ---------------------------------------------------------------------------

"""
    WindowSolverCache

Reusable solver state for the cost-optimal harmonic extension of a fixed window
sheaf, built once per window length and cached across MPC steps.  Each step then
reduces to `z[interior] = M * x_B` — a single dense matvec.
Construct with [`tracking_extension_operator`](@ref).

Fields: `M` (interior×boundary), `boundary`/`interior` (DOF index vectors),
`stalks`, `laplacian` (for energy), `null_dim`, `window`.
"""
struct WindowSolverCache
    M         :: Matrix{Float64}
    boundary  :: Vector{Int}
    interior  :: Vector{Int}
    stalks    :: Vector{Int}
    laplacian :: SparseMatrixCSC{Float64,Int}
    null_dim  :: Int
    window    :: Int
    scratch   :: Vector{Float64}
end

"""
    tracking_extension_operator(window_prob; cost=1.0) -> WindowSolverCache

Factorize the window Laplacian once and fold in the control-cost nullspace
projection to produce a dense operator `M` so each MPC step is `M * x_B`.
"""
function tracking_extension_operator(window_prob::TrackingProblem; cost::Union{Number,Function}=1.0)
    W      = window_prob.k
    sheaf  = build_time_expanded_tracking_sheaf(window_prob; state_only_initial=true)
    L      = sparse(sheaf_laplacian_matrix(sheaf))
    stalks = sheaf.vertex_stalks

    boundary = _window_boundary_dofs(window_prob, stalks)
    interior = setdiff(1:size(L, 1), boundary)

    # L is already sparse, so indexing it with index vectors returns sparse
    # blocks — no explicit `sparse(...)` needed (asserted by a unit test).
    L_II = L[interior, interior]
    L_IB = L[interior, boundary]
    F    = ldlt!(ChordalLDLt(L_II), RowMaximum(); check=false)
    tol  = sqrt(eps(Float64)) * max(1.0, maximum(i -> abs(F.D[i,i]), 1:size(F.D,1); init=0.0))
    _, N = ldlt_pseudoinverse_and_null(F, zeros(length(interior)); tol=tol)

    # Harmonic response to every boundary DOF: H = -L_II⁺ L_IB (pinv, many RHS).
    # `H` is dense `n_I × n_B`, so densifying `L_IB` here costs no extra order of
    # memory and avoids sparse-column scalar getindex in the per-column solve.
    H = ldlt_pinv_solve(F, Matrix(L_IB); tol=tol)
    H .*= -1

    if size(N, 2) > 0 && !(cost isa Number && iszero(cost))
        # Keep the interior cost block sparse — only the small `nd × n_I` product
        # `NtQ_II` is materialized dense, never the full `n_I × n_I` matrix.
        Q_II   = build_control_cost_matrix(window_prob, stalks, cost)[interior, interior]
        NtQ_II = N' * Q_II
        NtQN   = NtQ_II * N
        Fq     = ldlt!(DenseLDLtPivoted(NtQN), RowMaximum(); check=false)
        nd_q   = size(N, 2)
        tol_q  = nd_q * eps(Float64) * max(1.0, maximum(i -> abs(Fq.D[i,i]), 1:nd_q; init=0.0))
        coeff  = ldlt_pinv_solve(Fq, NtQ_II * H; tol=tol_q)
        M = H - N * coeff
    else
        M = H
    end

    return WindowSolverCache(M, boundary, interior, stalks, L, size(N, 2), W,
                             Vector{Float64}(undef, length(interior)))
end

"""
    mul!(z, op::WindowSolverCache, x_B) -> z

Apply the cached operator in place: scatter the boundary values `x_B` and write
the harmonic interior `M * x_B` directly into `z` (length `sum(op.stalks)`),
allocating nothing.
"""
function LinearAlgebra.mul!(z::AbstractVector, op::WindowSolverCache, x_B::AbstractVector)
    @argcheck length(x_B) == length(op.boundary)
    @argcheck length(z) == sum(op.stalks)
    fill!(z, 0.0)
    z[op.boundary] = x_B
    mul!(op.scratch, op.M, x_B)
    z[op.interior] = op.scratch
    return z
end

(op::WindowSolverCache)(x_B::AbstractVector) = mul!(zeros(sum(op.stalks)), op, x_B)

_uniform_activation(ts, k::Int) = isempty(ts) || sort(unique(ts)) == collect(0:k)

"""
    run_mpc_scenario(label, prob, x0_agents, target_trajs, times;
                     window, y_col=1, z_col=2, cost=1.0, solver=:cached) -> ScenarioResult

Receding-horizon MPC: at each step `t` solves the sheaf QP on `[t, min(t+window,k)]`,
applies the first control, and advances agent dynamics.  `solver=:cached` reuses a
[`WindowSolverCache`](@ref) for recurring window lengths; falls back to `:naive` on
inhomogeneous problems.
"""
function run_mpc_scenario(
    label::String,
    prob::TrackingProblem,
    x0_agents::Vector{Vector{Float64}},
    target_trajs::Vector{Vector{Vector{Float64}}},
    times::AbstractVector{<:Real};
    window::Int,
    y_col::Int = 1,
    z_col::Int = 2,
    cost::Union{Number,Function} = 1.0,
    solver::Symbol = :cached,
)
    k = prob.k
    @argcheck window >= 1 "window must be >= 1, got $window"
    @argcheck solver in (:cached, :naive) "solver must be :cached or :naive, got :$solver"
    @argcheck length(x0_agents) == prob.n_agents "x0_agents length must equal n_agents=$(prob.n_agents), got $(length(x0_agents))"
    @argcheck length(target_trajs) == prob.n_targets "target_trajs length must equal n_targets=$(prob.n_targets), got $(length(target_trajs))"
    @argcheck length(times) == k + 1 "times must have length k+1=$(k+1), got $(length(times))"
    for j in 1:prob.n_targets
        @argcheck length(target_trajs[j]) >= k + 1 "target_trajs[$j] must have >= k+1=$(k+1) elements, got $(length(target_trajs[j]))"
    end
    # Per-entity dimension checks (cover both the cached and naive paths, which
    # otherwise surface mismatches as opaque BoundsErrors deeper in the solve).
    for i in 1:prob.n_agents
        nx_i = size(prob.Ad[i], 1)
        @argcheck length(x0_agents[i]) == nx_i "x0_agents[$i] must have length nx=$nx_i, got $(length(x0_agents[i]))"
    end
    for j in 1:prob.n_targets
        nstalk = size(prob.target_Ad[j], 1) + size(prob.target_Bd[j], 2)
        for s in 1:(k + 1)
            @argcheck length(target_trajs[j][s]) == nstalk "target_trajs[$j][$s] must have length $nstalk (target stalk), got $(length(target_trajs[j][s]))"
        end
    end

    use_cached = solver === :cached &&
                 _uniform_activation(prob.consensus_timesteps, k) &&
                 _uniform_activation(prob.tracking_timesteps, k)

    nx          = [size(prob.Ad[i], 1) for i in 1:prob.n_agents]
    x_now       = copy.(x0_agents)
    agent_trajs = [Matrix{Float64}(undef, k+1, nx[i]) for i in 1:prob.n_agents]
    for i in 1:prob.n_agents; agent_trajs[i][1, :] = x_now[i]; end

    window_counts = Dict{Int,Int}()
    if use_cached
        for t in 0:k-1
            w = min(t + window, k) - t
            window_counts[w] = get(window_counts, w, 0) + 1
        end
    end
    ops      = Dict{Int,WindowSolverCache}()
    xbufs    = Dict{Int,Vector{Float64}}()   # reused boundary vectors, by window length
    zbufs    = Dict{Int,Vector{Float64}}()   # reused full sections, by window length
    nd       = 0
    residual = 0.0

    for t in 0:k-1
        t_end         = min(t + window, k)
        local_prob    = window_problem(prob, t, t_end)
        local_targets = window_targets(target_trajs, t, t_end)

        if use_cached && get(window_counts, t_end - t, 0) >= 2
            op = get!(ops, t_end - t) do; tracking_extension_operator(local_prob; cost=cost); end
            t == 0 && (nd = op.null_dim)
            x_B = get!(() -> Vector{Float64}(undef, length(op.boundary)), xbufs, t_end - t)
            _fill_boundary_vector!(x_B, local_prob, op.stalks, x_now, local_targets)
            z = get!(() -> Vector{Float64}(undef, sum(op.stalks)), zbufs, t_end - t)
            mul!(z, op, x_B)
            x_now, residual = _finish_step(local_prob, z, op.laplacian, op.stalks)
        else
            sheaf    = build_time_expanded_tracking_sheaf(local_prob; state_only_initial=true)
            boundary = _assemble_boundary(local_prob, x_now, local_targets)
            z, N     = harmonic_extension(sheaf, boundary)
            z_arr    = Array(z)
            if cost != 0.0
                z_arr = solve_quadratic_on_basis(z_arr, N, build_control_cost_matrix(local_prob, sheaf.vertex_stalks, cost))
            end
            t == 0 && (nd = size(N, 2))
            x_now, residual = _finish_step(local_prob, z_arr, sheaf_laplacian_matrix(sheaf), sheaf.vertex_stalks)
        end

        for i in 1:prob.n_agents; agent_trajs[i][t+2, :] = x_now[i]; end
    end

    return ScenarioResult(label, collect(times), agent_trajs, target_trajs, nd, residual, y_col, z_col)
end

function animate_tracking_xy end

end # module MultiAgentTracking
