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
                          sheaf_laplacian_matrix_direct,
                          ldlt_pseudoinverse_and_null, ldlt_pinv_solve

export TrackingEdge, TrackingProblem, BobbingTarget, ScenarioResult
export trajectory
export selector_matrix, state_projection_matrix
export agent_vertex, target_vertex
export build_time_expanded_tracking_sheaf
export generate_reference_trajectory, extract_state_trajectories, extract_target_trajectories, run_scenario
export animate_tracking_xy
export run_mpc_scenario
export TrackingExtension, tracking_extension_operator
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
    boundary::Union{Dict{Int,Vector{Float64}}, Dict{Int,Float64}},
    times::AbstractVector{<:Real};
    target_trajs::Union{Nothing,Vector{Vector{Vector{Float64}}}} = nothing,
    y_col::Int = 1,
    z_col::Int = 2,
    cost::Union{Number,Function}=1.0   # (agent, time) -> Q_control matrix or a number
)
    sheaf = build_time_expanded_tracking_sheaf(prob)
    z_harmonic, null_basis = harmonic_extension(sheaf, boundary)
    L  = sheaf_laplacian_matrix_direct(sheaf)
    z_harmonic_array = Array(z_harmonic)
    nd = size(null_basis, 2)

    if cost != 0.0
        Q = build_control_cost_matrix(prob, cost)
        z_opt = solve_quadratic_on_basis(z_harmonic_array, null_basis, Q)
    else
        z_opt = z_harmonic_array
    end
    residual = sqrt(max(0.0, dot(z_opt, L * z_opt)))

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
    _assemble_boundary(prob, sheaf, x_agents, target_window) -> Dict{Int,Float64}

Build the DOF-indexed boundary dictionary for `harmonic_extension`:
- Agent `i`'s state at the window's first timestep is pinned to `x_agents[i]`.
- Target `j`'s full stalk at each timestep is pinned to `target_window[j][t+1]`.

`sheaf` must be the sheaf built from `prob` (used for stalk offsets).
`target_window` must have `prob.k + 1` samples per target.
"""
function _assemble_boundary(
    prob::TrackingProblem,
    sheaf::EuclideanSheaf,
    x_agents::Vector{Vector{Float64}},
    target_window::Vector{Vector{Vector{Float64}}},
)
    offsets = [0; cumsum(sheaf.vertex_stalks)]
    pinned  = Dict{Int,Float64}()
    for i in 1:prob.n_agents
        va = agent_vertex(prob, i, 0)
        for c in 1:size(prob.Ad[i], 1)
            pinned[offsets[va] + c] = x_agents[i][c]
        end
    end
    for j in 1:prob.n_targets, ts in 0:prob.k
        vt  = target_vertex(prob, j, ts)
        val = target_window[j][ts + 1]
        for c in 1:length(val)
            pinned[offsets[vt] + c] = val[c]
        end
    end
    return pinned
end

"""
    _apply_first_control(window_prob, z_opt, stalks) -> Vector{Vector{Float64}}

Extract the control at the window's first timestep from `z_opt` and advance
each agent's dynamics: `x_{t+1} = Ad_i * x_t + Bd_i * u_t`.
The agent state `x_t` is read from the section (it was pinned when the window was solved).
"""
function _apply_first_control(
    window_prob::TrackingProblem,
    z_opt::AbstractVector,
    stalks::Vector{Int},
)
    z_block = BlockArray(z_opt, stalks)
    return map(1:window_prob.n_agents) do i
        stalk = Array(z_block[Block(agent_vertex(window_prob, i, 0))])
        nx_i  = size(window_prob.Ad[i], 1)
        nu_i  = size(window_prob.Bd[i], 2)
        x_t   = stalk[1:nx_i]
        u_t   = stalk[nx_i+1 : nx_i+nu_i]
        window_prob.Ad[i] * x_t + window_prob.Bd[i] * u_t
    end
end


# ---------------------------------------------------------------------------
# Cost-optimal harmonic extension operator (cached MPC core)
# ---------------------------------------------------------------------------

"""
    TrackingExtension

Precomputed operator for the cost-optimal harmonic extension of a fixed window sheaf.
Each MPC step reduces to `z[interior] = M * x_B` — a single dense matvec.
Construct with [`tracking_extension_operator`](@ref).

Fields: `M` (interior×boundary), `boundary`/`interior` (DOF index vectors),
`stalks`, `laplacian` (for energy), `null_dim`, `window`.
"""
struct TrackingExtension
    M         :: Matrix{Float64}
    boundary  :: Vector{Int}
    interior  :: Vector{Int}
    stalks    :: Vector{Int}
    laplacian :: SparseMatrixCSC{Float64,Int}
    null_dim  :: Int
    window    :: Int
end

"""
    tracking_extension_operator(window_prob; cost=1.0) -> TrackingExtension

Factorize the window Laplacian once and fold in the control-cost nullspace
projection to produce a dense operator `M` so each MPC step is `M * x_B`.
"""
function tracking_extension_operator(window_prob::TrackingProblem; cost::Union{Number,Function}=1.0)
    W       = window_prob.k
    sheaf   = build_time_expanded_tracking_sheaf(window_prob)
    L       = sheaf_laplacian_matrix_direct(sheaf)
    stalks  = sheaf.vertex_stalks
    offsets = [0; cumsum(stalks)]
    n       = size(L, 1)
    nx      = [size(window_prob.Ad[i], 1) for i in 1:window_prob.n_agents]

    boundary = Int[]
    for i in 1:window_prob.n_agents
        v = agent_vertex(window_prob, i, 0)
        for c in 1:nx[i]; push!(boundary, offsets[v] + c); end
    end
    for j in 1:window_prob.n_targets, t in 0:W
        v = target_vertex(window_prob, j, t)
        for c in 1:stalks[v]; push!(boundary, offsets[v] + c); end
    end
    sort!(boundary)
    interior = setdiff(1:n, boundary)

    L_II = sparse(L[interior, interior])
    L_IB = sparse(L[interior, boundary])
    F    = ldlt!(ChordalLDLt(L_II), RowMaximum(); check=false)
    tol  = sqrt(eps(Float64)) * max(1.0, maximum(i -> abs(F.D[i,i]), 1:size(F.D,1); init=0.0))
    _, N = ldlt_pseudoinverse_and_null(F, zeros(length(interior)); tol=tol)

    H = Matrix{Float64}(undef, length(interior), length(boundary))
    for col in axes(L_IB, 2)
        H[:, col] = ldlt_pinv_solve(F, -Vector(L_IB[:, col]); tol=tol)
    end

    if size(N, 2) > 0 && !(cost isa Number && iszero(cost))
        Q_II   = Matrix(build_control_cost_matrix(window_prob, cost)[interior, interior])
        NtQ_II = N' * Q_II
        Fq     = ldlt!(DenseLDLtPivoted(NtQ_II * N), RowMaximum(); check=false)
        rhs    = NtQ_II * H
        coeff  = similar(rhs)
        for col in axes(rhs, 2); coeff[:, col] = ldlt_pinv_solve(Fq, rhs[:, col]); end
        M = H - N * coeff
    else
        M = H
    end

    return TrackingExtension(M, boundary, interior, stalks, sparse(L), size(N, 2), W)
end

function (op::TrackingExtension)(x_B::AbstractVector)
    @argcheck length(x_B) == length(op.boundary)
    z = zeros(sum(op.stalks))
    z[op.boundary] = x_B
    z[op.interior] = op.M * x_B
    return z
end

_uniform_activation(ts, k::Int) = isempty(ts) || sort(unique(ts)) == collect(0:k)

"""
    run_mpc_scenario(label, prob, x0_agents, target_trajs, times;
                     window, y_col=1, z_col=2, cost=1.0, solver=:cached) -> ScenarioResult

Receding-horizon MPC: at each step `t` solves the sheaf QP on `[t, min(t+window,k)]`,
applies the first control, and advances agent dynamics.  `solver=:cached` reuses a
[`TrackingExtension`](@ref) for recurring window lengths; falls back to `:naive` on
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
    ops      = Dict{Int,TrackingExtension}()
    nd       = 0
    residual = 0.0

    for t in 0:k-1
        t_end         = min(t + window, k)
        local_prob    = window_problem(prob, t, t_end)
        local_targets = window_targets(target_trajs, t, t_end)

        if use_cached && get(window_counts, t_end - t, 0) >= 2
            op = get!(ops, t_end - t) do; tracking_extension_operator(local_prob; cost=cost); end
            t == 0 && (nd = op.null_dim)
            x_B = Vector{Float64}(undef, length(op.boundary))
            p = 1
            for i in 1:prob.n_agents, c in 1:nx[i]; x_B[p] = x_now[i][c]; p += 1; end
            for ts in 0:op.window, j in 1:prob.n_targets
                for c in eachindex(local_targets[j][ts+1]); x_B[p] = local_targets[j][ts+1][c]; p += 1; end
            end
            z        = op(x_B)
            residual = sqrt(max(0.0, dot(z, op.laplacian * z)))
            x_now    = _apply_first_control(local_prob, z, op.stalks)
        else
            sheaf    = build_time_expanded_tracking_sheaf(local_prob)
            boundary = _assemble_boundary(local_prob, sheaf, x_now, local_targets)
            z, N     = harmonic_extension(sheaf, boundary)
            z_arr    = Array(z)
            if cost != 0.0
                z_arr = solve_quadratic_on_basis(z_arr, N, build_control_cost_matrix(local_prob, cost))
            end
            L        = sheaf_laplacian_matrix_direct(sheaf)
            t == 0 && (nd = size(N, 2))
            residual = sqrt(max(0.0, dot(z_arr, L * z_arr)))
            x_now    = _apply_first_control(local_prob, z_arr, sheaf.vertex_stalks)
        end

        for i in 1:prob.n_agents; agent_trajs[i][t+2, :] = x_now[i]; end
    end

    return ScenarioResult(label, collect(times), agent_trajs, target_trajs, nd, residual, y_col, z_col)
end

function animate_tracking_xy end

end # module MultiAgentTracking
