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
export generate_reference_trajectory, extract_state_trajectories, extract_target_trajectories, run_scenario
export animate_tracking_xy

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
    L  = sheaf_laplacian_matrix_direct(sheaf)
    z_harmonic_array = Array(z_harmonic)
    nd = size(null_basis, 2)

    # -------------------------------------------------------------------
    # Optional quadratic control cost
    # -------------------------------------------------------------------
    # Build the global quadratic cost matrix (only on control components) and
    # solve the reduced problem on the nullspace of the harmonic‑extension
    # constraints.  The heavy lifting is delegated to `QuadraticCost` utilities.
    if cost != 0.0
        Q = build_control_cost_matrix(prob, cost)
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
    animate_tracking_xy(result; kwargs...)

Create a 2D animation for a `ScenarioResult` with agent and target trajectories.
This method is implemented by the optional plotting extension when `Plots` is
loaded.
"""
function animate_tracking_xy end

end # module MultiAgentTracking
