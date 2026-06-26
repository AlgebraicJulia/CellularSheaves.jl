using Test
using SparseArrays
using LinearAlgebra
using CellularSheaves.ControlSheaves.MultiAgentTracking: TrackingProblem, TrackingEdge, agent_vertex
using CellularSheaves.ControlSheaves.MultiAgentTracking.QuadraticCosts: build_control_cost_matrix, solve_quadratic_on_basis

# Helper to construct a minimal TrackingProblem with 1 agent, no targets
function make_simple_prob(; k=2, nx=2, nu=2, control_weight=1.0)
    Ad = [Matrix{Float64}(I, nx, nx)]
    Bd = [Matrix{Float64}(I, nx, nu)]   # simple control mapping
    prob = TrackingProblem(
        1,                # n_agents
        0,                # n_targets
        k,
        Ad,
        Bd,
        Vector{Matrix{Float64}}(),   # target_Ad
        Vector{Matrix{Float64}}(),   # target_Bd
        Tuple{Int,Int}[],            # agent_edges
        TrackingEdge[],              # tracking_edges (empty)
        Matrix{Float64}(I, nx+nu, nx+nu),   # consensus_restriction (unused)
        Int[],                      # consensus_timesteps
        Int[],                      # tracking_timesteps
        false,                      # include_target_dynamics
        1.0,                        # consensus_weight
        1.0                         # tracking_weight
    )
    return prob
end

@testset "QuadraticCost utilities" begin
    prob = make_simple_prob(k=2, nx=2, nu=2, control_weight=2.0)
    # Stalk layout of the (default) time-expanded sheaf: nx+nu per agent vertex.
    stalks = [size(prob.Ad[1],1) + size(prob.Bd[1],2) for _ in 0:prob.k]
    # Build Q with simple 2-norm weight (2.0 * I)
    Q = build_control_cost_matrix(prob, stalks, 2.0)
    total_dim = sum(stalks)
    @test size(Q) == (total_dim, total_dim)
    # Expected control indices for each time step
    control_idxs = Int[]
    for t in 0:prob.k
        v = agent_vertex(prob, 1, t)
        nx = size(prob.Ad[1],1)
        nu = size(prob.Bd[1],2)
        offset = (v-1)*(nx+nu)   # zero‑based offset as in original code
        append!(control_idxs, (offset+nx+1):(offset+nx+nu))
    end
    # Diagonal of Q should be 2.0 at control indices, zero elsewhere
    diagQ = diag(Q)
    @test all(diagQ[control_idxs] .== 2.0)
    @test all(diagQ[setdiff(1:total_dim, control_idxs)] .== 0.0)

    # Test solving on nullspace: use a trivial nullspace (vector of ones)
    N = ones(total_dim, 1)  # null_basis spanning constant offset
    z_harm = zeros(total_dim)
    z_opt = solve_quadratic_on_basis(z_harm, N, Q)
    @test length(z_opt) == total_dim
    @test norm(z_opt) ≈ 0.0
end
