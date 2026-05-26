using Test
using CellularSheaves
using CellularSheaves.ControlSheaves.MultiAgentTracking
using LinearAlgebra
using BlockArrays

@testset "MultiAgentTracking" begin

    # Minimal 2×2 system: nx=2, nu=1 (double-integrator-like).
    # Ad = I (2×2), Bd = [0; 1] (2×1).
    nx = 2; nu = 1
    Ad = Matrix{Float64}(I, nx, nx)
    Bd = [0.0; 1.0;;]   # nx × nu

    k = 2
    n_agents = 2
    n_targets = 2

    R_full = state_projection_matrix([1, 2], nx, nu)   # 2 × 3: projects both state coords
    R1     = state_projection_matrix([1],    nx, nu)   # 1 × 3: projects coord 1 only

    te_full = [TrackingEdge(1, 1, R_full, R_full), TrackingEdge(2, 2, R_full, R_full)]
    te_1    = [TrackingEdge(1, 1, R1, R1),          TrackingEdge(2, 2, R1, R1)]

    @testset "selector_matrix" begin
        S = selector_matrix([1, 3], 4)
        @test size(S) == (2, 4)
        @test S * [10.0, 20.0, 30.0, 40.0] ≈ [10.0, 30.0]
        @test_throws ArgumentError selector_matrix([0, 1], 4)
        @test_throws ArgumentError selector_matrix([1, 5], 4)
    end

    @testset "state_projection_matrix" begin
        P = state_projection_matrix([1], 2, 1)
        @test size(P) == (1, 3)
        @test P * [1.0, 2.0, 3.0] ≈ [1.0]
    end

    @testset "vertex indices" begin
        prob = TrackingProblem(
            n_agents, n_targets, k, [Ad, Ad], [Bd, Bd], [Ad, Ad], [Bd, Bd],
            [(1, 2)], te_full, R_full,
            collect(0:k), collect(0:k),
            true, 1.0, 1.0,
        )
        n_per_t = n_agents + n_targets
        @test agent_vertex(prob, 1, 0) == 1
        @test agent_vertex(prob, 2, 0) == 2
        @test target_vertex(prob, 1, 0) == 3
        @test target_vertex(prob, 2, 0) == 4
        # At t=1 the block is shifted by n_per_t
        @test agent_vertex(prob, 1, 1) == 1 + n_per_t
        @test target_vertex(prob, 1, 1) == 3 + n_per_t
    end

    @testset "build_time_expanded_tracking_sheaf — vertex count" begin
        prob = TrackingProblem(
            n_agents, n_targets, k, [Ad, Ad], [Bd, Bd], [Ad, Ad], [Bd, Bd],
            [(1, 2)], te_full, R_full,
            collect(0:k), collect(0:k),
            true, 1.0, 1.0,
        )
        sheaf = build_time_expanded_tracking_sheaf(prob)
        expected_verts = (k + 1) * (n_agents + n_targets)
        @test length(sheaf.vertex_stalks) == expected_verts
        @test all(==(nx + nu), sheaf.vertex_stalks)
    end

    @testset "build_time_expanded_tracking_sheaf — @argcheck validation" begin
        # Negative consensus weight
        bad_prob = TrackingProblem(
            n_agents, n_targets, k, [Ad, Ad], [Bd, Bd], [Ad, Ad], [Bd, Bd],
            [(1, 2)], te_full, R_full,
            collect(0:k), collect(0:k),
            false, -1.0, 1.0,
        )
        @test_throws Exception build_time_expanded_tracking_sheaf(bad_prob)

        # Negative tracking weight
        bad_prob2 = TrackingProblem(
            n_agents, n_targets, k, [Ad, Ad], [Bd, Bd], [Ad, Ad], [Bd, Bd],
            [(1, 2)], te_full, R_full,
            collect(0:k), collect(0:k),
            false, 1.0, -0.5,
        )
        @test_throws Exception build_time_expanded_tracking_sheaf(bad_prob2)

        # Out-of-range consensus timestep
        bad_prob3 = TrackingProblem(
            n_agents, n_targets, k, [Ad, Ad], [Bd, Bd], [Ad, Ad], [Bd, Bd],
            [(1, 2)], te_full, R_full,
            [k + 1], collect(0:k),
            false, 1.0, 1.0,
        )
        @test_throws Exception build_time_expanded_tracking_sheaf(bad_prob3)

        # Out-of-range agent index in agent_edges
        bad_prob4 = TrackingProblem(
            n_agents, n_targets, k, [Ad, Ad], [Bd, Bd], [Ad, Ad], [Bd, Bd],
            [(1, n_agents + 1)], te_full, R_full,
            [k], [k],
            false, 1.0, 1.0,
        )
        @test_throws Exception build_time_expanded_tracking_sheaf(bad_prob4)

        # Wrong length for Ad vector
        bad_prob5 = TrackingProblem(
            n_agents, n_targets, k, [Ad], [Bd, Bd], [Ad, Ad], [Bd, Bd],
            [(1, 2)], te_full, R_full,
            [k], [k],
            false, 1.0, 1.0,
        )
        @test_throws Exception build_time_expanded_tracking_sheaf(bad_prob5)

        # Wrong length for target_Ad vector
        bad_prob6 = TrackingProblem(
            n_agents, n_targets, k, [Ad, Ad], [Bd, Bd], [Ad], [Bd, Bd],
            [(1, 2)], te_full, R_full,
            [k], [k],
            false, 1.0, 1.0,
        )
        @test_throws Exception build_time_expanded_tracking_sheaf(bad_prob6)
    end

    @testset "run_scenario — zero residual when boundary is consistent" begin
        # A single agent, k=1, boundary pinning both endpoints exactly.
        # With consistent boundary data the harmonic extension is exact
        # and the Laplacian energy should be 0.
        prob_1a = TrackingProblem(
            1, 1, k, [Ad], [Bd],
            [Ad], [Bd],
            Tuple{Int,Int}[], [TrackingEdge(1, 1, R_full, R_full)],
            R_full,
            Int[], Int[],   # no consensus, no tracking edges
            false, 1.0, 1.0,
        )
        x0 = [0.0, 0.0]
        xk = [1.0, 0.0]
        bnd = Dict{Int,Vector{Float64}}()
        bnd[agent_vertex(prob_1a, 1, 0)] = vcat(x0, zeros(nu))
        bnd[agent_vertex(prob_1a, 1, k)] = vcat(xk, zeros(nu))
        # Pin target so it has no isolated vertices
        for t in 0:k
            bnd[target_vertex(prob_1a, 1, t)] = vcat(x0, zeros(nu))
        end
        times = Float64.(0:k)
        result = run_scenario("1-agent test", prob_1a, bnd, times)
        @test result.label == "1-agent test"
        @test result.null_dim >= 0
        @test result.residual >= 0.0
        # Boundary values are recovered exactly
        @test result.agent_trajs[1][1, :] ≈ x0
        @test result.agent_trajs[1][k + 1, :] ≈ xk
        # target_trajs auto-extracted
        @test length(result.target_trajs) == 1
        @test length(result.target_trajs[1]) == k + 1
    end

    @testset "ScenarioResult show" begin
        # Verify that the show method doesn't error and contains expected fields.
        result_dummy = ScenarioResult(
            "test", [0.0, 1.0],
            [zeros(2, nx)], Vector{Vector{Vector{Float64}}}(),
            3, 1.23, 1, 2,
        )
        plain = sprint(show, MIME("text/plain"), result_dummy)
        @test occursin("null dim", plain)
        @test occursin("test", plain)
        brief = sprint(show, result_dummy)
        @test occursin("null dim = 3", brief)
        @test occursin("1.23", brief)
    end

    @testset "BobbingTarget trajectory" begin
        bt = BobbingTarget(2.0, 5.0, 1.0, π)
        traj = trajectory(bt, 0:3, 0.5, 4, 2, 2, 3, 4)
        @test length(traj) == 4
        @test length(traj[1]) == 6   # nx + nu = 4 + 2
        # y component should be fixed
        @test all(v -> v[2] ≈ 2.0, traj)
        # z at t=0 is center (sin(0)=0)
        @test traj[1][3] ≈ 5.0
        # z at t=0.5 s: sin(π*0.5) = 1 → 5+1 = 6
        @test traj[2][3] ≈ 5.0 + 1.0 * sin(π * 0.5)
    end

end
