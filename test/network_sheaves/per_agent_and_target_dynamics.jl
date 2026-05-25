using Test
using CellularSheaves
using CellularSheaves.ControlSheaves.MultiAgentTracking
using LinearAlgebra

@testset "Heterogeneous per-agent and per-target dynamics" begin

    # ------------------------------------------------------------------
    # Heterogeneous agents: same stalk size (nx+nu = 3) but different A,B
    # ------------------------------------------------------------------
    Ad_a1 = [1.0 0.1; 0.0 1.0]          # 2-state double-integrator
    Bd_a1 = [0.0; 0.1;;]                 # 2×1

    Ad_a2 = [0.9 0.2; 0.0 0.9]          # 2-state, faster decay
    Bd_a2 = [0.05; 0.0;;]               # 2×1

    # ------------------------------------------------------------------
    # Heterogeneous targets: different stalk sizes
    # ------------------------------------------------------------------
    Ad_t1 = [1.0 0.05; 0.0 1.0]         # 2-state, stalk = 3
    Bd_t1 = [0.0; 0.0;;]                 # 2×1

    Ad_t2 = [0.95;;]                     # 1-state scalar, stalk = 1
    Bd_t2 = zeros(1, 0)                  # 1×0: no control input

    k = 3
    n_agents = 2
    n_targets = 2

    # Stalk sizes: agent 1 = 3, agent 2 = 3, target 1 = 3, target 2 = 1
    nx_a1 = size(Ad_a1, 1);  nu_a1 = size(Bd_a1, 2)   # 2, 1
    nx_a2 = size(Ad_a2, 1);  nu_a2 = size(Bd_a2, 2)   # 2, 1
    nx_t1 = size(Ad_t1, 1);  nu_t1 = size(Bd_t1, 2)   # 2, 1
    nx_t2 = size(Ad_t2, 1);  nu_t2 = size(Bd_t2, 2)   # 1, 0

    # Consensus restriction: project both state coordinates (3-column for stalk=3)
    R_cons = [1.0 0.0 0.0; 0.0 1.0 0.0]   # 2×3

    # Tracking edges: use consistent restriction matrix dimensions
    # Agent 1 (stalk 3) ↔ Target 1 (stalk 3): identity on full stalk
    Ra1 = Matrix{Float64}(I, 3, 3)
    Rt1 = Matrix{Float64}(I, 3, 3)

    # Agent 2 (stalk 3) ↔ Target 2 (stalk 1): select first state component
    Ra2 = [1.0 0.0 0.0]   # 1×3
    Rt2 = [1.0;;]          # 1×1

    te = [TrackingEdge(1, 1, Ra1, Rt1), TrackingEdge(2, 2, Ra2, Rt2)]

    prob = TrackingProblem(
        n_agents, n_targets, k,
        [Ad_a1, Ad_a2], [Bd_a1, Bd_a2],
        [Ad_t1, Ad_t2], [Bd_t1, Bd_t2],
        [(1, 2)], te, R_cons,
        collect(0:k), collect(0:k),
        true, 1.0, 1.0,
    )

    # ------------------------------------------------------------------
    # 1. Verify heterogeneous stalk sizes
    # ------------------------------------------------------------------
    @testset "heterogeneous stalk sizes" begin
        sheaf = build_time_expanded_tracking_sheaf(prob)

        # Agents (same stalk at every timestep)
        for t in 0:k
            @test sheaf.vertex_stalks[agent_vertex(prob, 1, t)] == nx_a1 + nu_a1   # 3
            @test sheaf.vertex_stalks[agent_vertex(prob, 2, t)] == nx_a2 + nu_a2   # 3
        end

        # Targets (different stalk sizes)
        for t in 0:k
            @test sheaf.vertex_stalks[target_vertex(prob, 1, t)] == nx_t1 + nu_t1  # 3
            @test sheaf.vertex_stalks[target_vertex(prob, 2, t)] == nx_t2 + nu_t2  # 1
        end

        # Total vertex count
        @test length(sheaf.vertex_stalks) == (k + 1) * (n_agents + n_targets)
    end

    # ------------------------------------------------------------------
    # 2. Solve with pinned targets and verify extraction helpers
    # ------------------------------------------------------------------
    @testset "solve and extract trajectories" begin
        # Pin target 1 as a constant stalk
        tgt1_stalk = [1.0, 0.5, 0.0]    # state + control for target 1 (stalk=3)
        # Pin target 2 as a constant scalar
        tgt2_stalk = [2.0]               # stalk=1

        boundary = Dict{Int,Vector{Float64}}()
        for t in 0:k
            boundary[target_vertex(prob, 1, t)] = tgt1_stalk
            boundary[target_vertex(prob, 2, t)] = tgt2_stalk
        end
        # Pin agent 1 at t=0
        boundary[agent_vertex(prob, 1, 0)] = [0.0, 0.0, 0.0]
        # Pin agent 2 at t=0
        boundary[agent_vertex(prob, 2, 0)] = [0.0, 0.0, 0.0]

        result = run_scenario("hetero-both", prob, boundary, 0.0:Float64(k))

        # Agent trajectories: one per agent, correct state dimension
        @test length(result.agent_trajs) == n_agents
        @test size(result.agent_trajs[1], 1) == k + 1
        @test size(result.agent_trajs[1], 2) == nx_a1   # 2
        @test size(result.agent_trajs[2], 1) == k + 1
        @test size(result.agent_trajs[2], 2) == nx_a2   # 2

        # Target trajectories auto-extracted: one per target, k+1 stalks each
        @test length(result.target_trajs) == n_targets
        @test length(result.target_trajs[1]) == k + 1
        @test length(result.target_trajs[2]) == k + 1

        # Target 1 is pinned so the harmonic extension must return the boundary value
        for t in 0:k
            @test result.target_trajs[1][t + 1] ≈ tgt1_stalk
        end

        # Target 2 is pinned to a scalar
        for t in 0:k
            @test result.target_trajs[2][t + 1] ≈ tgt2_stalk
        end

        # Pinned agent initial conditions are recovered
        @test result.agent_trajs[1][1, :] ≈ [0.0, 0.0]
        @test result.agent_trajs[2][1, :] ≈ [0.0, 0.0]

        # Residual and null dim are non-negative
        @test result.null_dim >= 0
        @test result.residual >= 0.0
    end

    # ------------------------------------------------------------------
    # 3. Explicit target_trajs override still works
    # ------------------------------------------------------------------
    @testset "explicit target_trajs override" begin
        tgt1_stalk = [1.0, 0.5, 0.0]
        tgt2_stalk = [2.0]

        boundary = Dict{Int,Vector{Float64}}()
        for t in 0:k
            boundary[target_vertex(prob, 1, t)] = tgt1_stalk
            boundary[target_vertex(prob, 2, t)] = tgt2_stalk
        end
        boundary[agent_vertex(prob, 1, 0)] = [0.0, 0.0, 0.0]
        boundary[agent_vertex(prob, 2, 0)] = [0.0, 0.0, 0.0]

        custom_tt = [fill([9.9, 9.9, 9.9], k + 1), fill([8.8], k + 1)]
        result = run_scenario("explicit-tt", prob, boundary, 0.0:Float64(k);
                              target_trajs = custom_tt)

        @test result.target_trajs[1][1] ≈ [9.9, 9.9, 9.9]
        @test result.target_trajs[2][1] ≈ [8.8]
    end

    # ------------------------------------------------------------------
    # 4. Validation: wrong dynamics vector lengths
    # ------------------------------------------------------------------
    @testset "validation — wrong vector lengths" begin
        # Too few agent Ad matrices
        bad1 = TrackingProblem(
            n_agents, n_targets, k,
            [Ad_a1], [Bd_a1, Bd_a2],
            [Ad_t1, Ad_t2], [Bd_t1, Bd_t2],
            Tuple{Int,Int}[], te, R_cons,
            Int[], Int[], false, 1.0, 1.0,
        )
        @test_throws Exception build_time_expanded_tracking_sheaf(bad1)

        # Too few target_Bd matrices
        bad2 = TrackingProblem(
            n_agents, n_targets, k,
            [Ad_a1, Ad_a2], [Bd_a1, Bd_a2],
            [Ad_t1, Ad_t2], [Bd_t1],
            Tuple{Int,Int}[], te, R_cons,
            Int[], Int[], false, 1.0, 1.0,
        )
        @test_throws Exception build_time_expanded_tracking_sheaf(bad2)
    end

    # ------------------------------------------------------------------
    # 5. Validation: non-square Ad matrix
    # ------------------------------------------------------------------
    @testset "validation — non-square Ad" begin
        bad_ad = [1.0 0.0 0.1; 0.0 1.0 0.0]   # 2×3, not square
        bad3 = TrackingProblem(
            n_agents, n_targets, k,
            [bad_ad, Ad_a2], [Bd_a1, Bd_a2],
            [Ad_t1, Ad_t2], [Bd_t1, Bd_t2],
            Tuple{Int,Int}[], te, R_cons,
            Int[], Int[], false, 1.0, 1.0,
        )
        @test_throws Exception build_time_expanded_tracking_sheaf(bad3)
    end

    # ------------------------------------------------------------------
    # 6. Regression: uniform dynamics (backward-compatibility)
    # ------------------------------------------------------------------
    @testset "regression — uniform dynamics" begin
        Ad_uni = [1.0 0.1; 0.0 1.0]
        Bd_uni = [0.0; 0.1;;]
        R_uni  = Matrix{Float64}(I, 3, 3)
        te_uni = [TrackingEdge(1, 1, R_uni, R_uni),
                  TrackingEdge(2, 2, R_uni, R_uni)]

        prob_uni = TrackingProblem(
            2, 2, k,
            [Ad_uni, Ad_uni], [Bd_uni, Bd_uni],
            [Ad_uni, Ad_uni], [Bd_uni, Bd_uni],
            [(1, 2)], te_uni, R_uni,
            collect(0:k), collect(0:k),
            true, 1.0, 1.0,
        )
        sheaf = build_time_expanded_tracking_sheaf(prob_uni)
        expected = (k + 1) * 4   # 2 agents + 2 targets
        @test length(sheaf.vertex_stalks) == expected
        @test all(==(3), sheaf.vertex_stalks)   # nx+nu = 2+1 = 3
    end

end
