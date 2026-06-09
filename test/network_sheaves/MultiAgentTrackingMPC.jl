using Test
using CellularSheaves
using CellularSheaves.ControlSheaves.MultiAgentTracking
using LinearAlgebra
using SparseArrays
using BlockArrays

@testset "MultiAgentTrackingMPC" begin

    # Fixture mirrors MultiAgentTracking.jl: same dynamics, same projection
    # matrices, same edge topology — horizon extended to k=8 for MPC tests.
    nx = 2; nu = 1
    Ad = Matrix{Float64}(I, nx, nx)
    Bd = [0.0; 1.0;;]   # nx × nu

    k         = 8
    n_agents  = 2
    n_targets = 2

    R_full = state_projection_matrix([1, 2], nx, nu)
    R1     = state_projection_matrix([1],    nx, nu)

    te_full = [TrackingEdge(1, 1, R_full, R_full), TrackingEdge(2, 2, R_full, R_full)]

    prob = TrackingProblem(
        n_agents, n_targets, k,
        [Ad, Ad], [Bd, Bd], [Ad, Ad], [Bd, Bd],
        [(1, 2)], te_full, R_full,
        collect(0:k), collect(0:k),
        false, 0.5, 1.0,
    )

    x0_agents    = [[0.0, 0.0], [1.0, 0.0]]
    target_trajs = [
        [Float64[0.0, Float64(t),     0.0] for t in 0:k],
        [Float64[1.0, Float64(k - t), 0.0] for t in 0:k],
    ]
    times = Float64.(0:k)

    @testset "run_mpc_scenario — returns ScenarioResult with correct dimensions" begin
        result = run_mpc_scenario("t1", prob, x0_agents, target_trajs, times; window=k)
        @test result isa ScenarioResult
        @test result.label == "t1"
        @test length(result.times) == k + 1
        @test length(result.agent_trajs) == n_agents
        @test all(i -> size(result.agent_trajs[i]) == (k + 1, nx), 1:n_agents)
        @test length(result.target_trajs) == n_targets
        @test result.residual >= 0.0
        # null_dim is the interior-nullspace dimension of the planning window
        # (free coordination DOFs after pinning agent initials + targets).
        @test result.null_dim == 2
    end

    @testset "run_mpc_scenario — initial states equal x0_agents" begin
        result = run_mpc_scenario("t2", prob, x0_agents, target_trajs, times; window=k)
        for i in 1:n_agents
            @test result.agent_trajs[i][1, :] ≈ x0_agents[i]
        end
    end

    @testset "run_mpc_scenario — all states remain finite" begin
        result = run_mpc_scenario("t3", prob, x0_agents, target_trajs, times; window=k, cost=1.0)
        for i in 1:n_agents
            @test all(isfinite, result.agent_trajs[i])
        end
    end

    @testset "run_mpc_scenario — window=1 greedy MPC runs without error" begin
        result = run_mpc_scenario("t4", prob, x0_agents, target_trajs, times; window=1)
        @test result isa ScenarioResult
        @test length(result.times) == k + 1
        @test all(i -> size(result.agent_trajs[i]) == (k + 1, nx), 1:n_agents)
    end

    @testset "run_mpc_scenario — window=k consistent with open-loop on constant target" begin
        # When the target is constant at x0 the optimal control is u=0 for all t.
        # Open-loop forces u_0=0 via the whole-stalk pin, which is also optimal here,
        # so both methods find the same section and must agree exactly.
        prob_const = TrackingProblem(
            1, 1, k, [Ad], [Bd], [Ad], [Bd],
            Tuple{Int,Int}[], [TrackingEdge(1, 1, R_full, R_full)],
            R_full, Int[], collect(0:k),
            false, 0.0, 1.0,
        )
        x0        = [0.0, 0.0]
        const_tgt = [Float64[0.0, 0.0, 0.0] for _ in 0:k]

        bnd = Dict{Int,Vector{Float64}}()
        bnd[agent_vertex(prob_const, 1, 0)] = vcat(x0, zeros(nu))
        for t in 0:k
            bnd[target_vertex(prob_const, 1, t)] = const_tgt[t + 1]
        end
        result_ol  = run_scenario("ol",  prob_const, bnd,       times; cost=1.0)
        result_mpc = run_mpc_scenario("mpc", prob_const, [x0], [const_tgt], times; window=k, cost=1.0)

        @test result_ol.agent_trajs[1] ≈ result_mpc.agent_trajs[1] atol=1e-10
        @test result_ol.residual       ≈ result_mpc.residual        atol=1e-10
    end

    @testset "run_mpc_scenario — bad arguments throw ArgumentError" begin
        @test_throws ArgumentError run_mpc_scenario(
            "bad", prob, x0_agents, target_trajs, times; window=0)
        @test_throws ArgumentError run_mpc_scenario(
            "bad", prob, [x0_agents[1]], target_trajs, times; window=k)
        @test_throws ArgumentError run_mpc_scenario(
            "bad", prob, x0_agents, [target_trajs[1]], times; window=k)
        @test_throws ArgumentError run_mpc_scenario(
            "bad", prob, x0_agents, target_trajs, Float64.(0:k+1); window=k)
        @test_throws ArgumentError run_mpc_scenario(
            "bad", prob, x0_agents, target_trajs, times; window=k, solver=:bogus)
    end

    @testset "run_mpc_scenario — cached and naive agree (homogeneous problem)" begin
        for W in (1, 3, k)
            res_c = run_mpc_scenario("c", prob, x0_agents, target_trajs, times;
                                     window=W, cost=1.0, solver=:cached)
            res_n = run_mpc_scenario("n", prob, x0_agents, target_trajs, times;
                                     window=W, cost=1.0, solver=:naive)
            @test res_c.null_dim == res_n.null_dim
            for i in 1:n_agents
                @test res_c.agent_trajs[i] ≈ res_n.agent_trajs[i] atol=1e-6
            end
        end
    end

    @testset "tracking_extension_operator — pinning and harmonic stationarity" begin
        op = tracking_extension_operator(prob; cost=1.0)
        n_tot = sum(op.stalks)
        @test length(op.boundary) + length(op.interior) == n_tot
        @test size(op.M) == (length(op.interior), length(op.boundary))
        @test op.window == k
        @test op.null_dim == 2

        x_B = randn(length(op.boundary))
        z = op(x_B)
        @test length(z) == n_tot
        # Boundary is pinned exactly.
        @test z[op.boundary] ≈ x_B
        # Interior is the harmonic extension: the interior energy gradient
        # g = L_II z_I + L_IB z_B lies in ker(L_II), so L_II * g ≈ 0.
        L_II = op.laplacian[op.interior, op.interior]
        L_IB = op.laplacian[op.interior, op.boundary]
        g = L_II * z[op.interior] + L_IB * x_B
        @test norm(L_II * g) < 1e-6
    end

    @testset "tracking_extension_operator — Laplacian blocks stay sparse" begin
        op = tracking_extension_operator(prob; cost=1.0)
        @test op.laplacian isa SparseMatrixCSC
        # Indexing the sparse Laplacian with index vectors must keep the blocks
        # sparse (the operator build relies on this instead of converting).
        @test op.laplacian[op.interior, op.interior] isa SparseMatrixCSC
        @test op.laplacian[op.interior, op.boundary] isa SparseMatrixCSC
    end

    @testset "run_mpc_scenario — solver=:naive default-equivalent on inhomogeneous problem" begin
        # Tracking active only at the terminal step → not time-homogeneous, so the
        # cached default must transparently fall back to the naive path.
        prob_inhom = TrackingProblem(
            n_agents, n_targets, k,
            [Ad, Ad], [Bd, Bd], [Ad, Ad], [Bd, Bd],
            [(1, 2)], te_full, R_full,
            Int[], [k],
            false, 0.5, 1.0,
        )
        res_default = run_mpc_scenario("d", prob_inhom, x0_agents, target_trajs, times;
                                       window=4, cost=1.0)
        res_naive   = run_mpc_scenario("n", prob_inhom, x0_agents, target_trajs, times;
                                       window=4, cost=1.0, solver=:naive)
        for i in 1:n_agents
            @test res_default.agent_trajs[i] ≈ res_naive.agent_trajs[i] atol=1e-10
        end
    end

end
