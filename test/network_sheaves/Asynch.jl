using Test
using CellularSheaves
using Graphs
using LinearAlgebra
using Distributions: Normal

@testset "AsynchSheaves" begin

    # shared test sheaf: 5-vertex cycle, stalk dim 2 everywhere
    g = cycle_graph(5)
    sheaf = constant_sheaf(g, 2)
    n = sum(vertex_stalks(sheaf))
    x0 = randn(n)

    L = sheaf_laplacian_matrix(sheaf)
    K = opnorm(L, 2)
    γ = 1.0 / K

    # ------------------------------------------------------------------
    # Random matrix generators
    # ------------------------------------------------------------------
    @testset "random_psd" begin
        M = random_psd(4)
        @test size(M) == (4, 4)
        @test issymmetric(M) || isapprox(M, M', atol=1e-12)
        @test all(eigvals(M) .>= -1e-10)
    end

    @testset "random_pd" begin
        M = random_pd(4)
        @test size(M) == (4, 4)
        @test all(eigvals(M) .> 0)
    end

    @testset "random_semi_orthogonal_matrix" begin
        Q = random_semi_orthogonal_matrix(5, 3)
        @test size(Q) == (5, 3)
        @test Q' * Q ≈ I(3) atol=1e-10
    end

    @testset "matrix_weighted_edge_generator" begin
        gen = matrix_weighted_edge_generator(pd_prob=1.0)
        M = gen(3)
        @test size(M, 2) == 3
    end

    # ------------------------------------------------------------------
    # Synchronous baseline
    # ------------------------------------------------------------------
    @testset "compute_trajectory convergence" begin
        traj = compute_trajectory(sheaf, x0, γ; max_iters=5000, tol=1e-6)
        @test length(traj) >= 2
        f = energy_function(sheaf_laplacian_matrix(sheaf))
        energies = f.(traj)
        @test energies[end] < energies[1]
    end

    # ------------------------------------------------------------------
    # Asynch — probabilistic scheduling
    # ------------------------------------------------------------------
    @testset "compute_trajectory_asynch ProbabilisticModelParams" begin
        params = ProbabilisticModelParams(0.8, 0.8)
        traj = compute_trajectory_asynch(sheaf, x0, γ, params; max_iters=200)
        @test length(traj) >= 2
        @test length(traj[1]) == n
    end

    # ------------------------------------------------------------------
    # Asynch — deterministic period scheduling
    # ------------------------------------------------------------------
    @testset "compute_trajectory_asynch deterministic periods" begin
        n_agents = nv(underlying_graph(sheaf))
        update_periods = fill(3, n_agents)
        broadcast_periods = fill(5, n_agents)
        traj = compute_trajectory_asynch(sheaf, x0, γ, update_periods, broadcast_periods; max_iters=200)
        @test length(traj) >= 2
        @test length(traj[1]) == n
    end

    # ------------------------------------------------------------------
    # Asynch — mixture model scheduling
    # ------------------------------------------------------------------
    @testset "compute_trajectory_asynch MixtureModelParams" begin
        B = 20
        update_model = MixtureModelParams(
            [Normal(0.1 * B, 0.01 * B), Normal(0.5 * B, 0.05 * B)],
            [0.5, 0.5],
        )
        broadcast_model = MixtureModelParams(
            [Normal(0.2 * B, 0.02 * B), Normal(0.7 * B, 0.07 * B)],
            [0.5, 0.5],
        )
        traj = compute_trajectory_asynch(sheaf, x0, γ, update_model, broadcast_model; max_iters=200, B=B)
        @test length(traj) >= 2
        @test length(traj[1]) == n
    end

    # ------------------------------------------------------------------
    # Batch over initial conditions
    # ------------------------------------------------------------------
    @testset "compute_trajectory_asynch batch x0s" begin
        B = 10
        update_model = MixtureModelParams([Normal(2.0, 0.2), Normal(5.0, 0.5)], [0.5, 0.5])
        broadcast_model = MixtureModelParams([Normal(3.0, 0.3), Normal(7.0, 0.7)], [0.5, 0.5])
        x0s = [randn(n) for _ in 1:3]
        trajs = compute_trajectory_asynch(sheaf, x0s, γ, update_model, broadcast_model; max_iters=100, B=B)
        @test length(trajs) == 3
        @test all(t -> length(t) >= 2, trajs)
    end

    # ------------------------------------------------------------------
    # Batch over step sizes
    # ------------------------------------------------------------------
    @testset "compute_trajectory_asynch batch γs" begin
        B = 10
        update_model = MixtureModelParams([Normal(2.0, 0.2), Normal(5.0, 0.5)], [0.5, 0.5])
        broadcast_model = MixtureModelParams([Normal(3.0, 0.3), Normal(7.0, 0.7)], [0.5, 0.5])
        γs = [γ, 0.5 * γ, 0.1 * γ]
        trajs = compute_trajectory_asynch(sheaf, x0, γs, update_model, broadcast_model; max_iters=100, B=B)
        @test length(trajs) == 3
        @test all(t -> length(t) >= 2, trajs)
    end

    # ------------------------------------------------------------------
    # I/O round-trip
    # ------------------------------------------------------------------
    @testset "trajectory I/O round-trip" begin
        traj = compute_trajectory(sheaf, x0, γ; max_iters=10)
        losses = energy_function(sheaf_laplacian_matrix(sheaf)).(traj)
        mktempdir() do dir
            f = joinpath(dir, "losses.csv")
            save_trajectory(f, losses)
            loaded = load_trajectory(f)
            @test vec(loaded) ≈ losses atol=1e-12
        end
    end

end
