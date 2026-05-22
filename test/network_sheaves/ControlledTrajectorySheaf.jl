using Test
using CellularSheaves
using LinearAlgebra
using BlockArrays

@testset "ControlledTrajectorySheaf" begin

    # -----------------------------------------------------------------------
    # Primary test case: scalar integrator ẋ = u, Ac = [0], Bc = [1]
    # With h = 1 this gives Ad = [1], Bd = [1].
    # -----------------------------------------------------------------------
    Ac_int = reshape([0.0], 1, 1)
    Bc_int = reshape([1.0], 1, 1)
    h_int  = 1.0
    F1     = EuclideanSheaf{Float64}(fill(1, 1))   # base sheaf: 1D stalk

    # -----------------------------------------------------------------------
    # 1. continuous_to_discrete_zoh returns the correct (Ad, Bd) for integrator
    # -----------------------------------------------------------------------
    @testset "ZOH discretization: scalar integrator" begin
        Ad, Bd = continuous_to_discrete_zoh(Ac_int, Bc_int, h_int)
        @test Ad ≈ [1.0;;]
        @test Bd ≈ [1.0;;]
    end

    # -----------------------------------------------------------------------
    # 2. ControlledTrajectorySheaf stores expected fields (k=3 for variety)
    # -----------------------------------------------------------------------
    @testset "constructor stores fields" begin
        k  = 3
        ts = ControlledTrajectorySheaf(F1, Ac_int, Bc_int, h_int, k)
        @test ts.k          == k
        @test ts.h          ≈  h_int
        @test ts.Ad         ≈  [1.0;;]
        @test ts.Bd         ≈  [1.0;;]
        @test ts.state_dim  == 1
        @test ts.control_dim == 1
    end

    # -----------------------------------------------------------------------
    # 3. vertex_stalks equals [1; fill(2, k); 1] for scalar integrator
    # -----------------------------------------------------------------------
    @testset "vertex stalks layout" begin
        k  = 4
        ts = ControlledTrajectorySheaf(F1, Ac_int, Bc_int, h_int, k)
        vs = vertex_stalks(ts.sheaf)
        @test vs == [1; fill(2, k); 1]
    end

    # -----------------------------------------------------------------------
    # 4. Every edge stalk has dimension 1 (= state dimension) for integrator
    # -----------------------------------------------------------------------
    @testset "all edge stalks have state dimension" begin
        k  = 3
        ts = ControlledTrajectorySheaf(F1, Ac_int, Bc_int, h_int, k)
        @test all(d -> d == 1, edge_stalk_dimensions(ts.sheaf))
    end

    # -----------------------------------------------------------------------
    # 5. Dummy initial edge restriction maps: [1] and [1 0]
    # -----------------------------------------------------------------------
    @testset "dummy initial edge restriction maps" begin
        k  = 3
        ts = ControlledTrajectorySheaf(F1, Ac_int, Bc_int, h_int, k)
        @test get_restriction_map(ts.sheaf, 1, 2) ≈ reshape([1.0], 1, 1)
        @test get_restriction_map(ts.sheaf, 2, 1) ≈ [1.0 0.0]
    end

    # -----------------------------------------------------------------------
    # 6. Dynamics edge restriction maps
    #    Non-terminal: [1 1] from left vertex, [1 0] from right vertex
    #    Final edge:   [1 1] from k+1-th vertex, [1] from k+2-th vertex
    # -----------------------------------------------------------------------
    @testset "dynamics edge restriction maps" begin
        k  = 3
        ts = ControlledTrajectorySheaf(F1, Ac_int, Bc_int, h_int, k)

        # Non-terminal dynamics edges: t in 1:k-1, edge (t+1, t+2)
        for t in 1:k-1
            @test get_restriction_map(ts.sheaf, t + 1, t + 2) ≈ [1.0 1.0]
            @test get_restriction_map(ts.sheaf, t + 2, t + 1) ≈ [1.0 0.0]
        end

        # Final edge: (k+1, k+2)
        @test get_restriction_map(ts.sheaf, k + 1, k + 2) ≈ [1.0 1.0]
        @test get_restriction_map(ts.sheaf, k + 2, k + 1) ≈ reshape([1.0], 1, 1)
    end

    # -----------------------------------------------------------------------
    # 7. feasible_control_trajectory_basis returns correct endpoint states
    #    in the public ordering z = (x₁, …, x_{k+1}, u₁, …, u_k)
    # -----------------------------------------------------------------------
    @testset "particular solution endpoint states preserved" begin
        k   = 3
        ts  = ControlledTrajectorySheaf(F1, Ac_int, Bc_int, h_int, k)
        x1  = [0.0]
        xk1 = [1.0]
        z_p, _ = feasible_control_trajectory_basis(ts, x1, xk1)

        n = ts.state_dim
        # First state
        @test z_p[1:n] ≈ x1
        # Last state (index (k+1)*n)
        @test z_p[k*n+1:(k+1)*n] ≈ xk1
    end

    # -----------------------------------------------------------------------
    # 8. Nullspace dimension is k-1 for scalar integrator with k steps
    #    (k+1 dofs for controls + state, minus 2 boundary states = k-1 free
    #     parameters when state constraints propagate)
    # For k=2 the scalar integrator has 1 free control degree of freedom.
    # -----------------------------------------------------------------------
    @testset "nullspace dimension for k=2" begin
        k   = 2
        ts  = ControlledTrajectorySheaf(F1, Ac_int, Bc_int, h_int, k)
        _, null_basis = feasible_control_trajectory_basis(ts, [0.0], [1.0])
        @test size(null_basis, 2) == 1
    end

    # -----------------------------------------------------------------------
    # 8b. Finite-horizon controllability helpers are part of the public API
    # -----------------------------------------------------------------------
    @testset "finite-horizon controllability helpers" begin
        k  = 3
        ts = ControlledTrajectorySheaf(F1, Ac_int, Bc_int, h_int, k)
        C  = finite_horizon_controllability(ts)

        @test size(C) == (1, 3)
        @test C ≈ [1.0 1.0 1.0]
        @test expected_feasible_dimension(ts) == 2
        @test_throws ArgumentError expected_feasible_dimension(ts; rtol=0.0)
    end

    # -----------------------------------------------------------------------
    # 9. Every null-basis column is a global section of the inner sheaf
    #    (after inverting the public→internal extraction, endpoints = 0)
    # -----------------------------------------------------------------------
    @testset "null-basis columns have zero endpoint states" begin
        k   = 3
        ts  = ControlledTrajectorySheaf(F1, Ac_int, Bc_int, h_int, k)
        _, null_basis = feasible_control_trajectory_basis(ts, [0.0], [Float64(k)])

        n = ts.state_dim
        for j in 1:size(null_basis, 2)
            col = null_basis[:, j]
            # Endpoint states must be zero for a null-basis column
            @test norm(col[1:n]) < 1e-10
            @test norm(col[k*n+1:(k+1)*n]) < 1e-10
        end
    end

    # -----------------------------------------------------------------------
    # 10. z_p + null_basis * α still hits the same endpoints for any α
    # -----------------------------------------------------------------------
    @testset "perturbed trajectory preserves endpoints" begin
        k   = 3
        ts  = ControlledTrajectorySheaf(F1, Ac_int, Bc_int, h_int, k)
        x1  = [0.0]
        xk1 = [1.0]
        z_p, null_basis = feasible_control_trajectory_basis(ts, x1, xk1)

        n = ts.state_dim
        for _ in 1:5
            α    = randn(size(null_basis, 2))
            z_per = z_p + null_basis * α
            @test z_per[1:n] ≈ x1
            @test z_per[k*n+1:(k+1)*n] ≈ xk1
        end
    end

    # -----------------------------------------------------------------------
    # 11. Wrong-sized Ac, Bc, x1, or xk1 throw ArgumentError
    # -----------------------------------------------------------------------
    @testset "argument validation" begin
        k  = 2
        ts = ControlledTrajectorySheaf(F1, Ac_int, Bc_int, h_int, k)

        # Non-square Ac
        @test_throws ArgumentError continuous_to_discrete_zoh(zeros(1, 2), Bc_int, h_int)
        # Bc row mismatch with Ac
        @test_throws ArgumentError continuous_to_discrete_zoh(Ac_int, zeros(2, 1), h_int)

        # Wrong-length x1
        @test_throws ArgumentError feasible_control_trajectory_basis(ts, [0.0, 0.0], [1.0])
        # Wrong-length xk1
        @test_throws ArgumentError feasible_control_trajectory_basis(ts, [0.0], [1.0, 0.0])
    end

    # -----------------------------------------------------------------------
    # 12. h <= 0 and k < 1 throw ArgumentError
    # -----------------------------------------------------------------------
    @testset "invalid h and k throw ArgumentError" begin
        @test_throws ArgumentError continuous_to_discrete_zoh(Ac_int, Bc_int, 0.0)
        @test_throws ArgumentError continuous_to_discrete_zoh(Ac_int, Bc_int, -1.0)
        @test_throws ArgumentError ControlledTrajectorySheaf(F1, Ac_int, Bc_int, h_int, 0)
        @test_throws ArgumentError ControlledTrajectorySheaf(F1, Ac_int, Bc_int, h_int, -1)
        @test_throws ArgumentError ControlledTrajectorySheaf(F1, Ac_int, Bc_int, 0.0, 2)
        @test_throws ArgumentError ControlledTrajectorySheaf(F1, Ac_int, Bc_int, -0.5, 2)
    end

    # -----------------------------------------------------------------------
    # Additional: 2D system for broader coverage
    # -----------------------------------------------------------------------
    @testset "2D double integrator" begin
        # ẍ = u → Ac = [0 1; 0 0], Bc = [0; 1]
        Ac2 = [0.0 1.0; 0.0 0.0]
        Bc2 = reshape([0.0, 1.0], 2, 1)
        F2  = EuclideanSheaf{Float64}(fill(2, 1))
        k   = 3
        h   = 0.1
        ts2 = ControlledTrajectorySheaf(F2, Ac2, Bc2, h, k)

        @test ts2.state_dim   == 2
        @test ts2.control_dim == 1
        @test vertex_stalks(ts2.sheaf) == [2; fill(3, k); 2]
        @test all(d -> d == 2, edge_stalk_dimensions(ts2.sheaf))

        x1_2  = [0.0, 0.0]
        xk1_2 = [1.0, 0.0]
        z_p2, _ = feasible_control_trajectory_basis(ts2, x1_2, xk1_2)

        n2 = ts2.state_dim
        @test z_p2[1:n2] ≈ x1_2
        @test z_p2[k*n2+1:(k+1)*n2] ≈ xk1_2
        @test size(finite_horizon_controllability(ts2)) == (2, 3)
        @test expected_feasible_dimension(ts2) == 1
    end

end
