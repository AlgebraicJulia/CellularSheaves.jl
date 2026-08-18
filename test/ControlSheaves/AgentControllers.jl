using Test
using CellularSheaves.AgentControllers
using CellularSheaves.TrajectorySheaves: continuous_to_discrete_zoh
using LinearAlgebra

@testset "AgentControllers" begin
    @testset "QuadrotorDynamics" begin
        dyn = QuadrotorDynamics(m=0.6, Ixx=0.02)
        @test dyn.m == 0.6
        @test dyn.Ixx == 0.02
        
        Ac, Bc = AgentControllers.continuous_matrices(dyn)
        @test size(Ac) == (10, 10)
        @test size(Bc) == (10, 3)
        @test Bc[8, 1] == 1.0 / 0.6
        
        dt = 0.05
        Ad, Bd = AgentControllers.discrete_matrices(dyn, dt)
        @test size(Ad) == (10, 10)
        @test size(Bd) == (10, 3)
    end

    @testset "DARE and LQRController" begin
        dyn = QuadrotorDynamics()
        dt = 0.05
        
        Q_diag = [150.0, 150.0, 150.0, 50.0, 50.0, 30.0, 30.0, 30.0, 1.0, 1.0]
        Q = Matrix(Diagonal(Q_diag))
        R = Matrix(Diagonal([0.01, 0.01, 0.01]))
        
        lqr = LQRController(dyn, dt, Q, R)
        @test size(lqr.K) == (3, 10)
        
        # Verify stabilizing properties
        Ad, Bd = AgentControllers.discrete_matrices(dyn, dt)
        A_cl = Ad - Bd * lqr.K
        rho = maximum(abs.(eigvals(A_cl)))
        @test rho < 1.0
    end

    @testset "position_indices / state_dim / initial_state" begin
        dyn_q = QuadrotorDynamics()
        dyn_p = PlanarQuadrotorDynamics()

        @test position_indices(dyn_q) == 1:3
        @test velocity_indices(dyn_q) == 6:8
        @test state_dim(dyn_q) == 10

        @test position_indices(dyn_p) == 1:2
        @test velocity_indices(dyn_p) == 4:5
        @test state_dim(dyn_p) == 6

        @testset "QuadrotorDynamics" begin
            # position only: zero velocity/acceleration (steady hover)
            x0 = initial_state(dyn_q, [1.0, 2.0, 3.0])
            @test length(x0) == 10
            @test x0[1:3] == [1.0, 2.0, 3.0]
            @test all(x0[6:8] .== 0.0)
            @test all(x0[4:5] .== 0.0) # level attitude, no trim

            # position + velocity, zero acceleration
            x1 = initial_state(dyn_q, [0.0, 0.0, 1.5], [0.5, -0.5, 0.0])
            @test x1[6:8] == [0.5, -0.5, 0.0]
            @test all(x1[4:5] .== 0.0)

            # position + velocity + nonzero lateral acceleration -> trimmed attitude
            g = dyn_q.g
            ax, ay = 1.0, 2.0
            x2 = initial_state(dyn_q, [0.0, 0.0, 1.5], [0.0, 0.0, 0.0], [ax, ay, 0.0])
            @test x2[4] ≈ -ay / g
            @test x2[5] ≈ ax / g

            @test_throws Exception initial_state(dyn_q, [1.0, 2.0]) # wrong position dim
        end

        @testset "PlanarQuadrotorDynamics" begin
            x0 = initial_state(dyn_p, [1.0, 1.5])
            @test length(x0) == 6
            @test x0[1:2] == [1.0, 1.5]
            @test all(x0[4:5] .== 0.0)
            @test x0[3] == 0.0

            g = dyn_p.g
            ay = 3.0
            x1 = initial_state(dyn_p, [0.0, 1.5], [0.0, 0.0], [ay])
            @test x1[3] ≈ -ay / g

            @test_throws Exception initial_state(dyn_p, [1.0, 2.0, 3.0]) # wrong position dim
        end
    end
end
