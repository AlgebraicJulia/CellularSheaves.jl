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
end
