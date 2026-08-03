using Test
using LinearAlgebra
using CellularSheaves
using CellularSheaves.ControlSheaves.Tikhonov
using CellularSheaves.ControlSheaves.AgentControllers
using CellularSheaves.ControlSheaves.DistributedLayeredControl

@testset "Feedforward Velocity Control Abstractions" begin
    @testset "JointTikhonovFilter" begin
        x0 = [1.0, 2.0]
        v0 = [0.1, -0.2]
        flt = JointTikhonovFilter(x0, v0; epsilon = 0.05)
        
        @test flt isa AbstractTikhonovFilter
        @test flt.x == x0
        @test flt.v == v0

        qstar = [1.5, 2.5]
        qstar_dot = [0.2, 0.0]
        x_ref, v_ref = tikhonov_step!(flt, qstar, qstar_dot, 0.01)
        
        @test length(x_ref) == 2
        @test length(v_ref) == 2
        @test x_ref[1] > x0[1] # Moving towards qstar
    end

    @testset "AgentState with Joint Velocity Tracking" begin
        dyn = PlanarQuadrotorDynamics()
        dt = 0.05
        Q = Matrix{Float64}(I, 6, 6) * 100.0
        R = Matrix{Float64}(I, 2, 2) * 0.01
        
        ctrl = LQRController(dyn, dt, Q, R)
        @test ctrl isa AbstractAgentController

        x0 = zeros(6)
        state = AgentState(x0, dyn, dt, ctrl.K, 0.02; use_velocity=true)
        @test state isa AbstractAgentState
        @test state.filter isa JointTikhonovFilter

        qstar = [1.0, 1.0]
        qstar_dot = [0.5, 0.0]
        x_act, x_ref = step_agent!(state, qstar, qstar_dot, dt)
        
        @test length(x_act) == 6
        @test length(x_ref) == 6
        @test norm(x_act) > 0.0 # Agent moved using LQR with velocity reference
    end
end
