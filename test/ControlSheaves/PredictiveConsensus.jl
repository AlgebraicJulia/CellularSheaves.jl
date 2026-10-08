using Test
using CellularSheaves
using CellularSheaves.ControlSheaves.CoordinationBenchmarks: coordination_scenario
using CellularSheaves.ControlSheaves.PredictiveConsensus
using CellularSheaves.NetworkSheaves.SchwarzMethods: MulticolorSweep, ParallelSweep
using CellularSheaves.ControlSheaves.PredictiveConsensus: init, solve, solve!
using LinearAlgebra
using Random

@testset "PredictiveConsensus" begin
    scenario = coordination_scenario(:grid; size_parameter = 3)       # 9 agents, 2 targets
    plant = PlanarDoubleIntegrator(0.2)
    targets = [0.0 4.0; 0.0 3.0]
    rng = Xoshiro(3)
    x0 = [3 .* randn(rng, 2, 9); randn(rng, 2, 9)]
    lq = ConsensusLQ(scenario, plant, 25; targets)
    prob = ConsensusProblem(lq, x0)

    @testset "plant and formation" begin
        @test plant.A ≈ exp([zeros(2, 2) I; zeros(2, 4)] .* 0.2)
        @test plant.B ≈ [0.02 * I(2); 0.2 * I(2)]
        qstar = formation(lq)
        @test size(qstar) == (2, 9)
        # q⋆ is the harmonic extension: zero sheaf disagreement.
        @test norm(scenario.H * vec(qstar) - scenario.Bmat * vec(targets)) < 1e-10
        @test_throws ArgumentError ConsensusLQ(coordination_scenario(:grid; size_parameter = 3, dim = 3), plant, 5;
            targets = zeros(3, 2))
        @test_throws ArgumentError ConsensusProblem(lq, zeros(4, 8))
    end

    riccati = solve(prob, JointRiccati())
    central = solve(prob, CentralizedQP())

    @testset "direct solvers agree" begin
        @test central.states ≈ riccati.states atol = 1e-5
        @test central.controls ≈ riccati.controls atol = 1e-5
        @test central.cost ≈ riccati.cost rtol = 1e-7
        @test riccati.states[:, 1, :] ≈ x0
        # The plan satisfies the dynamics.
        @test riccati.states[:, 2, 4] ≈ plant.A * x0[:, 4] + plant.B * riccati.controls[:, 1, 4]
        # Zero input is feasible, so the optimum is no worse.
        coast, idle = coasting_plan(prob)
        @test riccati.cost < team_cost(lq, coast, idle)
    end

    @testset "multicolor sweeps reach the team optimum" begin
        plan = solve(prob, PredictedTrajectorySweeps(tol = 1e-9))
        @test plan.converged
        @test plan.cost ≈ riccati.cost rtol = 1e-7
        @test plan.states ≈ riccati.states atol = 1e-4
        # Block coordinate descent on a convex objective: monotone.
        @test all(diff(plan.costs) .<= 1e-9 * riccati.cost)
    end

    @testset "damped parallel sweeps converge" begin
        plan = solve(prob, PredictedTrajectorySweeps(sweep = ParallelSweep(), tol = 1e-9, maxiter = 2000))
        @test plan.converged
        @test plan.cost ≈ riccati.cost rtol = 1e-6
        ws = init(prob, PredictedTrajectorySweeps(sweep = ParallelSweep()))
        @test ws.damping == 1 / length(ws.colors)
    end

    @testset "control bound (second-order cone)" begin
        bounded = ConsensusLQ(scenario, plant, 25; targets, control_bound = 1.0)
        bprob = ConsensusProblem(bounded, x0)
        @test_throws ArgumentError solve(bprob, JointRiccati())
        ref = solve(bprob, CentralizedQP())
        speeds = [norm(ref.controls[:, t, i]) for t in 1:25, i in 1:9]
        @test maximum(speeds) <= 1.0 + 1e-6
        @test maximum(speeds) >= 1.0 - 1e-3          # the bound is active
        @test ref.cost > riccati.cost
        plan = solve(bprob, PredictedTrajectorySweeps(tol = 1e-8))
        @test plan.converged
        @test plan.cost ≈ ref.cost rtol = 1e-5
        @test maximum(norm(plan.controls[:, t, i]) for t in 1:25, i in 1:9) <= 1.0 + 1e-6
    end

    @testset "closed loop" begin
        steps = 60
        lqr = rollout(lq, RiccatiFeedback(lq), x0, steps)
        diffusion = rollout(lq, SecondOrderDiffusion(lq, 1.0, 1.5), x0, steps)
        rhc = rollout(lq, RecedingHorizon(prob, PredictedTrajectorySweeps(); sweeps = 50), x0, steps)
        @test size(lqr.states) == (4, steps + 1, 9)
        @test lqr.cost < diffusion.cost
        # A converged long-horizon receding plan is close to the LQR optimum.
        @test rhc.cost ≈ lqr.cost rtol = 0.05
        # All three settle onto the formation.
        for run in (lqr, rhc)
            @test norm(run.states[1:2, end, :] - formation(lq)) < 0.1 * norm(x0[1:2, :] - formation(lq))
        end
        # A single sweep per step still stabilizes.
        one_sweep = rollout(lq, RecedingHorizon(prob, PredictedTrajectorySweeps(); sweeps = 1), x0, steps)
        @test one_sweep.cost < rollout(lq, SecondOrderDiffusion(lq, 0.0, 0.0), x0, steps).cost
    end
end
