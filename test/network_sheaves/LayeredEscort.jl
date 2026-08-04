using Test
using CellularSheaves
using CellularSheaves.AgentControllers
using CellularSheaves.ControlSheaves.Layered
using LinearAlgebra
using Graphs

@testset "Generic Layered Escort Architecture" begin
    # 1. Test Topology Specification
    rings = [
        RingSpec(1, 6, 0.3),
        RingSpec(2, 6, 0.36),
        RingSpec(3, 4, 0.25)
    ]
    supports = [
        SupportSpec(1, 2, 3),
        SupportSpec(2, 3, 3)
    ]

    spec = LayeredEscortSpec(rings, supports)
    @test spec.n_rings == 3
    @test spec.n_supports == 2
    @test spec.n_agents == 6 + 6 + 4 + 3 + 3  # 22 agents
    @test spec.n_targets == 3
    @test spec.total_nodes == 25

    # 2. Build Sheaf Topology and Homomorphism
    F = build_layered_escort_sheaf(spec)
    @test length(F.vertex_stalks) == 25

    f = build_layered_homomorphism(spec)
    @test f.n_target == 5 # 3 rings + 2 supports

    PfF = pushforward_sheaf(f, F)
    @test length(PfF.vertex_stalks) == 5

    # 3. Build Structured Fiber Bases
    bases = build_layered_fiber_bases(f, F, spec)
    @test length(bases.fiber_bases) == 5
    @test length(bases.target_subbases) == 3

    # 4. Test Coordinate Transformation (world_to_pf_stalk)
    p1 = [1.0, 0.0, 1.5, 1.0]
    q_pf_1 = world_to_pf_stalk(bases, 1, p1)
    @test length(q_pf_1) == size(bases.fiber_bases[1], 2)

    # Verify right-inverse reconstruction of Target 1
    B1_T1 = bases.target_subbases[1]
    p1_rec = B1_T1 * q_pf_1
    @test isapprox(p1_rec, p1, atol=1e-5)

    # 5. Test High-Level & Mid-Level Solvers
    p_targets = [
        [1.0, 0.0, 1.5, 1.0],
        [-1.0, 0.0, 1.5, 1.0],
        [0.0, 1.0, 1.5, 1.0]
    ]

    q_H = solve_high_level_harmonic(PfF, bases, p_targets)
    @test length(q_H) == 5

    qstar_agents = solve_mid_level_harmonic(q_H, bases)
    @test size(qstar_agents) == (22, 4)

    # 6. Test Direct Solve
    qstar_direct = solve_direct_harmonic(F, spec.target_nodes, p_targets)
    @test size(qstar_direct) == (22, 4)

    # 7. Test Simulation Execution
    target_trajs = [
        t -> [1.0 + 0.05*cos(t), 0.05*sin(t), 1.5, 1.0],
        t -> [-1.0 - 0.05*cos(t), -0.05*sin(t), 1.5, 1.0],
        t -> [0.0, 1.0 + 0.05*sin(t), 1.5, 1.0]
    ]
    v_trajs = [
        t -> [-0.05*sin(t), 0.05*cos(t), 0.0, 0.0],
        t -> [0.05*sin(t), -0.05*cos(t), 0.0, 0.0],
        t -> [0.0, 0.05*cos(t), 0.0, 0.0]
    ]
    a_trajs = [
        t -> [-0.05*cos(t), -0.05*sin(t), 0.0, 0.0],
        t -> [-0.05*cos(t), 0.05*sin(t), 0.0, 0.0],
        t -> [0.0, -0.05*sin(t), 0.0, 0.0]
    ]

    # Test standard simulation
    # Setup LQR gain for test
    dyn_test = QuadrotorDynamics()
    Ad, Bd = CellularSheaves.AgentControllers.discrete_matrices(dyn_test, 0.05)
    K_test = CellularSheaves.AgentControllers.solve_dare(Ad, Bd, Matrix{Float64}(I, 10, 10), Matrix{Float64}(I, 3, 3))

    prob = LayeredEscortProblem(
        spec, F, f, PfF, bases,
        HomogeneousDynamics(dyn_test, K_test),
        target_trajs, v_trajs, a_trajs,
        0.05, 10
    )

    res = run_layered_escort_simulation(prob)
    @test length(res.sim_data) == 10
    @test length(res.qstar_history) == 10
    @test length(res.target_history) == 10

    # Test feedforward simulation
    res_ff = run_layered_escort_simulation(prob; use_feedforward=true)
    @test length(res_ff.sim_data) == 10
end

@testset "Generalized Stalk Dimension Throws Exception" begin
    D = 6
    rings = [RingSpec(1, 4, 0.3)]
    supports = [SupportSpec(1, 1, 1)]
    
    spec = LayeredEscortSpec(rings, supports; D=D)
    @test_throws Exception F = build_layered_escort_sheaf(spec)
end
