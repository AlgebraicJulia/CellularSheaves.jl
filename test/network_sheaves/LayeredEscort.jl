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
    @test spec.D == 4
    @test spec.affine == true

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
        0.05, 10, nothing
    )

    res = run_layered_escort_simulation(prob)
    @test length(res.sim_data) == 10
    @test length(res.qstar_history) == 10
    @test length(res.target_history) == 10

    # By default (no initial_positions supplied) agents start from a simple
    # "airstrip" line, distinct from the target formation, so that convergence
    # into formation is visible: spread along x with fixed spacing, common
    # hover altitude on the last position coordinate (position_indices(dyn_test) == 1:3).
    for i in 1:spec.n_agents
        @test res.sim_data[1][i][1] ≈ (i - 1) * 0.5 atol=0.05
        @test res.sim_data[1][i][3] ≈ 1.5 atol=0.05
    end

    # Test feedforward simulation
    res_ff = run_layered_escort_simulation(prob; use_feedforward=true)
    @test length(res_ff.sim_data) == 10

    # Explicit initial_positions overrides the default airstrip layout.
    custom_positions = [zeros(3) for _ in 1:spec.n_agents]
    for i in 1:spec.n_agents
        custom_positions[i][3] = 3.0 # a different common altitude
    end
    prob_custom = LayeredEscortProblem(
        spec, F, f, PfF, bases,
        HomogeneousDynamics(dyn_test, K_test),
        target_trajs, nothing, nothing,
        0.05, 3, custom_positions
    )
    res_custom = run_layered_escort_simulation(prob_custom)
    for i in 1:spec.n_agents
        @test res_custom.sim_data[1][i][1:3] ≈ [0.0, 0.0, 3.0] atol=0.05
    end
end

@testset "Generalized Stalk Dimension" begin
    rings = [RingSpec(1, 4, 0.3)]
    supports = [SupportSpec(1, 1, 1)]

    # D != 4 now succeeds end-to-end for affine stalks.
    spec = LayeredEscortSpec(rings, supports; D=6)
    F = build_layered_escort_sheaf(spec)
    @test length(F.vertex_stalks) == spec.total_nodes
    @test F.vertex_stalks[1] == 6

    # affine=false with a nonzero ring radius is not representable.
    spec_bad = LayeredEscortSpec(rings, supports; D=6, affine=false)
    @test_throws Exception build_layered_escort_sheaf(spec_bad)
end

@testset "PlanarQuadrotorDynamics full pipeline (D=3, affine=true)" begin
    rings = [RingSpec(1, 4, 0.3), RingSpec(2, 4, 0.3)]
    supports = [SupportSpec(1, 2, 2)]
    spec = LayeredEscortSpec(rings, supports; D=3, affine=true)

    F = build_layered_escort_sheaf(spec)
    f = build_layered_homomorphism(spec)
    PfF = pushforward_sheaf(f, F)
    bases = build_layered_fiber_bases(f, F, spec)

    dyn_p = PlanarQuadrotorDynamics()
    Ad, Bd = CellularSheaves.AgentControllers.discrete_matrices(dyn_p, 0.05)
    K_p = CellularSheaves.AgentControllers.solve_dare(Ad, Bd, Matrix{Float64}(I, 6, 6), Matrix{Float64}(I, 2, 2))

    target_trajs = [
        t -> [1.0 + 0.05*cos(t), 1.5, 1.0],
        t -> [-1.0 - 0.05*cos(t), 1.5, 1.0]
    ]

    prob = LayeredEscortProblem(
        spec, F, f, PfF, bases,
        HomogeneousDynamics(dyn_p, K_p),
        target_trajs, nothing, nothing,
        0.05, 5, nothing
    )

    # Regression test: previously this crashed because agent state was
    # hardcoded to a 10D quadrotor pose regardless of dynamics_spec.
    res = run_layered_escort_simulation(prob)
    @test length(res.sim_data) == 5
    @test all(length(s) == 6 for s in res.sim_data[1])
end

@testset "RingSpec custom observers" begin
    rings = [RingSpec(1, 4, 0.3; observers=[1, 3])]
    supports = SupportSpec[]
    spec = LayeredEscortSpec(rings, supports)

    F = build_layered_escort_sheaf(spec)
    t_node = spec.target_nodes[1]

    @test !isnothing(CellularSheaves.get_restriction_map(F, 1, t_node))
    @test !isnothing(CellularSheaves.get_restriction_map(F, 3, t_node))
    @test_throws Exception CellularSheaves.get_restriction_map(F, 2, t_node)
end

@testset "SupportSpec asymmetric observers" begin
    rings = [RingSpec(1, 3, 0.3), RingSpec(2, 3, 0.3)]

    # Default: agent 1 sees both targets.
    supports_default = [SupportSpec(1, 2, 2)]
    spec_default = LayeredEscortSpec(rings, supports_default)
    F_default = build_layered_escort_sheaf(spec_default)
    s_range = spec_default.support_node_ranges[1]
    @test !isnothing(CellularSheaves.get_restriction_map(F_default, s_range[1], spec_default.target_nodes[1]))
    @test !isnothing(CellularSheaves.get_restriction_map(F_default, s_range[1], spec_default.target_nodes[2]))
    @test_throws Exception CellularSheaves.get_restriction_map(F_default, s_range[2], spec_default.target_nodes[1])
    @test_throws Exception CellularSheaves.get_restriction_map(F_default, s_range[2], spec_default.target_nodes[2])

    # Custom: replicate the old "agent 1 -> src, agent 2 -> tgt" behavior.
    supports_custom = [SupportSpec(1, 2, 2; src_observers=[1], tgt_observers=[2])]
    spec_custom = LayeredEscortSpec(rings, supports_custom)
    F_custom = build_layered_escort_sheaf(spec_custom)
    s_range_c = spec_custom.support_node_ranges[1]
    @test !isnothing(CellularSheaves.get_restriction_map(F_custom, s_range_c[1], spec_custom.target_nodes[1]))
    @test !isnothing(CellularSheaves.get_restriction_map(F_custom, s_range_c[2], spec_custom.target_nodes[2]))
    @test_throws Exception CellularSheaves.get_restriction_map(F_custom, s_range_c[1], spec_custom.target_nodes[2])
    @test_throws Exception CellularSheaves.get_restriction_map(F_custom, s_range_c[2], spec_custom.target_nodes[1])
end
