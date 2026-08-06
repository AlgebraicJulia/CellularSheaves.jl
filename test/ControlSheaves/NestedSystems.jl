using Test
using CellularSheaves
using CellularSheaves.ControlSheaves.NestedSystems
using CellularSheaves.ControlSheaves.AgentControllers
using CellularSheaves.ControlSheaves.Layered
using LinearAlgebra
using Statistics
using Graphs

include("nested_system_test_specs.jl")

@testset "NestedSystems" begin

@testset "tower — depth and level count" begin
    spec = two_level_spec()          # helper: 2 teams, 1 target each
    tower = build_sheaf_tower(spec)
    @test tower.depth == 2
    @test length(tower.levels) == 2
    @test length(tower.homs) == 1
    @test length(tower.bases) == 1
end

@testset "tower — targets are singleton fibres with identity bases at every level" begin
    tower = build_sheaf_tower(three_level_spec())
    for k in 1:length(tower.homs), t in tower.target_vertices
        @test length(fiber_vertices(tower.homs[k], t)) == 1
        @test tower.bases[k][t] ≈ I(tower.spec.D)
    end
end

@testset "tower — a rigid team collapses to exactly D dimensions" begin
    tower = build_sheaf_tower(two_level_spec())
    team_vertex = 1
    @test vertex_stalks(tower.levels[1])[team_vertex] == tower.spec.D
end

@testset "tower — one target observed by two teams" begin
    tower = build_sheaf_tower(shared_target_spec())
    t = tower.target_vertices[1]
    @test degree(tower.levels[1].underlying_graph, t) == 2
end

@testset "tower — irregular depth across siblings" begin
    # left child refined twice, right child a bare leaf
    tower = build_sheaf_tower(irregular_spec())
    @test tower.depth == 3
    @test length(tower.levels) == 3
end

@testset "tower — over-constrained fibre is rejected" begin
    @test_throws Exception build_sheaf_tower(degenerate_spec())
end

# ---------------------------------------------------------------------------
# Issue 010 — hierarchical solve, direct baseline, approximation gap
# ---------------------------------------------------------------------------

@testset "energy gap is nonnegative (feasibility theorem)" begin
    for spec in [two_level_spec(), three_level_spec(), irregular_spec(), shared_target_spec()]
        tower = build_sheaf_tower(spec)
        g = approximation_gap(tower, default_targets(spec))
        @test g.gap >= -1e-8
    end
end

@testset "rigid single-target teams: hierarchical == direct" begin
    tower = build_sheaf_tower(two_level_spec())   # rigid rings, one target each
    tv = default_targets(tower.spec)
    q_h = solve_hierarchical(tower, tv)[end]
    q_d = solve_direct(tower, tv)
    @test all(isapprox(a, b; atol=1e-8) for (a, b) in zip(q_h, q_d))
    @test approximation_gap(tower, tv).gap ≈ 0.0 atol=1e-8
end

@testset "team observing two separating targets: strict gap" begin
    tower = build_sheaf_tower(two_target_team_spec())
    near = [[0.0, 0.0, 1.5, 1.0], [0.2, 0.0, 1.5, 1.0]]
    far  = [[0.0, 0.0, 1.5, 1.0], [8.0, 0.0, 1.5, 1.0]]
    @test approximation_gap(tower, far).gap > approximation_gap(tower, near).gap
    @test approximation_gap(tower, far).gap > 1e-6

    # This configuration is exactly solvable, so pin the closed form rather than an
    # inequality. Each rigid ring translates for free, leaving a chain of unit springs
    # between the two pinned targets a distance Δ apart:
    #   direct       — t1 — repA — repB — t2, three springs in series: Δ²/3
    #   hierarchical — repA and repB locked together, two springs: Δ²/2
    # so the gap is Δ²/6 and the relative gap is 1/2, independent of Δ.
    for tv in (near, far)
        Δ = tv[2][1] - tv[1][1]
        g = approximation_gap(tower, tv)
        @test g.direct ≈ Δ^2 / 3 atol=1e-8
        @test g.hierarchical ≈ Δ^2 / 2 atol=1e-8
        @test g.gap ≈ Δ^2 / 6 atol=1e-8
        @test g.relative_gap ≈ 0.5 atol=1e-8
    end
end

@testset "targets are reproduced exactly at the finest level" begin
    tower = build_sheaf_tower(three_level_spec())
    tv = default_targets(tower.spec)
    levels = solve_hierarchical(tower, tv)
    @test length(levels) == tower.depth
    q = levels[end]
    for (t, v) in enumerate(tower.target_vertices)
        @test q[v] ≈ tv[t]
    end
    # ... and at every intermediate level too, since targets are singleton fibres throughout.
    for lvl in levels, (t, v) in enumerate(tower.target_vertices)
        @test lvl[v] ≈ tv[t]
    end
end

@testset "golden: a flat two-level tower hits the closed-form answer" begin
    # One rigid ring observing one target through an identity edge on its first agent.
    # The exact optimum is therefore zero-energy, with agent 1 sitting on the target and
    # the rest of the ring holding its formation offsets.
    tower = build_sheaf_tower(flat_equivalent_spec())
    @test tower.depth == 2
    F = tower.levels[end]

    p1 = [1.0, 0.0, 1.5, 1.0]
    q1 = solve_hierarchical(tower, [p1])[end]
    @test sheaf_energy(F, q1) ≈ 0.0 atol=1e-8
    @test q1[tower.target_vertices[1]] ≈ p1
    @test q1[tower.agent_vertices[1]] ≈ p1 atol=1e-8
    @test approximation_gap(tower, [p1]).gap ≈ 0.0 atol=1e-8

    # Rigidity: moving the target translates the whole formation without deforming it.
    p2 = [-2.0, 3.0, 1.5, 1.0]
    q2 = solve_hierarchical(tower, [p2])[end]
    for a in tower.agent_vertices
        @test (q1[a] - q1[tower.agent_vertices[1]]) ≈
              (q2[a] - q2[tower.agent_vertices[1]]) atol=1e-8
    end
    @test all(q2[a] - q1[a] ≈ p2 - p1 for a in tower.agent_vertices)
end

@testset "solver input validation" begin
    tower = build_sheaf_tower(two_level_spec())
    @test_throws Exception solve_hierarchical(tower, [zeros(3)])           # too few targets
    @test_throws Exception solve_direct(tower, [zeros(2), zeros(3)])       # wrong stalk dim
end

# ---------------------------------------------------------------------------
# Issue 011 — per-edge restriction maps
# ---------------------------------------------------------------------------

@testset "project(i) selects exactly child i's block" begin
    node = three_child_system()
    R = materialize_restriction(project(2), node, 4)
    @test size(R) == (4, 12)
    @test R[:, 5:8] ≈ I(4)
    @test all(iszero, R[:, 1:4]) && all(iszero, R[:, 9:12])
end

@testset "project(:name) resolves by child name" begin
    node = three_child_system()
    @test materialize_restriction(project(:bravo), node, 4) ≈
          materialize_restriction(project(2), node, 4)
    @test_throws Exception materialize_restriction(project(:nonexistent), node, 4)
end

@testset "centroid averages direct children" begin
    node = three_child_system()
    R = materialize_restriction(centroid(), node, 4)
    @test R[:, 1:4] ≈ I(4) / 3
    @test R * repeat([1.0, 2.0, 3.0, 1.0], 3) ≈ [1.0, 2.0, 3.0, 1.0]
end

@testset "centroid treats a refined child as one opaque unit" begin
    # `wide` refines into 6 agents, `narrow` into 2; centroid must still weight them 1/2 each.
    node = mixed_arity_system()
    R = materialize_restriction(centroid(), node, 4)
    @test R[:, 1:4] ≈ I(4) / 2
end

@testset "raw matrix escape hatch validates its shape" begin
    node = three_child_system()
    @test_throws Exception materialize_restriction(RawRestriction(zeros(4, 5)), node, 4)
end

@testset "default project(1) reproduces Issue 009 behaviour" begin
    tv = default_targets(spec_default())
    q_explicit = solve_hierarchical(build_sheaf_tower(spec_explicit_project1()), tv)[end]
    q_default = solve_hierarchical(build_sheaf_tower(spec_default()), tv)[end]
    @test q_explicit == q_default
end

@testset "centroid-wired tower still satisfies the energy-gap theorem" begin
    tower = build_sheaf_tower(centroid_wired_spec())
    @test approximation_gap(tower, default_targets(tower.spec)).gap >= -1e-8
end

@testset "centroid-wired tower reproduces a rigid translation" begin
    # Both rings are rigid and the single target is observed through only one of them (via the
    # default project(1) on the observation edge), so moving the target still translates the
    # *whole* two-ring assembly rigidly -- the centroid coupling adds no extra deformation for a
    # single, unopposed target. Gap must stay ~0, matching the project(1)-only rigid case.
    tower = build_sheaf_tower(centroid_wired_spec())
    tv = default_targets(tower.spec)
    q = solve_hierarchical(tower, tv)[end]
    @test sheaf_energy(tower.levels[end], q) ≈ 0.0 atol=1e-8
    @test approximation_gap(tower, tv).gap ≈ 0.0 atol=1e-8
end

# ---------------------------------------------------------------------------
# Issue 012 — nested dynamics context with most-specific-wins cascade
# ---------------------------------------------------------------------------

@testset "root default cascades to every leaf" begin
    spec = three_level_spec()
    ctx = SystemBinding(dynamics=QuadrotorDynamics())
    resolved = resolve_dynamics(spec, ctx)
    @test length(resolved) == n_agents(spec)
    @test all(r -> r.dynamics isa QuadrotorDynamics, resolved)
end

@testset "most specific wins: agent over team over ancestor" begin
    spec = two_level_spec()   # team1 (3 agents), team2 (3 agents)
    ctx = SystemBinding(
        dynamics = QuadrotorDynamics(),
        children = Dict(:team1 => SystemBinding(
            dynamics = PlanarQuadrotorDynamics(),
            agents = Dict(2 => AgentBinding(dynamics = QuadrotorDynamics())))))
    r = resolve_dynamics(spec, ctx)
    @test r[1].dynamics isa PlanarQuadrotorDynamics   # team default
    @test r[2].dynamics isa QuadrotorDynamics         # per-agent override
    @test last(r).dynamics isa QuadrotorDynamics      # root default, other subtree (team2)
end

@testset "fields resolve independently" begin
    # override only the initial position; dynamics must still be inherited
    spec = two_level_spec()
    ctx = SystemBinding(dynamics=QuadrotorDynamics(),
                        children=Dict(:team1 => SystemBinding(
                            agents=Dict(1 => AgentBinding(initial_position=[9.0, 9.0, 9.0])))))
    r = resolve_dynamics(spec, ctx)
    @test r[1].dynamics isa QuadrotorDynamics
    @test r[1].initial_position == [9.0, 9.0, 9.0]
end

@testset "unbound dynamics throws naming the path" begin
    spec = two_level_spec()
    err = try resolve_dynamics(spec, SystemBinding()); nothing catch e; e end
    @test err !== nothing
    @test occursin("team1", sprint(showerror, err))
end

@testset "typo'd child name is rejected" begin
    spec = two_level_spec()
    ctx = SystemBinding(dynamics=QuadrotorDynamics(),
                        children=Dict(:teamTypo => SystemBinding()))
    @test_throws Exception resolve_dynamics(spec, ctx)
end

@testset "out-of-range agent index is rejected" begin
    spec = two_level_spec()   # team1 has 3 agents
    ctx = SystemBinding(dynamics=QuadrotorDynamics(),
                        children=Dict(:team1 => SystemBinding(
                            agents=Dict(99 => AgentBinding()))))
    @test_throws Exception resolve_dynamics(spec, ctx)
end

@testset "agents declared on a RefinedSystem are rejected" begin
    spec = three_level_spec()   # root -> mid -> {teamA, teamB}
    ctx = SystemBinding(dynamics=QuadrotorDynamics(),
                        children=Dict(:mid => SystemBinding(agents=Dict(1 => AgentBinding()))))
    @test_throws Exception resolve_dynamics(spec, ctx)
end

@testset "root-only binding applies uniformly (flat trichotomy retired, Issue 013)" begin
    spec = two_level_spec()
    dyn = QuadrotorDynamics()
    r = resolve_dynamics(spec, SystemBinding(dynamics=dyn))
    @test all(a -> a.dynamics === dyn, r)
end

@testset "resolved order matches the tower's agent ordering" begin
    spec = three_level_spec()
    tower = build_sheaf_tower(spec)
    r = resolve_dynamics(spec, SystemBinding(dynamics=QuadrotorDynamics()))
    @test [x.agent_index for x in r] == collect(eachindex(tower.agent_vertices))
end

@testset "homogeneous_binding shorthand matches an explicit root binding" begin
    spec = two_level_spec()
    dyn = QuadrotorDynamics()
    K = ones(3, 10)
    r1 = resolve_dynamics(spec, homogeneous_binding(dyn; K_lqr=K))
    r2 = resolve_dynamics(spec, SystemBinding(dynamics=dyn, K_lqr=K))
    @test all(a -> a.dynamics === dyn && a.K_lqr === K, r1)
    @test [a.initial_position for a in r1] == [a.initial_position for a in r2]
end

# ---------------------------------------------------------------------------
# Simulation driver and multi-pin helpers (promoted from docs/literate/nested/)
#
# These also serve as quantitative regression tests for the three nested/ literate examples --
# the same specs, closed forms, and thresholds verified numerically while writing them, pinned
# here so their behavior doesn't silently drift and get caught only by eyeballing a plot.
# ---------------------------------------------------------------------------

@testset "translation_pin leaves the homogeneous row unconstrained" begin
    r = translation_pin(4, 4, 2)
    @test r isa RawRestriction
    @test size(r.M) == (4, 16)
    @test all(iszero, r.M[4, :])                       # homogeneous row entirely zero
    @test r.M[1:3, 5:7] ≈ I(3)                          # member 2's translation block selected
    @test all(iszero, r.M[1:3, setdiff(1:16, 5:7)])     # nothing else touched
end

@testset "redundant_pin: project(1) for k=1, translation_pin otherwise" begin
    @test redundant_pin(4, 4, 1) isa ProjectMember
    @test redundant_pin(4, 4, 1).member == 1
    @test redundant_pin(4, 4, 3) isa RawRestriction
end

@testset "agent_index_ranges reconstructs contiguous blocks" begin
    @test agent_index_ranges([3, 2, 4]) == [1:3, 4:5, 6:9]
end

@testset "run_nested_escort_simulation basic smoke test" begin
    spec = two_level_spec()
    tower = build_sheaf_tower(spec)
    dyn = QuadrotorDynamics()
    Ad, Bd = CellularSheaves.AgentControllers.discrete_matrices(dyn, 0.05)
    K = CellularSheaves.AgentControllers.solve_dare(Ad, Bd, Matrix{Float64}(I, 10, 10), Matrix{Float64}(I, 3, 3))
    bindings = homogeneous_binding(dyn; K_lqr=K)
    tv = default_targets(spec)
    trajs = [t -> tv[i] for i in eachindex(tv)]
    prob = NestedEscortProblem(tower, bindings, trajs; dt=0.05, steps=5)
    res = run_nested_escort_simulation(prob; use_feedforward=false)
    @test length(res.sim_data) == 5
    @test length(res.sim_data[1]) == length(tower.agent_vertices)
    @test all(x -> all(isfinite, x), res.sim_data[end])
end

@testset "quantitative Example A: redundant pin gives exact midpoint, gap grows with separation" begin
    ringA = LeafTeam(:ringA, :ring, 4, 1.0)
    ringB = LeafTeam(:ringB, :ring, 5, 1.0)
    mid = RefinedSystem(:mid, AbstractSystemNode[ringA, ringB],
                        [SystemEdge(1, 2; src_map=centroid(), dst_map=centroid())])
    root = RefinedSystem(:root, AbstractSystemNode[mid])
    targets = [TargetSpec(:t1), TargetSpec(:t2)]
    observations = vcat(
        [Observation([1, 1], 1; system_map=redundant_pin(ringA.n_agents, 4, k)) for k in 1:2:ringA.n_agents],
        [Observation([1, 2], 2; system_map=redundant_pin(ringB.n_agents, 4, k)) for k in 1:2:ringB.n_agents],
    )
    spec = NestedSystemSpec(root, targets, observations, 4, true)
    tower = build_sheaf_tower(spec)

    t1v, t2v = [3.0, 1.0, 1.5, 1.0], [-1.0, -2.0, 1.5, 1.0]
    q = solve_hierarchical(tower, [t1v, t2v])[end]
    cA = sum(q[v] for v in tower.agent_vertices[1:4]) / 4
    cB = sum(q[v] for v in tower.agent_vertices[5:9]) / 5
    midpoint = (t1v .+ t2v) ./ 2
    @test cA ≈ cB atol=1e-8
    @test cA[1:3] ≈ midpoint[1:3] atol=1e-8
    @test all(q[v][4] ≈ 1.0 for v in tower.agent_vertices)   # homogeneous coordinate untouched

    g_near = approximation_gap(tower, [[0.5, 0.0, 1.5, 1.0], [-0.5, 0.0, 1.5, 1.0]])
    g_far = approximation_gap(tower, [[2.0, 0.0, 1.5, 1.0], [-2.0, 0.0, 1.5, 1.0]])
    @test g_near.direct ≈ 0.0 atol=1e-8   # only the first, H_N-representable pin is visible to direct
    @test g_far.direct ≈ 0.0 atol=1e-8
    @test g_far.gap > g_near.gap
end

@testset "quantitative Example B: redundant pin measurably reduces cyclic-coupling tracking error" begin
    n, m, R_big, support_m = 3, [5, 5, 5], 3.0, 2
    escort_radius(m_i, m_max; safety=0.35) = safety * R_big * sin(pi / n) * (m_i / m_max)
    m_max = maximum(m)

    function build_cyclic_tower(observation_builder)
        rings = [LeafTeam(Symbol(:ring, i), :ring, m[i], escort_radius(m[i], m_max)) for i in 1:n]
        pods = [LeafTeam(Symbol(:pod, i), :path, support_m, 0.3 * escort_radius(m_max, m_max)) for i in 1:n]
        children = AbstractSystemNode[]
        for i in 1:n
            push!(children, rings[i]); push!(children, pods[i])
        end
        cyc_edges = SystemEdge[]
        for i in 1:n
            push!(cyc_edges, SystemEdge(2i - 1, 2i; src_map=centroid(), dst_map=centroid()))
        end
        for i in 1:n
            push!(cyc_edges, SystemEdge(2i, (2i % 2n) + 1; src_map=centroid(), dst_map=centroid()))
        end
        root = RefinedSystem(:root, children, cyc_edges)
        targets = [TargetSpec(Symbol(:t, i)) for i in 1:n]
        spec = NestedSystemSpec(root, targets, observation_builder(), 4, true)
        return build_sheaf_tower(spec)
    end

    single_pin() = [Observation([2i - 1], i) for i in 1:n]
    redundant() = [Observation([2i - 1], i; system_map=redundant_pin(m[i], 4, k))
                  for i in 1:n for k in 1:2:m[i]]

    tv = [[R_big * cos(2π * (i - 1) / n), R_big * sin(2π * (i - 1) / n), 1.5, 1.0] for i in 1:n]

    function mean_tracking_error(tower)
        q = solve_hierarchical(tower, tv)[end]
        off, errs = 0, Float64[]
        for i in 1:n
            rng = tower.agent_vertices[(off + 1):(off + m[i])]
            off += m[i] + support_m
            c = sum(q[v] for v in rng) / m[i]
            push!(errs, norm(c[1:2] .- tv[i][1:2]))
        end
        return sum(errs) / n
    end

    err_single = mean_tracking_error(build_cyclic_tower(single_pin))
    err_redundant = mean_tracking_error(build_cyclic_tower(redundant))
    @test err_redundant < err_single
    ## Exact values verified numerically while writing n_ring_formation.jl; pinned as a
    ## regression check, not just a directional one.
    @test isapprox(err_single, 2.01; atol=0.05)
    @test isapprox(err_redundant, 1.80; atol=0.05)

    ## Why the tracking error never reaches zero: because the pods have no target of their own,
    ## each ring's actual position is an exact linear blend of *every* target, not just its own --
    ## measured directly by perturbing one target at a time. This was caught by exactly this kind
    ## of direct, index-checked measurement after the docs page's own version of this computation
    ## initially used the wrong indexing convention (`agent_index_ranges`'s *local* agent-position
    ## ranges directly against `solve_hierarchical`'s *actual-vertex-numbered* output) and silently
    ## produced a plausible-looking but wrong number -- worth a regression test on its own.
    function target_response_weights(tower, ring_idx::Int)
        off = 0
        verts = Int[]
        for i in 1:n
            rng = tower.agent_vertices[(off + 1):(off + m[i])]
            i == ring_idx && (verts = rng)
            off += m[i] + support_m
        end
        base = fill([0.0, 0.0, 1.5, 1.0], n)
        ring_centroid(tv) = sum(solve_hierarchical(tower, tv)[end][v] for v in verts) / length(verts)
        c0 = ring_centroid(base)
        h = 1.0
        return [begin
                    perturbed = copy(base)
                    perturbed[j] = base[j] .+ [h, 0.0, 0.0, 0.0]
                    (ring_centroid(perturbed)[1] - c0[1]) / h
                end
                for j in 1:n]
    end

    w_single = target_response_weights(build_cyclic_tower(single_pin), 1)
    w_redundant = target_response_weights(build_cyclic_tower(redundant), 1)
    @test sum(w_single) ≈ 1.0 atol=1e-8      # translation-invariance: weights must blend to 1
    @test sum(w_redundant) ≈ 1.0 atol=1e-8
    @test isapprox(w_single[1], 0.5558; atol=1e-3)
    @test isapprox(w_redundant[1], 0.6; atol=1e-8)
    @test w_redundant[1] > w_single[1]        # redundant pin measurably raises the own-target weight
    @test w_redundant[2] ≈ w_redundant[3] atol=1e-8   # symmetric under the cycle's own symmetry
end

@testset "quantitative Example C: full redundant pins rescale a rigid team, never shear it" begin
    n = 6
    ring = LeafTeam(:formation, :ring, n, 1.0)
    root = RefinedSystem(:root, AbstractSystemNode[ring])
    targets = [TargetSpec(Symbol(:t, k)) for k in 1:n]
    ## Deliberately full, untruncated project(k) pins -- the point is to let w drift.
    observations = [Observation([1], k; system_map=project(k)) for k in 1:n]
    spec = NestedSystemSpec(root, targets, observations, 4, true)
    tower = build_sheaf_tower(spec)

    h_alt, base_rate = 1.5, 0.4
    scale_factor(k) = 0.6 + 0.8 * (k - 1) / (n - 1)
    angle(k) = 2π * (k - 1) / n
    target_traj(k, t) = [base_rate * scale_factor(k) * t * cos(angle(k)),
                         base_rate * scale_factor(k) * t * sin(angle(k)), h_alt, 1.0]

    for t in [1.0, 5.0, 10.0]
        tv = [target_traj(k, t) for k in 1:n]
        q = solve_hierarchical(tower, tv)[end]
        ws = [q[v][4] for v in tower.agent_vertices]
        @test all(w -> isapprox(w, ws[1]; atol=1e-8), ws)     # one shared homogeneous coordinate
        @test !isapprox(ws[1], 1.0; atol=1e-3)                # it genuinely drifted -- the mechanism firing
        positions = [q[v][1:3] for v in tower.agent_vertices]
        c = sum(positions) / n
        @test std([norm(p .- c) for p in positions]) < 1e-8   # perfectly regular: rescaled, not sheared
    end

    ## The same non-uniform divergence, tracked with no shape constraint at all, is measurably
    ## irregular -- the "shear" this rigid formation avoids.
    tv = [target_traj(k, 10.0) for k in 1:n]
    c = sum(p[1:3] for p in tv) / n
    @test std([norm(p[1:3] .- c) for p in tv]) > 0.5
end

end # testset NestedSystems
