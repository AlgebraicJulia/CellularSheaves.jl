using Test
using CellularSheaves
using CellularSheaves.ControlSheaves.NestedSystems
using LinearAlgebra
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

end # testset NestedSystems
