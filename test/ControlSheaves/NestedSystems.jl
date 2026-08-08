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

end # testset NestedSystems
