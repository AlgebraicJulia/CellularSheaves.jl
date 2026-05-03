using Test
using CellularSheaves
using LinearAlgebra
using Graphs
using SparseArrays

@testset "GraphHomomorphism" begin
    hom = GraphHomomorphism([1, 2, 1, 2])
    @test hom.n_target == 2
    @test fiber_vertices(hom, 1) == [1, 3]
    @test fiber_vertices(hom, 2) == [2, 4]

    # Explicit n_target
    hom2 = GraphHomomorphism([1, 2, 1, 2], 3)
    @test hom2.n_target == 3

    # Validation errors
    @test_throws ArgumentError GraphHomomorphism([1, 2, 0])
    @test_throws ArgumentError GraphHomomorphism([1, 2, 3], 2)
    @test_throws ArgumentError GraphHomomorphism(Int[], -1)
end

@testset "fiber_edges / cross_edges" begin
    g = Graph(4)
    add_edge!(g, 1, 2)
    add_edge!(g, 3, 4)
    add_edge!(g, 1, 3)  # cross edge
    add_edge!(g, 2, 4)  # cross edge

    hom = GraphHomomorphism([1, 1, 2, 2])

    fe1 = fiber_edges(hom, g, 1)
    @test (1, 2) in fe1
    @test length(fe1) == 1

    fe2 = fiber_edges(hom, g, 2)
    @test (3, 4) in fe2
    @test length(fe2) == 1

    ce = cross_edges(hom, g)
    @test length(ce) == 2
    for (e, pu, pv) in ce
        @test pu != pv
    end
end

@testset "pushforward_sheaf — injective edge map" begin
    # Simple path graph: 4 vertices, edges 1-2, 2-3, 3-4
    # Homomorphism: vertices 1,2 → target 1; vertices 3,4 → target 2
    # Cross-edges: edge 2-3 is the only cross-edge → injective edge map

    g = path_graph(4)
    F = EuclideanSheaf{Float64}(repeat([2], 4))
    # Identity restriction maps everywhere (R^2 → R^2)
    for e in edges(g)
        add_sheaf_edge!(F, src(e), dst(e), Matrix{Float64}(I, 2, 2), Matrix{Float64}(I, 2, 2))
    end

    hom = GraphHomomorphism([1, 1, 2, 2])
    PfF = pushforward_sheaf(hom, F)

    # Target graph has 2 vertices and 1 edge
    @test nv(underlying_graph(PfF)) == 2
    @test ne(underlying_graph(PfF)) == 1

    # Global sections of PfF ↔ global sections of F
    d_F   = sparse(coboundary_map(F))
    d_PfF = sparse(coboundary_map(PfF))
    ns_F   = nullspace_ldlt(d_F'   * d_F)
    ns_PfF = nullspace_ldlt(d_PfF' * d_PfF)
    @test size(ns_F, 2) == size(ns_PfF, 2)
end

@testset "pushforward_sheaf — non-injective edge map (two source edges → one target edge)" begin
    # 4-vertex graph; hom maps: 1→1, 2→2, 3→1, 4→2
    # Edges: 1-2 and 3-4 are both cross-edges mapping to target edge 1-2.
    # This is the non-injective case the bug was about.

    g = Graph(4)
    add_edge!(g, 1, 2)  # cross-edge → target (1,2)
    add_edge!(g, 3, 4)  # cross-edge → target (1,2)

    # 1-dimensional stalks everywhere, identity restriction maps
    F = EuclideanSheaf{Float64}(repeat([1], 4))
    for e in edges(g)
        add_sheaf_edge!(F, src(e), dst(e), ones(1, 1), ones(1, 1))
    end

    hom = GraphHomomorphism([1, 2, 1, 2])
    PfF = pushforward_sheaf(hom, F)

    # Target graph has 2 vertices and 1 edge
    @test nv(underlying_graph(PfF)) == 2
    @test ne(underlying_graph(PfF)) == 1

    # Edge stalk of PfF at (1,2) should be 2-dimensional (direct sum of two 1-dim stalks)
    @test get_edge_stalk(PfF, 1, 2) == 2
end

@testset "pushforward_sheaf — non-injective: global section count preserved" begin
    # Same setup: 4 vertices, hom: 1→1, 2→2, 3→1, 4→2
    # Two cross-edges (1-2) and (3-4) both mapping to target (1,2)
    # Plus fiber edges: 1-3 (within fiber 1) and 2-4 (within fiber 2)

    g = Graph(4)
    add_edge!(g, 1, 2)  # cross-edge
    add_edge!(g, 3, 4)  # cross-edge
    add_edge!(g, 1, 3)  # fiber edge for target vertex 1
    add_edge!(g, 2, 4)  # fiber edge for target vertex 2

    # 1-dimensional stalks, identity restriction maps
    F = EuclideanSheaf{Float64}(repeat([1], 4))
    for e in edges(g)
        add_sheaf_edge!(F, src(e), dst(e), ones(1, 1), ones(1, 1))
    end

    hom = GraphHomomorphism([1, 2, 1, 2])
    PfF = pushforward_sheaf(hom, F)

    # Target graph has 2 vertices and 1 edge
    @test nv(underlying_graph(PfF)) == 2
    @test ne(underlying_graph(PfF)) == 1

    # Edge stalk should be 2-dimensional (direct sum of the two cross-edge stalks)
    @test get_edge_stalk(PfF, 1, 2) == 2

    # Global section count: F and PfF should have the same number of global sections
    d_F   = sparse(coboundary_map(F))
    d_PfF = sparse(coboundary_map(PfF))
    ns_F   = nullspace_ldlt(d_F'   * d_F)
    ns_PfF = nullspace_ldlt(d_PfF' * d_PfF)
    @test size(ns_F, 2) == size(ns_PfF, 2)
end

@testset "pushforward_transfer_map — maps global sections correctly" begin
    # Path graph 1-2-3-4 with hom: 1,2→1; 3,4→2
    g = path_graph(4)
    F = EuclideanSheaf{Float64}(repeat([1], 4))
    for e in edges(g)
        add_sheaf_edge!(F, src(e), dst(e), ones(1, 1), ones(1, 1))
    end

    hom = GraphHomomorphism([1, 1, 2, 2])
    PfF = pushforward_sheaf(hom, F)
    T   = pushforward_transfer_map(hom, F)

    d_PfF = sparse(coboundary_map(PfF))

    # For every global section s of F, d_PfF * T * s ≈ 0
    d_F  = sparse(coboundary_map(F))
    ns_F = nullspace_ldlt(d_F' * d_F)

    for k in axes(ns_F, 2)
        s = ns_F[:, k]
        residual = d_PfF * T * s
        @test norm(residual) < 1e-10
    end
end
