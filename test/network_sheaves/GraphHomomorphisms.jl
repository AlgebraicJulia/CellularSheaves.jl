using Test
using CellularSheaves
using Graphs

@testset "Graph pushout" begin

    @testset "Path gluing — standard case" begin
        # K = 1-edge path (2 vertices)
        # G = 3-vertex path, φ maps {1,2} → {2,3}
        # H = 3-vertex path, ψ maps {1,2} → {1,2}
        # Expected: Q is a 4-vertex path; nv(Q) == 4, ne(Q) == 3
        G = path_graph(3)
        H = path_graph(3)
        φ = GraphHomomorphism([2, 3], 3)
        ψ = GraphHomomorphism([1, 2], 3)
        Q, jG, jH = graph_pushout(φ, ψ, G, H, 2)
        @test nv(Q) == 4
        @test ne(Q) == 3
        # Pushout property: jG ∘ φ == jH ∘ ψ
        @test compose(φ, jG).vertex_map == compose(ψ, jH).vertex_map
    end

    @testset "Vertex identification only (K = single vertex)" begin
        G = path_graph(3)
        H = path_graph(3)
        # K = 1 vertex, φ maps to vertex 1 of G, ψ maps to vertex 1 of H
        φ = GraphHomomorphism([1], 3)
        ψ = GraphHomomorphism([1], 3)
        Q, jG, jH = graph_pushout(φ, ψ, G, H, 1)
        @test nv(Q) == nv(G) + nv(H) - 1
        @test compose(φ, jG).vertex_map == compose(ψ, jH).vertex_map
    end

    @testset "Disjoint union (K = empty graph, nK = 0)" begin
        G = path_graph(3)
        H = cycle_graph(4)
        φ = GraphHomomorphism(Int[], 3)
        ψ = GraphHomomorphism(Int[], 4)
        Q, jG, jH = graph_pushout(φ, ψ, G, H, 0)
        @test nv(Q) == nv(G) + nv(H)
        @test ne(Q) == ne(G) + ne(H)
    end

    @testset "Full identification (φ and ψ both identity on same graph)" begin
        G = path_graph(4)
        φ = GraphHomomorphism([1, 2, 3, 4], 4)
        ψ = GraphHomomorphism([1, 2, 3, 4], 4)
        Q, jG, jH = graph_pushout(φ, ψ, G, G, 4)
        @test nv(Q) == nv(G)
        @test ne(Q) == ne(G)
        @test compose(φ, jG).vertex_map == compose(ψ, jH).vertex_map
    end

    @testset "Edge merging" begin
        # K has an edge that maps to edges in both G and H
        # G = single edge (2 vertices), H = single edge (2 vertices)
        # φ and ψ both map K's edge endpoints to G's and H's endpoints respectively
        G = SimpleGraph(2); add_edge!(G, 1, 2)
        H = SimpleGraph(2); add_edge!(H, 1, 2)
        φ = GraphHomomorphism([1, 2], 2)
        ψ = GraphHomomorphism([1, 2], 2)
        Q, jG, jH = graph_pushout(φ, ψ, G, H, 2)
        # Vertices fully identified, so Q ≅ a single edge
        @test nv(Q) == 2
        @test ne(Q) == 1   # merged edge appears exactly once
    end

end

@testset "compose" begin
    # For a known 3-hop chain A →f B →g C
    f = GraphHomomorphism([2, 3, 1], 3)
    g = GraphHomomorphism([3, 1, 2], 3)
    h = compose(f, g)
    @test h.vertex_map == g.vertex_map[f.vertex_map]
    @test h.n_target == g.n_target

    # compose with identity
    id = GraphHomomorphism([1, 2, 3], 3)
    @test compose(f, id).vertex_map == f.vertex_map
    @test compose(id, f).vertex_map == f.vertex_map
end
