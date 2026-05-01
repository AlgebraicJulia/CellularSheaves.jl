using Test
using CellularSheaves
using LinearAlgebra
using Graphs

# ---------------------------------------------------------------------------
# Small helper: build an EuclideanSheaf on a path graph 1–2–…–n with
# constant stalk dimension `d` and identity restriction maps.
# ---------------------------------------------------------------------------
function constant_path_sheaf(n::Int, d::Int)
    s = EuclideanSheaf{Float64}(Int[])
    for _ in 1:n
        add_vertex_stalk!(s, d)
    end
    for i in 1:(n-1)
        rm = Matrix{Float64}(I, d, d)
        add_sheaf_edge!(s, i, i+1, rm, rm)
    end
    return s
end

@testset "GraphHomomorphism construction" begin
    # Valid: fold path 1-2-3 onto edge 1-2, sending vertex 3 → 1
    G = path_graph(3)
    H = path_graph(2)
    φ = GraphHomomorphism(G, H, [1, 2, 1])
    @test φ.vertex_map == [1, 2, 1]

    # Invalid: map vertex to out-of-range index
    @test_throws AssertionError GraphHomomorphism(G, H, [1, 2, 3])

    # Invalid: wrong length vertex map
    @test_throws AssertionError GraphHomomorphism(G, H, [1, 2])

    # Invalid: non-collapsed edge maps to non-existent target edge
    # complete_graph(3) has edge 1-3; path_graph(3) has only 1-2 and 2-3 (not 1-3)
    G3 = complete_graph(3)
    H3 = path_graph(3)
    @test_throws AssertionError GraphHomomorphism(G3, H3, [1, 3, 2])
end

@testset "fiber helpers" begin
    G = path_graph(3)   # 1-2-3
    H = path_graph(2)   # 1-2
    φ = GraphHomomorphism(G, H, [1, 2, 1])

    @test sort(fiber_vertices(φ, 1)) == [1, 3]
    @test fiber_vertices(φ, 2) == [2]

    # No edges collapse to any vertex in a path folding
    @test fiber_vertex_edges(φ, 1) == []
    @test fiber_vertex_edges(φ, 2) == []

    # Both source edges map to the single target edge {1,2}
    fe = fiber_edge_edges(φ, 1, 2)
    @test length(fe) == 2
    # Each tuple (u_s, u_d) must have φ(u_s)=1 and φ(u_d)=2
    for (u_s, u_d) in fe
        @test φ.vertex_map[u_s] == 1
        @test φ.vertex_map[u_d] == 2
    end

    # Reversed target vertices → same fibers with flipped orientation
    fe_rev = fiber_edge_edges(φ, 2, 1)
    @test length(fe_rev) == 2
    for (u_s, u_d) in fe_rev
        @test φ.vertex_map[u_s] == 2
        @test φ.vertex_map[u_d] == 1
    end
end

@testset "fiber_vertex_edges — collapsing triangle onto vertex" begin
    # G = triangle (3-cycle), H = single vertex
    G = cycle_graph(3)
    H = SimpleGraph(1)
    φ = GraphHomomorphism(G, H, [1, 1, 1])

    @test sort(fiber_vertices(φ, 1)) == [1, 2, 3]
    ce = fiber_vertex_edges(φ, 1)
    @test length(ce) == 3           # all 3 edges collapse
end

@testset "global_sections_basis — no collapsed edges" begin
    F = constant_path_sheaf(3, 2)
    # Fiber over vertex 1 = {vertex 1}, no collapsed edges
    basis = global_sections_basis(F, [1], Tuple{Int,Int}[])
    @test size(basis) == (2, 2)
    @test basis ≈ I(2)

    # Fiber over two isolated vertices
    basis2 = global_sections_basis(F, [1, 3], Tuple{Int,Int}[])
    @test size(basis2) == (4, 4)
    @test basis2 ≈ I(4)
end

@testset "global_sections_basis — with collapsed edges (triangle)" begin
    # Sheaf on triangle with identity restriction maps: global sections = R^d
    # (constant sections), so basis should have 1 column per stalk-dim unit.
    n, d = 3, 2
    G = cycle_graph(n)
    F = EuclideanSheaf{Float64}(Int[])
    for _ in 1:n
        add_vertex_stalk!(F, d)
    end
    for i in 1:n
        v1 = i; v2 = mod1(i+1, n)
        rm = Matrix{Float64}(I, d, d)
        add_sheaf_edge!(F, v1, v2, rm, rm)
    end

    all_verts = [1, 2, 3]
    all_edges = [(1,2), (2,3), (1,3)]
    basis = global_sections_basis(F, all_verts, all_edges)

    # Global sections of the constant sheaf on a connected graph = R^d
    @test size(basis, 2) == d
    # Each column should be a "constant" section: equal stalk values on all vertices
    for j in 1:size(basis, 2)
        col = basis[:, j]
        # blocks for vertices 1, 2, 3 should all be proportional to the same d-vector
        b1 = col[1:d];  b2 = col[d+1:2d];  b3 = col[2d+1:3d]
        @test norm(b1 - b2) < 1e-10
        @test norm(b1 - b3) < 1e-10
    end
end

@testset "pushforward — identity map recovers original sheaf" begin
    n, d = 4, 2
    G = path_graph(n)
    F = constant_path_sheaf(n, d)

    # Identity graph homomorphism
    φ = GraphHomomorphism(G, G, collect(1:n))
    pF = pushforward(φ, F)

    # Stalk dimensions should equal those of F
    @test vertex_stalks(pF) == vertex_stalks(F)
    @test edge_stalk_dimensions(pF) == edge_stalk_dimensions(F)

    # Coboundary maps should be equal (up to numerical noise from nullspace)
    dF  = Array(coboundary_map(F))
    dpF = Array(coboundary_map(pF))
    @test size(dF) == size(dpF)
    # The pushforward coboundary should have the same rank / null space as F
    @test rank(dF) == rank(dpF)
end

@testset "pushforward — fold path 1-2-3 onto edge 1-2" begin
    # G = path 1-2-3 with constant stalk R^1 and identity restriction maps
    # H = single edge 1-2
    # φ: 1→1, 2→2, 3→1
    G = path_graph(3)
    H = path_graph(2)
    F = constant_path_sheaf(3, 1)
    φ = GraphHomomorphism(G, H, [1, 2, 1])
    pF = pushforward(φ, F)

    # vertex 1 fiber = {1, 3}, no collapsed edges → stalk = R^2
    # vertex 2 fiber = {2}                        → stalk = R^1
    @test get_vertex_stalk(pF, 1) == 2
    @test get_vertex_stalk(pF, 2) == 1

    # Two source edges map to the single target edge, so edge stalk = R^2
    @test get_edge_stalk(pF, 1, 2) == 2

    # Check restriction maps have correct shapes
    rm1 = get_restriction_map(pF, 1, 2)   # (φ_*F)(1) → (φ_*F)(e)  : 2×2
    rm2 = get_restriction_map(pF, 2, 1)   # (φ_*F)(2) → (φ_*F)(e)  : 2×1
    @test size(rm1) == (2, 2)
    @test size(rm2) == (2, 1)

    # Verify the pushforward sheaf satisfies the sheaf condition:
    # for each global section s of pF, the coboundary d(s) = 0.
    # Global sections of pF lie in ker(coboundary_map(pF)).
    B = Array(coboundary_map(pF))
    # There should be at least one non-trivial global section
    ns = nullspace(B)
    @test size(ns, 2) >= 1
end

@testset "pushforward — collapse a triangle to a point" begin
    # G = triangle, H = single vertex, F = constant sheaf R^2
    # φ_*(F) has a single vertex whose stalk = global sections of F on the triangle
    #   = R^2 (constant sections), and no edges.
    n, d = 3, 2
    G = cycle_graph(n)
    H = SimpleGraph(1)
    F = EuclideanSheaf{Float64}(Int[])
    for _ in 1:n
        add_vertex_stalk!(F, d)
    end
    for i in 1:n
        v1 = i; v2 = mod1(i+1, n)
        rm = Matrix{Float64}(I, d, d)
        add_sheaf_edge!(F, v1, v2, rm, rm)
    end

    φ = GraphHomomorphism(G, H, [1, 1, 1])
    pF = pushforward(φ, F)

    @test nv(underlying_graph(pF)) == 1
    @test ne(underlying_graph(pF)) == 0
    @test get_vertex_stalk(pF, 1) == d   # global sections of constant sheaf on triangle
end

@testset "pushforward — two separate edges folded to one" begin
    # G: two disjoint edges 1-2 and 3-4  (path graph on 4 vertices would have 1-2-3-4;
    # instead build a custom graph)
    # H: single edge 1-2
    # φ: 1→1, 2→2, 3→1, 4→2
    G = SimpleGraph(4)
    add_edge!(G, 1, 2); add_edge!(G, 3, 4)
    H = path_graph(2)
    F = EuclideanSheaf{Float64}(Int[])
    for _ in 1:4; add_vertex_stalk!(F, 1); end
    for (u, v) in [(1,2),(3,4)]
        rm = ones(1,1)
        add_sheaf_edge!(F, u, v, rm, rm)
    end

    φ = GraphHomomorphism(G, H, [1, 2, 1, 2])
    pF = pushforward(φ, F)

    # Each vertex fiber has 2 isolated vertices → stalk = R^2
    @test get_vertex_stalk(pF, 1) == 2
    @test get_vertex_stalk(pF, 2) == 2
    # Two source edges map to the single target edge → edge stalk = R^2
    @test get_edge_stalk(pF, 1, 2) == 2

    # Restriction maps
    rm1 = get_restriction_map(pF, 1, 2)
    rm2 = get_restriction_map(pF, 2, 1)
    @test size(rm1) == (2, 2)
    @test size(rm2) == (2, 2)
end
