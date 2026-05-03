using Test
using CellularSheaves
using LinearAlgebra
using Graphs

# ---------------------------------------------------------------------------
# Helper: build the 6-cycle sheaf used throughout the pushforward example
# ---------------------------------------------------------------------------

function make_6cycle_sheaf()
    n = 6; stalk_dim = 3
    F = EuclideanSheaf{Float64}(Int[])
    for _ in 1:n
        add_vertex_stalk!(F, stalk_dim)
    end
    for i in 1:n
        v1 = i; v2 = (i % n) + 1
        rm = Matrix{Float64}(I, stalk_dim - 1, stalk_dim)   # 2×3 projection
        add_sheaf_edge!(F, v1, v2, rm, rm)
    end
    return F
end

# ---------------------------------------------------------------------------
# GraphHomomorphisms
# ---------------------------------------------------------------------------

@testset "GraphHomomorphism construction" begin
    hom = GraphHomomorphism([1, 1, 2, 2, 3, 3])
    @test hom.n_target == 3
    @test length(hom.vertex_map) == 6

    # Explicit n_target
    hom2 = GraphHomomorphism([1, 1, 2, 2, 3, 3], 3)
    @test hom2.n_target == 3

    # Validation: zero index
    @test_throws ArgumentError GraphHomomorphism([0, 1, 2])
    # Validation: negative n_target
    @test_throws ArgumentError GraphHomomorphism(Int[], -1)
    # Validation: entry exceeds n_target
    @test_throws ArgumentError GraphHomomorphism([1, 2, 4], 3)
end

@testset "fiber_vertices / fiber_edges / cross_edges" begin
    hom = GraphHomomorphism([1, 1, 2, 2, 3, 3])
    @test fiber_vertices(hom, 1) == [1, 2]
    @test fiber_vertices(hom, 2) == [3, 4]
    @test fiber_vertices(hom, 3) == [5, 6]

    F = make_6cycle_sheaf()
    g = underlying_graph(F)

    # Each fiber of the 6-cycle partition has exactly one intra-fiber edge
    for tv in 1:3
        fverts = fiber_vertices(hom, tv)
        fedges = fiber_edges(hom, g, tv)
        @test length(fedges) == 1
        (a, b) = fedges[1]
        @test a in fverts
        @test b in fverts
    end

    # There are exactly 3 cross edges (the three edges of the triangle)
    cedges = cross_edges(hom, g)
    @test length(cedges) == 3
    for (_, tu, tv) in cedges
        @test tu != tv
    end
end

# ---------------------------------------------------------------------------
# nullspace_ldlt
# ---------------------------------------------------------------------------

@testset "nullspace_ldlt on matrices" begin
    # Zero matrix → full nullspace
    NS = nullspace_ldlt(zeros(3, 3))
    @test size(NS, 1) == 3
    @test size(NS, 2) == 3

    # Rank-1 matrix → 2-dimensional nullspace
    A = zeros(3, 3); A[1, 1] = 1.0
    NS2 = nullspace_ldlt(A)
    @test size(NS2, 2) == 2

    # Full-rank identity → trivial nullspace
    NS3 = nullspace_ldlt(Matrix{Float64}(I, 4, 4))
    @test size(NS3, 2) == 0
end

@testset "nullspace_ldlt on EuclideanSheaf" begin
    F = make_6cycle_sheaf()
    # 6-cycle with 2×3 projection maps:
    # global sections have shared (a,b)∈ℝ² and 6 free scalars → dim 8
    NS_F = nullspace_ldlt(F)
    @test size(NS_F, 2) == 8

    # Columns should be in the nullspace of d'*d
    d = coboundary_map(F)
    @test norm(d * NS_F) < 1e-10
end

# ---------------------------------------------------------------------------
# Pushforward sheaf
# ---------------------------------------------------------------------------

@testset "pushforward_sheaf stalk dimensions" begin
    F   = make_6cycle_sheaf()
    hom = GraphHomomorphism([1, 1, 2, 2, 3, 3])
    PfF = pushforward_sheaf(hom, F)

    # Each fiber: 2 vertices (ℝ³ each) joined by 1 edge (ℝ²) → dim = 3+3-2 = 4
    @test vertex_stalks(PfF) == [4, 4, 4]
    # Three triangle edges, each with the same 2-dimensional stalk as F's edges
    @test all(==(2), collect(values(edge_stalks(PfF))))
end

@testset "pushforward_sheaf nullspace dimension" begin
    F   = make_6cycle_sheaf()
    hom = GraphHomomorphism([1, 1, 2, 2, 3, 3])
    PfF = pushforward_sheaf(hom, F)

    NS_F   = nullspace_ldlt(F)
    NS_PfF = nullspace_ldlt(PfF)

    # Both nullspaces have the same dimension (isomorphic as vector spaces)
    @test size(NS_F, 2) == size(NS_PfF, 2)
    @test size(NS_PfF, 2) == 8
end

# ---------------------------------------------------------------------------
# Transfer map
# ---------------------------------------------------------------------------

@testset "pushforward_transfer_map residual" begin
    F   = make_6cycle_sheaf()
    hom = GraphHomomorphism([1, 1, 2, 2, 3, 3])
    PfF = pushforward_sheaf(hom, F)

    NS_F = nullspace_ldlt(F)
    T    = pushforward_transfer_map(hom, F)
    d_PfF = coboundary_map(PfF)

    # Every global section of F maps to a global section of φ_*F
    @test norm(d_PfF * T * NS_F) < 1e-10
end

# ---------------------------------------------------------------------------
# Input validation (O(1) check)
# ---------------------------------------------------------------------------

@testset "vertex_map length validation" begin
    F   = make_6cycle_sheaf()
    # Too-short vertex map (4 entries for a 6-vertex sheaf)
    bad_hom = GraphHomomorphism([1, 1, 2, 2])
    @test_throws ArgumentError pushforward_sheaf(bad_hom, F)
    @test_throws ArgumentError pushforward_transfer_map(bad_hom, F)
end
