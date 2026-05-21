using Test
using CellularSheaves
using LinearAlgebra
using Graphs
using BlockArrays

@testset "harmonic_extension" begin

    # Build a path graph with 5 vertices, 1D stalks, identity restriction maps.
    # The sheaf Laplacian equals the standard graph Laplacian, so harmonic
    # extension is linear interpolation — expected values are known exactly.
    g = path_graph(5)
    s = sheaf_from_graph(g, 1, d -> Matrix{Float64}(I, d, d); symmetric_edges=true)
    boundary = Dict(1 => [0.0], 5 => [4.0])
    x_p, null_basis = harmonic_extension(s, boundary)

    # 1. Interpolation is exact on a path (1-D stalks)
    @test x_p[Block(1)] ≈ [0.0]
    @test x_p[Block(2)] ≈ [1.0]
    @test x_p[Block(3)] ≈ [2.0]
    @test x_p[Block(4)] ≈ [3.0]
    @test x_p[Block(5)] ≈ [4.0]

    # 2. Interior satisfies Lx = 0
    L = sheaf_laplacian_matrix(s)
    offsets = [0; cumsum(s.vertex_stalks)]
    for v in [2, 3, 4]
        idx = offsets[v]+1 : offsets[v+1]
        @test norm((L * Array(x_p))[idx]) < 1e-10
    end

    # 3. Boundary values are exactly preserved
    @test Array(x_p[Block(1)]) ≈ [0.0]
    @test Array(x_p[Block(5)]) ≈ [4.0]

    # 4. Unique solution has empty null basis
    @test size(null_basis, 2) == 0

    # 5. Multi-dimensional stalks (2-D)
    s2 = sheaf_from_graph(g, 2, d -> Matrix{Float64}(I, d, d); symmetric_edges=true)
    boundary2 = Dict(1 => [0.0, 1.0], 5 => [4.0, 5.0])
    x_p2, nb2 = harmonic_extension(s2, boundary2)
    L2 = sheaf_laplacian_matrix(s2)
    offsets2 = [0; cumsum(s2.vertex_stalks)]
    for v in [2, 3, 4]
        idx = offsets2[v]+1 : offsets2[v+1]
        @test norm((L2 * Array(x_p2))[idx]) < 1e-10
    end
    @test size(nb2, 2) == 0
    # Boundary values preserved
    @test Array(x_p2[Block(1)]) ≈ [0.0, 1.0]
    @test Array(x_p2[Block(5)]) ≈ [4.0, 5.0]

    # 6. All-boundary vertex case: short-circuit returns boundary values verbatim
    all_bnd = Dict(v => rand(s.vertex_stalks[v]) for v in 1:nv(s.underlying_graph))
    x_all, nb_all = harmonic_extension(s, all_bnd)
    @test size(nb_all, 2) == 0
    for v in 1:nv(s.underlying_graph)
        @test Array(x_all[Block(v)]) ≈ all_bnd[v]
    end

    # 7. All-interior vertex case (empty boundary): particular solution is zero,
    #    null basis columns are global sections
    x_none, nb_none = harmonic_extension(s, Dict{Int,Vector{Float64}}())
    @test norm(Array(x_none)) < 1e-10
    d_map = coboundary_map(s)
    @test size(nb_none, 2) > 0
    for j in 1:size(nb_none, 2)
        @test norm(d_map * nb_none[:, j]) < 1e-10
    end

    # 8. Disconnected interior component gives non-trivial null basis
    # Two disjoint 3-vertex path components: 1-2-3 and 4-5-6.
    # Boundary only on component 1 (vertices 1 and 3); vertices 4, 5, 6 are
    # isolated from the boundary, giving a 1-D null space on that component.
    g2 = Graph(6)
    add_edge!(g2, 1, 2); add_edge!(g2, 2, 3)
    add_edge!(g2, 4, 5); add_edge!(g2, 5, 6)
    s_two = sheaf_from_graph(g2, 1, d -> Matrix{Float64}(I, d, d); symmetric_edges=true)
    boundary_c1 = Dict(1 => [0.0], 3 => [2.0])
    x_p2c, nb2c = harmonic_extension(s_two, boundary_c1)
    @test size(nb2c, 2) > 0
    d2 = coboundary_map(s_two)
    for j in 1:size(nb2c, 2)
        @test norm(d2 * nb2c[:, j]) < 1e-10
        @test isapprox(nb2c[1, j], 0.0; atol=1e-10)
        @test isapprox(nb2c[3, j], 0.0; atol=1e-10)
    end

    # 9. Global section round-trip: boundary consistent with a global section
    #    → particular solution recovers it and null basis is empty
    sec = nullspace_ldlt(s)[:, 1]
    bv = Dict(1 => sec[1:1], 5 => sec[end:end])
    x_sec, nb_sec = harmonic_extension(s, bv)
    @test energy_function(s)(Array(x_sec)) < 1e-10
    @test size(nb_sec, 2) == 0

    # 10. Mismatched boundary vector length throws ArgumentError
    @test_throws ArgumentError harmonic_extension(s, Dict(1 => [0.0, 1.0], 5 => [4.0]))

    # 11. Invalid boundary vertex index throws ArgumentError
    @test_throws ArgumentError harmonic_extension(s, Dict(1 => [0.0], 6 => [4.0]))

    # 12. LDL diagnostics describe the restricted harmonic-extension problem
    ldl_diag = harmonic_extension_ldl_diagnostics(s, boundary)
    @test ldl_diag isa HarmonicExtensionLDLDiagnostics
    @test ldl_diag.boundary_vertices == [1, 5]
    @test ldl_diag.interior_vertices == [2, 3, 4]
    @test ldl_diag.restricted_size == (3, 3)
    @test ldl_diag.nullity_estimate == 0
    @test ldl_diag.rank_estimate == 3
    @test ldl_diag.smallest_positive_pivot > 0.0
    @test isfinite(ldl_diag.pivot_gap)
    @test occursin("HarmonicExtensionLDLDiagnostics", sprint(show, MIME("text/plain"), ldl_diag))

    # 13. Empty-interior problems return empty diagnostics consistently
    all_boundary_diag_ldl = harmonic_extension_ldl_diagnostics(s, all_bnd)
    @test all_boundary_diag_ldl.restricted_size == (0, 0)
    @test isempty(all_boundary_diag_ldl.pivots)
    @test all_boundary_diag_ldl.nullity_estimate == 0

end
