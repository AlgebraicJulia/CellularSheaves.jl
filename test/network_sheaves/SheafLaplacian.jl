using Test
using CellularSheaves
using LinearAlgebra
using Graphs
using SparseArrays

@testset "sheaf_laplacian_matrix" begin
    # ── helper: path graph with identity restriction maps ──────────────────
    g = path_graph(5)
    s = sheaf_from_graph(g, 1, d -> Matrix{Float64}(I, d, d); symmetric_edges=true)

    L = sheaf_laplacian_matrix(s)
    L_sparse = sparse(L)

    # 1. sparse() returns a SparseMatrixCSC
    @test L_sparse isa SparseMatrixCSC

    # 2. Symmetric
    @test issymmetric(Array(L))

    # ── 2-D stalk sheaf ────────────────────────────────────────────────────
    s2 = sheaf_from_graph(g, 2, d -> Matrix{Float64}(I, d, d); symmetric_edges=true)
    @test issymmetric(Array(sheaf_laplacian_matrix(s2)))

    # ── random non-identity restriction maps ──────────────────────────────
    g3 = complete_graph(4)
    rand_map = d -> randn(d - 1, d)   # maps from R^d to R^(d-1)
    s3 = sheaf_from_graph(g3, 3, rand_map)
    @test issymmetric(Array(sheaf_laplacian_matrix(s3)))

    # ── single-vertex sheaf (no edges) ────────────────────────────────────
    s0 = EuclideanSheaf{Float64}([3])
    L0 = sparse(sheaf_laplacian_matrix(s0))
    @test size(L0) == (3, 3)
    @test iszero(L0)
end

@testset "restricted_laplacian_blocks" begin
    g = path_graph(5)
    s = sheaf_from_graph(g, 1, d -> Matrix{Float64}(I, d, d); symmetric_edges=true)

    L = Array(sheaf_laplacian_matrix(s))
    offsets = [0; cumsum(s.vertex_stalks)]

    interior = [2, 3, 4]
    boundary = [1, 5]

    I_idx = vcat([offsets[v]+1:offsets[v+1] for v in interior]...)
    B_idx = vcat([offsets[v]+1:offsets[v+1] for v in boundary]...)

    L_II_ref = L[I_idx, I_idx]
    L_IB_ref = L[I_idx, B_idx]

    L_II, L_IB = restricted_laplacian_blocks(s, interior, boundary)

    # 1. Return types are sparse
    @test L_II isa SparseMatrixCSC
    @test L_IB isa SparseMatrixCSC

    # 2. Sizes
    @test size(L_II) == (length(I_idx), length(I_idx))
    @test size(L_IB) == (length(I_idx), length(B_idx))

    # 3. Numerically correct
    @test Array(L_II) ≈ L_II_ref
    @test Array(L_IB) ≈ L_IB_ref

    # ── 2-D stalk sheaf ────────────────────────────────────────────────────
    s2 = sheaf_from_graph(g, 2, d -> Matrix{Float64}(I, d, d); symmetric_edges=true)
    L2 = Array(sheaf_laplacian_matrix(s2))
    offsets2 = [0; cumsum(s2.vertex_stalks)]
    I_idx2 = vcat([offsets2[v]+1:offsets2[v+1] for v in interior]...)
    B_idx2 = vcat([offsets2[v]+1:offsets2[v+1] for v in boundary]...)
    L2_II, L2_IB = restricted_laplacian_blocks(s2, interior, boundary)
    @test Array(L2_II) ≈ L2[I_idx2, I_idx2]
    @test Array(L2_IB) ≈ L2[I_idx2, B_idx2]

    # ── random restriction maps (complete graph, 4 vertices) ──────────────
    g3 = complete_graph(4)
    s3 = sheaf_from_graph(g3, 3, d -> randn(d - 1, d))
    L3 = Array(sheaf_laplacian_matrix(s3))
    off3 = [0; cumsum(s3.vertex_stalks)]
    int3 = [1, 2]; bnd3 = [3, 4]
    I3 = vcat([off3[v]+1:off3[v+1] for v in int3]...)
    B3 = vcat([off3[v]+1:off3[v+1] for v in bnd3]...)
    L3_II, L3_IB = restricted_laplacian_blocks(s3, int3, bnd3)
    @test Array(L3_II) ≈ L3[I3, I3]
    @test Array(L3_IB) ≈ L3[I3, B3]

    # ── empty boundary ────────────────────────────────────────────────────
    L_II_eb, L_IB_eb = restricted_laplacian_blocks(s, 1:5, Int[])
    @test Array(L_II_eb) ≈ L
    @test size(L_IB_eb) == (5, 0)

    # ── empty interior ────────────────────────────────────────────────────
    L_II_ei, L_IB_ei = restricted_laplacian_blocks(s, Int[], 1:5)
    @test size(L_II_ei) == (0, 0)
    @test size(L_IB_ei) == (0, 5)
end
