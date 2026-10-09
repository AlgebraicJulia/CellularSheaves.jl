using Test
using CellularSheaves
using LinearAlgebra
using SparseArrays
using Graphs
const LinearOperator = CellularSheaves.NetworkSheaves.RestrictionMaps.LinearOperator

@testset "RestrictionMaps" begin
    M = [1.0 2.0 0.0; 0.0 0.0 3.0]
    x = [1.0, -1.0, 2.0]
    y = [0.5, 2.0]

    dense = DenseRestriction(M)
    sparse_map = SparseRestriction(sparse(M))
    selection = SelectionRestriction{Float64}([3, 1], 3)
    matrix_free = FunctionRestriction{Float64}((out, v) -> mul!(out, M, v), (out, v) -> mul!(out, M', v), 2, 3)

    @testset "interface" begin
        for R in (dense, sparse_map, matrix_free)
            @test size(R) == (2, 3) && size(R, 1) == 2 && eltype(R) == Float64
            @test R * x ≈ M * x
            @test R' * y ≈ M' * y
            @test Matrix(R) ≈ M
            @test sparse(R) ≈ sparse(M)
            @test R * [x 2x] ≈ M * [x 2x]
            @test [y 2y]' * R ≈ [y 2y]' * M
            @test LinearOperator(R) * x ≈ M * x
            @test LinearOperator(R)' * y ≈ M' * y
        end
        @test selection * x == [2.0, 1.0]
        @test selection' * y == [2.0, 0.0, 0.5]
        @test Matrix(selection) == [0 0 1; 1 0 0]
        @test sparse(selection) == sparse([0.0 0 1; 1 0 0])
        @test SelectionRestriction{Float64}([3, 1], 3) == selection
        @test_throws ArgumentError SelectionRestriction{Float64}([4], 3)
        @test restriction_map(M) isa DenseRestriction{Float64}
        @test restriction_map(sparse(M)) isa SparseRestriction{Float64}
        @test restriction_map(selection) === selection
        @test occursin("2×3", sprint(show, dense))
    end

    # The same sheaf on a triangle, with maps of every kind.
    function triangle(::Type{M}, maps) where {M}
        s = EuclideanSheaf{Float64,M}([3, 3, 2])
        add_sheaf_edge!(s, 1, 2, maps[1], maps[2])
        add_sheaf_edge!(s, 2, 3, maps[3], maps[4])
        add_sheaf_edge!(s, 1, 3, maps[5], maps[6])
        return s
    end
    B = [1.0 0.0; 0.0 1.0]
    C = [0.0 1.0 0.0; 1.0 0.0 1.0]
    dense_sheaf = triangle(Matrix{Float64}, [M, C, C, B, M, 2B])
    mixed = triangle(AbstractRestrictionMap{Float64},
                     [matrix_free, SparseRestriction(sparse(C)), DenseRestriction(C), B, sparse_map, 2B])

    @testset "sheaves with restriction-map storage" begin
        @test dense_sheaf isa DenseEuclideanSheaf{Float64}
        @test !(mixed isa DenseEuclideanSheaf)
        @test get_restriction_map(mixed, 1, 2) === matrix_free
        @test Matrix(coboundary_map(mixed)) ≈ Matrix(coboundary_map(dense_sheaf))
        d = coboundary_operator(mixed)
        x0 = randn(8)
        y0 = randn(6)
        @test d * x0 ≈ Matrix(coboundary_map(dense_sheaf)) * x0
        @test d' * y0 ≈ Matrix(coboundary_map(dense_sheaf))' * y0
        @test coboundary_operator(dense_sheaf) * x0 ≈ d * x0
        @test Matrix(sheaf_laplacian_matrix(mixed)) ≈ Matrix(sheaf_laplacian_matrix(dense_sheaf))
        L_II, L_IB = restricted_laplacian_blocks(mixed, [1, 2], [3])
        D_II, D_IB = restricted_laplacian_blocks(dense_sheaf, [1, 2], [3])
        @test L_II ≈ D_II && L_IB ≈ D_IB
        @test size(nullspace_ldlt(mixed), 2) == size(nullspace_ldlt(dense_sheaf), 2)

        # A dense sheaf accepts restriction maps by materializing them.
        s = EuclideanSheaf{Float64}([3, 3])
        add_sheaf_edge!(s, 1, 2, selection, matrix_free)
        @test get_restriction_map(s, 1, 2) == Matrix(selection)
        @test get_restriction_map(s, 2, 1) ≈ M

        # Selection maps stay index lists.
        sel = EuclideanSheaf{Float64,SelectionRestriction{Float64}}([3, 3])
        add_sheaf_edge!(sel, 1, 2, selection, SelectionRestriction{Float64}([1, 2], 3))
        @test get_restriction_map(sel, 1, 2) isa SelectionRestriction{Float64}
        @test coboundary_map(sel) isa SparseMatrixCSC
        @test nnz(coboundary_map(sel)) == 4
    end
end
