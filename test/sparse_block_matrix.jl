using Test
using JET
using LinearAlgebra
using CellularSheaves.Sheaves
using SparseArrays


function random_block_sparse(nout, nvtx, rowrng, colrng, density)
    outsize = rand(rowrng, nout)
    vtxsize = rand(colrng, nvtx)

    I = Int[]
    J = Int[]
    V = Matrix{Float64}[]

    for j in 1:nvtx
        for i in 1:nout
            if rand() < density
                push!(I, i)
                push!(J, j)
                push!(V, randn(outsize[i], vtxsize[j]))
            end
        end
    end

    return sparseblock(I, J, V, nout, nvtx)
end

function random_block_sparse_tri(nvtx, blkrng, density, uplo)
    vtxsize = rand(blkrng, nvtx)

    I = Int[]
    J = Int[]
    V = Matrix{Float64}[]

    for j in 1:nvtx
        if uplo == :L
            rng = j:nvtx
        else
            rng = 1:j
        end

        for i in rng
            if i == j || rand() < density
                push!(I, i)
                push!(J, j)
                push!(V, randn(vtxsize[i], vtxsize[j]))
            end
        end
    end

    return sparseblock(I, J, V, nvtx, nvtx)
end

@testset "SparseBlockMatrix" begin
    @testset "gemv" begin
        A = random_block_sparse(75, 125, 1:10, 1:10, 0.1)
        B = sparse(A)

        x = randn(size(A, 2))
        y = randn(size(A, 1))

        @test           A  * x ≈           B  * x
        @test   adjoint(A) * y ≈   adjoint(B) * y
        @test transpose(A) * y ≈ transpose(B) * y

        @test_opt           A  * x
        @test_opt   adjoint(A) * y
        @test_opt transpose(A) * y

        @test_call           A  * x
        @test_call   adjoint(A) * y
        @test_call transpose(A) * y
    end

    @testset "gemm" begin
        A = random_block_sparse(75, 125, 1:10, 1:10, 0.1)
        B = sparse(A)

        X  = randn(size(A, 2), 10)
        Y  = randn(size(A, 1), 10)
        Xt = randn(10, size(A, 2))
        Yt = randn(10, size(A, 1))

        @test           A  *           X   ≈           B  *           X
        @test           A  *   adjoint(Xt) ≈           B  *   adjoint(Xt)
        @test           A  * transpose(Xt) ≈           B  * transpose(Xt)
        @test   adjoint(A) *           Y   ≈   adjoint(B) *           Y
        @test   adjoint(A) *   adjoint(Yt) ≈   adjoint(B) *   adjoint(Yt)
        @test   adjoint(A) * transpose(Yt) ≈   adjoint(B) * transpose(Yt)
        @test transpose(A) *           Y   ≈ transpose(B) *           Y
        @test transpose(A) *   adjoint(Yt) ≈ transpose(B) *   adjoint(Yt)
        @test transpose(A) * transpose(Yt) ≈ transpose(B) * transpose(Yt)

        @test_opt           A  *           X
        @test_opt           A  *   adjoint(Xt)
        @test_opt           A  * transpose(Xt)
        @test_opt   adjoint(A) *           Y
        @test_opt   adjoint(A) *   adjoint(Yt)
        @test_opt   adjoint(A) * transpose(Yt)
        @test_opt transpose(A) *           Y
        @test_opt transpose(A) *   adjoint(Yt)
        @test_opt transpose(A) * transpose(Yt)

        @test_call           A  *           X
        @test_call           A  *   adjoint(Xt)
        @test_call           A  * transpose(Xt)
        @test_call   adjoint(A) *           Y
        @test_call   adjoint(A) *   adjoint(Yt)
        @test_call   adjoint(A) * transpose(Yt)
        @test_call transpose(A) *           Y
        @test_call transpose(A) *   adjoint(Yt)
        @test_call transpose(A) * transpose(Yt)
    end

    @testset "symv" begin
        for uplo in (:L, :U)
            T = random_block_sparse_tri(100, 1:10, 0.1, uplo)
            A = Symmetric(       T,  uplo)
            B = Symmetric(sparse(T), uplo)

            x = randn(size(A, 2))

            @test A * x ≈ B * x

            @test_opt  A * x

            @test_call A * x
        end
    end

    @testset "symm" begin
        for uplo in (:L, :U)
            T = random_block_sparse_tri(100, 1:10, 0.1, uplo)
            A = Symmetric(       T,  uplo)
            B = Symmetric(sparse(T), uplo)

            X  = randn(size(A, 2), 10)
            Xt = randn(10, size(A, 2))

            @test A *           X   ≈ B *           X
            @test A *   adjoint(Xt) ≈ B *   adjoint(Xt)
            @test A * transpose(Xt) ≈ B * transpose(Xt)

            @test_opt  A *           X
            @test_opt  A *   adjoint(Xt)
            @test_opt  A * transpose(Xt)

            @test_call A *           X
            @test_call A *   adjoint(Xt)
            @test_call A * transpose(Xt)
        end
    end

    @testset "getindex" begin
        A = random_block_sparse(75, 125, 1:10, 1:10, 0.1)
        B = sparse(A)
        @test all(A .== B)
    end

    @testset "submatrix" begin
        A = random_block_sparse(75, 125, 1:10, 1:10, 0.1)
        B = sparse(A)

        ncol1 = ncols(A, 1); ncol2 = ncols(A, nvtxs(A))
        nrow1 = nrows(A, 1); nrow2 = nrows(A, nouts(A))

        outmask = trues(nouts(A)); outmask[begin] = outmask[end] = false
        vtxmask = trues(nvtxs(A)); vtxmask[begin] = vtxmask[end] = false

        C = submatrix(A, outmask, vtxmask)

        @test sparse(C) == A[nrow1 + 1:end - nrow2, ncol1 + 1:end - ncol2]

        @test_opt  submatrix(A, outmask, vtxmask)

        @test_call submatrix(A, outmask, vtxmask)
    end
end
