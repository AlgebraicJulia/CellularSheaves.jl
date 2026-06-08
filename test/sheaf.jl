using Test
using JET
using LinearAlgebra
using CellularSheaves.Sheaves

# Opinion Dynamics on Discourse Sheaves
# Hansen and Ghrist
# Figure 2
function test_sheaf()
    I = [2, 1, 3, 2, 4, 3, 1, 4]

    J = [1, 2, 2, 3, 3, 4, 4, 1]

    V = [
        [-2.0;;],
        [-1.0 2.0],
        [1.0 1.0],
        [1.0;;],
        [-1.0;;],
        [1.0 -1.0],
        [0.0 -1.0; 1.0 0.0],
        [0.0; 1.0;;],
    ]

    return sheaf(I, J, V, 4)
end

@testset "Sheaf" begin
    S = test_sheaf()

    L = [
         5  -2   4   0  -1   0
        -2   2  -1  -1   0   0
         4  -1   5  -1   0   0
         0  -1  -1   2   1  -1
        -1   0   0   1   2  -1
         0   0   0  -1  -1   2
    ]

    @testset "coboundary" begin
        δ = coboundary(S)

        # Verify δᵀδ = L (edge ordering may differ from paper)
        @test δ' * δ == L

        @test_opt  coboundary(S)
        @test_call coboundary(S)
    end

    @testset "laplacian" begin
        @test laplacian(S) == L

        @test_opt  laplacian(S)
        @test_call laplacian(S)
    end

    @testset "trilaplacian" begin
        @test Symmetric(trilaplacian(S, :U), :U) == L
        @test Symmetric(trilaplacian(S, :L), :L) == L

        @test_opt  trilaplacian(S, :U)
        @test_opt  trilaplacian(S, :L)

        @test_call trilaplacian(S, :U)
        @test_call trilaplacian(S, :L)
    end
end
