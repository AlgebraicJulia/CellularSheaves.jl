using Test
using CellularSheaves.Formations
using LinearAlgebra

@testset "Formations" begin
    @testset "SE(3) Matrices" begin
        # Translation
        d = [1.0, 2.0, 3.0]
        Td = se3_translation_matrix(d)
        @test size(Td) == (4, 4)
        @test Td[1:3, 1:3] == I(3)
        @test Td[1:3, 4] == -d
        @test Td[4, 4] == 1.0

        # Rotation
        θz = π/2
        Rz = se3_rotation_matrix(; θz=θz)
        @test Rz[1, 1] ≈ 0.0 atol=1e-10
        @test Rz[1, 2] ≈ -1.0
        @test Rz[2, 1] ≈ 1.0
        @test Rz[2, 2] ≈ 0.0 atol=1e-10
        
        # Rotations by angles
        R_angles = se3_rotation_matrix(; θz=θz)
        @test Rz ≈ R_angles
        
        # Affine
        T_aff = se3_affine_matrix(Rz[1:3, 1:3], d)
        @test T_aff[1:3, 1:3] ≈ Rz[1:3, 1:3]
        @test T_aff[1:3, 4] == -d
    end

    @testset "Escort Ring Builder" begin
        sheaf = build_escort_ring(6, 7, 1.2; observers=[1, 3])

        @test length(sheaf.vertex_stalks) == 7
        @test sheaf.vertex_stalks[1] == 4

        # Edges: 6 ring edges + 2 pinning edges = 8 edges
        @test length(sheaf.edge_stalks) == 8
    end

    @testset "affine_translation_matrix (arbitrary dimension)" begin
        d2 = [1.0, -2.0]
        T2 = affine_translation_matrix(d2)
        @test size(T2) == (3, 3)
        @test T2[1:2, 1:2] == I(2)
        @test T2[1:2, 3] == -d2
        @test T2[3, 3] == 1.0

        # n=3 case matches se3_translation_matrix
        d3 = [1.0, 2.0, 3.0]
        @test affine_translation_matrix(d3) == se3_translation_matrix(d3)
    end

    @testset "build_escort_ring generalized D/affine" begin
        # D=3, affine=true: 2D translation ring (e.g. planar agents)
        sheaf_planar = build_escort_ring(4, 5, 0.5; D=3, affine=true)
        @test length(sheaf_planar.vertex_stalks) == 5
        @test sheaf_planar.vertex_stalks[1] == 3

        # D=6, affine=true: higher-dimensional homogeneous stalks now succeed
        sheaf_hi = build_escort_ring(4, 5, 0.5; D=6, affine=true)
        @test sheaf_hi.vertex_stalks[1] == 6

        # affine=false: non-affine stalks can only represent a zero-radius (consensus) ring
        sheaf_consensus = build_escort_ring(4, 5, 0.0; D=2, affine=false)
        @test sheaf_consensus.vertex_stalks[1] == 2

        @test_throws Exception build_escort_ring(4, 5, 0.5; D=2, affine=false)

        # observers out of range are rejected
        @test_throws Exception build_escort_ring(4, 5, 0.5; observers=[1, 5])
    end
end
