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
end
