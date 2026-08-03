using Test
using CellularSheaves
using CellularSheaves.Formations
import CellularSheaves.NetworkSheaves.EuclideanSheaves: _harmonic_extension_restricted_laplacian
using CellularSheaves.NetworkSheaves.DistributedSolve
using CliqueTrees.Multifrontal
using SparseArrays
using LinearAlgebra
using Distributed

@testset "EscortTracking" begin

    # Stalk dimension D = 4 for SE(3) homogeneous coordinates [x, y, z, 1]
    D = 4
    NA = 6
    NT = 1
    TV1 = NA + NT
    r_ring = 1.2
    sheaf = build_escort_ring(NA, TV1, r_ring; observers=[1])

    target1_pos = [0.5, 0.2, 1.5, 1.0]
    boundary0 = Dict(TV1 => target1_pos)

    # 1. Verify restricted Laplacian properties
    _, _, Hraw, LIBraw = _harmonic_extension_restricted_laplacian(sheaf, boundary0)
    H = Matrix(Hraw)
    LIB = Matrix(LIBraw)

    @test size(H) == (24, 24)
    @test rank(H) == 24
    @test isapprox(H, H')

    # 2. Verify spatial escort ring geometry from harmonic extension
    qstar_full = H \ (-LIB * target1_pos)
    t1 = target1_pos[1:3]

    for i in 1:6
        pos_i = qstar_full[(i-1)*D+1:(i-1)*D+3]
        dist = norm(pos_i - t1)
        ang = atan(pos_i[2] - t1[2], pos_i[1] - t1[1])
        
        # Verify distance is exactly 1.2m
        @test dist ≈ 1.2
        
        # Verify angle matches expected hexagonal slot (within tolerance)
        expected_angle = (i - 1) * 2π / 6
        @test cos(ang) ≈ cos(expected_angle) atol=1e-5
        @test sin(ang) ≈ sin(expected_angle) atol=1e-5
        
        # Verify homogeneous normalization scale is exactly 1.0
        @test qstar_full[i*D] ≈ 1.0
    end

    # 3. Verify distributed clique-tree solve vs centralized solve
    F = cholesky!(ChordalCholesky(sparse(H)), NoPivot())
    Lfac = F.L
    partition = partition_tree(Lfac, 3)
    nchunk = length(partition.chunks)

    # Spawn workers to test distributed solving
    pids = addprocs(nchunk; exeflags = ["--project=$(pkgdir(CellularSheaves))"])
    try
        @everywhere pids using CellularSheaves
        
        rhs = Vector(-LIB * target1_pos)
        rhs_p = Vector(F.P' \ rhs)
        
        y_sol = distributed_tree_solve(Lfac, rhs_p, nchunk; pids = pids)
        qstar_distributed = F.P \ y_sol

        @test qstar_full ≈ qstar_distributed atol=1e-12
    finally
        rmprocs(pids)
    end

end
