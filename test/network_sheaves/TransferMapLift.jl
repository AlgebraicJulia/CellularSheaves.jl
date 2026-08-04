using Test
using CellularSheaves
using CellularSheaves.Formations
using CellularSheaves.Pushforwards
using CellularSheaves.GraphHomomorphisms
import CellularSheaves.NetworkSheaves.EuclideanSheaves: _harmonic_extension_restricted_laplacian
using Graphs
using LinearAlgebra

@testset "Transfer Map Pseudoinverse Pullback Verification (T⁺)" begin
    TOTAL_NODES = 17
    F = EuclideanSheaf{Float64}(fill(4, TOTAL_NODES))
    
    # Ring 1 (agents 1..6) around Target 1 (node 16)
    ring1 = build_escort_ring(6, 7, 0.3; observers=[1])
    for e in edges(ring1.underlying_graph)
        u = src(e) == 7 ? 16 : src(e)
        v = dst(e) == 7 ? 16 : dst(e)
        add_sheaf_edge!(F, u, v, get_restriction_map(ring1, src(e), dst(e)), get_restriction_map(ring1, dst(e), src(e)))
    end

    # Ring 2 (agents 7..12) around Target 2 (node 17)
    ring2 = build_escort_ring(6, 7, 0.36; observers=[1])
    for e in edges(ring2.underlying_graph)
        u = src(e) == 7 ? 17 : src(e) + 6
        v = dst(e) == 7 ? 17 : dst(e) + 6
        add_sheaf_edge!(F, u, v, get_restriction_map(ring2, src(e), dst(e)), get_restriction_map(ring2, dst(e), src(e)))
    end

    # Ring 3 (agents 13..15): Ring consensus edges
    for i in 1:3
        add_sheaf_edge!(F, 12+i, 12+(i%3)+1, Matrix{Float64}(I, 4, 4), Matrix{Float64}(I, 4, 4))
    end
    add_sheaf_edge!(F, 13, 16, Matrix{Float64}(I, 4, 4), Matrix{Float64}(I, 4, 4))
    add_sheaf_edge!(F, 14, 17, Matrix{Float64}(I, 4, 4), Matrix{Float64}(I, 4, 4))

    # Graph Homomorphism f
    v_map = zeros(Int, 17)
    v_map[1:6] .= 1; v_map[16] = 1
    v_map[7:12] .= 2; v_map[17] = 2
    v_map[13:15] .= 3
    f = GraphHomomorphism(v_map, 3)

    PfF = pushforward_sheaf(f, F)
    T = Matrix(pushforward_transfer_map(f, F))
    T_pinv = pinv(T)

    # Compute High-Level Harmonic Section q_H
    p1 = [1.0, 0.0, 1.5, 1.0]
    p2 = [-1.0, 0.0, 1.5, 1.0]
    boundary = Dict(1 => p1, 2 => p2)
    _, _, H_mat, B_mat = _harmonic_extension_restricted_laplacian(PfF, boundary)
    q_interior = vec(H_mat \ (-Matrix(B_mat) * [p1; p2]))
    q_H = [p1, p2, q_interior]
    q_H_flat = vcat(q_H...)

    # Lift to Agent Space via Pseudoinverse
    q_G_ref = T_pinv * q_H_flat

    # --- Test 1: Right Inverse Property (T * T⁺ * q_H = q_H) ---
    @test isapprox(T * q_G_ref, q_H_flat, atol=1e-5)

    # --- Test 2: Orthogonality to Transfer Map Nullspace ---
    N = nullspace(T)
    if size(N, 2) > 0
        @test norm(N' * q_G_ref) < 1e-5
    end

    # --- Test 3: Stalk Dimensionality & Agent Reference Matrix ---
    q_agents = Matrix(reshape(q_G_ref[1:15*4], 4, 15)')
    @test size(q_agents) == (15, 4)
end
