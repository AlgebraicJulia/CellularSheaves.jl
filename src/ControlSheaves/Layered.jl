"""
    Layered

Module for multi-agent layered escort sheaves, pushforward solvers, and hierarchical reference generators.
"""
module Layered

export build_multilayer_sheaf,
       build_homomorphism,
       solve_high_level_harmonic,
       solve_mid_level_harmonic,
       solve_direct_harmonic,
       compute_formation_radius_error,
       animate_comprehensive_escort

function animate_comprehensive_escort end

using LinearAlgebra
using Graphs
using ...NetworkSheaves: vertex_stalks, get_restriction_map, EuclideanSheaf
using ...NetworkSheaves.EuclideanSheaves: add_sheaf_edge!, _harmonic_extension_restricted_laplacian
using ...NetworkSheaves.Formations: build_escort_ring, se3_translation_matrix
using ...NetworkSheaves.GraphHomomorphisms: GraphHomomorphism
using ...NetworkSheaves.Pushforwards: pushforward_sheaf, all_fiber_bases

"""
    build_multilayer_sheaf(na1::Int, na2::Int, na3::Int, r1::Float64, r2::Float64, r3::Float64) -> EuclideanSheaf

Construct a 17-node cellular sheaf `F` for a dual-target multi-agent escort mission:
- Nodes 1..na1: Escort Ring 1 around Target 1 (node 16)
- Nodes (na1+1)..(na1+na2): Escort Ring 2 around Target 2 (node 17)
- Nodes (na1+na2+1)..(na1+na2+na3): Support Ring 3 observing Target 1 and Target 2 directly
- Nodes 16 (`T1_NODE`) & 17 (`T2_NODE`): Target nodes
"""
function build_multilayer_sheaf(na1::Int, na2::Int, na3::Int, r1::Float64, r2::Float64, r3::Float64)
    TOTAL_NODES = na1 + na2 + na3 + 2
    T1_NODE = na1 + na2 + na3 + 1
    T2_NODE = na1 + na2 + na3 + 2

    F = EuclideanSheaf{Float64}(fill(4, TOTAL_NODES))

    # Ring 1 around Target 1
    ring1 = build_escort_ring(na1, na1 + 1, r1; observers=[1])
    for e in edges(ring1.underlying_graph)
        u_local, v_local = src(e), dst(e)
        u = u_local == na1 + 1 ? T1_NODE : u_local
        v = v_local == na1 + 1 ? T1_NODE : v_local
        add_sheaf_edge!(F, u, v, get_restriction_map(ring1, u_local, v_local), get_restriction_map(ring1, v_local, u_local))
    end

    # Ring 2 around Target 2
    ring2 = build_escort_ring(na2, na2 + 1, r2; observers=[1])
    for e in edges(ring2.underlying_graph)
        u_local, v_local = src(e), dst(e)
        u = u_local == na2 + 1 ? T2_NODE : u_local + na1
        v = v_local == na2 + 1 ? T2_NODE : v_local + na1
        add_sheaf_edge!(F, u, v, get_restriction_map(ring2, u_local, v_local), get_restriction_map(ring2, v_local, u_local))
    end

    # Ring 3 (Support agents 13..15): Consensus ring with identity restriction maps
    r3_agents = (na1 + na2 + 1):(na1 + na2 + na3)
    for i in 1:na3
        u = r3_agents[i]
        v = r3_agents[i % na3 + 1]
        add_sheaf_edge!(F, u, v, Matrix{Float64}(I, 4, 4), Matrix{Float64}(I, 4, 4))
    end

    # Support tracking edges: Ring 3 agents observe Target 1 and Target 2 directly
    add_sheaf_edge!(F, 13, T1_NODE, Matrix{Float64}(I, 4, 4), Matrix{Float64}(I, 4, 4))
    add_sheaf_edge!(F, 14, T2_NODE, Matrix{Float64}(I, 4, 4), Matrix{Float64}(I, 4, 4))

    return F
end

"""
    build_homomorphism(na1::Int, na2::Int, na3::Int) -> GraphHomomorphism

Construct the graph homomorphism `f : G -> H` mapping:
- Ring 1 agents + Target 1 -> Vertex 1 of H
- Ring 2 agents + Target 2 -> Vertex 2 of H
- Support Ring 3 agents -> Vertex 3 of H
"""
function build_homomorphism(na1::Int, na2::Int, na3::Int)
    TOTAL_NODES = na1 + na2 + na3 + 2
    T1_NODE = na1 + na2 + na3 + 1
    T2_NODE = na1 + na2 + na3 + 2

    v_map = zeros(Int, TOTAL_NODES)
    v_map[1:na1] .= 1
    v_map[T1_NODE] = 1

    v_map[(na1+1):(na1+na2)] .= 2
    v_map[T2_NODE] = 2

    v_map[(na1+na2+1):(na1+na2+na3)] .= 3

    return GraphHomomorphism(v_map, 3)
end

"""
    solve_high_level_harmonic(PfF, p1, p2, B1_T1, B2_T2)

Solve the harmonic extension on the coarse pushforward sheaf `PfF`.
Converts world target coordinates `p1`, `p2` to `PfF` stalk basis coordinates using `B1_T1` and `B2_T2`.
"""
function solve_high_level_harmonic(PfF, p1, p2, B1_T1, B2_T2)
    q_pf_1 = B1_T1 \ p1
    q_pf_2 = B2_T2 \ p2

    boundary = Dict(1 => q_pf_1, 2 => q_pf_2)
    _, _, H_mat, B_mat = _harmonic_extension_restricted_laplacian(PfF, boundary)
    q_pf_3 = vec(H_mat \ (-Matrix(B_mat) * [q_pf_1; q_pf_2]))

    return [q_pf_1, q_pf_2, q_pf_3]
end

"""
    solve_mid_level_harmonic(q_H, fiber_bases)

Lift high-level pushforward section `q_H` to 4D agent reference states on `G` using `fiber_bases`.
"""
function solve_mid_level_harmonic(q_H, fiber_bases)
    q_G_1 = fiber_bases[1] * q_H[1]
    q_G_2 = fiber_bases[2] * q_H[2]
    q_G_3 = fiber_bases[3] * q_H[3]

    q_agents = zeros(15, 4)
    q_agents[1:6, :]   .= reshape(q_G_1[1:24], 4, 6)'
    q_agents[7:12, :]  .= reshape(q_G_2[1:24], 4, 6)'
    q_agents[13:15, :] .= reshape(q_G_3[1:12], 4, 3)'

    return q_agents
end

"""
    solve_direct_harmonic(F, p1, p2)

Solve harmonic extension directly on `F` with Target 1 pinned to `p1` (node 16) and Target 2 pinned to `p2` (node 17).
"""
function solve_direct_harmonic(F, p1, p2)
    boundary = Dict(16 => p1, 17 => p2)
    _, _, H_mat, B_mat = _harmonic_extension_restricted_laplacian(F, boundary)
    q_interior = H_mat \ (-Matrix(B_mat) * [p1; p2])

    q_full = zeros(17, 4)
    q_full[1:15, :] .= reshape(q_interior, 4, 15)'
    q_full[16, :] .= p1
    q_full[17, :] .= p2

    return q_full[1:15, :]
end

"""
    compute_formation_radius_error(agent_positions, centroid, r_expected)

Compute mean absolute radial error from expected radius `r_expected`.
"""
function compute_formation_radius_error(agent_positions, centroid, r_expected)
    N = size(agent_positions, 1)
    if N == 0
        return 0.0
    end
    err_sum = sum(abs(norm(agent_positions[i, 1:2] - centroid[1:2]) - r_expected) for i in 1:N)
    return err_sum / N
end

end # module Layered
