"""
Module for the categorical pushout of cellular sheaves over a fixed graph.

Given a span of sheaf morphisms

    K --f--> F
    |
    g
    |
    v
    G

the *pushout* `P = F ⊕_K G` is the unique sheaf (up to isomorphism) equipped
with morphisms `iF : F → P` and `iG : G → P` such that `iF ∘ f = iG ∘ g`,
universal among all such cones.

For network sheaves over a fixed graph the pushout is computed stalkwise as a
cokernel:

- At each vertex `v`:  `P(v) = coker((f_v, -g_v) : K(v) → F(v) ⊕ G(v))`
- At each edge `e`:    `P(e) = coker((f_e, -g_e) : K(e) → F(e) ⊕ G(e))`
- Restriction maps of `P` are induced from those of `F` and `G` simultaneously
  via the naturality conditions for `iF` and `iG`.

This module provides [`pushout_sheaf`](@ref).
"""
module Pushouts

export pushout_sheaf

using LinearAlgebra
using Graphs

using ..SheafInterface: vertex_stalks, edge_stalk_dimensions, underlying_graph,
                        get_restriction_map, add_vertex_stalk!, add_sheaf_edge!
using ..EuclideanSheaves: EuclideanSheaf
using ..SheafMorphisms: SheafMorphism

# ---------------------------------------------------------------------------
# Cokernel helper
# ---------------------------------------------------------------------------

# Compute the left null space (cokernel projection basis) of `A`.
# Returns Q :: (m × p) with orthonormal columns spanning ker(A').
# When A has zero columns (empty map), Q = I (full space).
function _cokernel_basis(A::AbstractMatrix{T}) where T
    m, n = size(A)
    if n == 0
        # Trivial: cokernel is the whole target space.
        return Matrix{T}(I, m, m)
    end
    return Matrix{T}(nullspace(A'))
end

# ---------------------------------------------------------------------------
# pushout_sheaf
# ---------------------------------------------------------------------------

"""
    pushout_sheaf(f::SheafMorphism, g::SheafMorphism, K, F, G) -> (EuclideanSheaf, SheafMorphism, SheafMorphism)

Compute the pushout of the span `K --f--> F`, `K --g--> G` in the category of
network sheaves over the shared underlying graph.

Returns `(P, iF, iG)` where:
- `P` is the pushout sheaf.
- `iF : F → P` is the canonical injection from `F`.
- `iG : G → P` is the canonical injection from `G`.

The universal property holds: `iF ∘ f ≈ iG ∘ g` (up to numerical tolerance).

All three sheaves `K`, `F`, `G` must share the same underlying graph.
"""
function pushout_sheaf(f::SheafMorphism, g::SheafMorphism,
                       K::EuclideanSheaf{T}, F::EuclideanSheaf{T},
                       G::EuclideanSheaf{T}) where T
    g_graph = underlying_graph(K)
    nv_g    = nv(g_graph)
    ne_g    = ne(g_graph)
    edge_list = collect(edges(g_graph))

    vdimsK = vertex_stalks(K)
    vdimsF = vertex_stalks(F)
    vdimsG = vertex_stalks(G)

    edimsK = edge_stalk_dimensions(K)
    edimsF = edge_stalk_dimensions(F)
    edimsG = edge_stalk_dimensions(G)

    # ------------------------------------------------------------------
    # Vertex stalks: compute cokernel at each vertex
    # ------------------------------------------------------------------
    # Q_vertex[v]  : (dimF_v + dimG_v) × p_v  orthonormal cokernel basis
    # iF_vertex[v] : p_v × dimF_v  injection  F(v) → P(v)
    # iG_vertex[v] : p_v × dimG_v  injection  G(v) → P(v)
    Q_vertex  = Vector{Matrix{T}}(undef, nv_g)
    iF_vertex = Vector{Matrix{T}}(undef, nv_g)
    iG_vertex = Vector{Matrix{T}}(undef, nv_g)
    p_vdims   = Vector{Int}(undef, nv_g)

    for v in 1:nv_g
        dimK_v = vdimsK[v]
        dimF_v = vdimsF[f.Vmap[v]]
        dimG_v = vdimsG[g.Vmap[v]]

        # A_v : K(v) → F(v) ⊕ G(v)
        A_v = dimK_v == 0 ?
              zeros(T, dimF_v + dimG_v, 0) :
              vcat(f.Vmaps[v], -g.Vmaps[v])

        # Q_v : orthonormal basis for cokernel of A_v
        Q_v = _cokernel_basis(A_v)
        p_v = size(Q_v, 2)

        # Injection maps:  iF_v = Q_v' * [I; 0],  iG_v = Q_v' * [0; I]
        iF_vertex[v] = Q_v' * vcat(Matrix{T}(I, dimF_v, dimF_v),
                                    zeros(T, dimG_v, dimF_v))
        iG_vertex[v] = Q_v' * vcat(zeros(T, dimF_v, dimG_v),
                                    Matrix{T}(I, dimG_v, dimG_v))
        Q_vertex[v]  = Q_v
        p_vdims[v]   = p_v
    end

    # ------------------------------------------------------------------
    # Edge stalks: compute cokernel at each edge
    # ------------------------------------------------------------------
    Q_edge  = Vector{Matrix{T}}(undef, ne_g)
    iF_edge = Vector{Matrix{T}}(undef, ne_g)
    iG_edge = Vector{Matrix{T}}(undef, ne_g)
    p_edims = Vector{Int}(undef, ne_g)

    for k in 1:ne_g
        dimK_e = edimsK[k]
        dimF_e = edimsF[f.Emap[k]]
        dimG_e = edimsG[g.Emap[k]]

        A_e = dimK_e == 0 ?
              zeros(T, dimF_e + dimG_e, 0) :
              vcat(f.Emaps[k], -g.Emaps[k])

        Q_e = _cokernel_basis(A_e)
        p_e = size(Q_e, 2)

        iF_edge[k] = Q_e' * vcat(Matrix{T}(I, dimF_e, dimF_e),
                                   zeros(T, dimG_e, dimF_e))
        iG_edge[k] = Q_e' * vcat(zeros(T, dimF_e, dimG_e),
                                   Matrix{T}(I, dimG_e, dimG_e))
        Q_edge[k]  = Q_e
        p_edims[k] = p_e
    end

    # ------------------------------------------------------------------
    # Build pushout sheaf P
    # ------------------------------------------------------------------
    P = EuclideanSheaf{T}(Int[])
    for v in 1:nv_g
        add_vertex_stalk!(P, p_vdims[v])
    end

    for (k, e) in enumerate(edge_list)
        v1 = src(e)
        v2 = dst(e)

        # Restriction maps of F and G from each endpoint to the edge
        rm_F_v1 = get_restriction_map(F, v1, v2)  # dimF_e × dimF_v1
        rm_F_v2 = get_restriction_map(F, v2, v1)  # dimF_e × dimF_v2
        rm_G_v1 = get_restriction_map(G, v1, v2)  # dimG_e × dimG_v1
        rm_G_v2 = get_restriction_map(G, v2, v1)  # dimG_e × dimG_v2

        # Induced restriction maps of P satisfy BOTH naturality conditions:
        #   ρ_P(v→e) ∘ iF_v = iF_e ∘ ρ_F(v→e)
        #   ρ_P(v→e) ∘ iG_v = iG_e ∘ ρ_G(v→e)
        #
        # Since [iF_v, iG_v] = Q_v' (orthonormal rows), the unique solution is:
        #   ρ_P(v→e) = [iF_e ∘ ρ_F(v→e), iG_e ∘ ρ_G(v→e)] * Q_v
        #            = [iF_e ∘ ρ_F(v→e), iG_e ∘ ρ_G(v→e)] * pinv([iF_v, iG_v])
        rm_P_v1 = Matrix{T}(hcat(iF_edge[k] * rm_F_v1, iG_edge[k] * rm_G_v1) *
                             Q_vertex[v1])
        rm_P_v2 = Matrix{T}(hcat(iF_edge[k] * rm_F_v2, iG_edge[k] * rm_G_v2) *
                             Q_vertex[v2])

        add_sheaf_edge!(P, v1, v2, rm_P_v1, rm_P_v2)
    end

    # ------------------------------------------------------------------
    # Return injection morphisms iF : F → P  and  iG : G → P
    # ------------------------------------------------------------------
    iF_morph = SheafMorphism(collect(1:nv_g), collect(1:ne_g),
                             iF_vertex, iF_edge)
    iG_morph = SheafMorphism(collect(1:nv_g), collect(1:ne_g),
                             iG_vertex, iG_edge)

    return (P, iF_morph, iG_morph)
end

end # module Pushouts
