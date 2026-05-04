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

This module provides [`SheafSpan`](@ref) and [`pushout_sheaf`](@ref).
"""
module Pushouts

export pushout_sheaf, SheafSpan

using LinearAlgebra
using Graphs

using ..SheafInterface: vertex_stalks, edge_stalk_dimensions, underlying_graph,
                        get_restriction_map, add_vertex_stalk!, add_sheaf_edge!
using ..EuclideanSheaves: EuclideanSheaf
using ..SheafMorphisms: SheafMorphism

# ---------------------------------------------------------------------------
# SheafSpan
# ---------------------------------------------------------------------------

"""
    SheafSpan

A *span* of sheaf morphisms

    F ← K → G

where `left : K → F` and `right : K → G` share the same domain `apex = K`.

The inner constructor validates that the morphism index maps have the correct
lengths relative to the apex's vertex and edge counts, and that all `Vmap`
entries are valid indices into the respective feet.  An `AssertionError` is
thrown if any check fails.

Fields
------
- `left  :: SheafMorphism` — morphism `K → F`.
- `right :: SheafMorphism` — morphism `K → G`.
- `apex  :: EuclideanSheaf` — the apex sheaf `K`.
- `F     :: EuclideanSheaf` — left foot.
- `G     :: EuclideanSheaf` — right foot.
"""
struct SheafSpan
    left::SheafMorphism
    right::SheafMorphism
    apex::EuclideanSheaf
    F::EuclideanSheaf
    G::EuclideanSheaf

    function SheafSpan(left::SheafMorphism, right::SheafMorphism,
                       apex::EuclideanSheaf, F::EuclideanSheaf,
                       G::EuclideanSheaf)
        nv_K = length(vertex_stalks(apex))
        ne_K = length(edge_stalk_dimensions(apex))
        nv_F = length(vertex_stalks(F))
        nv_G = length(vertex_stalks(G))
    k_graph = underlying_graph(K)
    f_graph = underlying_graph(F)
    g_graph = underlying_graph(G)

    if f_graph != k_graph || g_graph != k_graph
        throw(ArgumentError("pushout_sheaf requires K, F, and G to share the same underlying graph"))
    end

    nv_k = nv(k_graph)
    ne_k = ne(k_graph)
    if nv(f_graph) != nv_k || nv(g_graph) != nv_k || ne(f_graph) != ne_k || ne(g_graph) != ne_k
        throw(ArgumentError("pushout_sheaf requires K, F, and G to have matching numbers of vertices and edges"))
    end

    nv_g    = nv_k
    ne_g    = ne_k
    edge_list = collect(edges(k_graph))
        @assert all(1 .<= left.Vmap  .<= nv_F) "left.Vmap entries must be valid vertex indices of F"
        @assert all(1 .<= right.Vmap .<= nv_G) "right.Vmap entries must be valid vertex indices of G"
        new(left, right, apex, F, G)
    end
end

# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

# Compute the left null space (cokernel projection basis) of `A`.
# Returns an AbstractMatrix whose columns are an orthonormal basis for ker(A').
# When A has zero columns (trivial map), returns a lazy Diagonal identity to
# avoid materialising a dense m×m matrix.
function _cokernel_basis(A::AbstractMatrix{T}) where T
    m, n = size(A)
    if n == 0
        return Diagonal(ones(T, m))
    end
    return Matrix{T}(nullspace(A'))
end

# Compute cokernel data for a single stalk.
#
# Given linear maps  fmap : K -> F  and  gmap : K -> G  (both with the same
# domain K, i.e. same number of columns), this function computes:
#
#   Q   -- AbstractMatrix, size (dimF+dimG) x p, orthonormal basis for the
#          cokernel of [fmap; -gmap].
#   iF  -- Matrix{T}, size p x dimF, injection map F -> P.
#   iG  -- Matrix{T}, size p x dimG, injection map G -> P.
#   p   -- Int, pushout stalk dimension = dim cokernel of [fmap; -gmap].
#
# When dimK == 0 (fmap and gmap both have 0 columns), the cokernel is all of
# F (+) G, so p = dimF + dimG and the injections are the standard embeddings.
#
# Note: _cokernel_basis returns a lazy Diagonal for the trivial (dimK==0)
# case to avoid an O(m^2) allocation there.  The conversion to dense via
# Matrix{T}(Q') here is unavoidable for the subsequent matrix multiplications
# that build iF and iG, so the restriction maps are always dense.
function _cokernel_stalk(::Type{T}, dimF::Int, dimG::Int,
                          fmap::AbstractMatrix, gmap::AbstractMatrix) where T
    A  = vcat(fmap, -gmap)          # (dimF + dimG) x dimK
    Q  = _cokernel_basis(A)
    p  = size(Q, 2)                 # pushout stalk dimension
    Qt = Matrix{T}(Q')              # ensure dense for subsequent multiplications
    iF = Matrix{T}(Qt * vcat(Matrix{T}(I, dimF, dimF), zeros(T, dimG, dimF)))
    iG = Matrix{T}(Qt * vcat(zeros(T, dimF, dimG), Matrix{T}(I, dimG, dimG)))
    return Q, iF, iG, p
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
- `iF : F -> P` is the canonical injection from `F`.
- `iG : G -> P` is the canonical injection from `G`.

The universal property holds: `iF ∘ f ≈ iG ∘ g` (up to numerical tolerance).

All three sheaves `K`, `F`, `G` must share the same underlying graph.
"""
function pushout_sheaf(f::SheafMorphism, g::SheafMorphism,
                       K::EuclideanSheaf{T}, F::EuclideanSheaf{T},
                       G::EuclideanSheaf{T}) where T
    g_graph   = underlying_graph(K)
    nv_g      = nv(g_graph)
    ne_g      = ne(g_graph)
    edge_list = collect(edges(g_graph))

    vdimsF = vertex_stalks(F)
    vdimsG = vertex_stalks(G)
    edimsF = edge_stalk_dimensions(F)
    edimsG = edge_stalk_dimensions(G)

    # ------------------------------------------------------------------
    # Vertex stalks: compute cokernel at each vertex
    # ------------------------------------------------------------------
    # Q_vertex[v]  : AbstractMatrix, (dimF_v + dimG_v) x p_v, cokernel basis
    # iF_vertex[v] : p_v x dimF_v injection  F(v) -> P(v)
    # iG_vertex[v] : p_v x dimG_v injection  G(v) -> P(v)
    # p_vdims[v]   : pushout vertex stalk dimension at v
    Q_vertex  = Vector{AbstractMatrix{T}}(undef, nv_g)
    iF_vertex = Vector{Matrix{T}}(undef, nv_g)
    iG_vertex = Vector{Matrix{T}}(undef, nv_g)
    p_vdims   = Vector{Int}(undef, nv_g)

    for v in 1:nv_g
        Q_vertex[v], iF_vertex[v], iG_vertex[v], p_vdims[v] =
            _cokernel_stalk(T, vdimsF[f.Vmap[v]], vdimsG[g.Vmap[v]],
                            f.Vmaps[v], g.Vmaps[v])
    end

    # ------------------------------------------------------------------
    # Edge stalks: compute cokernel at each edge
    # ------------------------------------------------------------------
    Q_edge  = Vector{AbstractMatrix{T}}(undef, ne_g)
    iF_edge = Vector{Matrix{T}}(undef, ne_g)
    iG_edge = Vector{Matrix{T}}(undef, ne_g)
    p_edims = Vector{Int}(undef, ne_g)

    for k in 1:ne_g
        Q_edge[k], iF_edge[k], iG_edge[k], p_edims[k] =
            _cokernel_stalk(T, edimsF[f.Emap[k]], edimsG[g.Emap[k]],
                            f.Emaps[k], g.Emaps[k])
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

        rm_F_v1 = get_restriction_map(F, v1, v2)  # dimF_e x dimF_v1
        rm_F_v2 = get_restriction_map(F, v2, v1)  # dimF_e x dimF_v2
        rm_G_v1 = get_restriction_map(G, v1, v2)  # dimG_e x dimG_v1
        rm_G_v2 = get_restriction_map(G, v2, v1)  # dimG_e x dimG_v2

        # Induced restriction maps of P satisfy BOTH naturality conditions:
        #   rho_P(v->e) * iF_v = iF_e * rho_F(v->e)
        #   rho_P(v->e) * iG_v = iG_e * rho_G(v->e)
        #
        # Since [iF_v | iG_v] = Q_v' has orthonormal rows, the unique solution is:
        #   rho_P(v->e) = [iF_e * rho_F(v->e), iG_e * rho_G(v->e)] * Q_v
        rm_P_v1 = Matrix{T}(hcat(iF_edge[k] * rm_F_v1, iG_edge[k] * rm_G_v1) *
                             Q_vertex[v1])
        rm_P_v2 = Matrix{T}(hcat(iF_edge[k] * rm_F_v2, iG_edge[k] * rm_G_v2) *
                             Q_vertex[v2])

        add_sheaf_edge!(P, v1, v2, rm_P_v1, rm_P_v2)
    end

    # ------------------------------------------------------------------
    # Return injection morphisms iF : F -> P  and  iG : G -> P
    # ------------------------------------------------------------------
    iF_morph = SheafMorphism(collect(1:nv_g), collect(1:ne_g),
                             iF_vertex, iF_edge)
    iG_morph = SheafMorphism(collect(1:nv_g), collect(1:ne_g),
                             iG_vertex, iG_edge)

    return (P, iF_morph, iG_morph)
end

"""
    pushout_sheaf(span::SheafSpan) -> (EuclideanSheaf, SheafMorphism, SheafMorphism)

Convenience overload: compute the pushout of a [`SheafSpan`](@ref).
"""
function pushout_sheaf(span::SheafSpan)
    return pushout_sheaf(span.left, span.right, span.apex, span.F, span.G)
end

end # module Pushouts
