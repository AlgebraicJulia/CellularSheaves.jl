"""
Module for the pushforward of a cellular sheaf along a graph homomorphism.

Given a graph homomorphism φ : G → H and a cellular sheaf F on G, the
*pushforward* ``\\varphi_* F`` is a cellular sheaf on H whose space of global
sections is isomorphic to that of F.

This module provides:

- [`fiber_section_basis`](@ref) — basis for the global sections of the
  restriction of F to a single fiber.
- [`all_fiber_bases`](@ref) — collects the fiber-section bases for all target
  vertices.
- [`pushforward_sheaf`](@ref) — construct the pushforward sheaf ``\\varphi_* F``.
- [`pushforward_transfer_map`](@ref) — construct the linear map
  ``T : C^0(F) \\to C^0(\\varphi_* F)`` that sends every global section of F to
  the corresponding global section of ``\\varphi_* F``.

See the [`GraphHomomorphisms`](@ref) module for the underlying combinatorial
representation of graph homomorphisms.
"""
module Pushforwards

export fiber_section_basis, all_fiber_bases, pushforward_sheaf, pushforward_transfer_map

using Graphs
using LinearAlgebra
using SparseArrays

using ..SheafInterface: vertex_stalks, coboundary_map, underlying_graph,
                        get_vertex_stalk, get_restriction_map
using ..EuclideanSheaves: EuclideanSheaf, add_vertex_stalk!, add_sheaf_edge!, nullspace_ldlt
using ..GraphHomomorphisms: GraphHomomorphism, fiber_vertices, fiber_edges, cross_edges

# ---------------------------------------------------------------------------
# Fiber-section bases
# ---------------------------------------------------------------------------

"""
    fiber_section_basis(F, verts::Vector{Int}, fedges::Vector{Tuple{Int,Int}}) -> Matrix

Compute a basis for the space of global sections of the sub-sheaf of `F`
restricted to the sub-graph induced by `verts` (with `fedges` as edge set).

Returns an `(total_fiber_stalk_dim × k)` matrix whose `k` columns span the
fiber's global-section space.  The basis is computed via
[`nullspace_ldlt`](@ref) applied to the fiber coboundary operator.

The rows of the returned matrix are ordered to match the concatenation of
the vertex stalks `F(v)` for `v ∈ verts`, in the order given by `verts`.
"""
function fiber_section_basis(F, verts::Vector{Int}, fedges::Vector{Tuple{Int,Int}})
    stalk_dims = [get_vertex_stalk(F, v) for v in verts]
    # Local vertex indices within the fiber sub-sheaf.
    v_idx = Dict(v => i for (i, v) in enumerate(verts))
    Fib = EuclideanSheaf{Float64}(stalk_dims)
    for (a, b) in fedges
        add_sheaf_edge!(Fib, v_idx[a], v_idx[b],
                        get_restriction_map(F, a, b),
                        get_restriction_map(F, b, a))
    end
    d_fib = sparse(coboundary_map(Fib))
    return nullspace_ldlt(d_fib' * d_fib)
end

"""
    all_fiber_bases(hom::GraphHomomorphism, F) -> Vector{Matrix}

Return a `Vector` of fiber-section basis matrices, one per target vertex
of `hom`.  Entry `tv` is the result of
`fiber_section_basis(F, fiber_vertices(hom, tv), fiber_edges(hom, g, tv))`
where `g = underlying_graph(F)`.
"""
function all_fiber_bases(hom::GraphHomomorphism, F)
    g = underlying_graph(F)
    return [
        fiber_section_basis(F,
                            fiber_vertices(hom, tv),
                            fiber_edges(hom, g, tv))
        for tv in 1:hom.n_target
    ]
end

# ---------------------------------------------------------------------------
# Pushforward sheaf
# ---------------------------------------------------------------------------

"""
    pushforward_sheaf(hom::GraphHomomorphism, F) -> EuclideanSheaf

Construct the *pushforward sheaf* ``\\varphi_* F`` of the cellular sheaf `F`
along the graph homomorphism `hom : G → H`.

**Vertex stalks.** The stalk ``(\\varphi_* F)(v)`` at a target vertex `v` is
(isomorphic to) the space of global sections of `F` restricted to the fiber
``\\varphi^{-1}(v)``; a basis for this space is computed via
[`fiber_section_basis`](@ref) (which calls [`nullspace_ldlt`](@ref)).

**Edge stalks and restriction maps.** For each *cross edge* ``e = (u, w)``
(with ``\\varphi(u) = p``, ``\\varphi(w) = q``, ``p \\neq q``) the edge stalk of
``\\varphi_* F`` at ``(p,q)`` equals the edge stalk of `F` at `e`.  The
restriction map from ``(\\varphi_* F)(p)`` to the edge stalk is obtained by
composing `F`'s restriction map at `u` with the fiber-basis embedding matrix:
``\\rho_{p \\to e} = F_{u \\to e} \\cdot B_p[\\text{rows of }u, :]``
where ``B_p`` is the fiber-basis matrix whose columns are expressed in the
concatenated stalk space of ``\\varphi^{-1}(p)``.
"""
function pushforward_sheaf(hom::GraphHomomorphism, F)
    g   = underlying_graph(F)
    phi = hom.vertex_map

    # Compute a basis for the global sections of each fiber sub-sheaf.
    fiber_bases = all_fiber_bases(hom, F)

    # The stalk dimension at target vertex tv equals the nullspace dimension.
    pf_stalk_dims = [size(fiber_bases[tv], 2) for tv in 1:hom.n_target]
    PfF = EuclideanSheaf{Float64}(pf_stalk_dims)

    for (e, tri_u, tri_v) in cross_edges(hom, g)
        u, v = src(e), dst(e)
        pu, pv = phi[u], phi[v]
        # Align (u,v) with (tri_u, tri_v).
        src_u, src_v = pu == tri_u ? (u, v) : (v, u)

        fverts_u = fiber_vertices(hom, tri_u)
        fverts_v = fiber_vertices(hom, tri_v)

        # Row offset of src_u within the concatenated fiber-u stalk space.
        offset_u = sum(get_vertex_stalk(F, w)
                       for w in fverts_u if w < src_u; init=0)
        # Row offset of src_v within the concatenated fiber-v stalk space.
        offset_v = sum(get_vertex_stalk(F, w)
                       for w in fverts_v if w < src_v; init=0)

        d_u = get_vertex_stalk(F, src_u)
        d_v = get_vertex_stalk(F, src_v)

        # Original restriction maps (source vertex → edge).
        rm_u = get_restriction_map(F, src_u, src_v)
        rm_v = get_restriction_map(F, src_v, src_u)

        # Pushforward restriction maps = F's rm ∘ (rows of fiber basis for src_u/v).
        rm_pf_u = rm_u * fiber_bases[tri_u][offset_u+1 : offset_u+d_u, :]
        rm_pf_v = rm_v * fiber_bases[tri_v][offset_v+1 : offset_v+d_v, :]

        add_sheaf_edge!(PfF, tri_u, tri_v, rm_pf_u, rm_pf_v)
    end

    return PfF
end

# ---------------------------------------------------------------------------
# Transfer map
# ---------------------------------------------------------------------------

"""
    pushforward_transfer_map(hom::GraphHomomorphism, F) -> Matrix

Construct the linear transfer map ``T : C^0(F) \\to C^0(\\varphi_* F)`` that
sends every global section of `F` to the corresponding global section of the
pushforward sheaf ``\\varphi_* F``.

The map is built fiber-by-fiber: for each target vertex `tv`, `T` restricts a
0-cochain of `F` to the rows belonging to ``\\varphi^{-1}(\\text{tv})``, then
expresses those rows in the fiber-basis coordinates via a pseudoinverse.

Concretely, if ``B_{\\text{tv}}`` is the fiber-basis matrix (columns = global
sections of the restriction of `F` to fiber `tv`), then

```math
T[\\text{rows of }tv, \\text{cols of fiber}(tv)] = B_{\\text{tv}}^+
```

where ``B^+`` denotes the Moore–Penrose pseudoinverse.

The identity ``d_{\\varphi_* F} \\circ T \\circ s = 0`` holds for every global
section ``s`` of `F` (up to floating-point rounding).
"""
function pushforward_transfer_map(hom::GraphHomomorphism, F)
    # Re-use the same fiber bases that were used to build the pushforward sheaf.
    fiber_bases = all_fiber_bases(hom, F)

    pf_stalk_dims = [size(fiber_bases[tv], 2) for tv in 1:hom.n_target]

    # Cumulative offsets into F's 0-cochain space and the pushforward's 0-cochain space.
    F_cum  = [0; cumsum(vertex_stalks(F))]
    pf_cum = [0; cumsum(pf_stalk_dims)]

    T = zeros(sum(pf_stalk_dims), sum(vertex_stalks(F)))

    for tv in 1:hom.n_target
        fverts = fiber_vertices(hom, tv)
        # Gather the column indices of F's 0-cochain that belong to this fiber.
        f_rows = vcat([collect(F_cum[v]+1 : F_cum[v+1]) for v in fverts]...)
        B = fiber_bases[tv]
        # Project onto the fiber-basis coordinates via the pseudoinverse.
        T[pf_cum[tv]+1 : pf_cum[tv+1], f_rows] = pinv(B)
    end

    return T
end

end # module Pushforwards
