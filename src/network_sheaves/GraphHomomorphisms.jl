"""
Module for graph homomorphisms and their induced operations on cellular sheaves.

A *graph homomorphism* φ : G → H is a map on vertices that preserves edges.
Given such a map, the *pushforward* φ_*F of a cellular sheaf F on G is a
cellular sheaf on H whose global sections are isomorphic to those of F.

This module provides:

- [`GraphHomomorphism`](@ref) — a graph homomorphism specified by a vertex map.
- [`fiber_vertices`](@ref) — vertices of the source graph that map to a given target vertex.
- [`fiber_edges`](@ref) — edges within a single fiber (both endpoints in the same fiber).
- [`cross_edges`](@ref) — edges whose endpoints lie in different fibers.
- [`pushforward_sheaf`](@ref) — construct the pushforward sheaf φ_*F.
- [`pushforward_transfer_map`](@ref) — construct the linear map T : C⁰(F) → C⁰(φ_*F)
  that sends every global section of F to the corresponding global section of φ_*F.
"""
module GraphHomomorphisms

export GraphHomomorphism, fiber_vertices, fiber_edges, cross_edges,
       pushforward_sheaf, pushforward_transfer_map

using Graphs
using LinearAlgebra
using SparseArrays

using ..SheafInterface: vertex_stalks, coboundary_map, underlying_graph,
                        get_vertex_stalk, get_restriction_map
using ..EuclideanSheaves: EuclideanSheaf, add_vertex_stalk!, add_sheaf_edge!, nullspace_ldlt

# ---------------------------------------------------------------------------
# Type
# ---------------------------------------------------------------------------

"""
    GraphHomomorphism

A graph homomorphism φ : G → H specified by a vertex map.

Fields
------
- `vertex_map :: Vector{Int}` — `vertex_map[i]` is the image of source vertex `i`
  in the target graph.
- `n_target :: Int` — number of vertices in the target graph (inferred from the
  maximum of `vertex_map` if not provided explicitly).
"""
struct GraphHomomorphism
    vertex_map::Vector{Int}
    n_target::Int
end

"""
    GraphHomomorphism(vertex_map)

Construct a [`GraphHomomorphism`](@ref) from a vertex map vector, inferring the
number of target vertices from the maximum entry.
"""
GraphHomomorphism(vertex_map::Vector{Int}) =
    GraphHomomorphism(vertex_map, maximum(vertex_map))

# ---------------------------------------------------------------------------
# Printing
# ---------------------------------------------------------------------------

function Base.show(io::IO, hom::GraphHomomorphism)
    n_src = length(hom.vertex_map)
    n_tgt = hom.n_target
    println(io, "GraphHomomorphism: $n_src source vertices → $n_tgt target vertices")
    println(io, "  vertex_map: ", hom.vertex_map)
    # Pre-compute fiber sizes for display (no graph needed — just the vertex map).
    fiber_sizes = [count(==(tv), hom.vertex_map) for tv in 1:n_tgt]
    println(io, "  fiber sizes: ", fiber_sizes)
end

# ---------------------------------------------------------------------------
# Statistics
# ---------------------------------------------------------------------------

"""
    fiber_vertices(hom::GraphHomomorphism, tv::Int) -> Vector{Int}

Return the source vertices whose image under `hom` is `tv`.
"""
function fiber_vertices(hom::GraphHomomorphism, tv::Int)
    return findall(==(tv), hom.vertex_map)
end

"""
    fiber_edges(hom::GraphHomomorphism, g::AbstractGraph, tv::Int) -> Vector{Tuple{Int,Int}}

Return all edges of `g` that are *fiber edges* of `tv`: both endpoints map to `tv`.
"""
function fiber_edges(hom::GraphHomomorphism, g::AbstractGraph, tv::Int)
    result = Tuple{Int,Int}[]
    for e in edges(g)
        u, v = src(e), dst(e)
        if hom.vertex_map[u] == tv && hom.vertex_map[v] == tv
            push!(result, (u, v))
        end
    end
    return result
end

"""
    cross_edges(hom::GraphHomomorphism, g::AbstractGraph) -> Vector{Tuple{Any,Int,Int}}

Return all *cross edges* of `g` — edges whose endpoints lie in different fibers.

Each entry is `(edge, src_target_vertex, dst_target_vertex)` where `edge` is a
`Graphs.AbstractEdge` and `src_target_vertex`, `dst_target_vertex` are the images
of the edge's source and destination vertices.
"""
function cross_edges(hom::GraphHomomorphism, g::AbstractGraph)
    result = Tuple{Any,Int,Int}[]
    for e in edges(g)
        u, v = src(e), dst(e)
        pu, pv = hom.vertex_map[u], hom.vertex_map[v]
        if pu != pv
            push!(result, (e, pu, pv))
        end
    end
    return result
end

# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

# Build a basis for the global sections of the sub-sheaf of F restricted to
# the fiber of target vertex `tv`.  Returns an (total_fiber_dim × k) matrix
# whose columns span the fiber's global-section space.
function _fiber_section_basis(F, verts::Vector{Int}, fedges::Vector{Tuple{Int,Int}})
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

# Compute the fiber-section basis for every target vertex.
function _all_fiber_bases(hom::GraphHomomorphism, F)
    g = underlying_graph(F)
    return [
        _fiber_section_basis(F,
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
[`nullspace_ldlt`](@ref).

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
    fiber_bases = _all_fiber_bases(hom, F)

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
    fiber_bases = _all_fiber_bases(hom, F)

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

end # module GraphHomomorphisms
