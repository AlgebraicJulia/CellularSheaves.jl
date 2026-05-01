"""
Module for the pushforward of a EuclideanSheaf along a graph homomorphism.

A *graph homomorphism* φ: G → H sends each vertex of G to a vertex of H and is
compatible with incidence: every edge {u1, u2} of G either

- maps to the edge {φ(u1), φ(u2)} in H  (when φ(u1) ≠ φ(u2)), or
- *collapses* to the vertex φ(u1) = φ(u2) in H.

The *pushforward sheaf* φ_*(F) on H is defined as follows.

**Stalks**

- At vertex w ∈ V(H): the space of global sections of F restricted to the
  *fiber subgraph* over w, which consists of all vertices u ∈ V(G) with
  φ(u) = w together with all edges of G whose both endpoints map to w.
  The stalk is represented concretely as ℝ^k where k = dim ker(B_w), and
  `global_sections_basis` returns an orthonormal basis matrix Π_w: ℝ^k → C⁰(F|_w).

- At edge e′ = {v1, v2} ∈ E(H): the direct sum
  ⊕_{e ∈ φ^{-1}(e′)} F(e)
  of the edge stalks of F over all preimage edges of e′.

**Restriction maps**

For a fiber edge (u_s, u_d) with φ(u_s) = v1 and φ(u_d) = v2, the
component of the restriction map (φ_*F)(v1) → (φ_*F)(e′) is

    F.restriction_maps[u_s => u_d] * Π_{v1}[rows for u_s, :]

and similarly for the v2-side.  The full restriction maps are assembled by
stacking these blocks vertically (in the order returned by `fiber_edge_edges`).
"""
module Pushforward

export GraphHomomorphism, fiber_vertices, fiber_vertex_edges, fiber_edge_edges,
    global_sections_basis, pushforward

using Graphs
using LinearAlgebra

using ..SheafInterface
import ..SheafInterface: add_vertex_stalk!, add_sheaf_edge!
using ..EuclideanSheaves: EuclideanSheaf, UnorderedPair
import ..SheafInterface: get_vertex_stalk, get_edge_stalk, get_restriction_map

# ---------------------------------------------------------------------------
# Type
# ---------------------------------------------------------------------------

"""
    GraphHomomorphism

A graph homomorphism φ: G → H.

The vertex map `vertex_map[u] = w` means vertex `u` of the source graph G maps
to vertex `w` of the target graph H.  The edge map is derived: edge {u1, u2}
in G maps to edge {vertex_map[u1], vertex_map[u2]} in H when those are
distinct (and the edge must exist in H), or collapses to the shared vertex
when vertex_map[u1] == vertex_map[u2].

The constructor validates that every non-collapsed edge has a corresponding
edge in the target.
"""
struct GraphHomomorphism
    source::Graph
    target::Graph
    vertex_map::Vector{Int}

    function GraphHomomorphism(source::Graph, target::Graph, vertex_map::Vector{Int})
        @assert length(vertex_map) == nv(source) "vertex_map must have one entry per source vertex (got $(length(vertex_map)), expected $(nv(source)))"
        @assert all(1 .<= vertex_map .<= nv(target)) "all vertex_map values must be valid target-graph vertex indices in 1:$(nv(target))"
        for e in edges(source)
            u1, u2 = src(e), dst(e)
            w1, w2 = vertex_map[u1], vertex_map[u2]
            if w1 != w2
                @assert has_edge(target, w1, w2) "source edge ($u1,$u2) maps to non-existent target edge ($w1,$w2)"
            end
        end
        new(source, target, vertex_map)
    end
end

# ---------------------------------------------------------------------------
# Fiber helpers
# ---------------------------------------------------------------------------

"""
    fiber_vertices(φ::GraphHomomorphism, w::Int) -> Vector{Int}

Return all vertices of the source graph that map to vertex `w` of the target.
"""
function fiber_vertices(φ::GraphHomomorphism, w::Int)
    return [u for u in vertices(φ.source) if φ.vertex_map[u] == w]
end

"""
    fiber_vertex_edges(φ::GraphHomomorphism, w::Int) -> Vector{Tuple{Int,Int}}

Return all edges (u1, u2) of the source graph whose *both* endpoints map to
vertex `w` of the target (i.e., edges that collapse to `w`).
"""
function fiber_vertex_edges(φ::GraphHomomorphism, w::Int)
    result = Tuple{Int,Int}[]
    for e in edges(φ.source)
        u1, u2 = src(e), dst(e)
        if φ.vertex_map[u1] == w && φ.vertex_map[u2] == w
            push!(result, (u1, u2))
        end
    end
    return result
end

"""
    fiber_edge_edges(φ::GraphHomomorphism, v1::Int, v2::Int) -> Vector{Tuple{Int,Int}}

Return all edges of the source graph that map to the undirected edge {v1, v2}
of the target.  Each returned tuple `(u_s, u_d)` is *oriented* so that
`φ(u_s) = v1` and `φ(u_d) = v2`.
"""
function fiber_edge_edges(φ::GraphHomomorphism, v1::Int, v2::Int)
    result = Tuple{Int,Int}[]
    for e in edges(φ.source)
        u1, u2 = src(e), dst(e)
        w1, w2 = φ.vertex_map[u1], φ.vertex_map[u2]
        if w1 == v1 && w2 == v2
            push!(result, (u1, u2))
        elseif w1 == v2 && w2 == v1
            push!(result, (u2, u1))   # re-orient so first maps to v1
        end
    end
    return result
end

# ---------------------------------------------------------------------------
# Global sections of a restricted sheaf
# ---------------------------------------------------------------------------

"""
    global_sections_basis(F, vertices, collapsed_edges) -> Matrix

Compute a matrix whose columns form an orthonormal basis for the space of
global sections of `F` restricted to the subgraph defined by `vertices`
(a `Vector{Int}` of source-graph vertex indices) and `collapsed_edges`
(a `Vector{Tuple{Int,Int}}` of edge endpoints within that vertex set).

The rows of the returned matrix are ordered according to the direct sum
⊕_{v ∈ vertices} F(v), with blocks in the order given by `vertices`.

When there are no collapsed edges the global sections equal the full
0-cochain space and the identity matrix is returned.  When the fiber is
empty an empty 0×0 matrix is returned.
"""
function global_sections_basis(F::EuclideanSheaf{T},
                                fib_verts::Vector{Int},
                                collapsed_edges::Vector{Tuple{Int,Int}}) where T
    if isempty(fib_verts)
        return Matrix{T}(undef, 0, 0)
    end

    stalk_dims  = [get_vertex_stalk(F, u) for u in fib_verts]
    total_v_dim = sum(stalk_dims)
    offsets     = [0; cumsum(stalk_dims)]

    if isempty(collapsed_edges)
        # No intra-fiber edges: coboundary is zero, sections = entire 0-cochain space.
        return Matrix{T}(I, total_v_dim, total_v_dim)
    end

    vertex_index = Dict(u => i for (i, u) in enumerate(fib_verts))
    edge_dims    = [get_edge_stalk(F, u1, u2) for (u1, u2) in collapsed_edges]
    total_e_dim  = sum(edge_dims)

    # Build the coboundary map B: C⁰ → C¹ for the restricted sheaf.
    # Convention (matching EuclideanSheaves.coboundary_map):
    #   B[edge e, :] = rm(u1→e) on column block u1  −  rm(u2→e) on column block u2
    B = zeros(T, total_e_dim, total_v_dim)
    e_off = 0
    for (k, (u1, u2)) in enumerate(collapsed_edges)
        d_e   = edge_dims[k]
        i1    = vertex_index[u1]
        i2    = vertex_index[u2]
        v1cols = (offsets[i1]+1):offsets[i1+1]
        v2cols = (offsets[i2]+1):offsets[i2+1]
        erows  = (e_off+1):(e_off+d_e)

        B[erows, v1cols] =  get_restriction_map(F, u1, u2)
        B[erows, v2cols] = -get_restriction_map(F, u2, u1)
        e_off += d_e
    end

    return nullspace(B)
end

# ---------------------------------------------------------------------------
# Pushforward
# ---------------------------------------------------------------------------

"""
    pushforward(φ::GraphHomomorphism, F::EuclideanSheaf{T}) -> EuclideanSheaf{T}

Compute the pushforward sheaf φ_*(F) on the target graph of φ.

**Vertex stalks** — For each vertex `w` of the target, the stalk is the space
of global sections of `F` restricted to the fiber subgraph over `w` (all
source vertices mapping to `w`, plus any source edges that collapse to `w`).
The dimension equals the nullity of the coboundary map of that sub-sheaf.

**Edge stalks** — For each edge `e′ = {v1, v2}` of the target, the stalk is
⊕_{e ∈ φ^{-1}(e′)} F(e).  Edges of the target with an empty fiber are
omitted from the result.

**Restriction maps** — The map from (φ_*F)(v1) to (φ_*F)(e′) is assembled
block-by-block: for each fiber edge `(u_s, u_d)` with `φ(u_s) = v1`,
the block is `F.restriction_maps[u_s → u_d]` composed with the projection of
global sections onto the `F(u_s)` component of the 0-cochain space.
"""
function pushforward(φ::GraphHomomorphism, F::EuclideanSheaf{T}) where T
    H      = φ.target
    result = EuclideanSheaf{T}(Int[])

    # -----------------------------------------------------------------------
    # Step 1: vertex stalks
    # For each target vertex w compute the global-sections basis of the fiber
    # subgraph.  Store basis, fiber vertex list, and row-offset vector.
    # -----------------------------------------------------------------------
    n_target    = nv(H)
    fib_bases   = Vector{Matrix{T}}(undef, n_target)
    fib_verts   = Vector{Vector{Int}}(undef, n_target)
    fib_offsets = Vector{Vector{Int}}(undef, n_target)

    for w in 1:n_target
        fv      = fiber_vertices(φ, w)
        fe      = fiber_vertex_edges(φ, w)
        basis   = global_sections_basis(F, fv, fe)
        sdims   = [get_vertex_stalk(F, u) for u in fv]
        offs    = [0; cumsum(sdims)]

        fib_bases[w]   = basis
        fib_verts[w]   = fv
        fib_offsets[w] = offs
        add_vertex_stalk!(result, size(basis, 2))
    end

    # -----------------------------------------------------------------------
    # Step 2: edge stalks and restriction maps
    # For each target edge e′ = {v1, v2}, collect the preimage edges and
    # assemble the two restriction maps.
    # -----------------------------------------------------------------------
    for e′ in edges(H)
        v1, v2  = src(e′), dst(e′)
        fib_es  = fiber_edge_edges(φ, v1, v2)
        isempty(fib_es) && continue   # no preimage edges → omit

        basis1  = fib_bases[v1];   verts1  = fib_verts[v1];   offs1 = fib_offsets[v1]
        basis2  = fib_bases[v2];   verts2  = fib_verts[v2];   offs2 = fib_offsets[v2]

        vi1 = Dict(u => i for (i, u) in enumerate(verts1))
        vi2 = Dict(u => i for (i, u) in enumerate(verts2))

        rm1_blocks = Matrix{T}[]
        rm2_blocks = Matrix{T}[]

        for (u_s, u_d) in fib_es
            # u_s ∈ fiber(v1), u_d ∈ fiber(v2)
            is_  = vi1[u_s]
            id_  = vi2[u_d]
            rows_s = (offs1[is_]+1):offs1[is_+1]
            rows_d = (offs2[id_]+1):offs2[id_+1]

            # Restriction from (φ_*F)(v1):  F(u_s → edge) ∘ Π_1|_{u_s rows}
            push!(rm1_blocks, get_restriction_map(F, u_s, u_d) * basis1[rows_s, :])

            # Restriction from (φ_*F)(v2):  F(u_d → edge) ∘ Π_2|_{u_d rows}
            push!(rm2_blocks, get_restriction_map(F, u_d, u_s) * basis2[rows_d, :])
        end

        add_sheaf_edge!(result, v1, v2, vcat(rm1_blocks...), vcat(rm2_blocks...))
    end

    return result
end

end # module Pushforward
