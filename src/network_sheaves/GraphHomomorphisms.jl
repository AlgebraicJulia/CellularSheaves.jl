"""
Module for representing and interacting with graph homomorphisms.

A *graph homomorphism* φ : G → H is a map on vertices (and edges) that
respects the graph structure.  This module provides:

- [`GraphHomomorphism`](@ref) — a graph homomorphism specified by a vertex map.
- [`fiber_vertices`](@ref) — vertices of the source graph that map to a given target vertex.
- [`fiber_edges`](@ref) — edges within a single fiber (both endpoints in the same fiber).
- [`cross_edges`](@ref) — edges whose endpoints lie in different fibers.
- [`compose`](@ref) — composition of two graph homomorphisms.
- [`graph_pushout`](@ref) — categorical pushout of a span of graph homomorphisms.

For the action of a graph homomorphism on a cellular sheaf (pushforward), see
the [`Pushforwards`](@ref) module.
"""
module GraphHomomorphisms

export GraphHomomorphism, fiber_vertices, fiber_edges, cross_edges, compose, graph_pushout

using Graphs

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

    function GraphHomomorphism(vertex_map::Vector{Int}, n_target::Int)
        n_target < 0 && throw(ArgumentError("n_target must be nonnegative, got $n_target"))
        any(<(1), vertex_map) && throw(ArgumentError("vertex_map entries must be positive target vertex indices"))
        any(>(n_target), vertex_map) && throw(ArgumentError("vertex_map entries must lie in 1:n_target"))
        new(vertex_map, n_target)
    end
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
    cross_edges(hom::GraphHomomorphism, g::AbstractGraph) -> Vector{Tuple{E,Int,Int}} where E

Return all *cross edges* of `g` — edges whose endpoints lie in different fibers.

Each entry is `(edge, src_target_vertex, dst_target_vertex)` where `edge` has the
concrete edge type produced by `edges(g)`, and `src_target_vertex`,
`dst_target_vertex` are the images of the edge's source and destination vertices.
"""
function cross_edges(hom::GraphHomomorphism, g::AbstractGraph)
    result = Tuple{eltype(edges(g)),Int,Int}[]
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
# Composition
# ---------------------------------------------------------------------------

"""
    compose(f::GraphHomomorphism, g::GraphHomomorphism) -> GraphHomomorphism

Compute the composition `g ∘ f` (apply `f` first, then `g`).
"""
function compose(f::GraphHomomorphism, g::GraphHomomorphism)
    return GraphHomomorphism(g.vertex_map[f.vertex_map], g.n_target)
end

# ---------------------------------------------------------------------------
# Pushout
# ---------------------------------------------------------------------------

"""
    graph_pushout(φ::GraphHomomorphism, ψ::GraphHomomorphism,
                  G::SimpleGraph, H::SimpleGraph, nK::Int)
        -> (SimpleGraph, GraphHomomorphism, GraphHomomorphism)

Compute the pushout `Q = G ⊔_K H` of the span `G <--φ-- K --ψ--> H`.

`nK` is the number of vertices in `K` (i.e. `length(φ.vertex_map)`).

Returns `(Q, jG, jH)` where `jG : G → Q` and `jH : H → Q` are the canonical
graph homomorphisms satisfying `jG ∘ φ == jH ∘ ψ` (as vertex maps).
"""
function graph_pushout(φ::GraphHomomorphism, ψ::GraphHomomorphism,
                       G::SimpleGraph, H::SimpleGraph, nK::Int)
    nG, nH = nv(G), nv(H)
    parent = collect(1:nG+nH)

    function uf_find(x)
        while parent[x] != x
            parent[x] = parent[parent[x]]   # path compression
            x = parent[x]
        end
        return x
    end

    function uf_union!(x, y)
        rx, ry = uf_find(x), uf_find(y)
        rx != ry && (parent[rx] = ry)
    end

    # Identify φ(k) in G with ψ(k) in H for each vertex k of K
    for k in 1:nK
        uf_union!(φ.vertex_map[k], nG + ψ.vertex_map[k])
    end

    # Assign contiguous IDs 1..nQ to equivalence classes
    roots = [uf_find(i) for i in 1:nG+nH]
    unique_roots = unique(roots)
    id = Dict(r => i for (i, r) in enumerate(unique_roots))
    component = [id[uf_find(i)] for i in 1:nG+nH]

    nQ = length(unique_roots)
    Q = SimpleGraph(nQ)
    for e in edges(G)
        add_edge!(Q, component[src(e)], component[dst(e)])   # no-op for self-loops
    end
    for e in edges(H)
        add_edge!(Q, component[nG + src(e)], component[nG + dst(e)])
    end

    jG = GraphHomomorphism(component[1:nG], nQ)
    jH = GraphHomomorphism(component[nG+1:end], nQ)
    return Q, jG, jH
end

end # module GraphHomomorphisms

