"""
Module for representing and interacting with graph homomorphisms.

A *graph homomorphism* φ : G → H is a map on vertices (and edges) that
respects the graph structure.  This module provides:

- [`GraphHomomorphism`](@ref) — a graph homomorphism specified by a vertex map.
- [`fiber_vertices`](@ref) — vertices of the source graph that map to a given target vertex.
- [`fiber_edges`](@ref) — edges within a single fiber (both endpoints in the same fiber).
- [`cross_edges`](@ref) — edges whose endpoints lie in different fibers.

For the action of a graph homomorphism on a cellular sheaf (pushforward), see
the [`Pushforwards`](@ref) module.
"""
module GraphHomomorphisms

export GraphHomomorphism, fiber_vertices, fiber_edges, cross_edges

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

end # module GraphHomomorphisms

