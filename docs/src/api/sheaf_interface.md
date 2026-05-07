# Sheaf Interface

The **Sheaf Interface** defines the abstract API for all network sheaves in the package.  
All concrete sheaf types (e.g. `EuclideanSheaf`) implement this interface, providing
the combinatorial data (graph, stalks, restriction maps) and the linear‑algebra
operations needed to compute coboundaries, Laplacians, and perform linear‑algebraic
tasks such as finding global sections.

Below is a quick reference of the most frequently used methods.  For a deeper
walk‑through, see the tutorial pages linked from the navigation bar.

## Core abstract type

```julia
abstract type AbstractNetworkSheaf end
```

Represents a sheaf defined on a graph.  It stores:

* a `SimpleGraph` underlying structure,
* a vector of vertex‑stalk dimensions,
* a dictionary mapping unordered vertex pairs to edge‑stalk dimensions,
* restriction maps (`Matrix` objects) for each incident vertex–edge pair.

Concrete subtypes must implement the methods listed below.

## Public API

| Function | Brief description |
|----------|-------------------|
| `vertex_stalks(s::AbstractNetworkSheaf)` | Return a vector of dimensions of vertex stalks. |
| `edge_stalks(s::AbstractNetworkSheaf)` | Return a dictionary mapping each edge to its stalk dimension. |
| `edge_stalk_dimensions(s::AbstractNetworkSheaf)` | Return a vector of edge‑stalk dimensions ordered by graph edges. |
| `underlying_graph(s::AbstractNetworkSheaf)` | Return the `SimpleGraph` on which the sheaf is defined. |
| `get_vertex_stalk(s::AbstractNetworkSheaf, v::Int)` | Dimension of the stalk at vertex `v`. |
| `get_edge_stalk(s::AbstractNetworkSheaf, v1::Int, v2::Int)` | Dimension of the stalk on edge `(v1,v2)`. |
| `get_restriction_map(s::AbstractNetworkSheaf, v1::Int, v2::Int)` | Return the restriction matrix from vertex `v1` to the edge `(v1,v2)`. |
| `add_vertex_stalk!(s::AbstractNetworkSheaf, stalk_dim::Int)` | Append a new vertex with the given stalk dimension, updating the graph. |
| `add_sheaf_edge!(s::AbstractNetworkSheaf, v1::Int, v2::Int, rm1, rm2)` | Add an edge between `v1` and `v2` together with its two restriction maps. |
| `coboundary_map(s::AbstractNetworkSheaf)` | Construct the coboundary operator as a block‑sparse matrix. |
| `sheaf_laplacian(s::AbstractNetworkSheaf)` | Return the Laplacian operator (as a function). |

All of these functions are lightweight wrappers that delegate to the concrete
implementation (e.g. `EuclideanSheaf`).  They are exported from the package and
documented automatically via `@autodocs`, so the reference page below contains the
full signature list.

## Minimal example

```julia
using CellularSheaves
using Graphs

# a triangle graph
g = Graph(3)
add_edge!(g, 1, 2); add_edge!(g, 2, 3); add_edge!(g, 3, 1)

# 1‑dimensional stalks (R) at each vertex
stalks = [1, 1, 1]

# create an empty Euclidean sheaf
F = EuclideanSheaf(g, stalks)

# add an edge with identity restriction maps
add_sheaf_edge!(F, 1, 2, [1], [1])
```

See the [Core Sheaf Workflows](@ref) guide.

## API Reference

```@autodocs
Modules = [CellularSheaves.SheafInterface]
```