module SheafInterface

export AbstractNetworkSheaf, vertex_stalks, edge_stalks, edge_stalk_dimensions, coboundary_map, add_vertex_stalk!, add_sheaf_edge!, underlying_graph,
    get_vertex_stalk, get_edge_stalk, get_restriction_map, sheaf_laplacian

import Base: show
using DocStringExtensions

"""     AbstractNetworkSheaf

An abstract type for network sheaves.
A network sheaf is a cellular sheaf on a graph, i.e. a sheaf on a 1-dimensional cell complex.
    It consists of:
    - A graph G = (V, E)
    - A stalk S_v for each vertex v in V
    - A stalk S_e for each edge e in E
    - A restriction map r_{e->v} : S_e -> S_v for each incident vertex-edge pair (v, e)
"""
abstract type AbstractNetworkSheaf end

"""
    $(TYPEDSIGNATURES)

Get the vertex stalks of a network sheaf.

# Arguments
- `s`: The network sheaf

# Returns
- Vector{Int}: Dimensions of each vertex stalk, where the i-th element corresponds to the dimension of the stalk at vertex i

# See also
- [`edge_stalks`](@ref)
- [`AbstractNetworkSheaf`](@ref)
"""
function vertex_stalks(s::AbstractNetworkSheaf)
    error("vertex_stalks not implemented")
end

"""
    $(TYPEDSIGNATURES)

Get the edge stalks of a network sheaf.

# Arguments
- `s`: The network sheaf

# Returns
- Dict{UnorderedPair{Int}, Int}: Mapping from edge (as unordered vertex pair) to stalk dimension

# See also
- [`vertex_stalks`](@ref)
- [`AbstractNetworkSheaf`](@ref)
"""
function edge_stalks(s::AbstractNetworkSheaf)
    error("edge_stalks not implemented")
end

function edge_stalk_dimensions(s::AbstractNetworkSheaf)
    error("edge_stalk_dimensions not implemented")
end

"""
    $(TYPEDSIGNATURES)

Get the underlying graph of a network sheaf.

# Arguments
- `s`: The network sheaf

# Returns
- Graphs.SimpleGraph: The graph on which the sheaf is defined

# See also
- [`AbstractNetworkSheaf`](@ref)
"""
function underlying_graph(s::AbstractNetworkSheaf)
    error("underlying_graph not implemented")
end

"""
    $(TYPEDSIGNATURES)

Get the dimension of the stalk at a specific vertex.

# Arguments
- `s`: The network sheaf
- `v`: The vertex index

# Returns
- Int: Dimension of the vertex stalk at vertex `v`

# See also
- [`vertex_stalks`](@ref)
- [`get_edge_stalk`](@ref)
"""
function get_vertex_stalk(s::AbstractNetworkSheaf, v)
    error("get_vertex_stalk not implemented")
end

"""
    $(TYPEDSIGNATURES)

Get the dimension of the stalk at a specific edge.

# Arguments
- `s`: The network sheaf
- `v1`: The first vertex of the edge
- `v2`: The second vertex of the edge

# Returns
- Int: Dimension of the edge stalk on edge (v1, v2)

# See also
- [`edge_stalks`](@ref)
- [`get_vertex_stalk`](@ref)
"""
function get_edge_stalk(s::AbstractNetworkSheaf, v1, v2)
    error("get_edge_stalk not implemented")
end

"""
    $(TYPEDSIGNATURES)

Get the restriction map from vertex v1 to the edge (v1, v2).

# Arguments
- `s`: The network sheaf
- `v1`: The source vertex
- `v2`: The target vertex (defines the edge with v1)

# Returns
- Matrix: The restriction map matrix from vertex v1 stalk to edge (v1, v2) stalk

# See also
- [`get_edge_stalk`](@ref)
- [`get_vertex_stalk`](@ref)
"""
function get_restriction_map(s::AbstractNetworkSheaf, v1, v2)
    error("get_restriction_map not implemented")
end


"""
    $(TYPEDSIGNATURES)

Add a new vertex stalk to the network sheaf.

# Arguments
- `s`: The network sheaf to modify
- `stalk`: The dimension of the new vertex stalk to add

# Returns
- Nothing: The sheaf is modified in-place

# See also
- [`add_sheaf_edge!`](@ref)
- [`vertex_stalks`](@ref)
"""
function add_vertex_stalk!(s::AbstractNetworkSheaf, stalk)
    error("add_vertex_stalk! not implemented")
end

"""
    $(TYPEDSIGNATURES)

Add a new edge to the network sheaf with specified restriction maps.

# Arguments
- `s`: The network sheaf to modify
- `v1`: The first vertex of the edge
- `v2`: The second vertex of the edge
- `rm1`: The restriction map from vertex v1 to the edge
- `rm2`: The restriction map from vertex v2 to the edge

# Returns
- Int: The number of edges in the sheaf after addition

# See also
- [`add_vertex_stalk!`](@ref)
- [`get_restriction_map`](@ref)
"""
function add_sheaf_edge!(s::AbstractNetworkSheaf, v1, v2, rm1, rm2)
    error("add_sheaf_edge! not implemented")
end

"""
    $(TYPEDSIGNATURES)

Compute the coboundary map of a network sheaf.

# Arguments
- `s`: The network sheaf

# Returns
- AbstractMatrix: The coboundary map as a matrix-like linear operator. Implementations may return a block-sparse matrix representation supporting operations such as `*` and `adjoint`.

# See also
- [`sheaf_laplacian`](@ref)
- [`AbstractNetworkSheaf`](@ref)
"""
function coboundary_map(s::AbstractNetworkSheaf)
    error("coboundary_map not implemented")
end

"""
    $(TYPEDSIGNATURES)

Compute the sheaf Laplacian of a network sheaf.

# Arguments
- `s`: The network sheaf

# Returns
- Function: The sheaf Laplacian operator (as a function that takes a vector and returns a vector)

# See also
- [`coboundary_map`](@ref)
- `sheaf_laplacian_matrix` — computes the matrix representation of the sheaf Laplacian
"""
function sheaf_laplacian(s::AbstractNetworkSheaf)
    error("laplacian not implemented")
end

"""
    $(TYPEDSIGNATURES)

Construct a basis-indexed family of feasible controlled trajectories from an
affine feasible-space parameterization.

# Arguments
- `ts`: controlled trajectory object
- `z_p`: particular feasible trajectory in public coordinates
- `null_basis`: matrix whose columns span endpoint-preserving feasible perturbations

# Returns
- Vector of trajectory objects, one (or two, if signed) per basis column.
"""
function nullspace_trajectory_family(ts, z_p, null_basis; amplitude=1, include_negative=false)
    error("nullspace_trajectory_family not implemented")
end

function Base.show(io::IO, s::AbstractNetworkSheaf)
    println(io, "A network sheaf with ", length(vertex_stalks(s)), " vertex stalks and ", length(edge_stalks(s)), " edge stalks.")
end

end