# Module for network sheaves valued in Euclidean spaces with linear restriction maps
module EuclideanSheaves

export EuclideanSheaf, UnorderedPair, sheaf_laplacian_matrix, sheaf_from_graph, energy_function,
    nearest_global_section, edge_stalk_dimensions, nullspace_ldlt,
    zero_sheaf, constant_sheaf, cycle_sheaf

using Graphs
using AutoHashEquals: @auto_hash_equals
using LinearOperators
using Krylov
using LinearAlgebra
using BlockArrays
using CliqueTrees.Multifrontal
using SparseArrays
import Base: hash, ==, isequal

using ..SheafInterface
import ..SheafInterface: vertex_stalks, edge_stalks, edge_stalk_dimensions, coboundary_map, add_vertex_stalk!, add_sheaf_edge!, underlying_graph,
    get_vertex_stalk, get_edge_stalk, get_restriction_map, sheaf_laplacian
using ..BlockSparseArrays

struct UnorderedPair{T}
    first::T
    second::T
end

function Base.hash(up::UnorderedPair{T}, h::UInt) where T
    return hash((up.first, up.second), h) + hash((up.second, up.first), h)
end

function Base.:(==)(up1::UnorderedPair{T}, up2::UnorderedPair{T}) where T
    return (up1.first == up2.first && up1.second == up2.second) ||
           (up1.first == up2.second && up1.second == up2.first)
end

function Base.:(isequal)(up1::UnorderedPair{T}, up2::UnorderedPair{T}) where T
    return (up1.first == up2.first && up1.second == up2.second) ||
           (up1.first == up2.second && up1.second == up2.first)
end

"""     EuclideanSheaf{T}

A Euclidean sheaf is a network sheaf where each vertex stalk is a Euclidean space R^n for some n,
each edge stalk is a Euclidean space R^m for some m, and each restriction map is a linear map
R^n -> R^m represented by a matrix of type T.
"""
@auto_hash_equals struct EuclideanSheaf{T} <: AbstractNetworkSheaf
    vertex_stalks::Vector{Int}
    edge_stalks::Dict{UnorderedPair{Int},Int}
    underlying_graph::Graph
    restriction_maps::Dict{Pair{Int},Matrix{T}}
end

EuclideanSheaf{T}(vertex_stalks::Vector{Int}) where T = EuclideanSheaf{T}(vertex_stalks, Dict{UnorderedPair{Int},Int}(), Graph(length(vertex_stalks)), Dict{Pair{Int},Matrix{T}}())


function sheaf_from_graph(g::Graph, stalk_dim::Int, rm_generator::Function; symmetric_edges=false)
    n = nv(g)
    s = EuclideanSheaf{Float64}(repeat([stalk_dim], n))

    for e in edges(g)
        i, j = src(e), dst(e)
        if symmetric_edges
            rm = rm_generator(stalk_dim)
            add_sheaf_edge!(s, i, j, rm, rm)
        else
            rm1 = rm_generator(stalk_dim)
            rm2 = rm_generator(stalk_dim)
            add_sheaf_edge!(s, i, j, rm1, rm2)
        end
    end
    return s
end

function sheaf_from_graph(g::Graph, stalk_dim::Int, rm1_generator::Function, rm2_generator::Function)
    n = nv(g)
    s = EuclideanSheaf{Float64}(repeat([stalk_dim], n))

    for e in edges(g)
        i, j = src(e), dst(e)
        rm1 = rm1_generator(stalk_dim)
        rm2 = rm2_generator(stalk_dim)
        add_sheaf_edge!(s, i, j, rm1, rm2)
    end
    return s
end

function vertex_stalks(s::EuclideanSheaf)
    return s.vertex_stalks
end

function edge_stalks(s::EuclideanSheaf)
    return s.edge_stalks
end

function underlying_graph(s::EuclideanSheaf)
    return s.underlying_graph
end

function get_vertex_stalk(s::EuclideanSheaf, v::Int)
    @assert v <= length(s.vertex_stalks)
    return s.vertex_stalks[v]
end

function get_edge_stalk(s::EuclideanSheaf, v1::Int, v2::Int)
    edge_key = UnorderedPair(v1, v2)
    @assert haskey(s.edge_stalks, edge_key)
    return s.edge_stalks[edge_key]
end

"""    edge_stalk_dimensions(s::EuclideanSheaf)

Return a vector of edge stalk dimensions ordered by edges in the underlying graph.
"""
function edge_stalk_dimensions(s::EuclideanSheaf)
    return [get_edge_stalk(s, src(e), dst(e)) for e in edges(underlying_graph(s))]
end

"""    get_restriction_map(s::EuclideanSheaf, v1::Int, v2::Int)

Get the restriction map from vertex v1 to the edge (v1, v2).
"""
function get_restriction_map(s::EuclideanSheaf, v1::Int, v2::Int)
    @assert haskey(s.restriction_maps, v1 => v2)
    return s.restriction_maps[v1=>v2]
end

function add_vertex_stalk!(s::EuclideanSheaf, stalk_size::Int)
    push!(s.vertex_stalks, stalk_size)
    add_vertex!(s.underlying_graph)
end

"""   add_sheaf_edge!(s::EuclideanSheaf, v1::Int, v2::Int, rm1, rm2)

Add an edge between vertices `v1` and `v2` with restriction maps `rm1` and
`rm2`.  `rm1` is the restriction map from vertex `v1` to the edge and `rm2`
is the restriction map from vertex `v2` to the edge.

Both `rm1` and `rm2` may be any `AbstractMatrix` (e.g. `Matrix`, `Diagonal`,
`SparseMatrixCSC`); they are converted to `Matrix{T}` when stored.
"""
function add_sheaf_edge!(s::EuclideanSheaf{T}, v1::Int, v2::Int, rm1::AbstractMatrix, rm2::AbstractMatrix) where T
    @assert v1 <= length(s.vertex_stalks) && v2 <= length(s.vertex_stalks)
    @assert size(rm1, 1) == size(rm2, 1)
    @assert size(rm1, 2) == s.vertex_stalks[v1]
    @assert size(rm2, 2) == s.vertex_stalks[v2]

    stalk_size = size(rm1, 1)
    add_edge!(s.underlying_graph, v1, v2)
    edge_key = UnorderedPair(v1, v2)
    s.edge_stalks[edge_key] = stalk_size
    s.restriction_maps[v1=>v2] = Matrix{T}(rm1)
    s.restriction_maps[v2=>v1] = Matrix{T}(rm2)
    return ne(s.underlying_graph)
end

"""    zero_sheaf(g::Graph) -> EuclideanSheaf{Float64}

Construct the *zero sheaf* on graph `g`: all vertex and edge stalks are
0-dimensional (i.e. the trivial bundle).
"""
function zero_sheaf(g::Graph)
    s = EuclideanSheaf{Float64}(zeros(Int, nv(g)))
    for e in edges(g)
        rm = zeros(Float64, 0, 0)
        add_sheaf_edge!(s, src(e), dst(e), rm, rm)
    end
    return s
end

"""    constant_sheaf(g::Graph, d::Int=1) -> EuclideanSheaf{Float64}

Construct the *constant sheaf* ``\\underline{\\mathbb{R}^d}`` on graph `g`:
each vertex and edge stalk is ``\\mathbb{R}^d`` and every restriction map is
the ``d\\times d`` identity.
"""
function constant_sheaf(g::Graph, d::Int=1)
    s = EuclideanSheaf{Float64}(fill(d, nv(g)))
    rm = Matrix{Float64}(I, d, d)
    for e in edges(g)
        add_sheaf_edge!(s, src(e), dst(e), rm, rm)
    end
    return s
end

"""    cycle_sheaf(n::Int, stalk_dim::Int, rm_pair_fn::Function) -> EuclideanSheaf{Float64}

Construct an `EuclideanSheaf` on the `n`-vertex cycle graph where each vertex
stalk has dimension `stalk_dim`.  For edge `k` (in graph iteration order) the
pair of restriction maps `(rm_src, rm_dst)` is produced by calling
`rm_pair_fn(n, stalk_dim, k)`.
"""
function cycle_sheaf(n::Int, stalk_dim::Int, rm_pair_fn::Function)
    g = cycle_graph(n)
    s = EuclideanSheaf{Float64}(fill(stalk_dim, n))
    for (k, e) in enumerate(edges(g))
        rm1, rm2 = rm_pair_fn(n, stalk_dim, k)
        add_sheaf_edge!(s, src(e), dst(e), rm1, rm2)
    end
    return s
end

function coboundary_map(s::EuclideanSheaf{T}) where T
    I = Int64[]
    J = Int64[]
    V = Matrix{T}[]

    for (e_idx, e) in enumerate(edges(s.underlying_graph))
        i = src(e)
        j = dst(e)
        rm1 = s.restriction_maps[i=>j]
        rm2 = s.restriction_maps[j=>i]

        push!(J, i, j)
        push!(I, e_idx, e_idx)
        push!(V, rm1, -rm2)
    end
    return blocksparse(I, J, V)
end

function sheaf_laplacian(s::EuclideanSheaf)
    B = coboundary_map(s)
    return x -> B' * (B * x)
end

function sheaf_laplacian_matrix(s::EuclideanSheaf)
    B = coboundary_map(s)
    return B' * B
end


function energy_function(L::AbstractMatrix)
    return x -> 0.5 * x' * (L * x)
end

function energy_function(s::EuclideanSheaf)
    return energy_function(sheaf_laplacian_matrix(s))
end

"""nearest_global_section_cg

Iterative Krylov (CG) backend. Returns a BlockArray global section.
"""
function nearest_global_section_cg(s::EuclideanSheaf, x; verbose=false)
    d = coboundary_map(s)
    eL = LinearOperator(d) * LinearOperator(d')
    b = d * x

    y, stats = cg(eL, Array(b))
    if verbose
        println(stats)
    end

    return BlockArray(x - d' * y, s.vertex_stalks)
end

"""nearest_global_section_ldl

Deterministic direct backend using LDL factorization on the symmetric
edge-space operator E = d*d'. This handles symmetric indefinite matrices
without Tikhonov regularization. Falls back to pinv if LDL fails.
"""
function nearest_global_section_ldl(s::EuclideanSheaf, x; verbose=false)
    d = coboundary_map(s)
    E = d * d'
    b = d * x

    # Dense symmetric conversion for docs / small examples. Replace with a
    # sparse LDL solver for larger problems if needed.
    A = Symmetric(Array(E))

    try
        F = ldlt(A)
        y = F \ Array(b)
        if verbose
            println("nearest_global_section: solved with LDL factorization")
        end
        return BlockArray(x - d' * y, s.vertex_stalks)
    catch err
        if verbose
            println("nearest_global_section: LDL failed, falling back to pinv; error=", err)
        end
        y = pinv(Array(A)) * Array(b)
        return BlockArray(x - d' * y, s.vertex_stalks)
    end

end

"""nearest_global_section_pinv

Robust pseudoinverse backend (SVD). Deterministic but potentially expensive.
"""
function nearest_global_section_pinv(s::EuclideanSheaf, x; verbose=false)
    d = coboundary_map(s)
    E = d * d'
    b = d * x

    y = pinv(Array(E)) * Array(b)
    return BlockArray(x - d' * y, s.vertex_stalks)
end

"""nearest_global_section

Entrypoint solver. Keyword `method` selects backend: :cg (iterative, default),
:ldl (direct LDL), or :pinv (SVD pseudoinverse)."""
function nearest_global_section(s::EuclideanSheaf, x; method::Symbol=:cg, verbose=false)
    if method == :cg
        return nearest_global_section_cg(s, x; verbose=verbose)
    elseif method == :ldl
        return nearest_global_section_ldl(s, x; verbose=verbose)
    elseif method == :pinv
        return nearest_global_section_pinv(s, x; verbose=verbose)
    else
        error("Unknown method: $(method). Valid options are :cg, :ldl, :pinv")
    end
end


"""    nullspace_ldlt(X::AbstractMatrix; tol=nothing) -> Matrix

Compute a basis for the nullspace of the symmetric positive-semidefinite matrix `X`
using a sparse LDLt factorisation from `CliqueTrees.Multifrontal` (`ChordalLDLt`
with `RowMaximum` pivoting).

The factorisation is ``X = P^\\mathsf{T} L D L^\\mathsf{T} P``.  Columns of `D`
whose absolute diagonal value is at or below `tol` (default:
``\\varepsilon_\\text{Float64} \\times \\max(1, \\|D\\|_\\infty)``) identify null
directions; the corresponding columns of `P^{-1} (L^\\mathsf{T})^{-1}` form the
returned basis matrix.
"""
function nullspace_ldlt(X::AbstractMatrix; tol=nothing)
    M = ldlt!(ChordalLDLt(X), RowMaximum())
    D = M.D; L = M.L; P = M.P
    max_abs = maximum(i -> abs(D[i, i]), 1:size(D, 1); init=0.0)
    threshold = isnothing(tol) ? eps(Float64) * max(1.0, max_abs) : tol
    ind = findall(i -> abs(D[i, i]) <= threshold, 1:size(D, 1))
    U = zeros(size(D, 1), length(ind))
    for j in eachindex(ind)
        U[ind[j], j] = 1
    end
    return P \ (L' \ U)
end

"""    nullspace_ldlt(s::EuclideanSheaf; tol=nothing) -> Matrix

Convenience overload: computes the sheaf Laplacian ``L = d^\\mathsf{T} d`` (where
``d`` is the coboundary map of `s`) and delegates to `nullspace_ldlt(L)`.

The returned columns form a basis for the space of *global sections* of `s`.
"""
function nullspace_ldlt(s::EuclideanSheaf; tol=nothing)
    d = sparse(coboundary_map(s))
    return nullspace_ldlt(d' * d; tol=tol)
end


end # EuclideanSheaves