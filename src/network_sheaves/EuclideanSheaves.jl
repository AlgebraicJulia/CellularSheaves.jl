# Module for network sheaves valued in Euclidean spaces with linear restriction maps
module EuclideanSheaves

export EuclideanSheaf, UnorderedPair, sheaf_laplacian_matrix, sheaf_laplacian_matrix_direct,
    restricted_laplacian_blocks, sheaf_from_graph, energy_function,
    nearest_global_section, edge_stalk_dimensions, nullspace_ldlt, harmonic_extension

using ArgCheck: @argcheck
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

"""   add_sheaf_edge!(s::EuclideanSheaf, v1::Int, v2::Int, rm1::Matrix{Float64}, rm2::Matrix{Float64})

Add an edge between vertices v1 and v2 with restriction maps rm1 and rm2.
rm1 is the restriction map from vertex v1 to the edge, and rm2 is the restriction map from vertex v2 to the edge.
"""
function add_sheaf_edge!(s::EuclideanSheaf{T}, v1::Int, v2::Int, rm1::Matrix{T}, rm2::Matrix{T}) where T
    @assert v1 <= length(s.vertex_stalks) && v2 <= length(s.vertex_stalks)
    @assert size(rm1)[1] == size(rm2)[1]
    @assert size(rm1)[2] == s.vertex_stalks[v1]
    @assert size(rm2)[2] == s.vertex_stalks[v2]

    stalk_size = size(rm1)[1]
    add_edge!(s.underlying_graph, v1, v2)
    edge_key = UnorderedPair(v1, v2)
    s.edge_stalks[edge_key] = stalk_size
    s.restriction_maps[v1=>v2] = rm1
    s.restriction_maps[v2=>v1] = rm2
    return ne(s.underlying_graph)
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

"""    sheaf_laplacian_matrix_direct(s::EuclideanSheaf) -> SparseMatrixCSC

Assemble the sheaf Laplacian ``L = d^\\mathsf{T} d`` directly from the restriction
maps of `s`, **without forming the coboundary matrix** ``d``.

For each edge ``e = (u, v)`` the following block contributions are accumulated:

```
L_{u,u} += ρ_{u→e}^T ρ_{u→e}
L_{v,v} += ρ_{v→e}^T ρ_{v→e}
L_{u,v} += -ρ_{u→e}^T ρ_{v→e}
L_{v,u} += -ρ_{v→e}^T ρ_{u→e}
```

This avoids constructing the ``(\\sum d_e) \\times (\\sum d_v)`` coboundary matrix
and the subsequent matrix multiply, saving both time and memory for large sheaves.
The result is numerically identical to `sheaf_laplacian_matrix`.
"""
function sheaf_laplacian_matrix_direct(s::EuclideanSheaf{T}) where T
    n_total = sum(s.vertex_stalks)
    offsets = [0; cumsum(s.vertex_stalks)]

    II = Int[]
    JJ = Int[]
    VV = T[]

    for e in edges(s.underlying_graph)
        u, v = src(e), dst(e)
        ρ_u = s.restriction_maps[u => v]   # d_e × d_u
        ρ_v = s.restriction_maps[v => u]   # d_e × d_v
        du  = s.vertex_stalks[u]
        dv  = s.vertex_stalks[v]
        ou  = offsets[u]
        ov  = offsets[v]

        # Diagonal block for u: ρ_u' * ρ_u
        Duu = ρ_u' * ρ_u
        for j in 1:du, i in 1:du
            push!(II, ou + i); push!(JJ, ou + j); push!(VV, Duu[i, j])
        end

        # Diagonal block for v: ρ_v' * ρ_v
        Dvv = ρ_v' * ρ_v
        for j in 1:dv, i in 1:dv
            push!(II, ov + i); push!(JJ, ov + j); push!(VV, Dvv[i, j])
        end

        # Off-diagonal blocks: -ρ_u' * ρ_v  and its transpose -ρ_v' * ρ_u
        Duv = -(ρ_u' * ρ_v)
        for j in 1:dv, i in 1:du
            push!(II, ou + i); push!(JJ, ov + j); push!(VV,  Duv[i, j])
            push!(II, ov + j); push!(JJ, ou + i); push!(VV,  Duv[i, j])
        end
    end

    return sparse(II, JJ, VV, n_total, n_total)
end

"""    restricted_laplacian_blocks(s::EuclideanSheaf,
                                interior::AbstractVector{Int},
                                boundary::AbstractVector{Int})
        -> (L_II::SparseMatrixCSC, L_IB::SparseMatrixCSC)

Compute the interior-interior and interior-boundary blocks of the sheaf Laplacian
**without assembling the full Laplacian**.

`interior` and `boundary` are disjoint vectors of vertex indices (any order).
Only edges with at least one interior endpoint contribute.

- `L_II` is the ``(\\sum_{v \\in I} d_v) \\times (\\sum_{v \\in I} d_v)``
  symmetric positive-semidefinite block.
- `L_IB` is the ``(\\sum_{v \\in I} d_v) \\times (\\sum_{v \\in B} d_v)``
  interior-boundary coupling block.

The DOF ordering within each block follows the order of `interior` / `boundary`.
Edges entirely within the boundary subgraph are skipped, so assembly cost is
``O(|E_I| \\cdot d_{\\max}^2)`` where ``E_I`` is the set of edges with at least one
interior endpoint.
"""
function restricted_laplacian_blocks(s::EuclideanSheaf{T},
                                      interior::AbstractVector{Int},
                                      boundary::AbstractVector{Int}) where T
    interior_pos = Dict(v => i for (i, v) in enumerate(interior))
    boundary_pos = Dict(v => i for (i, v) in enumerate(boundary))

    n_I = isempty(interior) ? 0 : sum(s.vertex_stalks[v] for v in interior)
    n_B = isempty(boundary) ? 0 : sum(s.vertex_stalks[v] for v in boundary)

    I_offsets = isempty(interior) ? Int[] : [0; cumsum([s.vertex_stalks[v] for v in interior])]
    B_offsets = isempty(boundary) ? Int[] : [0; cumsum([s.vertex_stalks[v] for v in boundary])]

    II_rows = Int[]; II_cols = Int[]; II_vals = T[]
    IB_rows = Int[]; IB_cols = Int[]; IB_vals = T[]

    for e in edges(s.underlying_graph)
        u, v = src(e), dst(e)

        u_pos_I = get(interior_pos, u, 0)
        v_pos_I = get(interior_pos, v, 0)
        u_pos_B = get(boundary_pos, u, 0)
        v_pos_B = get(boundary_pos, v, 0)

        # Skip edges entirely in the boundary subgraph
        (u_pos_I == 0 && v_pos_I == 0) && continue

        ρ_u = s.restriction_maps[u => v]   # d_e × d_u
        ρ_v = s.restriction_maps[v => u]   # d_e × d_v
        du  = s.vertex_stalks[u]
        dv  = s.vertex_stalks[v]

        # Diagonal block L_II[u,u] += ρ_u' * ρ_u
        if u_pos_I > 0
            i_u = I_offsets[u_pos_I]
            Duu = ρ_u' * ρ_u
            for j in 1:du, i in 1:du
                push!(II_rows, i_u + i); push!(II_cols, i_u + j); push!(II_vals, Duu[i, j])
            end
        end

        # Diagonal block L_II[v,v] += ρ_v' * ρ_v
        if v_pos_I > 0
            i_v = I_offsets[v_pos_I]
            Dvv = ρ_v' * ρ_v
            for j in 1:dv, i in 1:dv
                push!(II_rows, i_v + i); push!(II_cols, i_v + j); push!(II_vals, Dvv[i, j])
            end
        end

        # Off-diagonal L_II[u,v] and L_II[v,u] when both endpoints are interior
        if u_pos_I > 0 && v_pos_I > 0
            i_u = I_offsets[u_pos_I]
            i_v = I_offsets[v_pos_I]
            Duv = -(ρ_u' * ρ_v)
            for j in 1:dv, i in 1:du
                push!(II_rows, i_u + i); push!(II_cols, i_v + j); push!(II_vals,  Duv[i, j])
                push!(II_rows, i_v + j); push!(II_cols, i_u + i); push!(II_vals,  Duv[i, j])
            end
        end

        # L_IB[u, v] += -ρ_u' * ρ_v  when u ∈ I, v ∈ B
        if u_pos_I > 0 && v_pos_B > 0
            i_u = I_offsets[u_pos_I]
            b_v = B_offsets[v_pos_B]
            Duv_IB = -(ρ_u' * ρ_v)
            for j in 1:dv, i in 1:du
                push!(IB_rows, i_u + i); push!(IB_cols, b_v + j); push!(IB_vals, Duv_IB[i, j])
            end
        end

        # L_IB[v, u] += -ρ_v' * ρ_u  when v ∈ I, u ∈ B
        if v_pos_I > 0 && u_pos_B > 0
            i_v = I_offsets[v_pos_I]
            b_u = B_offsets[u_pos_B]
            Dvu_IB = -(ρ_v' * ρ_u)
            for j in 1:du, i in 1:dv
                push!(IB_rows, i_v + i); push!(IB_cols, b_u + j); push!(IB_vals, Dvu_IB[i, j])
            end
        end
    end

    L_II = sparse(II_rows, II_cols, II_vals, n_I, n_I)
    L_IB = sparse(IB_rows, IB_cols, IB_vals, n_I, n_B)
    return L_II, L_IB
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
    return nullspace_ldlt(sheaf_laplacian_matrix_direct(s); tol=tol)
end

# Private helper: given a ChordalLDLt factorization `M` of a symmetric
# positive-semidefinite matrix and a right-hand side `b`, return both the
# minimum-norm particular solution `x_p` satisfying `M * x_p ≈ b` (in the
# pseudoinverse sense) and a matrix `null_vecs` whose columns span the
# nullspace of the original matrix.
function ldlt_pseudoinverse_and_null(M, b; tol=nothing)
    D = M.D; Lfac = M.L; P = M.P
    n = size(D, 1)
    max_abs = maximum(i -> abs(D[i, i]), 1:n; init=0.0)
    threshold = isnothing(tol) ? eps(Float64) * max(1.0, max_abs) : tol

    null_idx = findall(i -> abs(D[i, i]) <= threshold, 1:n)

    # Null basis — identical to nullspace_ldlt
    U = zeros(n, length(null_idx))
    for j in eachindex(null_idx)
        U[null_idx[j], j] = 1.0
    end
    null_vecs = P \ (Lfac' \ U)

    # Particular solution via pseudoinverse: suppress null directions in D
    c = P' \ b          # permute rhs
    z = Lfac \ c        # forward solve
    w = zeros(n)
    for i in setdiff(1:n, null_idx)
        w[i] = z[i] / D[i, i]  # D⁺ on free directions only
    end
    x_p = P \ (Lfac' \ w)  # back solve + undo permutation

    return x_p, null_vecs
end

"""
    harmonic_extension(s::EuclideanSheaf, boundary::Dict{Int,<:AbstractVector})
        -> (BlockVector, Matrix)

Compute the harmonic extension of `boundary` over the interior vertices of `s`.

`boundary` maps each boundary vertex index to a vector of length `s.vertex_stalks[v]`.

Returns `(x_p, null_basis)` where:
- `x_p` is the minimum-norm particular solution as a `BlockArray` over `s.vertex_stalks`
- `null_basis` is a matrix whose columns span the solution's indeterminate directions,
  embedded in the full cochain space (zero at boundary dofs). Has 0 columns when the
  boundary conditions uniquely determine the harmonic extension.

The complete solution set is `{ x_p + null_basis * c : c ∈ Rᵏ }` where
`k = size(null_basis, 2)`.

Throws `ArgumentError` if any boundary vector has the wrong length for its stalk.
"""
function harmonic_extension(s::EuclideanSheaf, boundary::Dict{Int,<:AbstractVector})
    nv_graph = nv(s.underlying_graph)
    for (v, val) in boundary
        @argcheck 1 <= v <= nv_graph
        @argcheck length(val) == s.vertex_stalks[v]
    end

    offsets        = [0; cumsum(s.vertex_stalks)]
    boundary_verts = sort(collect(keys(boundary)))
    interior_verts = setdiff(1:nv_graph, boundary_verts)
    n_total        = sum(s.vertex_stalks)

    # Short-circuit: no interior dofs — nothing to solve
    if isempty(interior_verts)
        x_p = BlockArray(vcat([boundary[v] for v in 1:nv_graph]...), s.vertex_stalks)
        return x_p, zeros(n_total, 0)
    end

    I_idx = vcat([offsets[v]+1:offsets[v+1] for v in interior_verts]...)

    B_idx = isempty(boundary_verts) ? Int[] :
        vcat([offsets[v]+1:offsets[v+1] for v in boundary_verts]...)

    # Assemble L_II and L_IB directly without materialising the full Laplacian
    L_II, L_IB = restricted_laplacian_blocks(s, interior_verts, boundary_verts)

    if isempty(boundary_verts)
        b = zeros(length(I_idx))
    else
        x_B = vcat([boundary[v] for v in boundary_verts]...)
        b   = -L_IB * x_B
    end

    x_interior, null_interior = ldlt_pseudoinverse_and_null(
        ldlt!(ChordalLDLt(L_II), RowMaximum()), b)

    x        = zeros(n_total)
    x[I_idx] = x_interior
    if !isempty(boundary_verts)
        x[B_idx] = vcat([boundary[v] for v in boundary_verts]...)
    end

    null_basis          = zeros(n_total, size(null_interior, 2))
    null_basis[I_idx,:] = null_interior

    return BlockArray(x, s.vertex_stalks), null_basis
end


end # EuclideanSheaves