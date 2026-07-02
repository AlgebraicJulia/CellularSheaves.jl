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

See the [`CellularSheaves.NetworkSheaves.GraphHomomorphisms`](@ref) module for the underlying combinatorial
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
    # Precondition: `verts` must not contain duplicate entries.
    # If duplicates are present, `v_idx` keeps only the last occurrence and the
    # constructed fiber sheaf no longer matches the advertised row ordering.
    # No check is performed here to avoid the allocation cost in tight loops.
    stalk_dims = [get_vertex_stalk(F, v) for v in verts]
    total_dim = sum(stalk_dims; init=0)

    # When the fiber has no edges the coboundary is zero and every assignment
    # of vertex values is a global section, so the basis is the full identity.
    if isempty(fedges)
        return I(total_dim)
    end

    # Local vertex indices within the fiber sub-sheaf.
    v_idx = Dict(v => i for (i, v) in enumerate(verts))
    Fib = EuclideanSheaf{Float64}(stalk_dims)
    for (a, b) in fedges
        add_sheaf_edge!(Fib, v_idx[a], v_idx[b],
                        get_restriction_map(F, a, b),
                        get_restriction_map(F, b, a))
    end
    d_fib = coboundary_map(Fib)
    return nullspace_ldlt(d_fib' * d_fib)
end

"""
    all_fiber_bases(hom::GraphHomomorphism, F) -> Vector{Matrix}

Return a `Vector` of fiber-section basis matrices, one per target vertex
of `hom`.  Entry `tv` is the result of
`fiber_section_basis(F, fiber_vertices(hom, tv), fiber_edges(hom, g, tv))`
where `g = underlying_graph(F)` is derived from `F` internally.
"""
function all_fiber_bases(hom::GraphHomomorphism, F)
    g = underlying_graph(F)
    n_src = nv(g)
    length(hom.vertex_map) == n_src || throw(ArgumentError(
        "vertex_map has length $(length(hom.vertex_map)) but the sheaf's underlying graph has $n_src vertices"))
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

**Edge stalks and restriction maps.** For each target edge ``(p, q)``, all
source cross-edges mapping to ``(p, q)`` contribute to a single combined edge
stalk via direct sum.  That is, if cross-edges ``e_1, e_2, \\ldots`` all satisfy
``\\varphi(\\text{src}(e_i)) = p`` and ``\\varphi(\\text{dst}(e_i)) = q``, the
edge stalk of ``\\varphi_* F`` at ``(p, q)`` is
``F(e_1) \\oplus F(e_2) \\oplus \\cdots``, and the restriction maps are the
vertical concatenations of the per-edge pushforward restriction maps.

This ensures correctness when the edge map induced by `hom` is not injective
(multiple source edges mapping to the same target edge pair).
"""
function pushforward_sheaf(hom::GraphHomomorphism, F)
    g   = underlying_graph(F)
    phi = hom.vertex_map

    # Compute a basis for the global sections of each fiber sub-sheaf.
    fiber_bases = all_fiber_bases(hom, F)

    # The stalk dimension at target vertex tv equals the nullspace dimension.
    pf_stalk_dims = [size(fiber_bases[tv], 2) for tv in 1:hom.n_target]
    PfF = EuclideanSheaf{Float64}(pf_stalk_dims)

    # Group cross-edges by their canonical target edge pair (min, max) to handle
    # non-injective edge maps (multiple source edges → same target edge).
    edge_groups = Dict{Tuple{Int,Int}, Vector{Tuple{Matrix{Float64},Matrix{Float64}}}}()

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

        # Use canonical key (min, max) and store maps in (tri_u, tri_v) orientation.
        key = (min(tri_u, tri_v), max(tri_u, tri_v))
        if !haskey(edge_groups, key)
            edge_groups[key] = Tuple{Matrix{Float64},Matrix{Float64}}[]
        end
        if key == (tri_u, tri_v)
            push!(edge_groups[key], (rm_pf_u, rm_pf_v))
        else
            push!(edge_groups[key], (rm_pf_v, rm_pf_u))
        end
    end

    # Add one sheaf edge per target edge, combining multiple source edges via direct sum.
    # Iterate in a deterministic order so the pushed-forward sheaf's edge order
    # (and any downstream coboundary/edge-stalk indexing derived from it) is stable.
    for (p, q) in sort!(collect(keys(edge_groups)))
        rms = edge_groups[(p, q)]
        combined_rm_p = vcat([rm[1] for rm in rms]...)
        combined_rm_q = vcat([rm[2] for rm in rms]...)
        add_sheaf_edge!(PfF, p, q, combined_rm_p, combined_rm_q)
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

    n_pf = sum(pf_stalk_dims)
    n_F = sum(vertex_stalks(F))

    # Assemble only the nonzero fiber blocks; off-fiber blocks are structurally zero.
    I = Int[]
    J = Int[]
    V = Float64[]

    for tv in 1:hom.n_target
        fverts = fiber_vertices(hom, tv)
        # Gather the column indices of F's 0-cochain that belong to this fiber.
        f_rows = vcat([collect(F_cum[v]+1 : F_cum[v+1]) for v in fverts]...)
        B = fiber_bases[tv]
        # Project onto the fiber-basis coordinates via the pseudoinverse.
        P = pinv(B)
        row_offset = pf_cum[tv]

        for local_j in axes(P, 2)
            global_j = f_rows[local_j]
            for local_i in axes(P, 1)
                val = P[local_i, local_j]
                if !iszero(val)
                    push!(I, row_offset + local_i)
                    push!(J, global_j)
                    push!(V, val)
                end
            end
        end
    end

    return sparse(I, J, V, n_pf, n_F)
end

end # module Pushforwards
