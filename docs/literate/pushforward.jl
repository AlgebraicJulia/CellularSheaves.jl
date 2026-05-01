# # Pushforward of a Sheaf along a Graph Homomorphism
#
# This example shows how to compute the *pushforward* of a cellular sheaf
# along a graph homomorphism, and demonstrates that the global sections
# (nullspace of the sheaf Laplacian) of the original sheaf and the
# pushed-forward sheaf are isomorphic as vector spaces.
#
# We start with sheaf `G` from the `sheaf_morphisms` example—a 6-cycle with
# 3-dimensional vertex stalks and 2-dimensional edge stalks—then collapse
# paired vertices (1,2)→1, (3,4)→2, (5,6)→3 to obtain a triangle.

using CellularSheaves
using CliqueTrees.Multifrontal
using Graphs
using LinearAlgebra
using SparseArrays

# ## Source sheaf G on a 6-cycle
#
# Each vertex stalk is ℝ³ and each edge stalk is ℝ².  The restriction map at
# every vertex-edge incidence is the 2×3 matrix `I[1:2, 1:3]`, which projects
# onto the first two coordinates.

n = 6; stalk_dim = 3

G = EuclideanSheaf{Float64}(Int[])
for i in 1:n
    add_vertex_stalk!(G, stalk_dim)
end
for i in 1:n
    v1 = i; v2 = (i % n) + 1
    rm = Matrix{Float64}(I, stalk_dim - 1, stalk_dim)   # 2×3 projection
    add_sheaf_edge!(G, v1, v2, rm, rm)
end

println(G)

# ## Nullspace of G via LDLt
#
# The global sections of a sheaf are the nullspace of the sheaf Laplacian
# `L = Bᵀ B` where `B` is the coboundary map.  We extract them using the
# sparse LDLt factorisation from CliqueTrees: columns of `D` with a zero
# (or near-zero) diagonal entry identify null directions.

function nullspace_ldlt(X::AbstractMatrix)
    M = ldlt!(ChordalLDLt(X), RowMaximum())
    D = M.D; L = M.L; P = M.P
    max_abs = maximum(i -> abs(D[i, i]), 1:size(D, 1); init=0.0)
    tol = eps(Float64) * max(1.0, max_abs)
    ind = findall(i -> abs(D[i, i]) <= tol, 1:size(D, 1))
    U = zeros(size(D, 1), length(ind))
    for j in eachindex(ind); U[ind[j], j] = 1; end
    return P \ (L' \ U)
end

CG = sparse(coboundary_map(G))
NS_G = nullspace_ldlt(CG' * CG)
println("Nullspace dimension of G: ", size(NS_G, 2))

# A global section `(s₁,…,s₆)` of `G` must satisfy `rm · sᵢ = rm · sⱼ` for
# every edge `(i,j)`, i.e. `sᵢ[1:2] = sⱼ[1:2]`.  Since the 6-cycle is
# connected, all vertices share a common `(a, b) ∈ ℝ²`, while each `sᵢ[3]`
# is free.  Nullspace dimension = 2 + 6 = **8**.

# ## Graph homomorphism φ: 6-cycle → triangle
#
# We define φ by a vertex map that merges pairs of adjacent vertices:
# `1,2 → 1`,  `3,4 → 2`,  `5,6 → 3`.

phi_v = [1, 1, 2, 2, 3, 3]   # phi_v[i] = image of vertex i in triangle

# Classify each source edge as either a *fiber edge* (both endpoints map to
# the same triangle vertex) or a *cross edge* (endpoints map to different
# triangle vertices).  Fiber edges stay "inside" a vertex stalk; cross edges
# become the edges of the triangle.

fiber_verts = [findall(==(tv), phi_v) for tv in 1:3]
fiber_edges = [Tuple{Int,Int}[] for _ in 1:3]
cross_edges = Tuple{Any,Int,Int}[]

for e in edges(underlying_graph(G))
    u, v = src(e), dst(e)
    pu, pv = phi_v[u], phi_v[v]
    if pu == pv
        push!(fiber_edges[pu], (u, v))
    else
        push!(cross_edges, (e, pu, pv))
    end
end

println("Fiber vertices per triangle vertex: ", fiber_verts)
println("Fiber edges per triangle vertex:    ", fiber_edges)
println("Cross edges (source edge, tri_src, tri_dst): ",
        [(string(e), tu, tv) for (e, tu, tv) in cross_edges])

# Each fiber consists of two vertices connected by one edge: the 6-cycle edges
# (1,2), (3,4), (5,6) are fiber edges; (1,6), (2,3), (4,5) are cross edges.

# ## Vertex stalks of φ_*G: global sections of each fiber
#
# The stalk `(φ_*G)(v)` equals the space of global sections of `G` restricted
# to the sub-sheaf over the fiber of `v`.  We build each fiber as a mini
# `EuclideanSheaf` and call `nullspace_ldlt` to recover a basis.

function fiber_global_sections_basis(G, verts, edges_list)
    stalk_dims = [get_vertex_stalk(G, v) for v in verts]
    v_idx = Dict(v => i for (i, v) in enumerate(verts))
    Fib = EuclideanSheaf{Float64}(stalk_dims)
    for (a, b) in edges_list
        add_sheaf_edge!(Fib, v_idx[a], v_idx[b],
                        get_restriction_map(G, a, b),
                        get_restriction_map(G, b, a))
    end
    C = sparse(coboundary_map(Fib))
    return nullspace_ldlt(C' * C)
end

fiber_bases = [fiber_global_sections_basis(G, fiber_verts[tv], fiber_edges[tv])
               for tv in 1:3]

for tv in 1:3
    println("Pushforward stalk dim at triangle vertex $tv: ",
            size(fiber_bases[tv], 2))
end

# Each fiber has two vertices (dim 3 each) joined by one edge (dim 2, with the
# 2×3 projection map): global section dimension = 3 + 3 − 2 = **4**.
# A global section `(s_a, s_b)` over a fiber satisfies `s_a[1:2] = s_b[1:2]`,
# so it is parameterised by `(a, b, s_a[3], s_b[3]) ∈ ℝ⁴`.

# ## Restriction maps of φ_*G
#
# For each cross edge `e = (u, v)` mapping to triangle edge `(tri_u, tri_v)`:
# - **edge stalk**: equals `G`'s stalk at `e` (ℝ²).
# - **restriction map** from `(φ_*G)(tri_u)` to the edge stalk: compose `G`'s
#   restriction map at the endpoint `u ∈ fiber(tri_u)` with the fiber basis
#   (projection of the global section onto vertex `u`'s rows).

# ## Construct φ_*G

pf_stalk_dims = [size(fiber_bases[tv], 2) for tv in 1:3]
PfG = EuclideanSheaf{Float64}(pf_stalk_dims)

for (e, tri_u, tri_v) in cross_edges
    u, v = src(e), dst(e)
    pu, pv = phi_v[u], phi_v[v]
    # Identify which endpoint belongs to which triangle vertex.
    src_u, src_v = pu == tri_u ? (u, v) : (v, u)

    # Row offsets of src_u and src_v within their respective fiber bases.
    offset_u = sum(get_vertex_stalk(G, w)
                   for w in fiber_verts[tri_u] if w < src_u; init=0)
    offset_v = sum(get_vertex_stalk(G, w)
                   for w in fiber_verts[tri_v] if w < src_v; init=0)
    d_u = get_vertex_stalk(G, src_u)
    d_v = get_vertex_stalk(G, src_v)

    rm_u = get_restriction_map(G, src_u, src_v)   # G: vertex src_u → edge
    rm_v = get_restriction_map(G, src_v, src_u)   # G: vertex src_v → edge

    # Pushforward restriction maps = G's rm composed with fiber basis rows.
    rm_pf_u = rm_u * fiber_bases[tri_u][offset_u+1 : offset_u+d_u, :]
    rm_pf_v = rm_v * fiber_bases[tri_v][offset_v+1 : offset_v+d_v, :]

    add_sheaf_edge!(PfG, tri_u, tri_v, rm_pf_u, rm_pf_v)
end

println(PfG)
println("Vertex stalk dims of φ_*G: ", vertex_stalks(PfG))
println("Edge stalk dims of φ_*G:   ", edge_stalk_dimensions(PfG))

# ## Nullspace of φ_*G

CPfG = sparse(coboundary_map(PfG))
NS_PfG = nullspace_ldlt(CPfG' * CPfG)
println("Nullspace dimension of φ_*G: ", size(NS_PfG, 2))

# A global section `(t₁, t₂, t₃)` of `φ_*G` must satisfy `tᵢ[1:2] = tⱼ[1:2]`
# for every triangle edge `(i,j)`.  Connectedness forces all three vertex
# stalks to share `(a, b) ∈ ℝ²`, with each `tᵢ[3:4]` free.
# Nullspace dimension = 2 + 3×2 = **8**.

# ## The two nullspaces are isomorphic
#
# The pushforward construction ensures that global sections of `G` biject with
# global sections of `φ_*G`.  We verify this by building the explicit linear
# transfer map `T` and checking that its image lands in the nullspace of `φ_*G`.
#
# `T` maps `G`'s 0-cochains (ℝ¹⁸) to `φ_*G`'s 0-cochains (ℝ¹²) by expressing
# each fiber cochain in the fiber-basis coordinates via a pseudo-inverse.

G_cum  = [0; cumsum(vertex_stalks(G))]
pf_cum = [0; cumsum(pf_stalk_dims)]
T = zeros(sum(pf_stalk_dims), sum(vertex_stalks(G)))

for tv in 1:3
    fverts = fiber_verts[tv]
    g_rows = vcat([collect(G_cum[v]+1 : G_cum[v+1]) for v in fverts]...)
    B = fiber_bases[tv]
    T[pf_cum[tv]+1 : pf_cum[tv+1], g_rows] = pinv(B)
end

# Apply T to the nullspace basis of G and check the coboundary residual.
mapped  = T * NS_G
residual = norm(CPfG * mapped)
println("Residual ‖d_{φ_*G} ∘ T ∘ NS_G‖ = ", residual)

# The near-zero residual confirms that every global section of `G` maps to a
# global section of `φ_*G` under `T`, establishing the isomorphism.

size(NS_G, 2) == size(NS_PfG, 2)
