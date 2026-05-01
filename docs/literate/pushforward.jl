# # Pushforward of a Sheaf along a Graph Homomorphism
#
# This example shows how to compute the *pushforward* of a cellular sheaf
# along a graph homomorphism, and demonstrates that the global sections
# (nullspace of the sheaf Laplacian) of the original sheaf and the
# pushed-forward sheaf are isomorphic as vector spaces.
#
# We start with sheaf `F` — a 6-cycle with 3-dimensional vertex stalks and
# 2-dimensional edge stalks — then collapse paired vertices
# ``1,2 \mapsto 1``, ``3,4 \mapsto 2``, ``5,6 \mapsto 3``
# to obtain a triangle.

using CellularSheaves
using Graphs
using LinearAlgebra
using SparseArrays

# ## Source sheaf F on a 6-cycle
#
# Each vertex stalk is ``\mathbb{R}^3`` and each edge stalk is ``\mathbb{R}^2``.
# The restriction map at every vertex-edge incidence is the ``2 \times 3``
# matrix ``I_{2 \times 3}``, which projects onto the first two coordinates.

n = 6; stalk_dim = 3

F = EuclideanSheaf{Float64}(Int[])
for i in 1:n
    add_vertex_stalk!(F, stalk_dim)
end
for i in 1:n
    v1 = i; v2 = (i % n) + 1
    rm = Matrix{Float64}(I, stalk_dim - 1, stalk_dim)   # 2×3 projection
    add_sheaf_edge!(F, v1, v2, rm, rm)
end

println(F)

# ## Nullspace of F via LDLt
#
# The global sections of a sheaf are the nullspace of the sheaf Laplacian
# ``L = d^\mathsf{T} d`` where ``d`` is the coboundary map.
# We extract them using `nullspace_ldlt`, which applies the sparse LDLt
# factorisation from CliqueTrees.

NS_F = nullspace_ldlt(F)
println("Nullspace dimension of F: ", size(NS_F, 2))

# A global section ``(s_1,\ldots,s_6)`` of `F` must satisfy
# ``\rho \, s_i = \rho \, s_j`` for every edge ``(i,j)``, i.e.
# ``s_i[1{:}2] = s_j[1{:}2]``.  Since the 6-cycle is connected, all vertices
# share a common ``(a, b) \in \mathbb{R}^2``, while each ``s_i[3]`` is free.
# Nullspace dimension ``= 2 + 6 = 8``.

# ## Graph homomorphism ``\varphi``: 6-cycle → triangle
#
# We define ``\varphi`` by the vertex map that merges pairs of adjacent vertices:
# ``1,2 \mapsto 1``, ``3,4 \mapsto 2``, ``5,6 \mapsto 3``.

hom = GraphHomomorphism([1, 1, 2, 2, 3, 3])
println(hom)

# The `GraphHomomorphism` type tracks which source edges are *fiber edges*
# (both endpoints collapse to the same target vertex) and which are *cross edges*
# (endpoints in different fibers).  Use `fiber_vertices`, `fiber_edges`, and
# `cross_edges` to inspect the structure of ``\varphi``.

g = underlying_graph(F)
println("Fiber vertices per triangle vertex:")
for tv in 1:hom.n_target
    println("  tv=$tv: ", fiber_vertices(hom, tv))
end

println("Fiber edges per triangle vertex:")
for tv in 1:hom.n_target
    println("  tv=$tv: ", fiber_edges(hom, g, tv))
end

println("Cross edges: ",
        [(string(e), tu, tv) for (e, tu, tv) in cross_edges(hom, g)])

# Each fiber consists of two source vertices joined by one fiber edge: the
# 6-cycle edges ``(1,2)``, ``(3,4)``, ``(5,6)`` are fiber edges while
# ``(1,6)``, ``(2,3)``, ``(4,5)`` are cross edges that become the three
# edges of the target triangle.

# ## Pushforward sheaf ``\varphi_* F``
#
# The pushforward is constructed by `pushforward_sheaf`:
#
# - **Vertex stalk** at target vertex `tv`:
#   the space of global sections of `F` restricted to fiber `tv`,
#   computed via `nullspace_ldlt`.  Each fiber has two vertices
#   (``\mathbb{R}^3`` each) joined by one edge (``\mathbb{R}^2``), so
#   ``\dim = 3 + 3 - 2 = 4``.
#
# - **Restriction maps** for each cross edge: compose `F`'s original
#   restriction map with the fiber-basis embedding.

PfF = pushforward_sheaf(hom, F)

println(PfF)
println("Vertex stalk dims of φ_*F: ", vertex_stalks(PfF))
println("Edge stalk dims of φ_*F:   ", edge_stalk_dimensions(PfF))

# ## Nullspace of ``\varphi_* F``

NS_PfF = nullspace_ldlt(PfF)
println("Nullspace dimension of φ_*F: ", size(NS_PfF, 2))

# A global section ``(t_1, t_2, t_3)`` of ``\varphi_* F`` must satisfy
# ``t_i[1{:}2] = t_j[1{:}2]`` for every triangle edge ``(i,j)``.
# Connectedness forces all three vertex stalks to share
# ``(a, b) \in \mathbb{R}^2``, with each ``t_i[3{:}4]`` free.
# Nullspace dimension ``= 2 + 3 \times 2 = 8``.

# ## The two nullspaces are isomorphic
#
# The pushforward construction ensures a bijection between global sections of
# `F` and global sections of ``\varphi_* F``.  We verify this with the
# explicit linear transfer map ``T : C^0(F) \to C^0(\varphi_* F)`` returned
# by `pushforward_transfer_map`.
#
# ``T`` expresses each fiber cochain in the fiber-basis coordinates via a
# pseudoinverse; the identity ``d_{\varphi_* F} \circ T \circ s = 0`` holds
# for every global section ``s`` of `F`.

T = pushforward_transfer_map(hom, F)

d_PfF = sparse(coboundary_map(PfF))
mapped   = T * NS_F
residual = norm(d_PfF * mapped)
println("Residual ‖d_{φ_*F} ∘ T ∘ NS_F‖ = ", residual)

# The near-zero residual confirms that every global section of `F` maps to a
# global section of ``\varphi_* F`` under `T`, establishing the isomorphism.

size(NS_F, 2) == size(NS_PfF, 2)

