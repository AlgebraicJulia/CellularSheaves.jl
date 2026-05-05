# # Computing the Nullspace of the Sheaf Laplacian

# ## Introduction
#
# This example demonstrates how to compute and visualize the nullspace of the Laplacian operator
# for a circular cellular sheaf. We construct a system of 6 agents, each with 3-dimensional stalks,
# arranged in a cycle topology. The nullspace of the Laplacian contains the global sections---the
# solutions that are preserved by the sheaf Laplacian operator. Understanding the structure of this
# nullspace is fundamental to understanding the harmonic analysis and signal processing on cellular sheaves.

using CellularSheaves
using LinearAlgebra
using Plots

# ## Visualization Helper Function
#
# We define a generalized function to visualize nullspace basis vectors for any sheaf with 3D stalks.
# This function creates a 3D scatter plot with each basis vector shown in a distinct color.

function plot_nullspace_basis(V, n, stalk_dim; title="Nullspace basis vectors")
    num_basis_vectors = size(V, 2)
    colors = palette(:tab10, num_basis_vectors)

    # Extract coordinate components for each agent
    # Coordinates are strided in V with period stalk_dim
    x_coords = V[1:stalk_dim:end, :]  # first row in each stalk_dim-sized block
    y_coords = V[2:stalk_dim:end, :]  # second row in each stalk_dim-sized block
    z_coords = V[3:stalk_dim:end, :]  # third row in each stalk_dim-sized block

    # Create 3D scatter plot
    p = scatter(legend=true, xlabel="x", ylabel="y", zlabel="z", title=title,
        markersize=8, markerstrokewidth=0.5, camera=(45, 30))

    # Plot each basis vector with color coding
    for b in 1:num_basis_vectors
        x_vals = x_coords[:, b]
        y_vals = y_coords[:, b]
        z_vals = z_coords[:, b]

        scatter!(p, x_vals, y_vals, z_vals, color=colors[b], label="basis $b",
            markersize=6, markerstrokewidth=0.5)
    end

    p
end

# # Example 1: Identity Restriction Maps
#
# In this first example, we demonstrate the nullspace of a sheaf where all edges are fully constrained
# by identity restriction maps. This means that adjacent agents must agree on all components of their
# stalk values. The resulting nullspace consists of *constant* sections---vectors that have the same
# value at every agent. For a 3-dimensional stalk, this gives a 3-dimensional nullspace.
#
# ## Building the Sheaf
#
# We start by creating a Euclidean sheaf with 6 vertices (agents), each equipped with a 3-dimensional
# vector space as its stalk. The stalks represent local data or state at each agent.

n = 6
stalk_dim = 3

F = EuclideanSheaf{Float64}(Int[])

for i in 1:n
    add_vertex_stalk!(F, stalk_dim)
end

# Next, we add edges connecting the agents in a cycle. Each edge is equipped with the identity
# restriction map, meaning that sections on adjacent agents must have the same value where they overlap.

for i in 1:n
    v1 = i
    v2 = (i % n) + 1
    rm = Matrix{Float64}(I, stalk_dim, stalk_dim)
    add_sheaf_edge!(F, v1, v2, rm, rm)
end

# ## Computing the Nullspace
#
# We call `nullspace_ldlt` directly on the sheaf. Internally it forms the
# Laplacian ``X = d^\mathsf{T} d`` from the coboundary map ``d``, factorises it as
# ``X = P^\mathsf{T} L D L^\mathsf{T} P`` via a sparse Chordal LDLt decomposition,
# and returns a basis for the nullspace from the zero-diagonal entries of ``D``.
# The nullspace of the Laplacian is exactly the space of *global sections* of the sheaf.

V = nullspace_ldlt(F)

# ## Visualizing the Nullspace
#
# We use our visualization helper function to show the three orthogonal projections of the nullspace
# basis vectors. Notice that the basis vectors are highly structured: they represent constant sections
# where each coordinate direction is independently constant across all agents.

plot_nullspace_basis(V, n, stalk_dim; title="Example 1: 3D Nullspace (fully constrained)")


# # Example 2: Partial Restriction Maps
#
# ## Effect of Partial Constraints
#
# In this second example, we modify the edge structure to allow one component of the stalk to vary
# freely across edges. We keep edge stalks at 3-dimensional but restrict them only in the x and y
# components using the restriction map ``[1\ 0\ 0;\ 0\ 1\ 0]``. This leaves the z component
# unconstrained on edges---agents can have different z values even on adjacent edges.
#
# This partial constraint has a dramatic effect on the nullspace: in addition to the 3 constant
# sections, we now get 3 additional basis vectors representing the free z degrees of freedom at
# each agent. The result is a 6-dimensional nullspace.
#
# ## Building the Modified Sheaf

F2 = EuclideanSheaf{Float64}(Int[])

for i in 1:n
    add_vertex_stalk!(F2, stalk_dim)
end

# Add edges with partial restriction maps: only x and y are constrained
for i in 1:n
    v1 = i
    v2 = (i % n) + 1 # Restriction map that selects only x and y components
    rm = [1.0 0.0 0.0; 0.0 1.0 0.0]
    add_sheaf_edge!(F2, v1, v2, rm, rm)
end

# ## Computing the Nullspace with Partial Constraints

V2 = nullspace_ldlt(F2)

# ## Visualizing the Extended Nullspace
#
# With partial constraints, the nullspace is now 6-dimensional. We can see this in the visualization:
# the basis vectors are more spread out, and we see additional degrees of freedom particularly in
# the z-component projections (y vs z and x vs z).

plot_nullspace_basis(V2, n, stalk_dim; title="Example 2: 6D Nullspace (z unconstrained)")
