# # Nearest Global Section: Iterative Method
#
# This example demonstrates the iterative conjugate-gradient (CG) backend for
# projecting a 0-cochain onto the nearest global section.  The CG solver is
# stochastic in its convergence behaviour and may occasionally fail to find a
# global section for ill-conditioned sheaves or small tolerances.  The example
# below wraps the call in a `try`/`catch` so that a docs build is never broken
# by a solver failure.
#
# For production use, prefer [`nearest_global_section`](@ref) with
# `method=:ldl` (the direct LDLt backend), which is deterministic and
# typically much faster for small-to-medium sheaves.

using CellularSheaves
using LinearAlgebra

# ## Building the Sheaf
#
# We build the same identity-restriction-map sheaf on a 4-cycle that appears
# in the main *Code Example*, but this time with 3-dimensional stalks so the
# sheaf Laplacian is slightly larger.

n = 4
stalk_dim = 3

s = EuclideanSheaf{Float64}(Int[])
for _ in 1:n
    add_vertex_stalk!(s, stalk_dim)
end
for i in 1:n
    v1 = i
    v2 = mod1(i + 1, n)
    rm = Matrix{Float64}(I, stalk_dim, stalk_dim)
    add_sheaf_edge!(s, v1, v2, rm, rm)
end

println(s)

# ## Running the Iterative Solver
#
# `nearest_global_section` with `method=:cg` uses the Krylov CG method from
# Krylov.jl.  Because the edge-space operator ``E = d d^\top`` is positive
# semi-definite (not strictly positive definite when the sheaf has non-trivial
# global sections), the iterative solver can sometimes stagnate.  We therefore
# wrap the call and print any error rather than letting it propagate.

x0 = rand(sum(vertex_stalks(s)))

gs = try
    nearest_global_section(s, x0; method=:cg, verbose=true)
catch e
    println("Iterative solver did not converge: ", e)
    nothing
end

if gs !== nothing
    d = coboundary_map(s)
    println("‖d · gs‖ = ", norm(d * gs))
else
    println("No global section returned by the iterative solver.")
end

# ## Comparison with the Direct Method
#
# For reference, the deterministic LDLt backend always succeeds:

gs_ldl = nearest_global_section(s, x0; method=:ldl)
d = coboundary_map(s)
println("‖d · gs_ldl‖ = ", norm(d * gs_ldl))
