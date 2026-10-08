# # Schwarz Domain Decomposition on a Cellular Sheaf
#
# Overlapping Schwarz methods solve a PDE on a large domain by repeatedly
# solving it on smaller, overlapping subdomains, each using its neighbours'
# current values as boundary data. This example sets the method up on a
# cellular sheaf:
#
# - the **vertices** are the subdomains ``\Omega_i``, each with stalk
#   ``\mathbb{R}^{\Omega_i}`` (the unknowns it owns a copy of);
# - the **edges** are the nonempty overlaps ``\Omega_i \cap \Omega_j``, with
#   stalk ``\mathbb{R}^{\Omega_i \cap \Omega_j}``;
# - the **restriction maps** read off a subdomain's values on an overlap.
#
# A Schwarz iterate is a 0-cochain of this sheaf. The coboundary
# ``\delta x`` measures how much overlapping subdomains disagree, and a
# converged iterate is a global section that glues into the solution.

using CellularSheaves
using LinearAlgebra
using SparseArrays
using Graphs
using Plots

# ## The PDE
#
# We solve the Poisson problem ``-\Delta u = f`` on the unit square with
# ``u = 0`` on the boundary, using the 5-point finite-difference Laplacian on an
# ``m \times m`` interior grid.

function poisson2d(m)
    h = 1 / (m + 1)
    T = spdiagm(-1 => fill(-1.0, m - 1), 0 => fill(2.0, m), 1 => fill(-1.0, m - 1))
    Id = sparse(1.0I, m, m)
    return (kron(Id, T) + kron(T, Id)) / h^2
end

m = 40
A = poisson2d(m)
xs = range(0, 1; length=m + 2)[2:end-1]
f = [8π^2 * sin(2π * x) * sin(2π * y) + 10 for y in xs for x in xs]
u_exact = A \ f;

# ## The decomposition
#
# We cut the grid into a ``3 \times 3`` array of boxes and grow each box by two
# layers of grid neighbours to get overlapping subdomains.

function box_partition(m, p)
    parts = Vector{Int}(undef, m * m)
    for jy in 1:m, jx in 1:m
        parts[(jy - 1) * m + jx] = (cld(jy * p, m) - 1) * p + cld(jx * p, m)
    end
    return parts
end

parts = box_partition(m, 3)
subdomains = overlapping_subdomains(A, parts; overlap=2)
dd = SchwarzDecomposition(A, subdomains; owner=parts)

# The overlap graph is the underlying graph of the sheaf. Diagonal boxes share
# a corner of overlap, so each interior subdomain has eight neighbours.

s = overlap_sheaf(dd)
(nv(underlying_graph(s)), ne(underlying_graph(s)))

# Global sections of this sheaf are exactly the functions on the whole grid:
# a cochain is a section precisely when all overlapping copies agree.

overlap_disagreement(dd, localize(dd, u_exact))

# ## Solving
#
# The **multiplicative** (alternating) method visits subdomains one at a time
# and pushes each new local solution through the restriction maps onto its
# neighbours, so the iterate stays a global section. The **parallel** (Lions)
# method solves every subdomain at once from the previous cochain. Its copies
# disagree on overlaps until convergence.

mult = schwarz_solve(dd, f; method=:multiplicative, tol=1e-10)
par = schwarz_solve(dd, f; method=:parallel, tol=1e-10)
(mult.iterations, par.iterations)

#-

norm(mult.u - u_exact) / norm(u_exact), norm(par.u - u_exact) / norm(u_exact)

# The residual of the glued iterate, and the disagreement ``\|\delta x\|`` on
# the overlaps, for both methods:

p1 = plot(mult.residuals; yscale=:log10, label="multiplicative", xlabel="sweep",
    ylabel="relative residual", lw=2)
plot!(p1, par.residuals; label="parallel", lw=2)
p2 = plot(par.disagreements[2:end] ./ norm(u_exact); yscale=:log10, label="parallel",
    xlabel="sweep", ylabel="‖δx‖ / ‖u‖", lw=2, color=2)
plot(p1, p2; layout=(1, 2), size=(900, 350))

# ## More overlap, fewer sweeps
#
# The convergence rate of one-level Schwarz improves with the overlap width
# ``\delta`` (roughly like ``1 - C\delta/H`` for subdomain size ``H``).

for δ in 1:4
    dδ = SchwarzDecomposition(A, overlapping_subdomains(A, parts; overlap=δ); owner=parts)
    println("overlap = $δ: ", schwarz_solve(dδ, f; method=:parallel, tol=1e-8).iterations, " sweeps")
end

# ## The solution

heatmap(xs, xs, reshape(par.u, m, m)'; aspect_ratio=1, title="u (parallel Schwarz)")
