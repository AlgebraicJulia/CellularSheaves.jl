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

# ## Multicolor sweeps
#
# Two subdomains *conflict* when they overlap or are coupled by ``A``.
# Subdomains of one color are independent, so a multiplicative sweep can solve
# a whole color class concurrently. The result is identical to a sequential
# sweep in color order. For the ``3 \times 3`` boxes four colors suffice, so a
# sweep takes four parallel steps instead of nine sequential ones.

dd.colors

#-

multicolor = schwarz_solve(dd, f; method=:multicolor, tol=1e-10)
multicolor.iterations

# ## Coarse spaces from the pushforward
#
# One-level methods only exchange information between neighbouring
# subdomains, so the iteration count grows with the number of subdomains. A
# coarse level fixes this. Both coarse spaces here start from a graph
# homomorphism ``\varphi : G \to H`` that groups subdomains into aggregates,
# and from the pushforward ``\varphi_* F`` of the overlap sheaf. Its stalk at an
# aggregate is the space of functions on the union of the subdomains in the
# fiber.
#
# - The **truncated pushforward** keeps one mode per aggregate: the constant
#   function, weighted by a partition of unity (with ``\varphi`` the identity
#   this is the Nicolaides coarse space). The coarse problem has one unknown per
#   aggregate.
# - The **exact pushforward** keeps the whole stalk. The coarse step is a
#   Schwarz sweep over the aggregates themselves.
#
# We group the boxes ``2 \times 2`` into aggregates.

function box_aggregation(p)
    q = cld(p, 2)
    return GraphHomomorphism([(cld(by, 2) - 1) * q + cld(bx, 2) for by in 1:p for bx in 1:p])
end

# The exact coarse space has the same stalks as `pushforward_sheaf`:

hom = box_aggregation(4)
small = SchwarzDecomposition(poisson2d(16), overlapping_subdomains(poisson2d(16), box_partition(16, 4); overlap=1);
    owner=box_partition(16, 4))
exact_small = ExactPushforwardCoarseSpace(small, hom)
vertex_stalks(pushforward_sheaf(hom, overlap_sheaf(small))) == length.(exact_small.decomposition.subdomains)

# ## Scalability
#
# We keep each box at ``8 \times 8`` grid points and add more boxes, so the grid
# grows with the number of subdomains. We compare sweeps to a relative residual
# of ``10^{-8}`` for the one-level method and for each coarse space, with their
# coarse dimensions (the number of unknowns the coarse level solves for).

function scaling_row(p)
    mp = 8p
    Ap = poisson2d(mp)
    pp = box_partition(mp, p)
    ddp = SchwarzDecomposition(Ap, overlapping_subdomains(Ap, pp; overlap=1); owner=pp)
    fp = ones(mp^2)
    coarse_spaces = (
        "one-level" => nothing,
        "truncated (φ = id)" => TruncatedPushforwardCoarseSpace(ddp),
        "truncated (2×2)" => TruncatedPushforwardCoarseSpace(ddp, box_aggregation(p)),
        "exact (2×2)" => ExactPushforwardCoarseSpace(ddp, box_aggregation(p)),
    )
    return map(coarse_spaces) do (name, c)
        its = schwarz_solve(ddp, fp; method=:multicolor, coarse=c, tol=1e-8, maxiter=5000).iterations
        dim = c === nothing ? 0 : coarse_dimension(c)
        (; p, name, its, dim)
    end
end

rows = reduce(vcat, [collect(scaling_row(p)) for p in (2, 4, 6, 8)])
println(rpad("boxes", 8), rpad("method", 22), rpad("sweeps", 8), "coarse dim")
for r in rows
    println(rpad("$(r.p)×$(r.p)", 8), rpad(r.name, 22), rpad(r.its, 8), r.dim)
end

# With the truncated coarse space and ``\varphi`` the identity, the sweep count
# levels off as boxes are added (23, 47, 56, 61 sweeps), while the one-level
# count grows roughly with the number of boxes. Its coarse problem has one
# unknown per subdomain. Coarser aggregates (``2 \times 2``) make the coarse
# problem four times smaller but less effective.
#
# The exact pushforward is the most accurate per sweep while there are few
# aggregates. With one aggregate it is a direct solve and converges in one
# sweep; on ``4 \times 4`` boxes it beats every truncated variant. But it
# carries no global information beyond its aggregates, so it is still a
# one-level method on larger subdomains. Its sweep count grows with the number
# of aggregates, overtaking the truncated space by ``6 \times 6`` boxes, and its
# coarse level costs as much as the whole problem plus overlaps. Higher accuracy
# per sweep, worse scalability.

ps = (2, 4, 6, 8)
plt = plot(; xlabel="boxes per side", ylabel="sweeps", yscale=:log10, legend=:topleft)
for name in unique(r.name for r in rows)
    plot!(plt, collect(ps), [r.its for r in rows if r.name == name]; label=name, marker=:circle, lw=2)
end
plt

# ## Robin transmission conditions
#
# Classical Schwarz passes *Dirichlet* data across each edge of the overlap
# graph. Optimized Schwarz methods pass *Robin* data
# ``(\partial_n + p)\,u`` instead, which damps low frequencies along the
# interface far better. Here this is done algebraically: the interface rows of
# each local matrix get a Neumann correction plus the Robin parameter ``p``. The
# exact solution stays the fixed point.
#
# Gander's optimized parameter for overlap width ``L`` is a continuous
# quantity. For our ``h^{-2}``-scaled matrix we divide it by ``h``. With overlap
# ``\delta`` grid layers on each side, the overlap width is ``L = (2\delta+1)h``.
#
# First on four vertical strips, where subdomain boundaries never meet:

h = 1 / (m + 1)
strips = [cld(jx * 4, m) for jy in 1:m for jx in 1:m]
strip_domains = overlapping_subdomains(A, strips; overlap=1)
pstar = optimized_robin_parameter(3h) / h
for scale in (nothing, 0.125, 0.5, 1, 2, 8)
    robin = scale === nothing ? nothing : scale * pstar
    dds = SchwarzDecomposition(A, strip_domains; owner=strips, robin)
    r = schwarz_solve(dds, f; method=:parallel, tol=1e-8, maxiter=3000)
    println(rpad(scale === nothing ? "Dirichlet" : "p = $(scale) p*", 14), r.iterations, " sweeps")
end

# The optimized parameter cuts the sweep count by an order of magnitude, and
# the formula's ``p^*`` is close to the best value.
#
# On boxes, four subdomains meet at each *cross point*. Discrete optimized
# Schwarz methods are known to be delicate there (Gander–Kwok 2013). On these
# ``3 \times 3`` boxes ``p^*`` converges in about half the Dirichlet sweep
# count, and ``2p^*`` does best. On a ``4 \times 4`` box grid with overlap 1,
# however, the stationary iteration diverged for ``p \le p^*`` and converged
# for ``p \ge 2p^*``. When in doubt, err towards larger ``p``, or use the CG
# solver below.

for scale in (1, 2, 4)
    ddr = SchwarzDecomposition(A, subdomains; owner=parts, robin=scale * optimized_robin_parameter(5h) / h)
    r = schwarz_solve(ddr, f; method=:multicolor, tol=1e-8, maxiter=3000)
    println("p = $(scale) p*: ", r.converged ? "$(r.iterations) sweeps" : "diverged")
end
println("Dirichlet: ", schwarz_solve(dd, f; method=:multicolor, tol=1e-8).iterations, " sweeps")

# ## Krylov acceleration
#
# The additive Schwarz operator ``\sum_i R_i^\mathsf{T} A_i^{-1} R_i`` is
# symmetric positive definite, so it can precondition conjugate gradients.
# Adding the truncated pushforward coarse space gives the classical two-level
# preconditioner. Its iteration count stays bounded as subdomains are added,
# while one-level CG grows.

for p in (4, 8, 16)
    mp = 8p
    Ap = poisson2d(mp)
    pp = box_partition(mp, p)
    ddp = SchwarzDecomposition(Ap, overlapping_subdomains(Ap, pp; overlap=1); owner=pp)
    fp = ones(mp^2)
    one_level = schwarz_cg(ddp, fp).iterations
    two_level = schwarz_cg(ddp, fp; coarse=TruncatedPushforwardCoarseSpace(ddp)).iterations
    println("$(p)×$(p) boxes: one-level CG $one_level, two-level CG $two_level")
end

# ## The solution

heatmap(xs, xs, reshape(par.u, m, m)'; aspect_ratio=1, title="u (parallel Schwarz)")
