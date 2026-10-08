# Schwarz Methods

`SchwarzMethods` implements overlapping Schwarz domain decomposition for a
sparse SPD system ``A u = f`` as message passing on a cellular sheaf. The
subdomains ``\Omega_i`` of the decomposition are the vertices, the nonempty
overlaps ``\Omega_i \cap \Omega_j`` are the edges, and the restriction maps
read off the values a subdomain holds on an overlap ([`overlap_sheaf`](@ref)).

A Schwarz iterate is a 0-cochain of this sheaf: every subdomain keeps its own
copy of the solution. Each sweep solves a local problem on every subdomain,
using data from overlapping neighbours ([`schwarz_step!`](@ref)). The
coboundary of the iterate measures how much the copies disagree on overlaps
([`overlap_disagreement`](@ref)). At convergence the iterate is a global
section, which glues into the solution ([`glue`](@ref)).

The data is bundled into a few structs:

- [`SchwarzDecomposition`](@ref) holds the matrix, the [`OverlapCover`](@ref),
  the [`Ownership`](@ref) of each dof, one [`LocalProblem`](@ref) per subdomain,
  and a coloring of the subdomains.
- A [`TransmissionCondition`](@ref) sets what neighbours exchange:
  [`DirichletTransmission`](@ref) (classical Schwarz) or
  [`RobinTransmission`](@ref) (optimized Schwarz, with
  [`optimized_robin_parameter`](@ref) as a starting value).
- A [`SchwarzSweep`](@ref) is one pass of local solves:
  [`MultiplicativeSweep`](@ref) keeps the iterate a section,
  [`MulticolorSweep`](@ref) solves each color class of non-conflicting
  subdomains concurrently, and [`ParallelSweep`](@ref) (Lions/RAS) lets the
  copies disagree until they converge.
- A coarse level comes from a graph homomorphism ``\varphi`` that groups
  subdomains into aggregates, together with the pushforward of the overlap sheaf
  along it. [`TruncatedPushforwardCoarseSpace`](@ref) keeps a few modes of each
  pushforward stalk; it is a small Galerkin problem that makes the iteration
  scalable. [`ExactPushforwardCoarseSpace`](@ref) keeps the full stalks; it
  solves on the aggregates themselves, which is more accurate per sweep while
  there are few aggregates, but scales worse and is costlier.
- A [`SchwarzProblem`](@ref) bundles the decomposition with the right-hand side
  and initial guess. [`solve`](@ref) runs it with a stationary
  [`SchwarzIteration`](@ref) or with conjugate gradients preconditioned by the
  additive two-level operator ([`SchwarzCG`](@ref),
  [`SchwarzPreconditioner`](@ref)).

```julia
dd = SchwarzDecomposition(A, overlapping_subdomains(A, parts; overlap=2); owner=parts)
prob = SchwarzProblem(dd, f)
solve(prob, SchwarzIteration(sweep=MulticolorSweep(), coarse=TruncatedPushforwardCoarseSpace(dd)))
solve(prob, SchwarzCG(coarse=TruncatedPushforwardCoarseSpace(dd)))
```

See the [Schwarz domain decomposition](../generated/schwarz_domain_decomposition.md)
example for a worked Poisson problem.

```@autodocs
Modules = [CellularSheaves.NetworkSheaves.SchwarzMethods]
```
