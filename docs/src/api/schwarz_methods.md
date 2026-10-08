# Schwarz Methods

`SchwarzMethods` implements overlapping Schwarz domain decomposition for a
sparse SPD system ``A u = f`` as message passing on a cellular sheaf. The
subdomains ``\Omega_i`` of the decomposition are the vertices, the nonempty
overlaps ``\Omega_i \cap \Omega_j`` are the edges, and the restriction maps
read off the values a subdomain holds on an overlap ([`overlap_sheaf`](@ref)).

A Schwarz iterate is a 0-cochain of this sheaf: every subdomain keeps its own
copy of the solution. Each sweep solves a Dirichlet problem on every subdomain
using boundary data from overlapping neighbours ([`schwarz_step!`](@ref)).
The coboundary of the iterate measures how much the copies disagree on overlaps
([`overlap_disagreement`](@ref)). At convergence the iterate is a global
section, which glues into the solution ([`glue`](@ref)). Two variants are
provided. The *multiplicative* (alternating) method keeps the iterate a
section at every step. The *parallel* (Lions/RAS) method lets the copies
disagree until they converge. The *multicolor* method is the multiplicative
method with each color class of non-conflicting subdomains solved
concurrently.

A coarse level comes from a graph homomorphism ``\varphi`` that groups
subdomains into aggregates, together with the pushforward of the overlap sheaf
along it. [`TruncatedPushforwardCoarseSpace`](@ref) keeps a few modes of each
pushforward stalk; it is a small Galerkin problem that makes the iteration
scalable. [`ExactPushforwardCoarseSpace`](@ref) keeps the full stalks; it
solves on the aggregates themselves, which is more accurate per sweep but
costlier.

See the [Schwarz domain decomposition](../generated/schwarz_domain_decomposition.md)
example for a worked Poisson problem.

```@autodocs
Modules = [CellularSheaves.NetworkSheaves.SchwarzMethods]
```
