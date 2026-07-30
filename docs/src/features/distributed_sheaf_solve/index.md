# Distributed Sheaf Solve

Coordinating a fleet on a cellular sheaf reduces, at every control step, to one
linear system, the **harmonic extension**

```math
H\, q^\star = -B\, p,
```

where ``H`` is the agent block of the sheaf Laplacian and ``p`` the current
target positions. This feature solves that system **exactly**, and in a
**distributed** way. The factor of ``H`` is computed once, cut into per-worker
slices along its elimination tree, and each solve is a short message-passing
protocol, so no machine ever holds the whole problem.

![The 12-agent two-ring formation and its elimination tree of nine supernodes](../../assets/figures/distributed_solve/fig_formation_tree.svg)

## The claims, and where each is established

| claim | where |
|---|---|
| Coordination *is* one SPD solve, and moving targets change only the right-hand side | [Coordination as a harmonic extension](theory_coordination.md) |
| The multifrontal factorization organizes that solve per supernode, with disjoint ownership | [The multifrontal factorization](theory_multifrontal.md) |
| The tree can be cut across machines, and the cross-cut messages are small and provably sufficient | [Distributing the solve](theory_distributed.md) |
| Under one fair radio model, it beats a centralized coordinator and accelerated diffusion by 10–36×, with measured caveats | [Three methods, one communication model](comparison.md) |
| Exact to ``10^{-16}``, ``O(\log n)`` slots under nested dissection, 12% worst-case memory per agent | [Benchmarks](benchmarks.md) |
| How to call all of it | [Usage and API](api.md) |

The worked example,
[Distributed Harmonic Extension: Tracking Without a Central Solver](@ref),
runs the full pipeline on twelve simulated agents across three real processes
and matches the central answer to ``8.9\times10^{-16}``.

## What it provides

- [`distributed_tree_solve`](@ref), the multi-process solve: one chunk of the
  factor per worker, boundary corrections over `RemoteChannel`s.
- [`partition_tree`](@ref) / [`worker_factorization`](@ref), the minimally-cut
  connected partition and each worker's self-contained slice.
- [`TreeWorkspace`](@ref) with three-argument `tree_forward_ldiv!` /
  `tree_backward_ldiv!`, the allocation-free cached re-solve for the per-step
  tracking loop.

## Why this is not "just call a sparse solver"

A monolithic sparse Cholesky stores its factor on one node and threads its solve
through one shared stack, and neither survives being split across machines. The
work here is doing the *identical* arithmetic with the tree cut into pieces:
same answer to machine precision, but the memory, the factor, and the
communication are all decentralized. The communication scales with the tree's
**depth**, which the elimination ordering keeps logarithmic, rather than with
the size of the fleet.
