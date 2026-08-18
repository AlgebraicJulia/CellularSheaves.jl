# Usage and API

This page shows the three ways the feature is used (a cached single-process
re-solve for the per-step tracking loop, a genuinely multi-process solve, and
choosing the elimination ordering), and flags the one permutation subtlety.
Full docstrings are in the [Distributed Solve](@ref) API reference.

## Factor once, re-solve every step

Build the Dirichlet block ``H`` once, factor it once, and reuse a
[`TreeWorkspace`](@ref) for every subsequent solve as the targets move:

```julia
using CellularSheaves, CliqueTrees.Multifrontal, LinearAlgebra, SparseArrays

F  = cholesky!(ChordalCholesky(sparse(H)), NoPivot())  # once per formation
L  = F.L
ws = TreeWorkspace(L, Float64)                         # once per factorization

# every control step, with a new right-hand side rhs = -B*p :
c = Vector(F.P' \ rhs)                                 # into the factor's space
tree_forward_ldiv!(L, c, ws)
tree_backward_ldiv!(L, c, ws)
qstar = F.P \ c                                        # back to problem space
```

After construction the two `tree_*_ldiv!` calls allocate almost nothing, so the
per-step cost is just the two triangular sweeps.

!!! note "The permutation is not optional"
    `ChordalCholesky` factors ``H = P^\top L\,L^\top P`` with a fill-reducing
    permutation `F.P`, so `L` lives in *permuted* coordinates. A right-hand side
    must be mapped in with `F.P' \ rhs` and the result mapped back out with
    `F.P \ y`. Skip this and the solve is exact, but for the wrong ordering of
    the unknowns. The `tree_*_ldiv!` and [`distributed_tree_solve`](@ref) routines
    all operate in `L`'s permuted space and leave this wrapping to the caller.

## A genuinely multi-process solve

[`partition_tree`](@ref) cuts the elimination tree into connected chunks and
[`distributed_tree_solve`](@ref) runs one chunk per worker process, exchanging
boundary corrections over `RemoteChannel`s. Provision the workers yourself. The
package does not call `addprocs`, matching how cluster provisioning is the
caller's decision:

```julia
using Distributed

partition = partition_tree(L, 6)
nchunk    = length(partition.chunks)          # may be < 6; see partition_tree

pids = addprocs(nchunk; exeflags = "--project=$(Base.active_project())")
@everywhere pids using CellularSheaves

rhs_p = Vector(F.P' \ rhs)
y     = distributed_tree_solve(L, rhs_p, 6; pids = pids)
qstar = F.P \ y
```

The result matches the single-process solve (and the monolithic `CliqueTrees`
solve) to machine precision. Note that the requested worker count is a target,
not a guarantee: a wide, shallow tree yields fewer, larger chunks (see
[`partition_tree`](@ref)), so size `pids` to `length(partition.chunks)` after
partitioning.

## Choosing the elimination ordering

The ordering controls the tree's depth and hence the communication makespan
(the [benchmarks](benchmarks.md) measure a 36× difference on a 510-agent
chain). `ChordalCholesky` accepts any `CliqueTrees` elimination algorithm, or
a plain permutation vector, through `alg`:

```julia
# default fill-reducing ordering
F = cholesky!(ChordalCholesky(sparse(H)), NoPivot())

# a hand-built nested-dissection order: eliminate interval middles last
function nd_order_path(n)
    out = Int[]
    rec(lo, hi) = begin
        lo > hi && return
        mid = (lo + hi) ÷ 2
        rec(lo, mid - 1); rec(mid + 1, hi)
        push!(out, mid)
    end
    rec(1, n)
    return out
end

F = cholesky!(ChordalCholesky(sparse(H); alg = dofs_of(nd_order_path(nagents))), NoPivot())
```

where `dofs_of` expands a vertex order to consecutive per-vertex dof indices.
For chain-, ring-, and grid-shaped formations this three-line construction
already achieves logarithmic tree depth. `CliqueTrees` also ships `ND()` for
general graphs (requires a dissection backend such as Metis).

## When to reach for which

- **Single process, one solve:** the two-argument `tree_forward_ldiv!` /
  `tree_backward_ldiv!`, or just `CliqueTrees`' own `ldiv!`.
- **Single process, many solves against one factor** (the tracking loop): build
  a [`TreeWorkspace`](@ref) and use the three-argument methods.
- **Across machines** (each agent holds only its slice): [`partition_tree`](@ref)
  + [`distributed_tree_solve`](@ref).

## API reference

See [Distributed Solve](@ref) for the complete docstrings of every type and
function referenced here.
