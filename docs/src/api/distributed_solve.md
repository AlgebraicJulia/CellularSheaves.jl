# Distributed Solve

`DistributedSolve` performs the multifrontal triangular solve of a
`CliqueTrees.Multifrontal` chordal factorization with each supernode's
contribution addressed by node id, rather than through the single shared
traversal-order stack the serial `CliqueTrees` solve uses — the change that
lets the solve be split across independent processes. It provides a
single-process solve (a recursive reference and a preallocated
[`TreeWorkspace`](@ref) fast path for repeated re-solves), a tree partitioner
([`partition_tree`](@ref) / [`worker_factorization`](@ref)), and a genuinely
multi-process backend ([`distributed_tree_solve`](@ref)) that exchanges
boundary corrections over `RemoteChannel`s. See the
[Distributed Sheaf Solve](@ref) feature guide for the theory and benchmarks.

```@autodocs
Modules = [CellularSheaves.NetworkSheaves.DistributedSolve]
```
