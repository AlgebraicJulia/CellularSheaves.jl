# Core Sheaf Workflows

This page highlights the main externally facing APIs for constructing sheaves, assembling operators, and solving section problems.

## When to use this page

Use this guide when you want to:

- Build a sheaf model on a fixed graph.
- Form the sheaf Laplacian and evaluate consistency energy.
- Compute or approximate globally consistent sections.

## Constructing and Inspecting Sheaves

- Create sheaves using `EuclideanSheaf`, `sheaf_from_graph`, and `@cellular_sheaf`.
- Inspect structure with `vertex_stalks`, `edge_stalks`, and `underlying_graph`.
- See: [Euclidean Sheaves](../api/euclidean_sheaves.md), [Sheaf Interface](../api/sheaf_interface.md), and [Cellular Sheaf DSL](../api/dsl.md).

```julia
using CellularSheaves, Graphs

g = path_graph(4)
s = sheaf_from_graph(g, 2, n -> Matrix{Float64}(I, n, n))

vertex_stalks(s)
edge_stalks(s)
underlying_graph(s)
```

## Laplacians and Section Computations

- Build operators with `coboundary_map` and `sheaf_laplacian_matrix`.
- Optimize and project with `energy_function`, `nearest_global_section`, `nullspace_ldlt`, and `harmonic_extension`.
- See: [Euclidean Sheaves](../api/euclidean_sheaves.md) and [Sheaf Interface](../api/sheaf_interface.md).

```julia
L = sheaf_laplacian_matrix(s)
b = rand(sum(vertex_stalks(s)))

x_star = nearest_global_section(s, b)
E = energy_function(s)(x_star)
```

Typical workflow:

1. Construct `s` from graph structure and restriction maps.
2. Evaluate `sheaf_laplacian_matrix(s)` to inspect conditioning and size.
3. Project data with `nearest_global_section` or solve constrained values with `harmonic_extension`.
