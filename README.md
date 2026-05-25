# CellularSheaves.jl

[![Stable Documentation](https://img.shields.io/badge/docs-stable-blue.svg)](https://AlgebraicJulia.github.io/CellularSheaves.jl/stable)
[![Development Documentation](https://img.shields.io/badge/docs-dev-blue.svg)](https://AlgebraicJulia.github.io/CellularSheaves.jl/dev)
[![Code Coverage](https://codecov.io/gh/AlgebraicJulia/CellularSheaves.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/AlgebraicJulia/CellularSheaves.jl)
[![CI/CD](https://github.com/AlgebraicJulia/CellularSheaves.jl/actions/workflows/julia_ci.yml/badge.svg)](https://github.com/AlgebraicJulia/CellularSheaves.jl/actions/workflows/julia_ci.yml)

**CellularSheaves.jl** is a Julia library for constructing, analyzing, and computing with cellular sheaves on graphs. It provides efficient sparse solvers for sheaf Laplacians, computation of global sections, harmonic extension, and sheaf-theoretic formulations of dynamical systems and control problems.

## What is a Cellular Sheaf?

A *cellular sheaf* on a graph G = (V, E) assigns a vector space — a *stalk* — to each vertex and each edge, together with a *restriction map* (a linear map) from each vertex stalk into the adjacent edge stalks. This data encodes consistency constraints across the graph: a *global section* is an assignment of a vector to each vertex stalk such that each pair of adjacent vertex values agrees when projected to the shared edge stalk.

Cellular sheaves generalize scalar-valued signal processing on graphs to vector-valued, heterogeneous data with structured consistency requirements. The *sheaf Laplacian* L = dᵀd (where d is the coboundary map) plays the role of a graph Laplacian adapted to the sheaf structure, and its nullspace is precisely the space of global sections.

For an accessible introduction to the theory, see:

- **Hansen & Ghrist (2019).** *Toward a Spectral Theory of Cellular Sheaves.* Journal of Applied and Computational Topology, 3(3), 315–358.
- **Curry (2014).** *Sheaves, Cosheaves and Applications.* PhD Thesis, University of Pennsylvania. [arXiv:1303.3405](https://arxiv.org/abs/1303.3405)

## Features

- **Euclidean sheaves** with sparse coboundary maps and sheaf Laplacians.
- **Global sections**: sparse LDLt-based nullspace computation (`nullspace_ldlt`) via [CliqueTrees.jl](https://github.com/AlgebraicJulia/CliqueTrees.jl).
- **Nearest global section**: project a 0-cochain onto the space of global sections.
- **Harmonic extension**: solve the Dirichlet problem L_II x_I = −L_IB x_B for sheaf-valued data.
- **Sheaf morphisms** and their composition.
- **Pushforward sheaves** along graph homomorphisms, with isomorphic global-section spaces.
- **Trajectory sheaves** for linear dynamical systems: encode multi-step dynamics x_{t+1} = A xₜ as a sheaf on a path graph and solve two-point boundary-value problems via harmonic extension.
- **Controlled trajectory sheaves** for zero-order-hold discretized LTI systems ẋ = Ac x + Bc u.
- **DSL** for specifying sheaves declaratively with the `@cellular_sheaf` macro.

## Installation

```julia
using Pkg
Pkg.add("CellularSheaves")
```

## Quick Example

```julia
using CellularSheaves
using LinearAlgebra

# Build a sheaf on a 4-cycle with 2-dimensional stalks
# and identity restriction maps
F = EuclideanSheaf{Float64}(Int[])
for i in 1:4
    add_vertex_stalk!(F, 2)
end
for i in 1:4
    add_sheaf_edge!(F, i, mod1(i+1, 4), I(2), I(2))
end

# Compute the space of global sections (nullspace of the sheaf Laplacian)
V = nullspace_ldlt(F)
println("Global section space dimension: ", size(V, 2))  # → 2

# Project a random 0-cochain onto the nearest global section
x0 = randn(sum(vertex_stalks(F)))
gs = nearest_global_section(F, x0)
println("Coboundary residual: ", norm(coboundary_map(F) * gs))  # → ≈ 0
```

## Applications

This package implements theory and algorithms from research on network sheaves, including:

- **Spectral analysis** of data on networks via the sheaf Laplacian (Hansen & Ghrist 2019).
- **Opinion dynamics** and multi-agent consensus modeled as sheaf diffusion (Hansen & Ghrist 2021).
- **Trajectory planning and colocation** for linear dynamical systems via harmonic extension of trajectory sheaves (Hanks, Riess, Hale, Dixon & Fairbanks).
- **Formation control and distributed estimation** using sheaf-theoretic encodings of multi-agent constraints.

See the [documentation](https://AlgebraicJulia.github.io/CellularSheaves.jl/dev) for worked examples, mathematical background, and the full API reference.

## Authors and Attribution

CellularSheaves.jl was developed by **Tyler Hanks** and **James Fairbanks** as part of the [AlgebraicJulia](https://algebraicjulia.org) ecosystem.

The package is an implementation of computational methods developed by **Tyler Hanks**, **Hans Riess**, **Matthew Hale**, **Warren Dixon**, and **James Fairbanks**, building on foundational mathematical research by **Robert Ghrist**, **Justin Curry**, and **Jakob Hansen**.

## References

- **Ghrist, R. (2014).** *Elementary Applied Topology.* Createspace. ISBN 978-1502880857.
- **Curry, J. (2014).** *Sheaves, Cosheaves and Applications.* PhD Thesis. [arXiv:1303.3405](https://arxiv.org/abs/1303.3405)
- **Hansen, J. & Ghrist, R. (2019).** *Toward a Spectral Theory of Cellular Sheaves.* Journal of Applied and Computational Topology, 3(3), 315–358. [DOI](https://doi.org/10.1007/s41468-019-00038-7)
- **Hansen, J. & Ghrist, R. (2021).** *Opinion Dynamics on Discourse Sheaves.* SIAM Journal on Applied Mathematics, 81(5), 2033–2060. [DOI](https://doi.org/10.1137/20M1341088)
- **Anwer, N., Riess, H., & Hale, M. (2026).** Multi-Agent System Identification with Nonlinear Sheaf Diffusion (arXiv:2605.11204). arXiv. https://doi.org/10.48550/arXiv.2605.11204
- **Ghrist, R., Lopez, M., North, P. R., & Riess, H. (2026).** Categorical Diffusion of Weighted Lattices (arXiv:2501.03890). arXiv. [DOI](https://doi.org/10.48550/arXiv.2501.03890)
- **Riess, H., & Ghrist, R. (2022). Diffusion of Information on Networked Lattices by Gossip. 2022** IEEE 61st Conference on Decision and Control (CDC), 5946–5952. [DOI](https://doi.org/10.1109/CDC51059.2022.9992539)
- **Riess, H., Munger, M., & Zavlanos, M. M. (2023).** Max-Plus Synchronization in Decentralized Trading Systems. 2023 62nd IEEE Conference on Decision and Control (CDC), 221–227. [DOI](https://doi.org/10.1109/CDC49753.2023.10383918)
- **Ghrist, R., & Riess, H. (2022).** Cellular sheaves of lattices and the Tarski Laplacian. Homology, Homotopy and Applications, 24(1), 325–345. [DOI](https://doi.org/10.4310/HHA.2022.v24.n1.a16)
