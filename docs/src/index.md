# CellularSheaves.jl

```@meta
CurrentModule = CellularSheaves
```

`CellularSheaves.jl` is a Julia library for constructing, analyzing, and computing with cellular sheaves on graphs. It provides efficient sparse solvers for sheaf Laplacians, global sections, harmonic extensions, and trajectory planning via sheaf-theoretic colocation methods.

The package is part of the [AlgebraicJulia](https://algebraicjulia.org) ecosystem and implements computational methods developed by Tyler Hanks, Hans Riess, Matthew Hale, Warren Dixon, and James Fairbanks, building on the mathematical foundations of Robert Ghrist, Justin Curry, and Jakob Hansen.

## Mathematical Background

### Cellular Sheaves

A *cellular sheaf* ``\mathcal{F}`` on a graph ``G = (V, E)`` assigns:

- a vector space ``\mathcal{F}(v)`` (the *stalk* at ``v``) to each vertex ``v \in V``,
- a vector space ``\mathcal{F}(e)`` (the *stalk* at ``e``) to each edge ``e \in E``, and
- a linear *restriction map* ``\rho_{v \trianglelefteq e} : \mathcal{F}(v) \to \mathcal{F}(e)`` for each vertex-edge incidence pair ``v \trianglelefteq e``.

For an oriented edge ``e = (u, v)``, the *coboundary map* ``d : C^0(\mathcal{F}) \to C^1(\mathcal{F})`` is defined entry-wise by

```math
(d\, x)_e = \rho_{u \trianglelefteq e}\, x_u \;-\; \rho_{v \trianglelefteq e}\, x_v.
```

The 0-cochain space ``C^0(\mathcal{F}) = \bigoplus_{v \in V} \mathcal{F}(v)`` collects all vertex-stalk assignments, and the 1-cochain space ``C^1(\mathcal{F}) = \bigoplus_{e \in E} \mathcal{F}(e)`` collects all edge-stalk assignments.

### Sheaf Laplacian

The *sheaf Laplacian* is

```math
L_{\mathcal{F}} = d^\top d : C^0(\mathcal{F}) \to C^0(\mathcal{F}).
```

It is symmetric positive semidefinite and generalizes the combinatorial graph Laplacian to vector-valued data with structured edge constraints. The *sheaf quadratic form* (Dirichlet energy) is

```math
h_{\mathcal{F}}(x) = x^\top L_{\mathcal{F}} x = \|d\, x\|^2 = \sum_{e=(u,v)\in E} \|\rho_{u \trianglelefteq e}\, x_u - \rho_{v \trianglelefteq e}\, x_v\|^2,
```

which measures the total disagreement in ``\mathcal{F}``-values across all edges.

### Global Sections

A *global section* of ``\mathcal{F}`` is a 0-cochain ``x \in C^0(\mathcal{F})`` satisfying the consistency condition

```math
\rho_{u \trianglelefteq e}\, x_u = \rho_{v \trianglelefteq e}\, x_v \quad \text{for every edge } e = (u,v) \in E.
```

Equivalently, ``x \in \ker(d) = \ker(L_{\mathcal{F}})``. The space of global sections ``H^0(\mathcal{F}) = \ker(L_{\mathcal{F}})`` captures all globally consistent assignments. Its dimension is the *zeroth sheaf cohomology* ``h^0(\mathcal{F})``.

Use [`nullspace_ldlt`](@ref) to compute a basis for ``H^0(\mathcal{F})`` via sparse LDLt factorization, and [`nearest_global_section`](@ref) to project any 0-cochain to the closest global section in the L² sense.

### Harmonic Extension

Given a partition of vertices into *boundary* ``B \subseteq V`` and *interior* ``I = V \setminus B``, and prescribed values ``x_B`` on the boundary stalks, the *harmonic extension* solves the Dirichlet problem

```math
L_{II}\, x_I = -L_{IB}\, x_B,
```

where ``L_{II}`` and ``L_{IB}`` are the corresponding blocks of the sheaf Laplacian. The solution ``x_I`` minimizes the sheaf energy over all 0-cochains agreeing with ``x_B`` on ``B``.

If ``L_{II}`` is singular (i.e. the boundary conditions do not uniquely determine the interior values), [`harmonic_extension`](@ref) returns both a minimum-norm particular solution and a basis for the null space of ``L_{II}``.

### Dynamics and Control

A *trajectory sheaf* encodes ``k`` steps of the linear discrete-time dynamics ``x_{t+1} = A x_t`` as a sheaf on a path graph of length ``k+1``. Vertex stalks are copies of the state space ``\mathbb{R}^d``; edge stalks also have dimension ``d``, with restriction maps

```math
\rho_{t \trianglelefteq e} = A, \quad \rho_{(t+1) \trianglelefteq e} = I_d.
```

The coboundary condition ``d\, x = 0`` is then exactly ``A x_t - x_{t+1} = 0``, so **global sections are valid trajectories**. Solving a two-point boundary-value problem (BVP) ``x_0 = \bar{x}_0,\; x_k = \bar{x}_k`` reduces to harmonic extension with boundary vertices ``\{0, k\}``.

For controlled systems ``\dot{x}(t) = A_c x(t) + B_c u(t)``, the [`ControlledTrajectorySheaf`](@ref) constructs the analogous sheaf using the zero-order-hold discretization ``x_{t+1} = A_d x_t + B_d u_t``.

See the [Spring Oscillator: Four Colocation Methods Compared](@ref) example for a worked comparison between the sheaf colocation method and the classical block-bidiagonal formulation.

## Getting Started

Install the package from the Julia general registry:

```julia
using Pkg
Pkg.add("CellularSheaves")
```

Then explore the Examples section, which includes:

- **[Code Example](@ref)**: constructing a sheaf with the `@cellular_sheaf` DSL and computing the nearest global section via the direct LDLt method.
- **[Computing the Nullspace of the Sheaf Laplacian](@ref)**: visualizing the nullspace of the sheaf Laplacian for fully and partially constrained sheaves.
- **[Sheaf Morphisms](generated/sheaf_morphisms.md)**: building morphisms between sheaves and verifying the naturality condition.
- **[Pushforward of a Sheaf along a Graph Homomorphism](@ref)**: computing the pushforward sheaf along a graph homomorphism and verifying the global-section isomorphism.
- **[Spring Oscillator: Four Colocation Methods Compared](@ref)**: solving a spring-oscillator BVP with four methods (classical QR vs. sheaf harmonic extension, Euler vs. exact dynamics).
- **[Nearest Global Section: Iterative Method](@ref)**: using the conjugate-gradient backend for `nearest_global_section`, with graceful error handling for convergence failures.

## References

The following papers provide the mathematical and computational foundations of this library:

- **Ghrist, R. (2014).** *Elementary Applied Topology.* Createspace. ISBN 978-1502880857.
- **Curry, J. (2014).** *Sheaves, Cosheaves and Applications.* PhD Thesis, University of Pennsylvania. [arXiv:1303.3405](https://arxiv.org/abs/1303.3405)
- **Hansen, J. & Ghrist, R. (2019).** *Toward a Spectral Theory of Cellular Sheaves.* Journal of Applied and Computational Topology, 3(3), 315–358. [DOI:10.1007/s41468-019-00038-7](https://doi.org/10.1007/s41468-019-00038-7)
- **Hansen, J. & Ghrist, R. (2021).** *Opinion Dynamics on Discourse Sheaves.* SIAM Journal on Applied Mathematics, 81(5), 2033–2060. [DOI:10.1137/20M1341088](https://doi.org/10.1137/20M1341088)
- **Riess, H. & Hansen, J. (2022).** *Multidimensional Persistence Module Classification via Lattice-Theoretic Measures.* Foundations of Computational Mathematics. [arXiv:2104.09036](https://arxiv.org/abs/2104.09036)
