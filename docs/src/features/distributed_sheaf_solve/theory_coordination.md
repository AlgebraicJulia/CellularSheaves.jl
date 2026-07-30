# Coordination as a harmonic extension

This page derives the one linear system that the whole feature exists to solve.
Everything later (the factorization, the partition, the message passing) is
machinery for this system, so we take the time to obtain it honestly, from the
energy it minimizes.

## The sheaf

Let ``\mathcal F`` be a cellular sheaf over a graph ``G = (V, E)``.

Each vertex ``v`` carries a stalk ``\mathcal F(v) = \mathbb R^{d_v}``, an agent's
local coordination variables. Each edge ``e`` carries a stalk
``\mathcal F(e) = \mathbb R^{d_e}``, and each incidence ``v \trianglelefteq e``
a restriction map

```math
\mathcal F_{v \trianglelefteq e} \colon \mathcal F(v) \longrightarrow \mathcal F(e).
```

*For an agent:* the edge stalk is the "language" the two endpoints compare
themselves in, and the restriction map is how each endpoint translates its own
state into that language.

## Disagreement, energy, Laplacian

The **coboundary** ``\delta`` measures translation mismatch on every edge. For a
0-cochain ``x = (x_v)_{v \in V}`` and an oriented edge ``e = (u, v)``,

```math
(\delta x)_e \;=\; \mathcal F_{v \trianglelefteq e}\, x_v \;-\; \mathcal F_{u \trianglelefteq e}\, x_u .
```

Half its squared norm is the **Dirichlet energy** of the configuration,

```math
\mathscr E(x) \;=\; \tfrac12 \lVert \delta x \rVert^2
\;=\; \tfrac12 \sum_{e = (u,v)} \bigl\lVert \mathcal F_{v \trianglelefteq e} x_v - \mathcal F_{u \trianglelefteq e} x_u \bigr\rVert^2 ,
```

and expanding the square exhibits the **sheaf Laplacian**:

```math
\mathscr E(x) \;=\; \tfrac12\, x^\top L_{\mathcal F}\, x,
\qquad
L_{\mathcal F} \;=\; \delta^\top \delta \;\succeq\; 0 .
```

*For an agent:* ``\mathscr E`` totals how badly every edge's two endpoints
disagree. Configurations with ``\mathscr E(x) = 0``, the kernel of
``L_{\mathcal F}`` and the **global sections**, are the perfectly coordinated
ones. For identity restriction maps on a connected graph that kernel is exactly
consensus. Richer maps encode offsets, rotations, or gauge relations.

## Pinning the targets

Coordination is not free consensus: some cells are driven. Split the vertices
into pinned cells ``P`` (targets, holding prescribed boundary values ``p``) and
free cells ``F`` (agents, to be solved for), and split ``L_{\mathcal F}``
conformally:

```math
L_{\mathcal F} \;=\;
\begin{bmatrix} H & B \\ B^\top & L_{PP} \end{bmatrix},
\qquad H = (L_{\mathcal F})_{FF}, \quad B = (L_{\mathcal F})_{FP}.
```

The coordination reference is the free configuration of least energy given the
pins. Substituting ``x = (q, p)`` and expanding,

```math
\mathscr E(q, p) \;=\; \tfrac12\, q^\top H q \;+\; q^\top B p \;+\; \tfrac12\, p^\top L_{PP}\, p .
```

Only the first two terms involve ``q``. Setting ``\nabla_q \mathscr E = 0``:

```math
\nabla_q \mathscr E \;=\; H q + B p \;=\; 0
\qquad \Longrightarrow \qquad
\boxed{\; H\, q^\star \;=\; -\,B\, p \;}
```

This ``q^\star`` is the **harmonic extension** of the boundary data. Equivalently
``(L_{\mathcal F}\, x)\big|_F = 0``: every free cell sits at the sheaf-weighted
average of its neighbours.

*For an agent:* ``q^\star_i`` is where the whole formation "wants" agent ``i`` to
be, consistent with its neighbours and, through them, with the targets. It is the
reference the control stack tracks.

## When ``H`` is invertible

``H`` inherits positive *semi*-definiteness from ``L_{\mathcal F}``; pinning is
what makes it definite.

!!! note "Dirichlet condition"
    ``H \succ 0`` if and only if every connected component of the free subgraph
    contains a vertex adjacent to a pinned cell whose restriction maps have full
    rank along that edge. In the identity-map case: every group of agents can
    "see" at least one target through the graph.

The proof is one line in each direction: ``q^\top H q = \lVert \delta(q, 0)
\rVert^2`` vanishes iff ``(q, 0)`` is a section, and a section that is zero on
the pins and connected to them must be zero.

*For an agent:* a component of agents with no path to any target has nothing
anchoring it, so its consensus value is genuinely undetermined, and ``H`` is
telling you so.

## Three consequences that shape everything downstream

**1. The solve is exact, not iterative.** ``H \succ 0`` has a Cholesky
factorization. The reference is obtained by triangular substitution, not by
running a diffusion to tolerance.

**2. Moving targets never touch ``H``.** The targets enter only through the
right-hand side ``-Bp``. Per control step the update is

```math
q^\star(t) \;=\; H^{-1} \bigl(-B\, p(t)\bigr)
\qquad \text{(same factor, new right-hand side),}
```

one forward and one backward substitution against a factorization computed
once. Re-factorizing per step would be pure waste, and the benchmarks never
charge it.

**3. Changing the *sheaf* itself is a separate regime.** Everything above holds
a *fixed* sheaf: the graph, the restriction maps, and hence ``H`` and its
factorization do not change while the targets move. When the sheaf *does* change
(an agent joins, a coupling edge appears or drops), ``H`` changes by
``\pm\, \delta_e^\top \delta_e`` for the affected edge ``e``, a low-rank
perturbation that a *dynamic* sheaf treatment can absorb by patching the factor
rather than rebuilding it. That is its own feature, out of scope here. This guide
solves the fixed system.

The remaining question is *how* to perform that factor-and-substitute when no
single machine is allowed to hold ``H``, its factor, or the whole of ``q^\star``.
That is the next two pages.
