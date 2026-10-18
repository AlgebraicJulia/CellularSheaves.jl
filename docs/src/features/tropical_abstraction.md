# Tropical Abstraction: Coarse Grids for Optimal Control

This page explains the min-plus ("tropical") side of the package
([`TropicalAbstraction`](@ref CellularSheaves.NetworkSheaves.TropicalAbstraction)):

- coarse grids for deterministic optimal control;
- a layered controller built on them;
- the certificates that come with that controller.

It covers three things: the idea, how it reads as sheaf theory, and how it relates
to the literature it borrows from. The studies that measured these constructions,
including where they do not pay off, accompany the paper on path lifting.

In one sentence: the multigrid ingredients we use for the HJB equation, the
pushforward along an aggregation and its Galerkin coarse operator, still make sense
when the semiring is ``(\min, +)``. There is no subtraction, so the coarse solution
cannot *correct* the fine one. Instead it *bounds* it. Lax and oplax morphisms of
min-plus sheaves say which way each bound goes.

## Optimal control as min-plus linear algebra

On a finite state graph whose edge ``x \to y`` costs ``W_{xy} \ge 0``, the cost-to-go
to a target set ``\Gamma`` solves the Bellman equation

```math
V(x) = \min\Big(b(x),\ \min_y W_{xy} + V(y)\Big), \qquad b = 0 \text{ on } \Gamma,\ \infty \text{ elsewhere}.
```

In the semiring ``\mathbb R_{\min} = (\mathbb R \cup \{\infty\}, \oplus = \min, \otimes = +)``
this is the *linear* fixed-point equation ``V = W \otimes V \oplus b``, solved by the
closure ``V = W^* \otimes b`` with ``W^* = I \oplus W \oplus W^{\otimes 2} \oplus \cdots``
([`minplus_closure`](@ref), [`cost_to_go`](@ref)). Dijkstra's algorithm computes one
column of ``W^*``. Other control problems are closures in other semirings:

| semiring | ``\oplus``, ``\otimes`` | closure computes | here |
|---|---|---|---|
| tropical | ``\min``, ``+`` | minimum cost / minimum time to a target | [`cost_to_go`](@ref), [`Abstraction`](@ref) |
| bottleneck | ``\max``, ``\min`` | reach-avoid value: best worst-case margin along a path | [`reach_avoid`](@ref), [`maxmin_closure`](@ref) |
| real | ``+``, ``\times`` | discounted policy evaluation ``V = P^* c``, ``P = D^{-1}N`` | `PathPackCoarsest` in [`Multigrid`](@ref CellularSheaves.ControlSheaves.DoubleIntegratorHJB.Multigrid) |

The first two are *idempotent* (``a \oplus a = a``), and that is what this page
exploits. The third is the ordinary linear algebra of the HJB solvers
([`DoubleIntegratorHJB`](../api/double_integrator_hjb.md)). It shows up here because
the same factorization code (PathPack's chordal semiring LU) solves all three.

The running example in the studies is a minimum-time double integrator on a
lattice:
- **States:** ``128^2`` positions times ``13^2`` velocities, 2.77M states, of which
  2.41M can reach the goal.
- **Moves:** one step changes each velocity component by at most one, at unit cost.
- **Obstacles:** two walls with gaps force the routes to snake.

## The sheaf picture

### The state graph as a min-plus sheaf

Let ``\psi : X \to H`` send each state to a *cell* of a coarse graph ``H``. In the
example ``H`` is the position grid and ``\psi(q, v) = q`` forgets the velocity
([`projection_homomorphism`](@ref)). Pushing the state graph forward along ``\psi``
gives a min-plus sheaf ``F`` on ``H``:

- the stalk over a cell ``q`` is ``\mathbb R_{\min}^{\psi^{-1}(q)}``, a value for every
  state in the fibre (for ``(q, v) \mapsto q``, a function of the velocity);
- the map on an edge ``q \to q'`` is the block ``W_{qq'}`` of ``W`` from the fibre over
  ``q'`` to the fibre over ``q``, and moves inside a cell are loops.

Two caveats on the word *sheaf*. The edges are directed (controls are not
reversible), so strictly ``F`` is a representation of the directed cell graph in
min-plus semimodules. And the maps are not restriction maps into an edge stalk but
transfer maps between vertex stalks. We keep the package's usage because
everything below is the min-plus counterpart of a construction on the real sheaves
elsewhere in the package, and the analogy guides the constructions.

### Sections become sub- and supersolutions

On a real sheaf a global section satisfies an *equation* on every edge. In
``\mathbb R_{\min}`` the natural notion is an *inequality*, and there are two of them:

- a **lax section** (subsolution): ``V \le W \otimes V \oplus b``;
- a **colax section** (supersolution): ``W \otimes V \oplus b \le V``. This is exactly a
  *control Lyapunov function*: from every state off the targets some move lowers ``V``
  by at least its cost ([`is_control_lyapunov`](@ref)).

The optimal value is the greatest subsolution and the least supersolution, and a
fixed point (the only one when there are no zero-cost cycles). Every subsolution lies below it and every supersolution above
it:

```math
\underbrace{V_{\text{sub}}}_{\text{lower bound, A* heuristic}} \;\le\; V^\star \;\le\; \underbrace{V_{\text{super}}}_{\text{upper bound, certified policy}} .
```

The vectors ``\overline{\mathbb R}^X`` form a complete lattice and the Bellman
operator is monotone, so this is Tarski's fixed-point theorem. That puts the
construction in the setting of Ghrist and Riess's sheaves of lattices, whose global
sections are fixed points of the Tarski Laplacian.

### Pushforward and pullback: a Galois connection

Aggregate the states along a map ``\mu`` (into blocks, or into modes of each
stalk, below). The *refinement* ``E = \mu^*`` copies a block's value to its states
([`tropical_pullback`](@ref), [`refine`](@ref)). Its residual, the *coarsening*
``P``, takes the minimum over each block ([`tropical_pushforward`](@ref),
[`coarsen`](@ref)). They form a Galois connection:

```math
E g \le f \iff g \le P f, \qquad P E = I, \qquad E P \le I .
```

The real block average needs a division; here ``P E = I`` holds with no
normalization because ``\min`` is idempotent. ``EP`` is the max-plus projection of
Cohen, Gaubert and Quadrat, the one used by max-plus finite elements. The Galerkin
product ``W_c = P \otimes W \otimes E`` is the **quotient graph**
([`quotient_graph`](@ref)): one vertex per block, the cheapest edge between two
blocks, and edges inside a block contracted.

### Lax and oplax morphisms give bounds in opposite directions

A min-plus matrix ``\Phi`` between two systems commutes with their Bellman
operators up to the order ([`intertwines`](@ref)):

- **oplax** when ``\Phi \otimes W_s \le W_t \otimes \Phi``;
- **lax** when ``W_t \otimes \Phi \le \Phi \otimes W_s``;
- **strict** when both hold.

Each kind of morphism carries one kind of section to the same kind:

| coarse graph | quotient rule | ``E`` is | ``E`` carries | so ``E V_c`` is |
|---|---|---|---|---|
| optimistic ("may", ∃∃) | ``\min_{x \in B} \min_{y \in C} W_{xy}`` ([`quotient_graph`](@ref)) | oplax | subsolutions to subsolutions | a lower bound, an admissible A* heuristic ([`lower_bound`](@ref)) |
| pessimistic ("must", ∀∃) | ``\max_{x \in B} \min_{y \in C} W_{xy}`` ([`robust_quotient`](@ref)) | lax | supersolutions to supersolutions | an upper bound and a control Lyapunov function ([`upper_bound`](@ref)) |

The proofs are one line each.
- **Optimistic:** if ``V_c \le W_c V_c \oplus P b``, then
  ``E V_c \le E W_c V_c \oplus E P b \le W E V_c \oplus b``. This uses the oplax
  square and ``EP \le I``.
- **Pessimistic:** if ``W_c V_c \oplus b_c \le V_c`` with ``b_c`` the *largest*
  terminal cost in each block, then ``W E V_c \oplus b \le E W_c V_c \oplus E b_c \le E V_c``.

Neither square commutes strictly unless the moves are constant on blocks. The gap
between the two quotients measures how much the aggregation hides.

### The lifting defect: when the optimistic bound is exact

[`lifting_defect`](@ref) makes that gap precise, arc by arc. For an arc ``B → C`` of
the optimistic quotient, the defect is the worst excess, over the states ``x`` of
``B``, of the cost of reaching ``C`` from ``x`` over the coarse arc's cost:

```math
ε(B, C) = \max_{x ∈ B} \big( d(x, C) − W_c(B, C) \big)^+ .
```

The arc *lifts* when its defect is zero: every state of ``B`` can follow it at no
extra cost. The path-lifting theorem says:

- **Exactness.** ``E V_c = V`` for every target that is a union of blocks if and only
  if every arc lifts.
- **Sandwich.** Adding the defects to the coarse arc costs gives an upper bound,
  ``E V_c ≤ V ≤ E V_c^{+ε} ≤ E V_c^{∀∃}``, never worse than the pessimistic quotient.
- **Two-sided error.** The worst error of the lower bound at a block is at least
  its largest defect.
- **Optimality.** The least sound upper weighting is the one-sided Hausdorff
  distance ``\max_{x ∈ B} d(x, C)`` between blocks; the defect-inflated costs
  approximate it locally.

Every value [`lifting_defect`](@ref) returns is an upper bound on the true defect
(a local scope or a search horizon only enlarges it), so these bounds stay sound.

### Truncated pushforwards: keeping modes of the stalk

Forgetting the velocity entirely (one value per stalk) gives a weak bound, because
a cell's cheapest state may be moving the wrong way. [`TruncatedPushforward`](@ref)
keeps a partition of each stalk into ``K`` *modes* (heading × speed bins of the
velocity). The coarse sheaf ``G`` then has stalk ``\mathbb R_{\min}^K`` over every cell, and
``E``, ``P`` act stalk by stalk over the identity of ``H``. With one mode it is the
velocity-blind quotient; with one mode per state it is an isomorphism. The modes
interpolate between the two, and the bound tightens monotonically.

### Changing the gauge

A change of gauge relabels the zero point of cost at every state. It leaves the
fine problem alone and changes what a coarse problem can see.

**The operation.** Pick a number ``\varphi(x)`` for every state, the *potential*,
and replace each edge cost by its reduced cost

```math
W^\varphi_{xy} = W_{xy} + \varphi(y) - \varphi(x),
\qquad\text{i.e.}\qquad
W^\varphi = D^{-1} \otimes W \otimes D,\quad D = \operatorname{diag}(\varphi).
```

In min-plus, ``D`` is an invertible diagonal matrix, so this is a similarity: an
isomorphism of min-plus sheaves that shifts each stalk's values by ``\varphi``. In
cochain terms, the edge costs are a 1-cochain ``c`` and the change is
``c \mapsto c + \delta\varphi``: it adds a coboundary.

**Why the fine problem does not change.** Along any path ``x_0 \to \cdots \to x_n``
the potentials telescope:

```math
\textstyle\sum W^\varphi(\text{path}) = \sum W(\text{path}) + \varphi(x_n) - \varphi(x_0).
```

Every path between the same endpoints shifts by the same amount, and every cycle
by nothing. So the optimal paths, the optimal policy and the reachable set are
unchanged. With ``\varphi = 0`` on the targets, the cost-to-go shifts by exactly
the potential: ``V^\varphi = V - \varphi``. Only differences of ``\varphi``
matter, which is why the name *gauge* fits.

The same move is familiar under other names:

- **A\* is Dijkstra in a gauge.** Dijkstra on the reduced costs with
  ``\varphi = h`` is exactly A\* with heuristic ``h``. The heuristic must be
  consistent, ``h(x) \le W_{xy} + h(y)``, precisely so that ``W^\varphi \ge 0``.
- **Johnson's algorithm** (1977) reweights a graph this way to remove negative
  edges.
- **Linear programming:** ``\varphi`` is the dual vector and ``W^\varphi`` the
  reduced costs.
- **Reinforcement learning:** potential-based reward shaping (Ng, Harada &
  Russell 1999), proved not to change the optimal policy.

**Why the coarse problem does change.** The pushforward takes a *minimum over each
block*, and a minimum does not commute with adding different amounts to different
states. So the same fine problem in two gauges has two different quotients.

A line of states ``0 - 1 - 2 - 3`` with unit steps and the goal at 3 has
``V = (3, 2, 1, 0)``. Aggregate it into blocks ``A = \{0, 1\}`` and
``B = \{2, 3\}``.

- **Original gauge.** The quotient has one edge ``A \to B`` of cost 1, and the
  steps inside each block are contracted away. The bound at state 0 is 1, against
  a true cost of 3.
- **Gauge ``\varphi = V``.** A step toward the goal costs 0 and a step away costs
  2. The coarse cost-to-go is 0 everywhere, and the bound is
  ``\varphi + 0 = 3``, which is exact.

The gauge moves information *out of* the edge costs, which the minimum over blocks
destroys, *into* the potential, which is added back state by state.

[`Abstraction`](@ref)`(W, ψ; potential = φ)` does this. The quotient is built
from ``W^\varphi``, and [`lower_bound`](@ref) returns ``\varphi + E V^\varphi_c``,
a refinement that copies block values *on top of* ``\varphi`` rather than
flat. The potential must be consistent, for example the cost-to-go of a
relaxation, and then the bound is again a consistent heuristic.

This is the min-plus counterpart of aggregation with a near-nullspace mode in
algebraic multigrid. There too the choice of basis decides what the coarse space
must represent: `TruncatedPushforwardCoarseSpace` aggregates ``D_h z`` rather than
constants. Here the "mode" is ``\varphi``, and the coarse space only represents
corrections to it.

### Elimination is min-plus harmonic extension

On a real sheaf, [`harmonic_extension`](@ref) fixes the values on boundary cells
and solves for the interior. The min-plus counterpart is the Schur complement
onto a separator ``S`` ([`schur_complement`](@ref)):

```math
W_S = W_{SS} \oplus W_{SI} \otimes W_{II}^* \otimes W_{IS}, \qquad (W_S)^* = (W^*)_{SS}.
```

The identity on the right says that eliminating the interior loses nothing on
``S``; it holds in any closed semiring (Lehmann 1977). Cutting the state space into
cells and eliminating every interior onto the *portals* gives an [`Overlay`](@ref).
The cells play the role of the cover in the Schwarz methods, and each shortcut is
the harmonic extension of boundary data through a cell.

Two commutation facts organize the layered controller:

- **Restriction commutes with pushforward** (base change). If ``U`` is a union of
  fibres of ``\mu``, restricting the fine graph to ``\mu^{-1}(U)`` and then pushing
  forward equals pushing forward and then restricting to ``U``. So the coarse bound
  stays admissible inside a corridor.
- **Aggregation and elimination commute only laxly.** The two orders give
  ``AE \le EA \le \text{exact}``, with equality when the aggregation is injective on the
  eliminated interior. Eliminating before aggregating is the tighter order. On the
  ``64^2 \times 5^2`` lattice the two squares agreed on 46% of the entries, and the
  cost-to-go reached 0.961 (aggregate first) and 0.987 (eliminate first) of exact
  on average.

### The same pushforward in the real semiring

The coarse operator ``E^{\mathsf T} W E`` is the Galerkin operator of the HJB
multigrid ([`GridMultigrid`](../api/grid_schwarz.md)). What changes is how the
coarse solution is used:

| | real ``(+, \times)``: HJB policy evaluation | tropical ``(\min, +)``: minimum time |
|---|---|---|
| pullback ``E`` | copy a block's value | copy a block's value |
| pushforward | block average ``T`` (divides) | block minimum ``P``, ``PE = I`` for free |
| Galerkin operator | ``T A E`` | ``P W E``, the quotient graph |
| coarse solution used as | a correction of the residual | a bound (heuristic or Lyapunov function) |
| because | subtraction exists | no subtraction; the order replaces it |
| exactness comes from | iterating to convergence | optimality of A* with an admissible heuristic |
| local solves | Schwarz: pinned harmonic extension | overlay: Schur complement onto portals |

## The layered controller and its certificates

[`LayeredPlanner`](@ref) splits control into two layers:

- a **planner** routes over the coarse sheaf: cells, each with a mode of its stalk
  ([`plan_route`](@ref));
- a **low-level controller** searches the full 4D state graph exactly, but only
  over the pullback ``\mu^{-1}(U)`` of a *corridor* ``U`` around the route
  ([`corridor`](@ref), [`layered_path`](@ref)). ``U`` contains the cells within
  Chebyshev distance `width` of the route, with modes within `mode_width` steps.

``U`` is a union of fibres, so base change keeps the coarse lower bound admissible
inside it. If the corridor holds no path, the controller widens it.

[`certified_path`](@ref) adds two certificates computed per query, neither of which
needs the optimum:

- **Stability (the lax side).** The low level follows the shortest-path tree of
  ``V_U``, the cost-to-go restricted to the corridor. [`is_lyapunov`](@ref) checks
  that ``V_U`` is a Lyapunov morphism for the closed loop, in the discrete setting
  of categorical Lyapunov theory (Ames, Moeller & Tabuada; Ames, Mattenet &
  Moeller). Concretely, it is positive definite and decrescent, its decrease is
  strict, and the domain is closed under the closed loop. Every state above the
  corridor that can reach the target does so within ``V_U`` steps.
- **Performance (the oplax side).** With the lower bound ``E V_c`` from an
  optimistic abstraction, the ratio ``V_U(x) / E V_c(x)`` is a guaranteed bound on
  the closed loop's cost relative to the optimum:

```math
E\,V_c^{\text{opt}} \;\le\; V^\star \;\le\; V_U .
```

The controller is sandwiched between an oplax and a lax morphism, and the width of
the sandwich is its suboptimality certificate.

We also tried transporting a Lyapunov function from a *pessimistic* coarse sheaf
([`upper_bound`](@ref)). It is sound but conservative. A ∀∃ coarse edge needs
every state in a mode to make the move, and wide speed bins rarely satisfy that.
The per-query corridor certificate is the practical form of the lax side.

## Relation to the literature

The ingredients are classical. What the sheaf view adds is a single vocabulary in
which the heuristic, the coarse grid, the abstraction, the Lyapunov function and
the elimination are the same few constructions:

- pushforwards and their residuals;
- lax and oplax morphisms;
- harmonic extension;
- base change.

**Path problems and dioids.** Shortest paths, bottleneck paths and transitive
closure are one algebraic path problem over a closed semiring (Carré 1971; Lehmann
1977; Tarjan 1981; Gondran & Minoux 2008). Gaussian elimination, and with it the
Schur complement identity used above, works in that generality. PathPack's chordal
semiring LU is a sparse, supernodal instance.

**Max-plus linear systems and residuation.** The Galois connection ``E \dashv P``,
the projection ``EP`` and the order-theoretic replacement for subtraction come from
residuation theory in idempotent semimodules (Baccelli, Cohen, Olsder & Quadrat
1992; Cohen, Gaubert & Quadrat 2004). The closest relative of our coarse problems
is the max-plus finite element method (Akian, Gaubert & Lakhoua 2008). It projects
the Bellman operator onto a max-plus subsemimodule with exactly this residuated
projection and gets one-sided error bounds. McEneaney's max-plus methods (2006)
use the same linearity to avoid grids altogether.

**Aggregation in dynamic programming.** Aggregating states and solving the
aggregate problem is an old idea (Bertsekas & Castañón 1989). So are multigrid
methods for discrete-time stochastic control (Chow & Tsitsiklis 1991). In the real
semiring these give corrections or approximations. In the min-plus semiring hard
aggregation gives *bounds*, which is why the coarse solution can drive an exact
search instead of an iteration.

**Abstraction heuristics and route planning.** An optimistic quotient used as an
A* heuristic is Hierarchical A* (Holte, Perez, Zimmer & MacDonald 1996). The same
idea underlies pattern databases (Culberson & Schaeffer 1998), the all-pairs table,
and merge-and-shrink abstractions in classical planning (Helmert, Haslum &
Hoffmann 2007). The overlay is customizable route planning (Delling, Goldberg,
Pajor & Werneck 2011), a cousin of contraction hierarchies (Geisberger, Sanders,
Schultes & Delling 2008). Its preprocessing is the Schur complement.

**Abstraction-based control.** The two quotients are the *may* and *must*
transitions of modal transition systems (Larsen & Thomsen 1988). The must
abstraction is what symbolic controller synthesis refines (Tabuada 2009). Feedback
refinement relations (Reissig, Weber & Rungger 2017) give conditions under which an
abstract controller refines soundly. The lax square ``W \otimes E \le E \otimes W_c``
plays that role here, stated as a morphism rather than a relation.

**Lyapunov certificates and layered control.** Verifying a controller by a
Lyapunov function along a planned route echoes LQR-trees and funnel libraries
(Tedrake, Manchester, Tobenkin & Roberts 2010). The planner/low-level split is a
layered control architecture in the sense of Matni, Ames & Doyle (2024). The
certificate is a Lyapunov morphism in categorical Lyapunov theory (Ames, Moeller &
Tabuada 2025; Ames, Mattenet & Moeller 2025). There a Lyapunov function is a lax
morphism into a stability object, and the lax/oplax sandwich above is its
two-sided form for optimal control.

**Sheaves.** The real side of the package follows Curry (2014) and Hansen & Ghrist
(2019): global sections, Laplacians, harmonic extension. Sheaves valued in
lattices and their fixed-point "Laplacian" (Ghrist & Riess 2022) are the natural
home for the min-plus constructions. There, sections become sub- and
supersolutions and Tarski's theorem replaces linear algebra.

## Limitations and open directions

- **Deterministic only.** Expectations bring back ``(+, \times)`` and with it
  subtraction, so stochastic problems belong to the real-semiring solvers
  ([`GridPolicyIteration`](@ref CellularSheaves.ControlSheaves.DoubleIntegratorHJB.GridPolicyIteration) with [`Multigrid`](@ref CellularSheaves.ControlSheaves.DoubleIntegratorHJB.Multigrid)).
- **Lattice dynamics.** The state graph is a discretization with finitely many
  controls and a unit time step. The bounds are exact *for the lattice*, not for
  the continuous system.
- **Fibres are the bottleneck.** Wide velocity fibres make the pessimistic quotient
  conservative and the overlay dense. Choosing modes adaptively, per cell and per
  query, is the obvious next step.
- **Choosing the abstraction** (block sizes, ``k``, modes) is manual. In planning,
  merge-and-shrink and counterexample-guided abstraction refinement automate it.
  The lax/oplax gap is a natural refinement signal: refine where the two quotients
  disagree.
