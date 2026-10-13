# Tropical Abstraction: Coarse Grids for Optimal Control

This page explains the min-plus ("tropical") side of the package
([`TropicalAbstraction`](@ref CellularSheaves.NetworkSheaves.TropicalAbstraction)):

- coarse grids for deterministic optimal control;
- a layered controller built on them;
- the certificates that come with that controller.

It covers three things: the idea, how it reads as sheaf theory, and how it relates
to the literature it borrows from. The last section reports what the studies in
`docs/scripts/tropical_*.jl` measured, including where the approach does not pay off.

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

## Cooperative multi-agent planning

Several agents share a map, one planner controls them all, and each pays for its
own steps, waiting included, until it rests at its goal. The joint state graph is
the product of the agents' graphs with the collisions removed (vertex and swap
conflicts), so it grows exponentially with the number of agents. The question is
how much of the *coordination cost*, the waiting and detours that the agents
impose on one another, a polynomial-size abstraction can see.

**Work in the decoupled gauge.** Let ``d_i`` be agent ``i``'s single-agent
cost-to-go and ``\varphi = \sum_i d_i``. In the gauge of ``\varphi``:

- an agent moving along a shortest path costs 0;
- waiting costs 1;
- stepping away from its goal costs 2.

So the joint cost-to-go is ``\varphi`` plus the total **delay**, and a heuristic
only has to bound the delay. Without the gauge, coarsening the other agent by
position throws away its travel cost, as in the line example above, and the bound
cannot beat the plain sum of solo distances. The summary below keeps the other
agent's remaining distance anyway, and with it the gauge turns out to be
redundant: it matters exactly for summaries that discard ``d_j``.

**Summarize the other agent by the Reeb graph of its value function.** In agent
``i``'s coarse problem, ``i`` keeps its exact state and agent ``j`` is replaced by
a mode, and the joint graph is pushed forward along that map. Which mode decides
everything:

- **Position blocks** of agent ``j`` leak *timing*. A block forgets where in it
  ``j`` is, so consecutive coarse steps can put ``j`` anywhere in the block. ``j``
  can then reach a corridor early and slip through before ``i`` arrives.
- **Level sets** of ``d_j`` (agent ``j``'s remaining distance) fix the timing,
  because ``j``'s delay-free moves go from level ``k`` to ``k - 1``. But a level
  set contains equidistant cells on opposite sides of an obstacle, so ``j`` can
  hop across a wall for free.
- **Arcs of the Reeb graph of ``d_j``** fix both. A level-``k`` state's arc is its
  component of the band ``k \le d_j \le k + 1``, the discrete level set
  ``d_j = k + \tfrac12``, so cells on opposite sides of an obstacle fall in
  different arcs. Reeb graphs are themselves cosheaves (de Silva, Munch & Patel
  2016), so the summary is a pushforward along a map with topological content,
  not a hand-made partition. The merge tree, which uses components of the
  superlevel sets ``\{d_j \ge k\}``, does *not* work: those sets reconnect through
  distant routes.

**The interaction sheaf.** For ``N`` agents, the coarse problems become edge
stalks on the interaction graph:

- **Vertex stalk:** an agent's own state graph in the gauge of its ``d_i``.
- **Edge stalk of a pair ``(i, j)``:** the pushforward of the *pair graph* (both
  agents moving at once, collisions removed) along a map that keeps ``i``'s state
  and summarizes ``j`` relative to it. Cover the map with blocks:
  - **Near** (``j``'s block next to ``i``'s): ``j``'s cell and velocity mode are
    kept. Every collision happens here.
  - **Far:** ``j`` becomes the pair (Reeb arc of ``d_j``, block). That records its
    progress, and where it is relative to ``i`` at the resolution of the blocks.

  The cost-to-go ``\Delta_{ij}`` of the edge stalk bounds the pair's delay.
- **Gluing:** an ``N``-agent plan restricts to a feasible plan of every pair, and
  delays are nonnegative. So for any matching of the interaction graph,
  ``\varphi + \sum_{(i,j)} \Delta_{ij}`` is admissible, and the heuristic takes
  the best matching. This is the pairwise heuristic of conflict-based search (Li,
  Felner et al. 2019), with the pair bounds coming from pushforwards instead of
  pair solves. The interaction graph keeps only pairs whose *tubes* meet: the
  cells each agent can use within a small delay. Dropping a pair never breaks
  admissibility, only tightness.

**Base change makes the edge stalk exact without the pair graph.** The far summary
depends on ``j``'s state alone, and the near/far split depends only on blocks. So
the far region is a union of fibres of a *product* map ``\mathrm{id} \times \pi_j``.
On a product, pushforward commutes with the tensor product (simultaneous moves),
[`tensor_product`](@ref). Over the far region the edge stalk is therefore
``W_i \boxtimes Q_{\pi_j}(W_j)``: agent ``i``'s graph times the quotient of
``j``'s own graph, with its loops ([`quotient_graph`](@ref) with
`loops = true`). Only the near region and its boundary are enumerated, and their
size does not grow with the map. The result *equals* the pushforward of the pair
graph, which the studies assert state by state; it is not an over-approximation.
[`reeb_arcs`](@ref) computes the arcs as a `GraphHomomorphism`, using bands between
consecutive distinct levels, so costs need not be unit.

An earlier version kept only ``j``'s Reeb arc when far, with no relative
coordinate. It had to be patched with an occupancy check and still lost a third of
the coordination cost on double integrators. The missing piece was the relation,
which is what an edge stalk is supposed to carry.

**Changing the gauge on the cover.** The matching uses each agent in at most one
edge. A gauge transformation on the cover uses every edge instead.

- **The split.** Each agent's delay is split, state by state, among the edges at
  it: weights ``w_{ie}(s) \ge 0`` with ``\sum_e w_{ie}(s) \le 1``. Edge
  ``e = \{i, j\}`` then carries the pair problem with arc costs
  ``w_{ie}\,\mathrm{delay}_i + w_{je}\,\mathrm{delay}_j``.
- **Why it is a gauge.** The split moves cost between the cells of the cover
  without changing the cost of any global section. This mirrors how ``\varphi``
  adds a coboundary on one agent's graph without changing the cost of any path.
- **The two gauges compose.** In the delay gauge every arc cost is nonnegative and
  an agent alone costs 0. So the split is well posed, and the whole bound comes
  from the edges.
- **Validity.** For every split, the sum over edges of the weighted pair optima is
  an admissible and consistent heuristic, because a joint move restricts to a move
  of every pair. A matching is the split with ``w \equiv 1`` on matched edges.
- **Optimization.** `docs/scripts/tropical_cover_gauge.jl` maximizes the split at
  the start by projected supergradient ascent from the matching. Its weighted edge
  stalk is the same pushforward as before: agent ``j``'s arc cost on a far coarse
  arc is the least weighted delay over the fibre.

**Edges are not enough: the cover needs 2-cells.** Weights are at most 1, so every
edge factor is at most its pair's own optimum. No split, however clever, can see
an interaction that no pair sees, and three agents meeting at a corridor is
exactly that. In one four-agent instance:
- every pair except ``(1, 4)`` has optimal delay 0, and ``(1, 4)`` has delay 3;
- the triple ``\{1, 3, 4\}`` has delay 5, which is also the optimum for all four.

The pairwise data is exhausted at 3. The gap is an obstruction that lives on a
triangle of the interaction graph, not on any edge.

So extend the cover to the clique complex of the interaction graph, with a stalk
on every triangle: the triple problem.
- **Bound.** The delays ``x_i`` of a joint plan restrict to a plan of every simplex
  ``\sigma``, so they satisfy ``\sum_{i \in \sigma} x_i \ge \Delta_\sigma`` for
  every edge and triangle. The smallest such integer 0-cochain, the *simplicial
  cover*, bounds the total delay. With edges alone it is the weighted vertex cover
  of Li, Felner et al. (2019).
- **Its relation to the gauge.** Its linear relaxation is the linear-programming
  dual of the scalar split over simplices,
  ``\max \sum_\sigma y_\sigma \Delta_\sigma`` with ``\sum_{\sigma \ni i} y_\sigma \le 1``.
  So the gauge on the cover and the cover by 0-cochains are two sides of one
  duality, and integrality strengthens the cover side.

**The cover as a sheaf on a 2-complex.** `TropicalCovers` builds this structure in
the package, using these pieces:

- **The complex** is the clique complex of the interaction graph
  ([`clique_complex`](@ref)): agents are vertices, interacting pairs are edges,
  and mutually interacting triples are triangles.
- **Stalks** ([`ProductStalk`](@ref)). Every simplex carries the pushforward of
  its agents' joint graph along a summary map that is *symmetric* in the agents.
  Each agent is replaced by its own summary (Reeb arc of its cost-to-go, block)
  unless its summary is in contact with that of another agent of the simplex, and
  then it keeps its exact state. No agent of a simplex is privileged, on edges or
  on triangles.
- **Exactness.** Contact is decided on summaries, so the coarse states are unions
  of fibres, and the fibres are products. As long as agents out of contact cannot
  collide, the stalk is enumerated move by move and *equals* the pushforward of
  the joint graph. The tests and the study check this against the pair graph.
- **Restriction maps** ([`restrict`](@ref)) forget agents. An agent kept exact
  only because of a contact outside the face falls back to its summary. The square
  with the summary maps commutes, and every arc restricts to an arc that costs no
  more. So a face's cost-to-go, pulled back, is a subsolution on the simplex: a
  cellular sheaf of min-plus modules with pullbacks as restrictions.
- **Contact may grow with the dimension.** A simplex must resolve its agents at
  least as finely as its faces, or the restriction is undefined. So triangles may
  keep agents exact over a wider range than edges.
- **Building and gluing.** Edge stalks are built and solved once
  ([`stalk_graph`](@ref)). Triangle stalks are searched on demand from the query
  state, with the pulled-back edge values as the A* heuristic
  ([`section_value`](@ref)). [`cover_bound`](@ref) glues the values by the
  simplicial cover.

## What the studies found

All timings are on one HiPerGator CPU node, averaged over 30 random
queries on the ``128^2 \times 13^2`` lattice. They come from
`docs/scripts/tropical_*.jl`, and the data are in `docs/figures/tropical_*`.

**Coarsen time with space.** Aggregating positions alone degrades the bound
quickly: each coarse step can cross a whole block. Using the ``k``-step graph
``W^{\otimes k}`` (`steps = k`) restores it, with an admissible but inconsistent
heuristic, so A* reopens states. With 1024 coarse states, the space–time
abstraction gives a better bound than the position-only one with 43k.

| abstraction | coarse states | bound / cost | A* expansions | time / query |
|---|---|---|---|---|
| exact Dijkstra | – | – | 1.73M | 0.333 s |
| positions 2×2 | 692k | 0.95 | 11k | 0.042 s |
| positions 8×8 | 43k | 0.57 | 248k | 0.073 s |
| positions 32×32 | 2.7k | 0.16 | 992k | 0.283 s |
| space–time 4×4, ``k = 2`` | 50k | 0.87 | 23k | 0.055 s |
| space–time 16×16, ``k = 4`` | 1024 | 0.66 | 185k | 0.059 s |

**Velocity modes make the planner useful.** Exact 4D A* with the truncated
pushforward's bound expands 27–30 times fewer states than Dijkstra:

| bound | expansions | time / query |
|---|---|---|
| none (Dijkstra) | 1.19M | 0.245 s |
| velocity-blind (1 mode) | 142k | 0.052 s |
| 8 headings × 3 speeds (25 modes) | 44k | 0.021 s |
| 16 headings × 3 speeds (49 modes) | 39k | 0.021 s |

Restricting the search to a corridor trades optimality for commitment. The
velocity-blind corridor of width 4 finds a path on the first try 57% of the time,
at 1.06 times the optimum on average and 1.55 in the worst case. With 49 modes and
a mode width of 1 the same width gives 1.006 on average and 1.087 in the worst
case. Widened to 16 it is optimal on every query.

**Certificates.** Every query was certified stable. The guaranteed suboptimality
bound is looser than the true ratio but tightens with the modes, because the lower
bound does:

| planner | true cost / optimum | guaranteed bound (mean, worst) |
|---|---|---|
| velocity-blind | 1.062 | 2.18, 7.50 |
| 8 headings × 3 speeds | 1.013 | 1.47, 2.50 |
| 16 headings × 3 speeds | 1.006 | 1.30, 1.79 |

The pessimistic quotient certifies 0% of the reachable states with one mode, 9
states with 25 or 49 heading-major modes, and 59% with 8 headings × 6 speeds. At
that resolution its upper bound averages 2.7 times the optimum. Only the
isomorphism (169 modes) certifies everything exactly.

**Elimination and tables.**
- **Overlay.** On ``64^2 \times 5^2`` it cuts expansions from 86k (Dijkstra) to
  6.7k (overlay with A*) but not time: 0.005 s against 0.0025 s for plain A*. The
  velocity fibres make portals numerous and the shortcuts dense (7M of them).
- **All-pairs tables.** Building a full coarse table, a pattern database and in
  control terms the min-plus fundamental solution of the coarse Bellman equation,
  takes 0.15–2 s. After that the bound for any goal is a column lookup.
- **Many goals at once.** PathPack's semiring LU factorizes the coarse closure once
  for every goal. On these sparse, nearly planar quotients, one Dijkstra per goal
  was faster in most configurations. PathPack won only narrowly in a few, such as
  1024 goals on the 16×16 quotient (0.28 s against 0.31 s). Fill-in costs about as
  much as repeated near-linear searches.

**Bottleneck semiring.** In the max-min semiring, [`reach_avoid`](@ref) finds the
route with the most clearance from the walls. The fastest routes take 23–29 steps
with 0.5 cells of clearance, while the safest take 63–72 steps with up to 9.7 cells.

**Multiple agents.** The coordination cost is the gap between the exact joint cost
and the sum of solo distances. "Recovered" is the share of that gap a heuristic's
bound closes at the start state, averaged over queries where the agents must
cross. All heuristics are admissible and every search returns the optimum, so the
expansions measure heuristic quality alone.

*Two grid agents, 32² map with a wall and two one-cell corridors (712k joint
states):*

| other agent summarized by | coarse states (both agents) | recovered | joint A* expansions |
|---|---|---|---|
| nothing (sum of solo distances) | – | 0% | 1243 |
| position blocks 4×4 + near field | 62k | 0% | 1237 |
| level sets of ``d_j`` + near field (+ exact corridor cells) | 105k (125k) | 0% | 1237 (1174) |
| **Reeb arcs of ``d_j``** | 173k | **89%** | **398** |
| exact joint cost-to-go | 712k | 100% | 380 |

With Reeb arcs, adding the near field or dropping the gauge changes nothing.

*Three grid agents, 10² map (593k joint states), pair bounds glued over the
interaction graph:*

| pair bound | coarse states | recovered | joint A* expansions |
|---|---|---|---|
| none | – | 0% | 306 |
| edge stalks (block 2), glued | 34k | 100% | 110 |
| exact pair delays, glued | 21k (pair graphs) | 100% | 109 |
| exact joint cost-to-go | 593k | 100% | 103 |

*Two double integrators, 16² map, velocities in ``\{-1, 0, 1\}^2`` (2052 states
each, 4.2M pair states). The near field keeps the other agent's velocity exactly
or as a mode; each row is the exact pushforward along its own map:*

| near field keeps | coarse states (both directions) | recovered | pair A* expansions | build |
|---|---|---|---|---|
| nothing (decoupled) | – | 0% | 1289 | – |
| exact velocity (9) | 1.97M | 100% | 571 | 7 s |
| rest + 4 headings (5) | 1.54M | 67% | 720 | 7 s |
| rest / moving (2) | 1.22M | 33% | 994 | 6 s |
| (exact pair delay) | 4.2M | 100% | 571 | 8 s |

*Four to eight grid agents, 20² map with two corridors, half crossing each way
(4 instances each).* Here the joint graph cannot be built. The **sandwich**
compares:
- the sum of solo distances;
- the sheaf bound at the start;
- the optimum, found by independence detection (agents merge into jointly searched
  groups only when their plans collide);
- a feasible plan from prioritized planning (each agent plans around the
  reservations of those before it).

| N | interacting pairs | recovered: matching / vertex cover / cover gauge / edges + triangles | certified ratio (feasible / bound): matching → edges + triangles | joint expansions, decoupled → sheaf heuristic |
|---|---|---|---|---|
| 4 | 5.5 of 6 | 60% / 60% / 60% / **100%** | 1.007 → **1.000** | 6444 → 888 |
| 6 | 12 of 15 | 85% / 85% / 85% / **100%** | 1.012 → **1.000** | 5265 → 1946 |
| 8 | 26 of 28 | – | 1.046 → 1.034 | groups of 5–7: out of reach |

Observations:

- **The gauge on the edges finds nothing the matching misses.** The ascent
  reproduces the matching on every instance, and so does the vertex cover (one
  N = 8 instance aside, by 1). Their ceiling is the sum of the pair optima, and
  the pairs had already given what they had.
- **Triangles close the gap.** On every solved instance with a gap, the
  remaining coordination cost came from triangles whose triple delay exceeds every
  one of their edges: 1 of 4 triangles at N = 4, 2 of 2 at N = 6, and 1–4 per
  instance at N = 8. With the triangles, the bound equals the optimum on all 8
  solved instances. At N = 8 it raises the bound on every instance, and certifies
  the prioritized plans within 1–8% of optimal (previously 2–9%), still without
  the optimum. Solving the 36–56 triangles exactly takes 1–2 minutes per
  instance.

**The 2-complex with symmetric stalks** (`docs/scripts/tropical_cover_triangles.jl`,
same instances, block 2). Edges keep agents exact within one block. The table
varies the agents' summary and the triangles' contact radius:

| summary | triangle contact radius | edges exact (of 175) | three-agent interactions seen (of 14) | recovered, N = 4 / 6 | bound at N = 8 vs. exact triangles | edge stalk states |
|---|---|---|---|---|---|---|
| (Reeb arc, block) | 1 block (same as edges) | 171 | 9 | 60% / 89% | below on 2 of 4 | 84k |
| (Reeb arc, block) | 5 blocks | 171 | 13 | 100% / 100% | equal on 4 of 4 | 84k |
| refined by future arcs | **1 block** | **174** | **13** | **100% / 100%** | below on 1 of 4 | 86k |
| refined by future arcs | 5 blocks | **174** | **14** | **100% / 100%** | **equal on 4 of 4** | 86k |
| (exact triple search) | – | – | 14 | 100% / 100% | – | (pair graph 132k) |

Each triangle takes on average 750–1900 expansions (at most about 12k) and at
most 7 s per instance for all of them, against up to 109 s for exact triple
searches.

Findings:

- **Symmetric edges are exact pushforwards and nearly as tight.** The edge stalk
  has 84–86k coarse states, against 132k for the pair graph. With refined
  summaries its value equals the exact pair delay on 174 of 175 edges (171 with
  the plain summary), against 175 for the asymmetric stalks. Each asymmetric edge
  takes the larger of two directions, in each of which one agent is never
  summarized.
- **The fix belongs to the vertices, not to the triangles.** With the plain
  summary, the triangles needed a contact radius of 5 blocks, half the map. With
  summaries refined by their future, radius 1, the same as the edges, already
  sees 13 of 14 interactions and closes every solved gap. Radius 5 sees all 14,
  and the bound then equals the one from exactly solved triangles on every
  instance. That includes all four N = 8 instances, where it certifies the
  prioritized plans within 1–8% of optimal.
- **Search on demand is cheap because of the restriction maps.** Pulled back, the
  edge values guide each triangle's search to a few thousand expansions. That is
  about 15× faster than the exact triple searches.
- **What leaks: sliding along a level set.** A summarized agent is known by its
  (Reeb arc, block). Each step may start from any cell of that fibre, so it can
  slide sideways along its level set by one cell per step at no cost. In open
  space a level set can touch two routes, and sliding switches routes for free.
  In one example, two agents meet head-on in a corridor, but in the stalk one of
  them slides onto the other corridor's route and avoids the conflict entirely.
  Widening the edges' contact to 2 blocks did not help: the slide happens while
  the agents are still far apart.
- **The cure: summaries whose moves lift.** The summary map should let every coarse
  zero-cost walk lift to a real one from each state of its fibre. Asking for that
  exactly ([`stable_refinement`](@ref), the coarsest partition stable for
  zero-delay moves; Paige & Tarjan 1987, regular equivalence in Everett & Borgatti
  1994) removes the leak but, on a grid, splits the summary down to single cells:
  131k states per edge stalk, the pair graph again. What the stalks need is
  weaker: the states of a fibre must agree on which *Reeb arcs* lie ahead along
  zero-delay moves ([`future_refinement`](@ref)). Then a summarized agent can only
  reach, for free, parts of its future its actual state could also reach. It
  cannot slide onto another corridor's route. The refinement costs about 2k
  coarse states per edge stalk.
- **What still leaks.** With refined summaries, one edge in 175 and, at radius 1,
  one triangle in 14 are still below the exact value. Sharing a set of reachable
  arcs does not fix *when* they are reached, so some lateral freedom remains.

**Parallel construction and larger maps** (`docs/scripts/tropical_cover_scaling.jl`:
eight agents, two instances per map, future-refined summaries, triangle radius 1).
The cover's work is a set of independent tasks:
- one summary per agent, which is vertex data;
- one explicit stalk per edge, built from its two agents alone;
- one search per triangle, which only reads finished edge stalks.

`CoverSheaf` and [`cover_bound`](@ref) run these tasks over Julia's threads. Only the
gluing, a small integer program on the agents, is sequential.

| map (block) | edge stalk vs. pair graph | certified ratio (prioritized / bound) | build, 1 thread | build, 16 threads | triangles, 1 → 16 threads |
|---|---|---|---|---|---|
| 40² (2) | 1.39M vs. 2.32M states (60%) | 1.017, 1.000 | 45–54 s | 13–21 s | 3.5–6.3 s → 0.6–1.2 s |
| 64² (4) | 3.5M vs. 15.8M states (22–23%) | 1.000, 1.009 | 177–241 s | 53–71 s | 3.3–10.2 s → 0.7–2.3 s |

- **Compression grows with the map, if the blocks do.** With blocks of 2 cells,
  (Reeb arc, block) has nearly one summary per cell, and the edge stalk is 60% of
  the pair graph. On 64² with blocks of 4 cells it is 22%. The certificate stays
  within 2%, and on one 64² instance it proves the prioritized plan optimal.
- **Triangles stay cheap.** At every size, the triangle searches guided by the
  pulled-back edge values take a few thousand expansions each, and seconds for
  all of them together.
- **Threads give 3–4× on 16 cores, limited by memory traffic.** The compute
  itself scales: a pure-arithmetic loop scales perfectly on these nodes, and the
  tasks share nothing. But one 40² edge stalk streams through about 32 million
  arcs. With 16 built at once, each takes about twice as long as alone (6.6 s
  against 3.1 s), so 16 concurrent edges finish 6.3× faster than in sequence. A
  full cover runs its 20–28 edges in two waves on 16 threads, hence 3–4×.
- **Allocation.** Allocation was the first bottleneck and is now small. Each
  explicit stalk used to allocate 1.7 GB of short-lived arc arrays, and garbage
  collection paused every thread.

  Four changes removed it:
  - precomputed move lists;
  - building the arcs in compressed form, one state at a time;
  - a dense index over the box of coarse states, in place of hashing;
  - reusing arc buffers across stalks when the arcs are not kept.

  An edge stalk now allocates 272 MB, mostly its own states, index and
  cost-to-go, and a 40² instance builds on one thread in 54 s instead of 117 s. Faster thread
  scaling would have to cut memory traffic, not allocation.
- **Across nodes, over MPI** (`docs/scripts/tropical_cover_mpi.jl`, `hpg-milan`,
  AMD EPYC 75F3, two instances per map). The decomposition follows the sheaf:
  - every rank builds the shared structure itself (summaries, complex, product
    stalks), deterministically;
  - edges are dealt to ranks, longest first, and built over each rank's threads;
  - each edge's cost-to-go is broadcast, and every rank enumerates its states and
    index again (`StalkGraph(P, Δ)`);
  - triangles are dealt to ranks, and only their values are reduced;
  - rank 0 glues the cover.

  Every configuration gives exactly the same bound as one process.

  | layout | 40² total (edges) | 64² total (edges) |
  |---|---|---|
  | 1 node, 1 rank × 16 threads | 10–12 s (9.1–9.4) | 43–49 s (41–44) |
  | 1 node, 4 ranks × 4 threads | 9.5–13 s (7.4–8.9) | 33–37 s (28–32) |
  | 2 nodes, 1 rank × 16 threads each | 6.4–9.0 s (4.9–7.0) | 25–29 s (22–23) |
  | 4 nodes, 1 rank × 16 threads each | 6.6–10 s (4.9–6.6) | 23–24 s (19) |
  | 4 nodes, 4 ranks × 4 threads each | 10–12 s (4.9–5.1) | 27–31 s (15) |

  - **More nodes do help.** On 64², two nodes halve the edge phase, because each
    node brings its own memory bandwidth.
  - **Four nodes reach the floor set by the longest edge.** Its stalk takes about
    12 s alone, and 20–25 edges over 64 threads leave little for a fifth node to
    do. Going further means splitting a single stalk's states across threads.
  - **NUMA matters within a node.** These nodes have 8 NUMA domains. Four ranks of
    four threads build the edges 1.3–1.5× faster than one rank of 16, since each
    process keeps its memory near its cores.
  - **The exchange is the cost of many ranks.** Broadcasting every table to every
    rank, one edge at a time, takes 4–11 s with 16 ranks, against under 1 s with
    four. Only the ranks that search a triangle need the tables of its edges, so
    targeted sends (or one all-gather) would remove most of that.
  - **Each 4-node job waited 31–78 minutes for nodes** on a nearly full
    partition, under the group's investment QOS and with a 15-minute limit.
- **The edge stalk is exact and tight.** It equals the pushforward of the pair
  graph, and that pushforward along (Reeb arc, block) loses nothing against the
  exact pair delay on these maps. With exact pair bounds the gluing reproduced the
  three-agent optimum on every query.
- **The certificate scales; the search does not yet.** At eight agents no joint
  search finishes, but with triangles the sandwich certifies the prioritized plans
  within 1–8% of optimal without knowing the optimum. That is the lax/oplax
  sandwich of the single-agent controller, now over the interaction sheaf.
- **The heuristic cuts joint expansions 3–7×,** but each evaluation costs more, so
  on these small groups wall time is similar.
- **Size.** At these map sizes the pair graph is still solvable and comparable in
  cost. Each direction of the double-integrator edge stalk uses about a quarter
  of the pair graph's states (0.98M against 4.2M). Its near field is the part that
  does not shrink with the map, and velocity modes trade its size for tightness
  gradually.

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

**Multi-agent path finding.** Planning in the joint graph with a decoupled
heuristic, and coupling agents only where they conflict, is the core of
independence detection (Standley 2010), subdimensional expansion and M\* (Wagner &
Choset 2015), and conflict-based search (Sharon, Stern, Felner & Sturtevant 2015).
Pairwise bounds combined over a matching or vertex cover of the conflict graph are
the weighted dependency-graph heuristics of Li, Felner, Boyarski, Ma & Koenig
(2019). The interaction sheaf replaces their exact pair solves by pushforwards along
a summary of the other agent. The gauge is potential-based shaping (Ng, Harada &
Russell 1999) and Johnson's reweighting (1977). The summary is the Reeb graph of the
other agent's value function (Reeb 1946); contour trees (Carr, Snoeyink & Axen
2003) and Reeb graphs as cosheaves (de Silva, Munch & Patel 2016) are the
computational and categorical forms of it.

Splitting cost among the factors of a decomposition without changing the total is
*reparametrization* in MAP inference. Maximizing the sum of local minima over
reparametrizations is the dual of the local-polytope relaxation (Wainwright &
Jordan 2008), and it is what tree-reweighted message passing optimizes (Kolmogorov
2006). When pairwise factors leave a gap, that literature adds triplet clusters,
found where the gap concentrates (Sontag, Meltzer, Globerson, Jaakkola & Weiss
2008). That is the same move as the 2-cells of the interaction sheaf above.

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
- **Joint search beyond four agents.** Searching a group of ``m`` agents moves
  all of them at once: ``5^m`` successors per step on a grid. Operator
  decomposition (one agent moves per step) is the standard fix, but the pair
  bounds are only proved admissible at synchronized states. A half-moved
  configuration needs its own argument. Until then, large groups get the
  certificate but not the optimum.
- **Summaries between "future" and "stable".** Refining by future Reeb arcs leaves
  a small leak (one edge, one triangle at radius 1); full stability on zero-delay
  moves removes it but keeps almost nothing summarized on a grid. Intermediate
  refinements, futures seen at a few steps' horizon or through a finer partition
  than the Reeb arcs, are untested.
- **Compression needs larger blocks.** The edge stalks are 60–65% of the pair
  graph with blocks of 2 cells, and 22% on a 64² map with blocks of 4 cells.
  Larger blocks widen the contact region, so how block size, contact radius and
  tightness trade off on larger maps is open. The 64² certificates stayed
  within 1%.
- **Scaling stops at the longest stalk.** On one node 16 threads give 3–4×, limited
  by memory bandwidth. Over MPI, two nodes halve the 64² edge phase, and four
  reach the floor of the longest single edge stalk. Splitting one stalk's states
  across threads would lower that floor. Targeted exchange, sending each table
  only to the ranks that use it, would let more ranks per node pay off.
- **Bounds at the start only.** The triangle values are searched from one joint
  state. A heuristic at every state of a joint search needs them cached or built
  across queries.
- **Dense interaction graphs.** On maps that funnel every agent through the same
  corridors nearly all pairs interact, and the gluing gains little from sparsity.
- **Scale of the edge stalks.** The near field is fixed in size but large (a
  ``3 \times 3`` block window times the velocity modes). The edge stalks pay off
  only where pair graphs are too large to solve, which is also where they can no
  longer be checked against them.

## Reproducing the results

- **Abstractions:** `docs/scripts/tropical_abstraction_study.jl`, plus
  `tropical_multi_rhs.jl` for many goals at once.
- **Tables and overlay:** `docs/scripts/tropical_overlay_study.jl`, which also
  checks the commutation of aggregation and elimination.
- **Layered planner and certificates:** `docs/scripts/tropical_layered_study.jl`.
- **Pessimistic quotient:** `docs/scripts/tropical_lyapunov_study.jl`.
- **Multiple agents:**
  - `docs/scripts/tropical_multiagent_study.jl`: two agents, the ablation of
    gauge, summary and near field;
  - `docs/scripts/tropical_pair_sheaf.jl`: the interaction sheaf (edge stalks by
    base change, tubes, matching), shared by the next three;
  - `docs/scripts/tropical_multiagent_sheaf.jl`: three grid agents;
  - `docs/scripts/tropical_multiagent_cover.jl`: four to eight agents, the
    sandwich, independence detection, the gauge on the cover
    (`tropical_cover_gauge.jl`) and the simplicial cover with exact triangles;
  - `docs/scripts/tropical_cover_mpi.jl`: the same study over MPI, run with
    `srun --mpi=pmix` (one or more ranks per node, each with threads);
  - `docs/scripts/tropical_cover_scaling.jl`: 20² to 64² maps, eight agents,
    the cover built and evaluated over threads (arguments: map size, agents,
    instances, `profile`; `COVER_BLOCK` sets the block size), run with
    `julia -t T`;
  - `docs/scripts/tropical_cover_triangles.jl`: the same instances on the
    2-complex of `TropicalCovers`, with symmetric edge and triangle stalks
    (arguments: map size, instances, triangle and edge contact radius, and the
    summary: `reeb`, `future` or `stable`);
  - `docs/scripts/tropical_multiagent_double_integrator.jl`: double-integrator
    pairs.

An interactive visualization is in `docs/viz/tropical_demo.html`, built from
`docs/scripts/tropical_viz_export.jl`. It shows the value function, the coarse
route and corridor, the closed-loop phase plot and the Lyapunov sandwich.
