# 019 — Path lifting: the lifting defect of an abstraction and the bound it gives

**Status:** implemented. `lifting_defect(W, ψ; scope, horizon)` is in
`TropicalAbstraction`, with tests in `test/network_sheaves/TropicalAbstraction.jl`.
The global scope runs over threads and accepts a `horizon` beyond which defects are
reported as ∞ (still upper bounds), which the double-integrator studies need,
because lifts pass through neighbouring blocks. The theorem is stated and proved in
`docs/papers/path-lifting/paper.tex` (Section 4) and checked by
`docs/scripts/path_lifting_check.jl`.

## Mathematical Background

**Setting.** A finite *cost graph* ``(X, W)``: states ``X``, arc costs
``W(x, y) ∈ [0, ∞]`` (``∞`` = no arc), stored as `W[x, y]` exactly as in
`TropicalAbstraction`. Its path category, enriched in the Lawvere quantale
``([0, ∞], ≥, +, 0)`` (Lawvere 1973), has hom-objects
``d(x, y) = \min_{\text{paths } x \to y} \text{cost}``: the shortest-path
quasi-metric, which `minplus_closure` computes. A terminal cost ``b`` is an arc
``x → ⋆`` of cost ``b(x)`` into a sink ``⋆``, and the cost-to-go is
``V(x) = d(x, ⋆)``.

**Abstraction.** A surjection ``ψ : X → Y`` onto blocks. The optimistic quotient
(`quotient_graph`) is

```math
W_ψ(B, C) = \min_{x ∈ B,\ y ∈ C} W(x, y) \quad (B ≠ C), \qquad b_ψ(B) = \min_{x ∈ B} b(x),
```

the greatest arc weighting on ``Y`` for which ``ψ`` does not increase the cost of
any arc. Its cost-to-go ``V_ψ`` pulls back to a lower bound,
``V_ψ ∘ ψ ≤ V`` (`lower_bound`). The pessimistic quotient
``W^{\text{must}}(B, C) = \max_{x ∈ B} \min_{y ∈ C} W(x, y)`` (`robust_quotient`)
gives an upper bound.

**The lifting defect.** A coarse arc ``B → C`` *lifts* from ``x ∈ B`` if ``x``
can reach ``C`` at no more than the arc's cost. The defect measures by how much
it fails:

```math
ε_ψ(B, C) = \max_{x ∈ B} \big( d(x, C) − W_ψ(B, C) \big)^+, \qquad
ε_ψ(B, ⋆) = \max_{x ∈ B} \big( V(x) − b_ψ(B) \big)^+ ,
```

with ``d(x, C) = \min_{y ∈ C} d(x, y)``. The path from ``x`` may take any number
of steps, inside ``B`` or elsewhere.

**Theorem (path lifting; paper Theorem 4.2 and Corollary 4.3).**

1. *Sandwich.* For every ``x ∈ X``,
   ``V_ψ(ψx) ≤ V(x) ≤ V_ψ^{+ε}(ψx) ≤ V^{\text{must}}(ψx)``, where ``V_ψ^{+ε}`` is
   the cost-to-go of the quotient with arc costs ``W_ψ + ε_ψ`` and terminal cost
   ``b_ψ + ε_ψ(·, ⋆)``.
2. *Exactness.* The following are equivalent:
   - (i) ``V_ψ ∘ ψ = V`` for every target that is a union of blocks;
   - (ii) ``ε_ψ(B, C) = 0`` for every arc ``B → C`` of ``W_ψ`` between distinct
     blocks.
   
   When they hold, ``V_ψ ∘ ψ = V`` for every terminal cost constant on blocks.
3. *Locality.* Any upper bound on ``d(x, C)`` gives an upper bound on
   ``ε_ψ(B, C)``, and item 1 still holds with it. In particular, distances
   computed inside the subgraph on ``B ∪ C`` suffice.

The proof of item 1 lifts an optimal walk of ``W_ψ + ε`` arc by arc: from any
state of ``B_i``, some state of ``B_{i+1}`` is reachable within
``W_ψ(B_i, B_{i+1}) + ε_ψ(B_i, B_{i+1})``. For (i) ⇒ (ii), take the target ``C``:
then ``d(x, C) = V(x) = V_ψ(B) ≤ W_ψ(B, C)``.

**Why it matters.** The defect turns "this abstraction is too coarse" into a
number per coarse arc. It explains the studies:
- forgetting velocity has defects that grow with speed (a state moving away from
  ``C`` must brake first);
- level sets of a cost-to-go have defects equal to detours around walls;
- Reeb arcs, velocity modes and `future_refinement` all reduce the defect.

Boolean special cases are already in the package. `is_stable(W, ψ; arcs)` says
that every selected coarse arc lifts in one step from every state, so on
unit-cost or zero-cost arcs it implies zero defect. `stable_refinement` is the
coarsest such refinement (Paige & Tarjan 1987).

## Codebase Orientation

| File | Why it matters |
|---|---|
| `src/network_sheaves/TropicalAbstraction.jl` | `quotient_graph` (``W_ψ``), `robust_quotient` (``W^{\text{must}}``), `cost_to_go` and `minplus_closure`, `lower_bound`; the Dijkstra and heap helpers (`_astar`, `_DijkstraClosure`) to reuse for the per-arc searches. The new function goes here. |
| `src/network_sheaves/TropicalCovers.jl` | `is_stable`, `stable_refinement`, `future_refinement`: the Boolean lifting conditions the defect quantifies. Tests compare against them. |
| `src/network_sheaves/GraphHomomorphisms.jl` | `GraphHomomorphism` (`vertex_map`, `n_target`) for ``ψ``. |
| `test/network_sheaves/TropicalAbstraction.jl` | Where the tests go. |
| `docs/scripts/path_lifting_check.jl` | Dense reference implementation of the defect and of every claim of the theorem; port its checks into the tests. |
| `docs/src/features/tropical_abstraction.md` | Add a short "lifting defect" paragraph after the may/must table. |

## Requested Implementation

```julia
"""
    lifting_defect(W, ψ; scope = :pair) -> SparseMatrixCSC{Float64,Int}

The lifting defect of the abstraction `ψ` of the cost graph `W` (`W[x, y]` the cost
of the arc ``x → y``): for every arc ``B → C`` of the optimistic quotient
[`quotient_graph`](@ref)`(W, ψ)` between distinct blocks, the worst excess, over the
states ``x`` of ``B``, of the cost of reaching ``C`` from ``x`` over the cost of
the coarse arc,

```math
ε(B, C) = \\max_{x ∈ B} \\big( d(x, C) − W_ψ(B, C) \\big)^+ .
```

A coarse arc *lifts* when its defect is zero: every state of its source block can
follow it at no extra cost. If every arc lifts, the pulled-back coarse cost-to-go
is exact for every target that is a union of blocks. In general, adding the
defect to the quotient's arc costs gives an upper bound: the cost-to-go of
`quotient_graph(W, ψ) + ε` pulled back is at least the fine cost-to-go, and never
worse than the pessimistic bound of [`robust_quotient`](@ref).

`scope = :pair` measures ``d(x, C)`` with paths inside the blocks ``B ∪ C`` only.
This is an upper bound on the defect, cheap and local, and enough for the bound.
`scope = :global` uses paths through the whole graph and gives the exact defect.
An infinite entry means some state of ``B`` cannot reach ``C`` (within the scope).

The result has the sparsity pattern of the quotient's arcs between distinct
blocks; stored zeros are arcs that lift.
"""
function lifting_defect(W::SparseMatrixCSC, ψ::GraphHomomorphism; scope::Symbol = :pair)
```

**Validation** (`@argcheck`): `W` square with nonnegative costs (`_check_costs`),
`size(W, 1) == length(ψ.vertex_map)`, `scope in (:pair, :global)`.

**Algorithm.**
1. `Wq = quotient_graph(W, ψ)`.
2. Group the states by block from `ψ.vertex_map`.
3. For every stored arc `(B, C)` of `Wq`, run a reverse Dijkstra from the states of
   ``C`` (multi-source, distance 0), along in-arcs (columns of `W`):
   - for `:pair`, restricted to states in ``B ∪ C``;
   - for `:global`, unrestricted. Group arcs by ``C`` so that one reverse search
     serves every ``B``.
4. Set ``ε = \max_{x ∈ B} \max(\text{dist}(x) − Wq[B, C], 0)``.

Cost:
- `:pair`: the sum over coarse arcs of the Dijkstra cost on ``|B| + |C|`` states.
- `:global`: one reverse Dijkstra per block ``C`` that is the head of an arc.

Store the result with `sparse(I, J, V, nc, nc)`, keeping explicit zeros. Build it
from triplets and pass `combine = max` so zeros are not dropped.

## Tests to Write

In `test/network_sheaves/TropicalAbstraction.jl`, `@testset "lifting_defect"`:

```julia
# The line 1 - 2 - 3 - 4 (unit steps toward 4), blocks {1,2}, {3,4}
W = sparse([1, 2, 3], [2, 3, 4], [1.0, 1.0, 1.0], 4, 4)
ψ = GraphHomomorphism([1, 1, 2, 2], 2)
ε = lifting_defect(W, ψ)
@test ε[1, 2] == 1.0                              # state 1 needs 2 steps, the arc costs 1
@test lifting_defect(W, ψ; scope = :global)[1, 2] == 1.0

# In the gauge of the exact cost-to-go every arc toward the goal costs 0, and the defect vanishes
V = cost_to_go(W, [4])
Wg = sparse([1, 2, 3], [2, 3, 4], [1.0 + V[2] - V[1], 1.0 + V[3] - V[2], 1.0 + V[4] - V[3]], 4, 4)
@test all(iszero, nonzeros(lifting_defect(Wg, ψ)))

# Zero-cost cycles inside blocks make every arc lift; the bound is then exact
rng = MersenneTwister(19)
W, ψ = random_cost_graph_with_zero_cycles(rng)    # helper as in docs/scripts/path_lifting_check.jl
@test all(iszero, nonzeros(lifting_defect(W, ψ)))
for C in 1:ψ.n_target
    targets = findall(==(C), ψ.vertex_map)
    @test lower_bound(Abstraction(W, ψ), targets) ≈ cost_to_go(W, targets)
end

# Against the dense reference on random graphs: :global equals the reference, :pair bounds it
for _ in 1:50
    W, ψ = random_cost_graph(rng)
    εg, εp, εref = lifting_defect(W, ψ; scope = :global), lifting_defect(W, ψ), dense_defect(W, ψ)
    @test findnz(εg)[1:2] == findnz(quotient_graph(W, ψ))[1:2]   # the quotient's arcs (it drops arcs inside blocks)
    @test all(isapprox.(nonzeros(εg), nonzeros(εref)))
    @test all(nonzeros(εp) .>= nonzeros(εg) .- 1e-12)
end

# The sandwich: lower bound ≤ cost-to-go ≤ defect-inflated quotient ≤ pessimistic quotient
for _ in 1:50
    W, ψ = random_cost_graph(rng); targets = rand_targets(rng, ψ)       # targets: a union of blocks
    V = cost_to_go(W, targets)
    Wq = quotient_graph(W, ψ); ε = lifting_defect(W, ψ)
    tc = unique(ψ.vertex_map[targets])
    up = cost_to_go(Wq .+ ε, tc)[ψ.vertex_map]                        # same pattern, so .+ is safe
    @test all(lower_bound(Abstraction(W, ψ), targets) .<= V .+ 1e-9)
    @test all(V .<= up .+ 1e-9)
    must = cost_to_go(robust_quotient(W, ψ), tc)[ψ.vertex_map]
    @test all(up .<= must .+ 1e-9)
end

# Boolean special case: on a unit-cost graph a stable partition has zero defect
W = unit_cost_grid(6, 6); ψ = stable_refinement(W, coarse_blocks(6, 6, 2))
@test all(iszero, nonzeros(lifting_defect(W, ψ)))
@test is_stable(W, ψ)

# Forgetting velocity on a small double-integrator lattice: positive defects that grow with speed
W, ψ = small_double_integrator_with_projection()  # e.g. 8² positions × 3² velocities, ψ = projection_homomorphism
ε = lifting_defect(W, ψ)
@test maximum(nonzeros(ε)) > 0
```

Helpers (`random_cost_graph`, `dense_defect`, and so on) go at the top of the
testset, ported from `docs/scripts/path_lifting_check.jl`. Do not add test
dependencies.

## Verification Checklist

- [x] `lifting_defect` exported from `TropicalAbstraction`; the docstring renders on the TropicalAbstraction API page (`docs/src/api/grid_schwarz.md`) and its four cross-references resolve (full docs build, all literate examples executed). `stable_refinement` lives in `TropicalCovers`, so its reference is module-qualified.
- [x] The `:global` defect matches the dense reference on random graphs; `:pair` is never below it.
- [x] Zero defect ⇔ `lower_bound` exact on every single-block target (both directions tested).
- [x] The sandwich `lower_bound ≤ cost_to_go ≤ inflated ≤ robust` holds on random graphs.
- [x] Line example: defect 1 in the plain gauge, 0 in the gauge of ``V``.
- [x] `stable_refinement` on a unit-cost graph gives zero defect.
- [x] `docs/scripts/path_lifting_check.jl` still passes (`julia docs/scripts/path_lifting_check.jl`).
- [x] Full test suite on a compute node: 5818 passed, 3 failed. The `TropicalAbstraction` tests pass (837 of 837, with 1 and with 4 threads). The 3 failures are the `GridSchwarz boxes` MPI launches, whose subprocess cannot load MPI in the sandboxed test environment. This change does not touch that test or the code it runs. Two test-file name clashes that hid tests in the full run (`tensor_product`, `constant_sheaf`) are fixed.
- [x] A short paragraph and the theorem statement added to `docs/src/features/tropical_abstraction.md`.

## Out of Scope

- An `upper_bound` method for `Abstraction` built on the inflated quotient (a
  separate issue once `lifting_defect` exists).
- Refining ``ψ`` where the defect is large (counterexample-guided abstraction
  refinement).
- Defects of the space–time quotient (`steps > 1`): the theorem is stated for one
  step.
- The terminal defect ``ε_ψ(B, ⋆)``. It is zero for targets that are unions of
  blocks, the only case `lower_bound` supports.
- Measuring the defect on the 128² × 13² lattice and on the multi-agent summaries.
  That is a study script for the paper (`docs/scripts/`), not library code.
- Stochastic (real-semiring) aggregation, where the analogous quantity is a residual,
  not a defect.
