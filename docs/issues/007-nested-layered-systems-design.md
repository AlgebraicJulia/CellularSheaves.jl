# Nested Layered Systems: Hierarchical Sheaf Compilation with a Direct-Solve Baseline

**Status**: design spec. Supersedes the earlier `002_nested_layered_systems_dsl.md` draft, which
was written before several load-bearing design questions were resolved. Where this document
disagrees with that draft, this document wins; Section 6 records what changed and why.

**Scope**: this is the umbrella design + rationale. Individual implementable issues are split
out per Section 9 and follow the repo's issue convention (`CLAUDE.md`).

---

## 1. Motivation

`ControlSheaves.Layered` (`src/ControlSheaves/Layered.jl`) implements a flat two-layer escort
architecture: `RingSpec`/`SupportSpec` describe escort rings bridged by support pods, compiled
via one `pushforward_sheaf` call into a coarse solve, then lifted back to per-agent references.

We want to generalize from "flat" to a **tree of nested systems** of arbitrary, irregular depth,
specified through an `MLStyle`-based DSL. The payoff is not merely expressiveness — it is
**computational**. The hierarchy replaces one large harmonic-extension problem with a small
harmonic extension at the top plus a sequence of block-diagonal linear lifts.

That speedup is bought with an approximation, and quantifying it is a first-class requirement
(Section 3.4), not an afterthought.

---

## 2. Codebase Orientation

| File | Why it matters |
|---|---|
| `src/ControlSheaves/Layered.jl` | The flat 2-layer architecture being generalized. Note `build_layered_homomorphism` (line 279) currently folds each target **into its ring's fiber** (`v_map[spec.target_nodes[r_idx]] = r_idx`) — Section 4.2 changes this. `build_layered_fiber_bases` (line 181) computes `pinv` of a target sub-basis; Section 4.2 removes the need for it. |
| `src/network_sheaves/Formations.jl` | `build_escort_ring` is already `D`/`affine`-parameterized. `build_escort_clique` (line 134) **already exists but is `D=4`-hardcoded with no `affine` flag** — so Issue 1 is "unify and generalize", not "purely additive". `affine_translation_matrix` generalizes SE(3) offsets to any dimension. |
| `src/network_sheaves/Pushforwards.jl` | `pushforward_sheaf`, `all_fiber_bases`, `fiber_section_basis`. Coarse stalk dimension = **dimension of the fiber's exact global-section space** (via `nullspace_ldlt`). This fact drives Sections 3.2–3.4. Never chained today; this feature chains it. |
| `src/network_sheaves/EuclideanSheaves.jl` | `EuclideanSheaf` (`vertex_stalks::Vector{Int}` is already per-vertex, and restriction maps need not be square). `_harmonic_extension_restricted_laplacian` (line 626) takes `boundary::Dict{Int,<:AbstractVector}` — it pins **whole vertices only**, which is why targets must be vertices (Section 4.2). |
| `src/network_sheaves/GraphHomomorphisms.jl` | `GraphHomomorphism`, `fiber_vertices`, `compose`. |
| `src/ControlSheaves/AgentControllers.jl` | Per-dynamics dispatch (`position_indices`, `state_dim`, `initial_state`), `AgentState`/`step_agent!`. Leaf agents reuse this unchanged. |
| `src/ControlSheaves/TrackingDSL/` | The 5-stage DSL template (ADT → Parser → Validator → Resolver → Lowering) this DSL mirrors. |
| `src/network_sheaves/{ADT,Parser}.jl` | The *other*, deliberately fused DSL (`@cellular_sheaf`). Read for contrast; not the template here. |
| `test/network_sheaves/LayeredEscort.jl` | Existing test patterns. **Note lines 61–69**: it computes both the hierarchical and direct solutions but asserts only their `size` — the two paths have **never been numerically compared**. Section 8 makes that comparison a real test. |

### 2.1 Blocking defect (fix before starting)

The package does not precompile on `jpf/multilayer-control`. `src/IPM/kkt/uzawa.jl:25` imports
`FactorizationWorkspace` / `DivisionWorkspace` from `CliqueTrees.Multifrontal`, which the
installed `CliqueTrees` no longer defines. Unrelated to this feature, but no test can run until
it is resolved (either re-pin `CliqueTrees` in `Project.toml` or update the import to the
current API).

---

## 3. Mathematical Background

### 3.1 The tower

A nested specification compiles to a tower of sheaves connected by graph homomorphisms:

```
H₀  (coarsest)      ← the only harmonic extension happens here
 ↑ f₁
H₁                  ← q₁ = B_v · q₀[v]   (block-diagonal lift, no solve)
 ↑ f₂
H₂                  ← q₂ = B_v · q₁[v]
 ⋮
H_N (raw agents)    ← per-agent LQR tracks the lifted reference
```

Each `H_{k-1}` is the **pushforward** of `H_k` along `f_k`. Chaining pushforwards is new to this
codebase; nothing today pushes a pushforward forward.

### 3.2 What the hierarchy actually computes

`pushforward_sheaf` sets each coarse vertex's stalk to the dimension of its fiber's **exact
global-section space**. So the columns of each fiber basis `B_v` span precisely the
configurations in which that fiber is internally consistent — for an escort ring, exactly the
rigid-formation configurations.

The lift `x = B·q` therefore has `B` block-diagonal over fibers, and every configuration the
tower can produce satisfies all fiber-internal constraints **exactly**. Internal Dirichlet
energy is identically zero by construction.

Stated as an optimization problem, the hierarchical solve is:

> Constrain every team to lie on its space of global sections (i.e. *in formation*).
> Subject to that, place the formations so as to minimize total sheaf energy.

The direct solve on the fully expanded sheaf `H_N` instead minimizes
`‖δx‖² = Σ_internal + Σ_cross` over **all** configurations.

### 3.3 The approximation, precisely

The hierarchical solution is a **feasible point** of the direct problem — it satisfies every
constraint the direct problem imposes, plus the extra constraint that each fiber be internally
exact. Therefore:

```
E_hierarchical  ≥  E_direct
```

with equality **iff** the direct optimum is already fiberwise-exact. Rigid escort rings around a
single target are the equality case, which is why the current flat code is exact. As soon as a
team observes multiple targets that pull apart, the direct optimum wants the team to *deform*,
the tower cannot represent that, and the gap becomes strictly positive.

The upside is computational: one small harmonic extension plus block-diagonal lifts, instead of
one large harmonic extension. The downside is that formations translate as rigid wholes rather
than deforming. This trade is the entire point of the architecture, so **the gap must be
measurable** — hence the dual compiler (Section 4.5).

### 3.4 Support-pod isolation invariant

Positions flow strictly root → leaves. Targets are boundary conditions at the top; everything
else is either an interior unknown of the single top-level solve or a lift of it. A support pod
cannot pull a team out of position, because it has no path by which to do so.

---

## 4. Design Specification

### 4.1 Tree structure

- A **system** is a graph. Vertices are either **leaf** teams (raw agents in a ring/star/path/
  clique formation) or **refined** into a recursive sub-specification.
- Depth is arbitrary and **irregular**: siblings may differ in depth, topology, and arity.
- Four primitive topologies: **ring, star, path, clique**, each with an `observers` field
  carrying today's `RingSpec.observers` semantics.

### 4.2 Targets are top-level vertices

**Decision.** Targets are modeled as agents over which we have no control authority: they have
autonomous dynamics and act purely as boundary conditions supplied at solve time. They are
**their own vertices at the top level `H₀`**, never folded into a team's fiber.

This is a change from `build_layered_homomorphism` (`Layered.jl:284`), which currently maps each
target into its ring's fiber. Three consequences, all improvements:

1. **No pseudo-inverse.** Pinning a target becomes a genuine vertex pin, which is exactly what
   `_harmonic_extension_restricted_laplacian` accepts. `world_to_pf_stalk` /
   `target_subbases_inv` (and their hidden rigidity assumption) disappear entirely.
2. **Targets are not owned by teams.** A team may observe several targets; a target may be
   observed by several teams. Arbitrary team↔target incidence is expressible.
3. **`D` need not be uniform.** A vertex's stalk is whatever its fiber's section dimension turns
   out to be; only edge stalks must agree at their endpoints.

This **supersedes** the earlier draft's "exactly one target per vertex" rule and the
uniform-`D`-per-branch constraint.

### 4.3 `R` is a per-edge restriction map

The earlier draft made `R` a single per-vertex map defining "the value a refined vertex presents
to its parent". With targets as top-level vertices and multi-target teams allowed, that
collapses into something simpler and more standard: **`R` is just the restriction map on an
edge**, from a team's joint state to that edge's stalk.

`R` is declared symbolically against the subsystem's **raw joint state** (dimension = sum of its
own vertices' stalks — known structurally, no rank computation), and composed with the fiber
basis `B` at lowering time to give the actual coarse-level restriction map `R·B`. Declaration
stays symbolic; only lowering is numeric.

Constructors:

- `project(i)` / `project(:name)` — selection matrix picking member `i`'s block.
- `centroid()` — averaging matrix over the subsystem's **direct** vertices (opaque to their
  internal composition).
- A raw `Matrix{Float64}` escape hatch.

`EuclideanSheaf` already permits non-square restriction maps, so no changes are needed there.

### 4.4 Support pods

Support pods are ordinary unpinned vertices at whatever level they appear. Their harmonic
character emerges from the single top-level energy minimization rather than being constructed
per level.

Because pods are just unpinned vertices, a pod incident to `k > 2` systems needs no special
handling — the earlier draft's open question about "star-shaped" pods dissolves.

### 4.5 Dual compiler: hierarchical and direct

**Requirement.** One specification compiles two ways:

- `solve_hierarchical` — build the tower, one harmonic extension at `H₀`, lifts down.
- `solve_direct` — build the fully expanded sheaf `H_N` and solve the harmonic extension
  directly. This is the **baseline**.

The direct path is not merely a test fixture; it is how the approximation of Section 3.3 is
measured. Both must be reachable from the same spec so the comparison is apples-to-apples.

### 4.6 DSL surface syntax

- Named component systems referenced by variable-like names.
- Dotted/indexed selection: `systemA.agent[2]`, `systemA.agents[1:2:end]`. **Selection only** —
  no general cross-tree equation sublanguage.
- `end` resolves **locally** to a graph's own declared arity; a refined vertex counts as one
  unit for its parent's indexing. This keeps validation purely symbolic.

### 4.7 Dynamics and initial conditions

Topology and dynamics assignment stay orthogonal. The DSL produces a structural tree; a separate
resolve step binds `AbstractAgentDynamics`, LQR gains, and initial positions to leaf agents.

**Decision — cascading defaults, most specific wins.** Dynamics may be declared at any node of
the tree and are inherited by every leaf beneath it. Where both an inherited value and a local
one apply, **the local one wins**; precedence runs

```
per-agent  >  leaf team  >  nearer ancestor  >  … >  root default
```

**Decision — the context is a nested structure, not a flat dict.** Binding data mirrors the
system tree rather than using flat dotted keys like `"systemA.child1.agent[2]"`. It is built from
custom types with keyword-argument constructors, so inputs are validated at construction and
typos in child names are caught eagerly rather than silently binding nothing.

Target trajectories are the exception: they are inherently per-target and do **not** cascade.
They are supplied separately at solve time, as the boundary values of Section 3.4.

---

## 5. DSL Architecture

Mirrors `src/ControlSheaves/TrackingDSL/` module-for-module, as `NestedLayeredDSL` under
`src/ControlSheaves/`:

- **`ADT.jl`** — immutable AST node types + a `NestedLayeredDSLError` exception hierarchy.
- **`Parser.jl`** — `@nested_system begin ... end` macro, MLStyle `@match` per statement shape,
  plus a functional `parse_nested_system(block::Expr)` entry point.
- **`Validator.jl`** — purely symbolic checks: duplicate names, unresolved references, arity and
  index-range checks. No numeric work.
- **`Resolver.jl`** — binds dynamics / initial positions / trajectories from a late-supplied `ctx`.
- **`Lowering.jl`** — builds the tower (and the direct baseline) from the resolved tree.
- **`NestedLayeredDSL.jl`** — thin aggregator with `@reexport`.

Staged rather than fused (`@cellular_sheaf`'s model) because the tree has real depth and
cross-level name references needing a symbol table resolved in dependency order.

---

## 6. What Changed From the Earlier Draft

| Earlier draft | Now | Why |
|---|---|---|
| Harmonic solve at **every** level (§4.4) | **One** solve at `H₀`, pure lifts below | Owner's design intent; simpler and faster. Reduces "recursive solve/lift cascade" to a fold over fiber bases. |
| `R` = one per-vertex outward value | `R` = per-edge restriction map | Falls out of targets-as-vertices; more general, less machinery. |
| Exactly one target per vertex | Arbitrary team↔target incidence | Owner requires shared/multi-target tracking. |
| Uniform `D` per branch | Per-vertex stalks, `D` only on edges | Enabled by targets-as-vertices; also dissolves the old §7.1. |
| Support pods = rings around a synthetic harmonic midpoint | Ordinary unpinned vertices | The midpoint property emerges from the top-level solve; no synthetic vertex needed. Also dissolves the old §7.2 (`k>2`). |
| "Don't route through `fiber_section_basis`" (§3.2/§6.4) | Fiber bases **are** the lift mechanism | The draft's *warning* was right (never assume a fiber's section space is `D`-dimensional); its *prescription* was wrong. `R` handles the dimension mismatch. |
| Direct solve as a test cross-check | Direct solve as a **first-class compiler target** | It is the baseline against which the hierarchy's approximation is measured. |

Still rejected, unchanged: a full cross-tree equation sublanguage; adding Catlab.jl or
AlgebraicDynamics.jl (borrow the UWD *shape*, not the dependency).

---

## 7. Resolved Questions

Both questions that previously blocked the DSL resolver are now settled; see Section 4.7 for the
normative statement.

**7.1 Dynamics-assignment cascade — RESOLVED.** Dynamics cascade down the tree, and the most
specific specification wins: per-agent overrides its leaf team, which overrides any ancestor
default. Target trajectories do not cascade.

**7.2 `ctx` ergonomics — RESOLVED.** The binding context is a **nested structure mirroring the
system tree**, built from custom types with keyword-argument constructors, chosen over flat
dotted-string keys specifically so that inputs can be validated at construction time.

No open questions remain. Issues 5–6 are unblocked.

---

## 8. Tests to Write

- **Energy-gap invariant** (correctness-critical): for any spec, assert
  `E_hierarchical ≥ E_direct - tol`. This is the theorem of Section 3.3 and must hold for every
  test case.
- **Equality in the rigid case**: for single-target rigid escort rings, assert
  `q_hierarchical ≈ q_direct`. This closes the gap left by `LayeredEscort.jl:61–69`, which
  computes both and compares neither.
- **Strict gap in the non-rigid case**: a team observing two separating targets must show
  `E_hierarchical > E_direct`, demonstrating the rigidity restriction is real and measured.
- **Golden regression**: a 2-level tree with no refinement must reproduce today's flat
  `LayeredEscortSpec` pipeline numerically.
- **`project`/`centroid` unit tests**, including that `centroid()` treats a refined member as one
  opaque unit rather than expanding it.
- **Shared-target test**: one target observed by two teams; both teams' references respond to it.
- **Irregular tree test**: only some vertices refined at a given depth, siblings differing in
  topology, arity, and depth.
- **DSL stage tests**: parse (AST shape), validate (`@test_throws` for duplicates/unbound
  references/arity), resolve (`ctx` binding), lower (tower shapes).

---

## 9. Issue Breakdown

The chain is **strictly sequential** — each issue's tests exercise the API the previous one
introduces. Do not parallelize.

| # | File | Summary | Prereq |
|---|---|---|---|
| 0 | — | Fix the `CliqueTrees` precompile break (§2.1). **Done.** | — |
| 1 | `008-unified-formation-topologies.md` | One `build_escort_topology(kind, ...)` covering ring/star/path/clique, all `D`/`affine`-parameterized; folds in the `D=4`-only `build_escort_clique`. | 0 |
| 2 | `009-nested-spec-and-tower-compiler.md` | Spec types + tower compiler: targets as singleton fibres, chained pushforwards building `H₀…H_N`. | 1 |
| 3 | `010-dual-compiler-and-energy-gap.md` | Hierarchical solve/lift, direct baseline, and the `E_hier ≥ E_direct` gap diagnostics of §3.3. | 2 |
| 4 | `011-per-edge-restriction-maps.md` | Per-edge `R` with `project`/`centroid`, replacing Issue 2's `TODO(011)` placeholder. | 2, 3 |
| 4b | `012-nested-dynamics-context-and-cascade.md` | Nested dynamics context with most-specific-wins cascade (§4.7). **Independent of 3 and 4** — may run in parallel once 2 lands. | 2 |
| 4c | `013-retire-dynamics-trichotomy.md` | Delete `HomogeneousDynamics`/`TeamHomogeneousDynamics`/`IndividualizedDynamics`, porting the flat `Layered` path onto `SystemBinding`. | 4b |
| — | `014-align-docs-module-coverage.md` | Add the missing `Formations` API page and align `makedocs(modules=...)` with the modules that actually have `@autodocs` blocks. Independent of this chain, but **touches `docs/make.jl`**, so do not run it concurrently with an issue that also edits that file. | none |
| 5 | — | DSL `ADT.jl` + `Parser.jl`, targeting the API from 2–4b. **Done** (`src/ControlSheaves/NestedDSL/`). | 4, 4b |
| 6 | — | DSL `Validator.jl` + `Lowering.jl`. **Done** (`src/ControlSheaves/NestedDSL/`). | 5 |

Issues 5–6 were deliberately left unwritten until the compilation target stabilized; they have
now landed together as the `NestedDSL` module. Two decisions there depart from §5's sketch, both
because the problems being specified turned out to be combinatorial rather than declarative:

- **No `Resolver.jl`.** §5 planned a resolve stage binding values from a late-supplied `ctx`,
  mirroring `TrackingDSL`. `NestedDSL` instead *executes* its block, so every numeric value is
  evaluated by ordinary Julia in the caller's own scope at the point the declaration runs and
  there is nothing left to resolve. `@bind` still expresses §4.7's cascade in full; what
  disappeared is the dict of names to look values up in.
- **Fragments are first-class.** `@nested_system` returns a `SystemFragment` describing one node,
  with relative paths that lowering rewrites against wherever it lands. Composition —
  `merge`, `@include`, `@system name = fragment` — replaces the DSL-level iteration and
  abstraction §4.6 would otherwise have needed, and lets specification text interleave freely
  with the Julia that computes what it refers to.

---

## 10. Future Work: Soft (Spectral) Fiber Sections

Recorded so it is not lost. The rigidity limitation of Section 3.3 is entirely a consequence of
`B_v` spanning only the fiber Laplacian's **nullspace** (`λ = 0`, exact sections).

**Idea**: span instead the eigenvectors of the `m` smallest eigenvalues of the fiber Laplacian.
Teams then deform along their *softest* modes, `m` becomes a per-fiber compliance budget, and the
one-solve tower survives unchanged — only the basis-construction function changes. Setting
`m = dim ker` recovers exactly the behavior specified in this document, so it is a strict
generalization with a continuous knob between "rigid formations, cheap" and "fully deformable,
expensive".

`pushforward_sheaf`'s exact section-preservation weakens to an approximation, which is the point:
it trades a controlled amount of the Section 3.3 gap for a slightly larger top-level solve.
`Krylov` is already a dependency, so the eigensolve needs no new one.

**Deferred by owner decision** — implement Sections 1–9 first.
