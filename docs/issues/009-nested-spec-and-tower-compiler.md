# Issue 009 — Nested System Specification and Sheaf Tower Compiler

**Parent design**: `007-nested-layered-systems-design.md` (Issue 2)
**Prerequisites**: 008
**Blocks**: 010, 011

---

## Mathematical Background

A nested specification compiles to a **tower** of sheaves linked by graph homomorphisms:

```
H₀  (coarsest)
 ↑ f₁
H₁
 ↑ f₂
 ⋮
H_N (finest — raw agents and targets)
```

where each `H_{k-1} = (f_k)_* H_k` is the pushforward of the level below. `pushforward_sheaf`
sets each coarse stalk to the dimension of that fibre's **exact global-section space**, so a
rigid escort team (Issue 008) collapses to a single `D`-dimensional vertex representing its
centre.

### Targets are singleton fibres

Design doc §4.2 requires targets to remain top-level vertices — they are agents we have no
control authority over, they are not owned by any team, and they may be observed by several
teams at once.

The tower realizes this cleanly: **every homomorphism `f_k` maps each target vertex to its own
vertex at the coarser level**, i.e. targets are singleton fibres at every level. A singleton
fibre has no internal edges, so its coboundary is zero, `fiber_section_basis` returns an
identity basis, and the target's stalk stays `D` all the way to `H₀` with the identity as its
lift. No pseudo-inverse is needed anywhere, and `_harmonic_extension_restricted_laplacian` — which
pins whole vertices only — can pin targets directly at `H₀`.

This is what removes `world_to_pf_stalk` / `target_subbases_inv` (`Layered.jl:181–216`) and the
hidden rigidity assumption they carry.

---

## Codebase Orientation

| File | Why |
|---|---|
| `src/network_sheaves/Pushforwards.jl` | `pushforward_sheaf(hom, F)` (line 125), `all_fiber_bases` (line 86), `fiber_section_basis` (line 52). Chained here for the first time in this repo. |
| `src/network_sheaves/GraphHomomorphisms.jl` | `GraphHomomorphism(vertex_map, n_target)` (line 43), `fiber_vertices`. |
| `src/network_sheaves/EuclideanSheaves.jl` | `EuclideanSheaf{T}(vertex_stalks)`, `add_sheaf_edge!` (line 193 — asserts restriction-map column counts match endpoint stalks). |
| `src/ControlSheaves/Layered.jl` | `LayeredEscortSpec` (line 120) and `build_layered_homomorphism` (line 279) are the flat precedent. Note line 284 folds targets **into** ring fibres — this issue does the opposite. |
| `src/network_sheaves/Formations.jl` | `build_escort_topology` from Issue 008. |

---

## Requested Implementation

New module `src/ControlSheaves/NestedSystems.jl`, included from `ControlSheaves.jl`.

### Specification types

```julia
abstract type AbstractSystemNode end

"""
    LeafTeam(name, kind, n_agents, radius; observers=[1])

A team of raw agents in a `kind` formation (see `build_escort_topology`). Leaves are where
actual agents live; every other node in the tree is structural.
"""
struct LeafTeam <: AbstractSystemNode
    name::Symbol
    kind::Symbol
    n_agents::Int
    radius::Float64
    observers::Vector{Int}
end

"""
    RefinedSystem(name, children, internal_edges)

A subsystem whose vertices are themselves systems. `internal_edges` are consensus edges among
`children`, given as `(i, j)` index pairs into `children`.
"""
struct RefinedSystem <: AbstractSystemNode
    name::Symbol
    children::Vector{AbstractSystemNode}
    internal_edges::Vector{Tuple{Int,Int}}
end

"""
    TargetSpec(name)

An uncontrolled agent supplying a boundary condition at solve time. Targets are top-level
vertices and are never owned by a team.
"""
struct TargetSpec
    name::Symbol
end

"""
    Observation(system_path, target_index)

Declares that the system at `system_path` (a path of child indices from the root) observes
target `target_index`. Arbitrary many-to-many incidence is allowed: one system may observe
several targets and one target may be observed by several systems.
"""
struct Observation
    system_path::Vector{Int}
    target_index::Int
end

struct NestedSystemSpec
    root::RefinedSystem
    targets::Vector{TargetSpec}
    observations::Vector{Observation}
    D::Int
    affine::Bool
end
```

### Tower type and compiler

```julia
"""
    SheafTower

Compiled tower of sheaves. `levels[1]` is the coarsest sheaf `H₀`; `levels[end]` is the finest
`H_N` containing one vertex per raw agent plus one per target. `homs[k] : levels[k+1] → levels[k]`,
and `bases[k][v]` is the fibre-section basis over vertex `v` of `levels[k]`, used to lift a value
at `v` down into `levels[k+1]`.

`target_vertices[t]` is target `t`'s vertex index, which by construction is valid at **every**
level (targets are singleton fibres). `agent_vertices[a]` is agent `a`'s vertex index in
`levels[end]`.
"""
struct SheafTower
    spec::NestedSystemSpec
    levels::Vector{EuclideanSheaf{Float64}}
    homs::Vector{GraphHomomorphism}
    bases::Vector{Vector{Matrix{Float64}}}
    target_vertices::Vector{Int}
    agent_vertices::Vector{Int}
    depth::Int
end

"""
    build_sheaf_tower(spec::NestedSystemSpec) -> SheafTower

Compile a nested specification into a tower of sheaves connected by pushforwards.
"""
function build_sheaf_tower(spec::NestedSystemSpec)
```

### Algorithm sketch

1. **Flatten to the finest level.** Walk the tree depth-first, assigning consecutive vertex
   indices to every leaf agent. Append one vertex per target after all agents. Record each
   subtree's vertex range and each leaf's agent range.
2. **Build `H_N`.** For each `LeafTeam`, call `build_escort_topology` and transplant its edges
   into `H_N` at the leaf's vertex offsets (mirroring `Layered.jl:238–243`). For each
   `RefinedSystem`'s `internal_edges`, add a consensus edge between representative vertices of
   the two children — see "Open detail" below. For each `Observation`, add an edge from the
   observing system to its target vertex.
3. **Build each `f_k` bottom-up.** At depth `k`, every subtree rooted at depth `k` collapses to
   one coarse vertex; every target maps to its own singleton fibre. Construct the
   `vertex_map::Vector{Int}` and `GraphHomomorphism(vertex_map, n_coarse)`.
4. **Push forward.** `levels[k] = pushforward_sheaf(homs[k], levels[k+1])`, and
   `bases[k] = all_fiber_bases(homs[k], levels[k+1])`. Cache the bases rather than recomputing.
5. **Assert the target invariant.** After each pushforward, check every target vertex still has
   stalk `D` and an identity basis. Fail loudly with a clear message if not — this is the
   invariant the whole architecture rests on.
6. **Assert non-degenerate fibres.** If any fibre's section space has dimension 0, throw: the
   team is over-constrained and the tower cannot represent it (design doc §3.2).

### Open detail for the implementer

Step 2's "representative vertex" for an internal edge between two child systems is the one place
this issue touches Issue 011's territory. For **this** issue, use the child's first agent as its
representative (equivalent to `project(1)`). Issue 011 replaces that with a declared per-edge
restriction map. Leave a clearly-marked `TODO(011)` at that call site rather than inventing a
mechanism here.

---

## Tests to Write

New file `test/ControlSheaves/NestedSystems.jl`, included from `test/runtests.jl`:

```julia
@testset "tower — depth and level count" begin
    spec = two_level_spec()          # helper: 2 teams, 1 target each
    tower = build_sheaf_tower(spec)
    @test tower.depth == 2
    @test length(tower.levels) == 2
    @test length(tower.homs) == 1
    @test length(tower.bases) == 1
end

@testset "tower — targets are singleton fibres with identity bases at every level" begin
    tower = build_sheaf_tower(three_level_spec())
    for k in 1:length(tower.homs), t in tower.target_vertices
        @test length(fiber_vertices(tower.homs[k], t)) == 1
        @test tower.bases[k][t] ≈ I(tower.spec.D)
    end
end

@testset "tower — a rigid team collapses to exactly D dimensions" begin
    tower = build_sheaf_tower(two_level_spec())
    team_vertex = 1
    @test vertex_stalks(tower.levels[1])[team_vertex] == tower.spec.D
end

@testset "tower — one target observed by two teams" begin
    tower = build_sheaf_tower(shared_target_spec())
    t = tower.target_vertices[1]
    @test degree(tower.levels[1].underlying_graph, t) == 2
end

@testset "tower — irregular depth across siblings" begin
    # left child refined twice, right child a bare leaf
    tower = build_sheaf_tower(irregular_spec())
    @test tower.depth == 3
    @test length(tower.levels) == 3
end

@testset "tower — over-constrained fibre is rejected" begin
    @test_throws Exception build_sheaf_tower(degenerate_spec())
end
```

---

## Verification Checklist

- [ ] `NestedSystems` module created, included from `ControlSheaves.jl`, exports added.
- [ ] Targets verified as singleton fibres with identity bases at every level.
- [ ] No `pinv` anywhere in the new code path.
- [ ] Fibre bases computed once and cached in the tower, not recomputed per solve.
- [ ] Matrices kept sparse where the existing code does; no `Matrix(A)` in assembly loops.
- [ ] `@argcheck` at the public boundary (`build_sheaf_tower`, spec constructors).
- [ ] `TODO(011)` marker at the representative-vertex call site.
- [ ] `docs/make.jl` and `docs/src/api/` updated for the new module, per `CLAUDE.md`.
- [ ] Test helpers (`two_level_spec` etc.) defined once and shared with Issue 010's tests.

---

## Out of Scope

- Solving anything — Issue 010.
- Declared per-edge restriction maps (`project`/`centroid`) — Issue 011.
- Any DSL surface syntax — Issues 012+.
- Dynamics and initial-condition binding.
