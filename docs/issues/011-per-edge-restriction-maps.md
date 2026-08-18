# Issue 011 — Per-Edge Restriction Maps (`project` / `centroid`)

**Parent design**: `007-nested-layered-systems-design.md` (Issue 4)
**Prerequisites**: 009, 010
**Blocks**: DSL issues (012+)

---

## Mathematical Background

When a refined subsystem is wired to something else — a sibling system, a support pod, or a
target — that edge needs a restriction map from the subsystem's state to the edge's stalk.

Earlier design drafts made this a single per-vertex map defining "the value a refined vertex
presents to its parent". With targets promoted to top-level vertices and teams permitted to
observe several targets (design doc §4.2), that idea collapses into ordinary sheaf structure:
**it is just the restriction map on an edge**, and a subsystem may carry a different one on each
of its edges.

### Declaration is symbolic, lowering is numeric

A restriction map is declared against the subsystem's **raw joint state** — the plain
concatenation of its own vertices' values, of dimension `Σ vertex_stalks(children)`. That
dimension is known **structurally**, from declared arities alone, with no rank or nullspace
computation. This keeps validation purely symbolic, which is what lets the DSL's validator stage
stay numeric-free.

At lowering time the declared map `R : R^total → R^D` is composed with the fibre basis to give
the actual coarse-level restriction map:

```
R_coarse = R · B_v        (D × total) · (total × k)  =  D × k
```

where `k = size(B_v, 2)` is the fibre's section dimension, discovered by the pushforward.
`EuclideanSheaf` permits non-square restriction maps (`add_sheaf_edge!` only requires the column
count to match the endpoint stalk), so `R_coarse` needs no special handling.

### Why `centroid` is over direct vertices only

`centroid()` averages a subsystem's **direct** children, not its recursively flattened leaves. A
refined child counts as exactly one unit regardless of how many agents it eventually expands to.
This keeps refinement opaque to the parent — the same principle that makes `end` resolve locally
in the DSL — and keeps the map's dimensions computable from one level of declared arity.

---

## Codebase Orientation

| File | Why |
|---|---|
| `src/ControlSheaves/NestedSystems.jl` | Home of the change. Issue 009 left a `TODO(011)` at the representative-vertex call site — that is exactly what this issue replaces. |
| `src/network_sheaves/EuclideanSheaves.jl:193` | `add_sheaf_edge!` asserts `size(rm, 2) == vertex_stalks[v]`; `R_coarse` must satisfy this. |
| `src/network_sheaves/Pushforwards.jl:52` | `fiber_section_basis` documents its **row ordering** — rows follow the concatenation of `verts`' stalks in `verts` order. `R` is declared against that same ordering, so composition is valid only if the implementer matches it exactly. |
| `src/ControlSheaves/Layered.jl:181` | `build_layered_fiber_bases`'s `pinv`-based `target_subbases` is the mechanism this replaces. Do not reintroduce it. |

---

## Requested Implementation

```julia
"""
    RestrictionSpec

Declarative description of a restriction map from a subsystem's raw joint state to an edge
stalk. Materialized into a matrix at lowering time by [`materialize_restriction`](@ref).
"""
abstract type RestrictionSpec end

"""
    project(i::Int)  -> RestrictionSpec
    project(name::Symbol) -> RestrictionSpec

Select direct child `i` (or the child named `name`) and expose its own block unchanged. The
materialized matrix is the selection matrix `[0 … I_D … 0]`.
"""
project(i::Union{Int,Symbol}) = ProjectMember(i)

"""
    centroid() -> RestrictionSpec

Expose the unweighted average of the subsystem's **direct** children, `(1/N) Σ`. A refined child
counts as one unit; the map is opaque to its internal composition.
"""
centroid() = Centroid()

"""
    materialize_restriction(r::RestrictionSpec, node::RefinedSystem, D::Int) -> Matrix{Float64}

Build the `D × total_dim` matrix for `r`, where `total_dim = Σ vertex_stalks(node.children)`.
Purely structural — depends only on declared arities, never on a rank or nullspace computation.
"""
function materialize_restriction(r::RestrictionSpec, node::RefinedSystem, D::Int)
```

Concrete types: `ProjectMember(member::Union{Int,Symbol})`, `Centroid()`, and
`RawRestriction(M::Matrix{Float64})` as the escape hatch, constructed directly from a
`Matrix{Float64}`.

### Wiring into the spec

Extend Issue 009's edge and observation declarations to carry an optional `RestrictionSpec`,
defaulting to `project(1)` so Issue 009's behaviour is preserved exactly:

```julia
struct SystemEdge
    src::Int
    dst::Int
    src_map::RestrictionSpec   # default project(1)
    dst_map::RestrictionSpec   # default project(1)
end
```

`Observation` gains a `system_map::RestrictionSpec` with the same default.

### Algorithm sketch

1. `materialize_restriction(::ProjectMember, node, D)` — compute the child's row offset as
   `sum(child_stalk_dims[1:i-1])`, then place `I_D` at those columns. Resolve a `Symbol` member
   against `child.name` and throw a clear error if unmatched.
2. `materialize_restriction(::Centroid, node, D)` — place `(1/N)·I_D` at every child's block.
3. `materialize_restriction(::RawRestriction, node, D)` — `@argcheck size(M) == (D, total_dim)`
   and return `M`.
4. At tower assembly, for an edge incident to refined vertex `v` at level `k`, compute
   `R_coarse = materialize_restriction(spec_map, node, D) * tower.bases[k][v]` and pass it to
   `add_sheaf_edge!`.
5. **Row-ordering check.** Assert that the child ordering used to build `R` matches
   `fiber_vertices(homs[k], v)`'s ordering, which is what `fiber_section_basis` used for `B_v`'s
   rows. A mismatch here silently produces a wrong answer rather than an error, so make it an
   explicit assertion, not a comment.

---

## Tests to Write

```julia
@testset "project(i) selects exactly child i's block" begin
    node = three_child_system()      # children of stalk D each
    R = materialize_restriction(project(2), node, 4)
    @test size(R) == (4, 12)
    @test R[:, 5:8] ≈ I(4)
    @test all(iszero, R[:, 1:4]) && all(iszero, R[:, 9:12])
end

@testset "project(:name) resolves by child name" begin
    node = three_child_system()
    @test materialize_restriction(project(:bravo), node, 4) ≈
          materialize_restriction(project(2), node, 4)
    @test_throws Exception materialize_restriction(project(:nonexistent), node, 4)
end

@testset "centroid averages direct children" begin
    node = three_child_system()
    R = materialize_restriction(centroid(), node, 4)
    @test R[:, 1:4] ≈ I(4) / 3
    @test R * repeat([1.0, 2.0, 3.0, 1.0], 3) ≈ [1.0, 2.0, 3.0, 1.0]
end

@testset "centroid treats a refined child as one opaque unit" begin
    # `wide` refines into 6 agents, `narrow` into 2; centroid must still weight them 1/2 each.
    node = mixed_arity_system()
    R = materialize_restriction(centroid(), node, 4)
    @test R[:, 1:4] ≈ I(4) / 2
end

@testset "raw matrix escape hatch validates its shape" begin
    node = three_child_system()
    @test_throws Exception materialize_restriction(RawRestriction(zeros(4, 5)), node, 4)
end

@testset "default project(1) reproduces Issue 009 behaviour" begin
    @test solve_hierarchical(build_sheaf_tower(spec_explicit_project1()), tv)[end] ≈
          solve_hierarchical(build_sheaf_tower(spec_default()), tv)[end]
end

@testset "centroid-wired tower still satisfies the energy-gap theorem" begin
    tower = build_sheaf_tower(centroid_wired_spec())
    @test approximation_gap(tower, default_targets(tower.spec)).gap >= -1e-8
end
```

---

## Verification Checklist

- [ ] `project`, `centroid`, `RestrictionSpec` exported; raw `Matrix{Float64}` accepted.
- [ ] `materialize_restriction` performs **no** rank, nullspace, or `pinv` computation.
- [ ] Row ordering asserted against `fiber_vertices`, not assumed.
- [ ] Default `project(1)` leaves Issue 009's results bit-for-bit unchanged.
- [ ] The `TODO(011)` marker from Issue 009 is removed.
- [ ] Energy-gap theorem still holds for centroid-wired towers.
- [ ] `@argcheck` on all public constructors; clear errors for unresolved symbol members.
- [ ] Docstrings explain why declaration is symbolic and composition happens at lowering.

---

## Out of Scope

- Any DSL syntax for declaring restriction maps — Issues 012+. This issue provides the Julia API
  those will target.
- Dimension-changing restriction maps between branches with different `D`. Deferred until a real
  use case appears.
- Soft/spectral fibre bases — design doc §10.
