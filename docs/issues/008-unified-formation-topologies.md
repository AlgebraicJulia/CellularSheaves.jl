# Issue 008 — Unified Formation Topology Primitives

**Parent design**: `007-nested-layered-systems-design.md` (Issue 1)
**Prerequisites**: none
**Blocks**: 009, 010, 011

---

## Mathematical Background

An escort formation is a cellular sheaf on a small graph whose vertices are agents plus one
target node. Each agent `i` is assigned a geometric offset `d_i` from the formation centre;
the restriction map on an incident edge is the homogeneous affine translation
`affine_translation_matrix(d_i)`, so an edge `(i,j)` imposes `x_i - d_i = x_j - d_j`.

Two concerns are **independent** and this issue keeps them so:

- **Geometry** — where agents sit (the offsets `d_i`), currently angular placement on a circle.
- **Consensus topology** — which agents are wired to which (ring, star, path, clique).

Because every edge constraint has the form "agent `i` = centre + `d_i`", the constraint system is
globally realizable for *any connected* consensus topology. Hence for all four topologies the
fibre's exact global-section space is exactly `D`-dimensional, parameterized by the centre. This
is the property the hierarchical architecture depends on (design doc §3.2), so it must hold for
every topology this issue adds.

---

## Codebase Orientation

| File | Why |
|---|---|
| `src/network_sheaves/Formations.jl` | Home of the change. `build_escort_ring` (line 88) is already `D`/`affine`-parameterized and is the model to follow. `build_escort_clique` (line 134) does the same job but is **`D=4`-hardcoded with no `affine` flag** — it must be folded into the unified function. |
| `src/ControlSheaves/Layered.jl:237` | Existing caller of `build_escort_ring`; must keep working unchanged. |
| `test/network_sheaves/TransferMapLift.jl:16` | Another `build_escort_ring` caller. |
| `test/network_sheaves/Formations.jl` | Existing test file to extend. |
| `src/network_sheaves/NetworkSheaves.jl` | Export list to update. |

---

## Requested Implementation

```julia
"""
    build_escort_topology(kind::Symbol, n_agents::Int, target_node::Int, radius::Float64;
                          observers=1:n_agents, D::Int=4, affine::Bool=true) -> EuclideanSheaf{Float64}

Construct a `EuclideanSheaf` with stalk dimension `D` for an `n_agents` escort formation around
`target_node`.

`kind` selects the **consensus topology** — which agents share a constraint edge:

- `:ring`   — cycle, `i — (i mod n)+1`
- `:path`   — open chain, `i — i+1` for `i in 1:n-1`
- `:star`   — hub-and-spoke, agent `1 — i` for `i in 2:n`
- `:clique` — all pairs `i < j`

`kind` does **not** change the formation geometry: in every case agent `i` is placed at angle
`2π(i-1)/n_agents` and distance `radius` in the plane spanned by the first two translation
coordinates, matching `build_escort_ring`. Topology is who agrees with whom; geometry is where
they sit.

When `affine=true` (default) stalks carry `D-1` translation coordinates plus a homogeneous row
and restriction maps are `affine_translation_matrix` offsets. When `affine=false` stalks are `D`
plain Euclidean coordinates, every restriction map is the identity, and `radius` must be `0.0`.

`observers` names which agents (local indices `1:n_agents`) are pinned to the target.

For every `kind` the consensus graph is connected, so the resulting sheaf's space of global
sections is exactly `D`-dimensional — the formation is rigid and parameterized by its centre.
"""
function build_escort_topology(kind::Symbol, n_agents::Int, target_node::Int, radius::Float64;
                               observers=1:n_agents, D::Int=4, affine::Bool=true)
```

### Algorithm sketch

1. `@argcheck kind in (:ring, :path, :star, :clique)`.
2. `@argcheck affine || radius == 0.0` and `@argcheck all(1 .<= o .<= n_agents for o in observers)`
   — reuse the existing wording from `build_escort_ring`.
3. `@argcheck n_agents >= 2` for `:ring` (a 2-cycle degenerates; document the choice) and
   `n_agents >= 1` otherwise.
4. Build `restriction_matrix(i)` exactly as `build_escort_ring` does today (extract the existing
   closure verbatim — angular offsets when `affine`, identity otherwise).
5. Generate the consensus edge list from `kind`, then `add_sheaf_edge!(sheaf, i, j,
   restriction_matrix(i), restriction_matrix(j))` for each.
6. Pin observers to `target_node` with `restriction_matrix(i)` on the agent side and `I_D` on the
   target side.

### Backwards compatibility

Keep both existing names as thin wrappers so no caller changes:

```julia
build_escort_ring(n, t, r; observers=1:n, D=4, affine=true) =
    build_escort_topology(:ring, n, t, r; observers, D, affine)

build_escort_clique(n, t, r; observers=1:n, D=4, affine=true) =
    build_escort_topology(:clique, n, t, r; observers, D, affine)
```

Note this *widens* `build_escort_clique`, which previously accepted no `D`/`affine` keywords and
always produced `D=4`. Its default behaviour is unchanged.

---

## Tests to Write

Add to `test/network_sheaves/Formations.jl`:

```julia
@testset "build_escort_topology — edge counts per kind" begin
    n = 5
    for (kind, n_consensus) in [(:ring, 5), (:path, 4), (:star, 4), (:clique, 10)]
        s = build_escort_topology(kind, n, n+1, 0.3; observers=[1])
        @test ne(s.underlying_graph) == n_consensus + 1   # +1 observer edge
    end
end

@testset "build_escort_topology — every kind is rigid (section space is D-dimensional)" begin
    for kind in (:ring, :path, :star, :clique)
        s = build_escort_topology(kind, 4, 5, 0.3; observers=[1])
        B = fiber_section_basis(s, collect(1:5),
                                [(src(e), dst(e)) for e in edges(s.underlying_graph)])
        @test size(B, 2) == 4
    end
end

@testset "build_escort_topology — ring wrapper reproduces build_escort_ring" begin
    a = build_escort_topology(:ring, 6, 7, 0.3; observers=[1])
    b = build_escort_ring(6, 7, 0.3; observers=[1])
    @test a == b
end

@testset "build_escort_topology — clique honours D and affine" begin
    s = build_escort_topology(:clique, 3, 4, 0.0; D=3, affine=false)
    @test all(==(3), s.vertex_stalks)
    @test_throws Exception build_escort_topology(:clique, 3, 4, 0.5; affine=false)
end

@testset "build_escort_topology — rejects unknown kind" begin
    @test_throws Exception build_escort_topology(:hexagon, 4, 5, 0.3)
end
```

---

## Verification Checklist

- [ ] `build_escort_topology` exported from `Formations` and re-exported via `NetworkSheaves.jl`.
- [ ] `build_escort_ring` / `build_escort_clique` retained as wrappers; all existing callers
      (`Layered.jl:237`, `TransferMapLift.jl`) pass unchanged.
- [ ] Section space is `D`-dimensional for all four kinds (rigidity property).
- [ ] `@argcheck` used for all input validation, per `CLAUDE.md`.
- [ ] Docstring explains topology-vs-geometry for non-specialists.
- [ ] `docs/src/api/` autodocs updated if `Formations` gains a new public name.
- [ ] Full test suite passes.

---

## Out of Scope

- Custom (non-circular) geometric placement — offsets stay angular for every kind. Arbitrary
  geometry belongs to the raw-matrix escape hatch in Issue 011.
- Per-agent radii or non-planar formations.
- Support-pod construction (Issue 009).
- Anything DSL-related.
