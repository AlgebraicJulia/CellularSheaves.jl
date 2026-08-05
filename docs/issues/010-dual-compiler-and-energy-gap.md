# Issue 010 — Hierarchical Solve, Direct Baseline, and Approximation-Gap Diagnostics

**Parent design**: `007-nested-layered-systems-design.md` (Issue 3)
**Prerequisites**: 009
**Blocks**: 011

---

## Mathematical Background

### The two compilations

A single `NestedSystemSpec` must compile two ways, so the hierarchy's approximation can be
measured rather than assumed.

**Hierarchical.** One harmonic extension at the coarsest level `H₀` with targets pinned, then
pure block-diagonal lifts down the tower:

```
q₀ = harmonic_extension(H₀, targets pinned)
q_k = B_v · q_{k-1}[v]     for each fibre v,  k = 1 … N
```

**Direct (baseline).** One harmonic extension on the fully expanded finest sheaf `H_N`, with the
same targets pinned. This is what `solve_direct_harmonic` (`Layered.jl:369`) does for the flat
case.

### Why the hierarchical answer differs

`pushforward_sheaf` builds each fibre basis `B_v` from the fibre's **exact** global-section space,
so `B` is block-diagonal over fibres and every configuration the tower can express satisfies all
fibre-internal constraints exactly. Internal Dirichlet energy is identically zero.

The hierarchical solve therefore answers a **constrained** version of the direct problem:

> Constrain every team to lie on its space of global sections (in formation).
> Subject to that, place the formations to minimize total sheaf energy.

### The invariant this issue must pin

The hierarchical solution is a **feasible point** of the direct problem — it satisfies every
constraint the direct problem imposes, plus extra ones. Therefore

```
E_hierarchical ≥ E_direct
```

with equality **iff** the direct optimum is already fibrewise-exact. Rigid single-target escort
teams are the equality case, which is why today's flat code is exact. A team observing two
targets that pull apart is the strict case: the direct optimum deforms the team, the tower
cannot, and the gap opens.

This is a theorem, not a hypothesis. The tests below pin it rather than investigate it.

---

## Codebase Orientation

| File | Why |
|---|---|
| `src/ControlSheaves/NestedSystems.jl` | `SheafTower` from Issue 009 — the input to both solvers. |
| `src/network_sheaves/EuclideanSheaves.jl` | `_harmonic_extension_restricted_laplacian(s, boundary)` (line 626) returns `(boundary_verts, interior_verts, L_II, L_IB)`; solve `L_II \ (-L_IB * p_boundary)`. `coboundary_map` (line 254) gives `δ` for the energy computation. |
| `src/ControlSheaves/Layered.jl` | `solve_high_level_harmonic` (line 304), `solve_mid_level_harmonic` (line 339), `solve_direct_harmonic` (line 369) — the flat precedents. Note the boundary-ordering handling at lines 312–329, which the general version must not get wrong. |
| `test/network_sheaves/LayeredEscort.jl:61–69` | Computes both solutions today and compares **only their sizes**. This issue supplies the comparison that has never been made. |
| `CLAUDE.md` | Solver conventions — `ChordalLDLt` for new sparse symmetric PSD systems. |

---

## Requested Implementation

Extend `src/ControlSheaves/NestedSystems.jl`.

```julia
"""
    solve_hierarchical(tower::SheafTower, target_values::Vector{Vector{Float64}})
        -> Vector{Vector{Vector{Float64}}}

Solve the tower top-down: one harmonic extension at `tower.levels[1]` with `target_values`
pinned at the target vertices, then successive fibre-basis lifts to each finer level.

Returns one cochain per level, outermost first, each as a vector of per-vertex values. The last
entry holds the per-agent reference states.
"""
function solve_hierarchical(tower::SheafTower, target_values::Vector{Vector{Float64}})

"""
    solve_direct(tower::SheafTower, target_values::Vector{Vector{Float64}})
        -> Vector{Vector{Float64}}

Baseline: solve the harmonic extension directly on the fully expanded finest sheaf
`tower.levels[end]`, pinning the same targets. Returns per-vertex values on the finest level.

This is the unconstrained optimum against which `solve_hierarchical` is measured; see
`approximation_gap`.
"""
function solve_direct(tower::SheafTower, target_values::Vector{Vector{Float64}})

"""
    sheaf_energy(F::EuclideanSheaf, x::AbstractVector) -> Float64

Dirichlet energy `‖δx‖²` of a flat cochain `x` on `F`.
"""
function sheaf_energy(F::EuclideanSheaf, x::AbstractVector)

"""
    approximation_gap(tower::SheafTower, target_values)
        -> (; hierarchical, direct, gap, relative_gap)

Energy of both solutions on the finest sheaf, and their difference. `gap` is guaranteed
nonnegative up to floating-point tolerance: the hierarchical solution is feasible for the direct
problem, so it can never achieve lower energy. A gap of zero means the direct optimum was already
fibrewise-exact (rigid formations suffice); a positive gap quantifies what rigidity costs.
"""
function approximation_gap(tower::SheafTower, target_values)
```

### Algorithm sketch — `solve_hierarchical`

1. Build `boundary::Dict{Int,Vector{Float64}}` mapping `tower.target_vertices[t] =>
   target_values[t]`.
2. `bverts, iverts, L_II, L_IB = _harmonic_extension_restricted_laplacian(tower.levels[1], boundary)`.
3. Assemble `p_boundary` **in `bverts` order** — do not assume targets are `1:n_targets`. This is
   the ordering trap in `Layered.jl:312`.
4. Solve `q_interior = L_II \ (-L_IB * p_boundary)` via `ChordalLDLt` per `CLAUDE.md`, keeping
   matrices sparse.
5. Scatter boundary and interior values into a full cochain `q₀` using each vertex's own stalk
   dimension (stalks are **not** uniform — see design doc §4.2).
6. For `k = 1 … N`: for each coarse vertex `v`, compute `tower.bases[k][v] * q_{k-1}[v]` and
   scatter the result across `fiber_vertices(tower.homs[k], v)` in that function's ordering, which
   is the ordering `fiber_section_basis` documents for its rows.

### Algorithm sketch — `solve_direct`

Same as steps 1–5 against `tower.levels[end]`, with targets pinned at their finest-level vertex
indices. No lifting.

### Algorithm sketch — `sheaf_energy`

`d = coboundary_map(F); return norm(d * x)^2`. Keep `d` sparse.

---

## Tests to Write

Add to `test/ControlSheaves/NestedSystems.jl`, reusing Issue 009's spec helpers:

```julia
@testset "energy gap is nonnegative (feasibility theorem)" begin
    for spec in [two_level_spec(), three_level_spec(), irregular_spec(), shared_target_spec()]
        tower = build_sheaf_tower(spec)
        g = approximation_gap(tower, default_targets(spec))
        @test g.gap >= -1e-8
    end
end

@testset "rigid single-target teams: hierarchical == direct" begin
    tower = build_sheaf_tower(two_level_spec())   # rigid rings, one target each
    tv = default_targets(tower.spec)
    q_h = solve_hierarchical(tower, tv)[end]
    q_d = solve_direct(tower, tv)
    @test all(isapprox(a, b; atol=1e-8) for (a, b) in zip(q_h, q_d))
    @test approximation_gap(tower, tv).gap ≈ 0.0 atol=1e-8
end

@testset "team observing two separating targets: strict gap" begin
    tower = build_sheaf_tower(two_target_team_spec())
    near = [[0.0, 0.0, 1.5, 1.0], [0.2, 0.0, 1.5, 1.0]]
    far  = [[0.0, 0.0, 1.5, 1.0], [8.0, 0.0, 1.5, 1.0]]
    @test approximation_gap(tower, far).gap > approximation_gap(tower, near).gap
    @test approximation_gap(tower, far).gap > 1e-6
end

@testset "targets are reproduced exactly at the finest level" begin
    tower = build_sheaf_tower(three_level_spec())
    tv = default_targets(tower.spec)
    q = solve_hierarchical(tower, tv)[end]
    for (t, v) in enumerate(tower.target_vertices)
        @test q[v] ≈ tv[t]
    end
end

@testset "golden regression against the flat Layered pipeline" begin
    # A 2-level tree with no refinement must reproduce today's flat result.
    tower = build_sheaf_tower(flat_equivalent_spec())
    q_new = solve_hierarchical(tower, flat_targets())[end]
    q_old = solve_mid_level_harmonic(
        solve_high_level_harmonic(flat_PfF(), flat_bases(), flat_targets()), flat_bases(), 4)
    @test agents_of(q_new) ≈ q_old atol=1e-8
end
```

Also close the long-standing hole in the existing flat tests:

```julia
# test/network_sheaves/LayeredEscort.jl, replacing the size-only assertion at line 69
@test qstar_agents ≈ qstar_direct atol=1e-8
```

---

## Verification Checklist

- [ ] Boundary values assembled in `bverts` order, not assumed to be `1:n_targets`.
- [ ] Per-vertex stalk dimensions respected everywhere; no `fill(D, n)` assumption.
- [ ] `ChordalLDLt` used for the harmonic solves, per `CLAUDE.md`; no dense conversions.
- [ ] Fibre bases read from the tower, never recomputed inside a solve.
- [ ] `gap >= 0` holds across every spec in the test suite.
- [ ] The rigid case shows equality, the non-rigid case shows a strict gap — both directions
      tested, since a one-sided test would pass trivially on a broken implementation.
- [ ] `LayeredEscort.jl:69` upgraded from a size assertion to a value comparison.
- [ ] Docstrings explain the feasibility argument for readers who are not specialists.

---

## Out of Scope

- Soft/spectral fibre bases — design doc §10, explicitly deferred.
- Declared per-edge restriction maps — Issue 011.
- Agent dynamics, LQR tracking, and simulation loops; this issue produces reference states only.
- Performance tuning of the tower build.
