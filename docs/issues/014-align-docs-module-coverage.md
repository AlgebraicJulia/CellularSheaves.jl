# Issue 014 — Add `Formations` API Page and Align Documented Module Coverage

**Parent design**: `007-nested-layered-systems-design.md` (§2, docs gap noted during Issue 008)
**Prerequisites**: none — independent of the whole nested-systems chain
**Blocks**: nothing

---

## Background

While implementing Issue 008 it emerged that `Formations` has **no API page at all**: none of
`se3_translation_matrix`, `se3_rotation_matrix`, `se3_affine_matrix`, `affine_translation_matrix`,
`build_escort_ring`, `build_escort_clique`, or the new `build_escort_topology` render anywhere in
the docs, despite all carrying docstrings.

Investigating that surfaced a second, larger inconsistency. `CLAUDE.md` states:

> Keep `makedocs(modules=...)` in `docs/make.jl` aligned with all documented modules.

It currently is not. `docs/make.jl:67` lists nine modules:

```
CellularSheaves, ControlSheaves, ControlSheaves.Tikhonov, ControlSheaves.AgentControllers,
ControlSheaves.DistributedLayeredControl, ControlSheaves.Layered,
ControlSheaves.MultiAgentTracking, ControlSheaves.MultiAgentTracking.QuadraticCosts,
AsynchSheaves
```

but `docs/src/api/*.md` contains `@autodocs` blocks for **eleven modules absent from that list**:

```
NetworkSheaves.EuclideanSheaves      NetworkSheaves.Pushforwards
NetworkSheaves.GraphHomomorphisms    NetworkSheaves.Pushouts
NetworkSheaves.PotentialSheaves      NetworkSheaves.SheafMorphisms
NetworkSheaves.TrajectorySheaves     NetworkSheaves.DistributedSolve
NetworkSheaves.CellularSheafParser   BlockSparseArrays
SheafInterface
```

These pages still render — `@autodocs Modules=[X]` pulls docstrings whether or not `X` appears in
`modules=`. What the omission actually costs is **coverage checking**: Documenter's `checkdocs`
pass only inspects modules named in `modules=`, so undocumented exports in eleven of the repo's
core modules are invisible. Compounding this, `docs/make.jl:75` sets `checkdocs=:none`, which
disables the pass entirely even for the nine listed modules.

Two further modules are loaded but documented nowhere: the `TrackingDSL` submodules
(`TrackingDSLParser`, `TrackingDSLValidator`, `TrackingDSLResolver`, `TrackingDSLLowering`,
`TrackingDSLTerm`) and `IPM`. Note `dsl.md` documents `CellularSheafParser` — the `@cellular_sheaf`
DSL — not `TrackingDSL`, so the two are easy to conflate.

---

## Codebase Orientation

| File | Why |
|---|---|
| `docs/make.jl:67` | The `modules=[...]` list to extend. |
| `docs/make.jl:75–76` | `checkdocs=:none` and `warnonly=[:cross_references]`. Relaxing `checkdocs` is the optional part of this issue — see the staging note below. |
| `docs/make.jl:77+` | The `pages=Any[...]` nav tree; a new page must be registered here or it renders orphaned. |
| `docs/src/api/pushforwards.md` | A minimal, well-formed example page to copy for the new `formations.md`. |
| `src/network_sheaves/Formations.jl` | The module to document. All public names already carry docstrings, so no docstring writing should be needed — verify rather than assume. |
| `src/network_sheaves/NetworkSheaves.jl:33` | `@reexport using .Formations` — confirms the module is loaded and its names are public. |

---

## Requested Implementation

### 1. Add the `Formations` API page

Create `docs/src/api/formations.md` mirroring the structure of an existing page:

```markdown
# Formations

Sheaf constructions for rigid multi-agent formations: homogeneous affine transforms and
the escort-topology primitives used by the layered control architectures.

```@autodocs
Modules = [CellularSheaves.NetworkSheaves.Formations]
```
```

Register it in `docs/make.jl`'s `pages=Any[...]` under the same API section as the other sheaf
pages.

### 2. Align `modules=` with reality

Extend `docs/make.jl:67` to include every module that has an `@autodocs` block — the eleven listed
above plus `Formations`. Verify by grepping `Modules = [...]` across `docs/src/api/*.md` and
diffing against the `modules=` list; the two sets should match exactly afterwards.

### 3. Document `TrackingDSL` and `IPM` (optional, only if cheap)

If these modules' public names already carry docstrings, add pages for them the same way. If they
do not, **stop and leave them out** — writing missing docstrings for two substantial modules is a
separate piece of work, not part of a docs-wiring fix. Report which case applies.

### 4. Tighten `checkdocs` (staged, optional)

Once `modules=` is aligned, try `checkdocs=:exports`. This is the setting that makes the coverage
guarantee real rather than nominal.

**Expect this to surface a backlog of undocumented exports.** If the count is small, fix them in
this issue. If it is large, revert to `checkdocs=:none`, record the count in the issue's closing
report, and open a follow-up. Do not let an incidental docs-wiring fix balloon into writing dozens
of docstrings.

---

## Tests to Write

Documentation has no unit tests; verification is the build itself.

```bash
julia --project=docs docs/make.jl
```

Additionally add a lightweight consistency check to the test suite so this cannot silently drift
again:

```julia
@testset "docs: every @autodocs module is declared in makedocs(modules=...)" begin
    make_src = read(joinpath(@__DIR__, "..", "docs", "make.jl"), String)
    api_dir  = joinpath(@__DIR__, "..", "docs", "src", "api")
    autodoc_mods = Set{String}()
    for f in readdir(api_dir; join=true)
        endswith(f, ".md") || continue
        for m in eachmatch(r"Modules\s*=\s*\[([^\]]+)\]", read(f, String))
            for mod in split(m.captures[1], ",")
                push!(autodoc_mods, strip(mod))
            end
        end
    end
    for mod in autodoc_mods
        @test occursin(mod, make_src)
    end
end
```

Place it in a new `test/docs_consistency.jl` included from `test/runtests.jl`. It reads files
only — no Documenter dependency, so it will not slow the suite or require `docs/Project.toml`.

---

## Verification Checklist

- [ ] `docs/src/api/formations.md` created and registered in `pages=`.
- [ ] All `Formations` public names render, including `build_escort_topology` from Issue 008.
- [ ] `modules=` in `docs/make.jl` matches the set of modules referenced by `@autodocs` blocks
      exactly — verified by the new consistency test, not by eye.
- [ ] `julia --project=docs docs/make.jl` completes; no new warnings introduced.
- [ ] Cross-reference warnings checked, per `CLAUDE.md` ("they appear as warnings, not errors").
- [ ] `TrackingDSL`/`IPM` either documented or explicitly deferred with a reason reported.
- [ ] `checkdocs` outcome reported: tightened, or left at `:none` with the backlog count recorded.
- [ ] Full test suite still passes.

---

## Out of Scope

- Writing missing docstrings in bulk. If `checkdocs=:exports` reveals a large backlog, record it
  and open a follow-up issue.
- Restructuring the docs nav tree or rewriting existing pages.
- `ConicSheaves` — it exists at `src/network_sheaves/ConicSheaves.jl` but is **not** `include`d by
  `NetworkSheaves.jl`, so it is not loaded and cannot be documented. Its status is a separate
  question worth raising, not resolving here.
- Any literate example changes.
