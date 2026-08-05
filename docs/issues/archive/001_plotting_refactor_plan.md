# Plotting Refactor Plan for CellularSheaves.jl

**Purpose**: Centralise all plotting logic currently scattered across examples, literate notebooks, and scripts into the optional `CellularSheavesPlots` extension using **Plots.jl recipes**.  This provides a single source‑of‑truth for colours, layout, axis handling, and enables downstream users to customise plots easily.

---

## Extension Skeleton
- `ext/CellularSheavesPlots/Project.toml` – already created; declares dependencies (`Plots`, `CSV`, `DataFrames`, `JSON3`) and a weak dependency on the core package.
- `ext/CellularSheavesPlots/src/`
  - `ScenarioPlots.jl` – **completed** (recipe for `ScenarioResult`).
  - `ErrorNormPlots.jl` – recipe for error‑norm time‑series with RMS annotation.
  - `Trajectory3D.jl` – recipe for 3‑D agent/target trajectories.
  - `SnapshotPlots.jl` – recipe for 2‑D formation snapshots (edges, pinning, triangle groups, square group).
  - `SpyPlots.jl` – simple `spy` view and side‑by‑side `spy` recipe.
  - `Utils.jl` – shared helpers (`_extract_num`, `_list_state_files`, `_rms`, `_edge_pairs`, `_agent_neighbors`, `_pos_at`, `_xy_limits`, `uniform_rows`, `exp_decay_rows`, etc.).
  - `Style.jl` – optional helper to set project‑wide default plot attributes (`set_default_style!`, `set_large_style!`).

## Workflow for Adding a New Plot
1. **Identify the visual pattern** (e.g., two‑panel time series, 3‑D trajectories).  If it matches an existing recipe, reuse it; otherwise add a new module.
2. **Implement the recipe** in the appropriate module, re‑using utilities from `Utils.jl`.
3. **Export the recipe** from the module and, if desired, add it to the extension’s `Project.toml` under `[extensions]`.
4. **Replace duplicated code** in examples/literate files with a single call to the recipe, e.g. `plot(sr)` or `plot(d)`.
5. **Add a small unit test** under `test/ext/CellularSheavesPlots/` that builds a synthetic object and verifies `isa(plot(...), Plot)`.
6. **Update documentation** (`docs/make.jl`) to include `CellularSheavesPlots` in the `modules=` list so the new recipes appear in the API docs.

---

## Immediate Tasks
- **Create `Utils.jl`** with the helper functions copied from `examples/multi_agent_target_tracking/src/visualization/plotter.jl` and add appropriate docstrings.
- **Add `ErrorNormPlots.jl`** and `Trajectory3D.jl` modules (use the code snippets from the planning notes).
- **Implement `SnapshotPlots.jl`** – wrap the existing `_snapshot_context`/`draw_single_snapshot` logic into a recipe.
- **Add `SpyPlots.jl`** for `spy` visualisations.
- **Write tests** in `test/ext/CellularSheavesPlots/` for each recipe.
- **Run the docs build** to ensure no missing cross‑references.

---

## Long‑Term Benefits
- Consistent visual style across the whole project.
- Easy customization for researchers: they can call `plot!(p, …)` after the recipe or overload style keywords.
- Optional dependency – core package stays lightweight.
- Centralised maintenance: fixing a bug or tweaking colours needs to be done only once.
- Ready foundation for future visualisation features (animations, heat‑maps, etc.).

---

*Prepared on $(date)*
