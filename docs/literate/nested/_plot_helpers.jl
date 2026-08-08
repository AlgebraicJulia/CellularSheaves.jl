# Shared plotting helpers for the `nested/` literate examples -- purely about how these pages
# compose static composite plots, so kept docs-side rather than promoted into `NestedSystems`
# alongside the simulation driver and multi-pin helpers (see `src/ControlSheaves/NestedSystems.jl`).
#
# Not a standalone documentation page: `docs/make.jl`'s literate walk skips files whose name
# starts with `_`, so this produces no generated markdown/notebook.

"""
    snapshot_steps(steps::Int, n_snapshots::Int) -> Vector{Int}

`n_snapshots` step indices evenly spaced across `1:steps` (endpoints included). The temporal
resolution a simulation needs for accurate dynamics (every `dt`) is far finer than a static plot
needs for a legible composite of a formation's shape over time -- this picks the sparse subset
worth actually drawing.
"""
snapshot_steps(steps::Int, n_snapshots::Int) = unique(round.(Int, range(1, steps, length=n_snapshots)))

"""
    fade_alpha(i::Int, n::Int; lo=0.15, hi=0.85) -> Float64

Linear alpha ramp from `lo` (first of `n` items, `i == 1`) to `hi` (last, `i == n`) -- used so a
sequence of overlaid snapshots reads as "earlier is fainter, later is more solid" rather than a
uniform, hard-to-parse overlay.
"""
fade_alpha(i::Int, n::Int; lo::Real=0.15, hi::Real=0.85) = n <= 1 ? hi : lo + (hi - lo) * (i - 1) / (n - 1)
