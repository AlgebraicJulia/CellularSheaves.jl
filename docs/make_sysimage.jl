# Build a precompiled system image for the documentation workflow.
#
# This script is called by .github/workflows/docs.yml when the sysimage cache
# is stale. It compiles a sysimage from *all* dependencies of the docs
# environment **except** the package being documented (CellularSheaves).
# The list is generated automatically, so you no longer need to maintain a
# manual hard‑coded list.
#
# Usage (from repository root):
#   julia --project=docs docs/make_sysimage.jl

using Pkg

# ️1. Activate the docs environment (which already has CellularSheaves
#    develop‑ed and all docs dependencies listed in docs/Project.toml).
Pkg.activate(@__DIR__)
Pkg.instantiate()

# 2️⃣ Query the **direct** dependencies listed in the docs Project.toml.
#    `Pkg.project()` returns the active project's manifest where `deps`
#    is a Dict mapping package name → UUID for *direct* dependencies only.
proj = Pkg.project()
# `proj.dependencies` contains the direct dependencies (including CellularSheaves).
# Extract the keys (the package names) and drop CellularSheaves.
package_names = filter(name -> name != "CellularSheaves", collect(keys(proj.dependencies)))

println("Creating sysimage with direct docs dependencies (excluding CellularSheaves):")
println(package_names)

using PackageCompiler

PackageCompiler.create_sysimage(
    package_names;                     # all deps except CellularSheaves
    sysimage_path = joinpath(@__DIR__, "..", "sysimage.so"),
    cpu_target = "generic",
)

