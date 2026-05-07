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

# 2. Query the resolved dependency graph. `Pkg.dependencies()` returns a
#    dictionary keyed by UUID where each entry contains a `name` field.
#    We collect every package name and then filter out the target package
#    itself, ensuring the sysimage contains only the *dependencies*.
all_deps = collect(values(Pkg.dependencies()))
package_names = map(d -> d.name, all_deps)
filter!(n -> n != "CellularSheaves", package_names)  # exclude the core package

println("Creating sysimage with the following packages (excluding CellularSheaves):")
println(package_names)

using PackageCompiler

PackageCompiler.create_sysimage(
    package_names;                     # all deps except CellularSheaves
    sysimage_path = joinpath(@__DIR__, "..", "sysimage.so"),
    cpu_target = "generic",
)

