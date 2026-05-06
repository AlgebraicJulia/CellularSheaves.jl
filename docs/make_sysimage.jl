# Build a precompiled system image for the documentation workflow.
#
# This script is called by .github/workflows/docs.yml when the sysimage cache
# is stale. It compiles a sysimage from the manually maintained package list
# below, which is intended to cover dependencies from Project.toml,
# test/Project.toml, and docs/Project.toml (except CellularSheaves itself,
# which changes with every source edit). Keep that list in sync as
# dependencies change. Heavy packages like Plots.jl and Documenter.jl
# benefit most from precompilation.
#
# Usage (from repo root):
#   julia --project=docs docs/make_sysimage.jl

using Pkg
Pkg.develop(PackageSpec(path=joinpath(@__DIR__, "..")))
Pkg.instantiate()

using PackageCompiler

create_sysimage(
    [
        # CellularSheaves and its transitive deps are included automatically
        "CellularSheaves",
        # Documentation dependencies (docs/Project.toml)
        "Documenter", "Literate", "Plots",
    ],
    sysimage_path=joinpath(@__DIR__, "..", "sysimage.so"),
    cpu_target="generic",
)
