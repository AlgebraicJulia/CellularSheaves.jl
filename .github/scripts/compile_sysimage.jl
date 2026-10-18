using Pkg
using TOML

println("📥 Instantiating project and tracking artifacts...")
Pkg.activate(".")
Pkg.instantiate()
Pkg.precompile()

using PackageCompiler

# Every library dependency from Project.toml goes into the image, except those
# pinned by [sources] (git URLs or paths): Pkg refuses to add a package by URL
# when it is already in the system image, so Pkg.test could not set up its
# environment.
project = TOML.parsefile("Project.toml")
from_sources = keys(get(project, "sources", Dict{String,Any}()))
pkgs = [Symbol(k) for k in keys(project["deps"]) if k != "PackageCompiler" && !(k in from_sources)]

println("📦 Building system image blob for packages: $pkgs")
create_sysimage(
    pkgs;
    project=".",
    sysimage_path=ARGS[1],
    cpu_target="generic;x86-64-v3;x86-64-v4"
)
println("🎯 System image compiled successfully!")
