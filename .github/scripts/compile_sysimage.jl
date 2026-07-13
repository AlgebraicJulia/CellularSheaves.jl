using Pkg

println("📥 Instantiating project and tracking artifacts...")
Pkg.activate(".")
Pkg.instantiate()
Pkg.precompile()

using PackageCompiler

# Automatically extract all your library dependencies from Project.toml
pkgs = [Symbol(k) for k in keys(Pkg.project().dependencies) if k != "PackageCompiler"]

println("📦 Building system image blob for packages: $pkgs")
create_sysimage(
    pkgs;
    project=".",
    sysimage_path=ARGS[1],
    cpu_target="generic;x86-64-v3;x86-64-v4"
)
println("🎯 System image compiled successfully!")
