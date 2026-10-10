# Build a Julia system image holding every dependency of CellularSheaves.jl,
# MPI.jl and CUDA.jl, and the test dependencies, but not the packages under
# development: CellularSheaves itself, CliqueTrees and Mumblebee (nor anything
# that depends on them). Those load as ordinary packages on top of the image,
# so they can be edited (`Pkg.develop`) without rebuilding it.
#
#   julia sysimage/build.jl
#
# Environment variables:
#   CS_SYSIMAGE_DIR   where the environment, builder and image go
#                     (default $HOME/sysimage)
#   CS_DEVELOP        comma-separated paths of local checkouts to develop into
#                     the environment first (e.g. a CliqueTrees or Mumblebee clone)
#   CS_MPI_LIBRARY    the system libmpi to bake in (default: HiPerGator's
#                     OpenMPI 5.0.7); "jll" keeps MPI.jl's bundled MPICH
#   CS_CUDA_VERSION   CUDA toolkit version of the cluster module (default 12.9);
#                     "artifact" lets CUDA.jl download its own runtime
#   CS_CPU_TARGET     JULIA_CPU_TARGET of the image (default
#                     "x86-64-v3;x86-64-v4,clone_all": AVX2 everywhere, AVX-512
#                     clones for CPUs that have it)
#
# The image is written to $CS_SYSIMAGE_DIR/cellularsheaves-deps-<julia>-<hash>.so,
# keyed by a hash of the environment's manifest, and linked as current.so. Run
# Julia through sysimage/julia.sh, which uses the image, that environment, its
# preferences and the matching CPU target.
#
# A package compiled into the image cannot be replaced by another version: if
# a change to CliqueTrees or Mumblebee needs a different version of one of
# their dependencies, rebuild. New dependencies are fine; they load normally
# until the next build.
using Pkg
using SHA: sha1

const ROOT = dirname(@__DIR__)
const DIR = get(ENV, "CS_SYSIMAGE_DIR", joinpath(homedir(), "sysimage"))
const ENVDIR = joinpath(DIR, "env")
const BUILDER = joinpath(DIR, "builder")
const UNDER_DEVELOPMENT = Set(["CellularSheaves", "CliqueTrees", "Mumblebee"])
# Loaded alongside CellularSheaves: MPI and GPU support, and the test suite's dependencies.
const EXTRA = ["MPI", "MPIPreferences", "CUDA", "KernelAbstractions", "Adapt", "Statistics", "Test", "Aqua",
               "JET", "Distributed", "Distributions", "Graphs", "BlockArrays", "Random", "LinearAlgebra",
               "SparseArrays"]
const MPI_LIBRARY = get(ENV, "CS_MPI_LIBRARY", "/apps/mpi/gcc/14.2.0/openmpi/5.0.7_el97/lib/libmpi.so")
const CUDA_VERSION = get(ENV, "CS_CUDA_VERSION", "12.9")
const CPU_TARGET = get(ENV, "CS_CPU_TARGET", "x86-64-v3;x86-64-v4,clone_all")

mkpath(ENVDIR)
mkpath(BUILDER)

# PackageCompiler lives in its own environment, outside the image.
Pkg.activate(BUILDER)
haskey(Pkg.project().dependencies, "PackageCompiler") || Pkg.add("PackageCompiler")

# The environment the image is built from and run with.
Pkg.activate(ENVDIR)
Pkg.develop(path=ROOT)
for path in split(get(ENV, "CS_DEVELOP", ""), ","; keepempty=false)
    Pkg.develop(path=String(strip(path)))
end
missing_extra = filter(p -> !haskey(Pkg.project().dependencies, p), EXTRA)
isempty(missing_extra) || Pkg.add(missing_extra)
Pkg.instantiate()

# Compile-time preferences, set in a fresh process so they are written to the
# environment's LocalPreferences.toml before anything is compiled.
julia = Base.julia_cmd()
if MPI_LIBRARY != "jll"
    isfile(MPI_LIBRARY) || error("CS_MPI_LIBRARY $MPI_LIBRARY not found (load the MPI module, or set CS_MPI_LIBRARY=jll)")
    run(`$julia --project=$ENVDIR --startup-file=no -e "using MPIPreferences; MPIPreferences.use_system_binary(;
        library_names = [\"$MPI_LIBRARY\"], mpiexec = \"srun\")"`)
end
if CUDA_VERSION != "artifact"
    run(`$julia --project=$ENVDIR --startup-file=no -e "using CUDA; CUDA.set_runtime_version!(v\"$CUDA_VERSION\"; local_toolkit = true)"`)
end

# Every package in the manifest except the ones under development and anything
# that (transitively) depends on them.
deps = Pkg.dependencies()
name_of = Dict(uuid => info.name for (uuid, info) in deps)
tainted = Dict{Base.UUID,Bool}()
function depends_on_development(uuid)
    get!(tainted, uuid) do
        info = deps[uuid]
        info.name in UNDER_DEVELOPMENT && return true
        tainted[uuid] = false                       # break cycles
        any(depends_on_development, values(info.dependencies))
    end
end
packages = sort!([name_of[u] for u in keys(deps) if !depends_on_development(u)])
excluded = sort!([name_of[u] for u in keys(deps) if depends_on_development(u)])
println("image: ", length(packages), " packages; loaded on top: ", join(excluded, ", "))

manifest = joinpath(ENVDIR, "Manifest.toml")
key = bytes2hex(sha1(read(manifest)))[1:12]
image = joinpath(DIR, "cellularsheaves-deps-$(VERSION)-$key.so")
if isfile(image)
    println("up to date: ", image)
else
    # PackageCompiler compiles only direct dependencies of the project it builds
    # from: give it a project listing every image package directly, with the
    # environment's exact manifest and preferences.
    image_project = joinpath(DIR, "image-project")
    rm(image_project; recursive=true, force=true)
    mkpath(image_project)
    uuid_of = Dict(info.name => uuid for (uuid, info) in deps)
    open(joinpath(image_project, "Project.toml"), "w") do io
        println(io, "[deps]")
        for p in packages
            println(io, p, " = \"", uuid_of[p], "\"")
        end
    end
    cp(manifest, joinpath(image_project, "Manifest.toml"))
    preferences = joinpath(ENVDIR, "LocalPreferences.toml")
    isfile(preferences) && cp(preferences, joinpath(image_project, "LocalPreferences.toml"))
    Pkg.activate(image_project)
    Pkg.instantiate()
    pushfirst!(LOAD_PATH, BUILDER)
    @eval using PackageCompiler
    # The workload loads CellularSheaves, which is not an image package: run it
    # with the full environment on the load path.
    push!(LOAD_PATH, ENVDIR)
    PackageCompiler.create_sysimage(Symbol.(packages); sysimage_path=image, project=image_project,
        precompile_execution_file=joinpath(@__DIR__, "precompile_workload.jl"), cpu_target=CPU_TARGET)
    Pkg.activate(ENVDIR)
end
current = joinpath(DIR, "current.so")
rm(current; force=true)
symlink(image, current)
open(joinpath(DIR, "cpu_target"), "w") do io
    println(io, CPU_TARGET)
end
println("system image: ", image, "\nrun Julia with: ", joinpath(@__DIR__, "julia.sh"))
