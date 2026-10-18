# The MPI tests, run apart from Pkg.test (GitHub CI has no MPI):
#
#   julia --project=<env> test/mpi/runtests.jl
#
# <env> is an environment with CellularSheaves (developed from this repository), MPI
# and Test. Each test script in this directory is launched with MPI.jl's `mpiexec` on
# several rank counts, with that same environment. When MPI is not available (not in
# the environment, or its library or launcher does not work) the tests are skipped
# with a message, and the script exits successfully.
#
# Environment variables:
#   CELLULARSHEAVES_MPI_RANKS   rank counts to run, comma separated (default 1,3,4)
#   CELLULARSHEAVES_MPI_PREFIX  a command put before each rank's julia, for instance
#                               `in-sandbox` on HiPerGator (default none)

skip(reason) = (@info "Skipping the MPI tests: $reason"; exit(0))

Base.find_package("MPI") === nothing && skip("MPI is not in the active environment ($(Base.active_project())).")
Base.find_package("CellularSheaves") === nothing &&
    skip("CellularSheaves is not in the active environment ($(Base.active_project())).")

using Test
try
    @eval using MPI
catch err
    skip("MPI.jl does not load: $(sprint(showerror, err))")
end

const RANKS = parse.(Int, split(get(ENV, "CELLULARSHEAVES_MPI_RANKS", "1,3,4"), ','))
const PREFIX = split(get(ENV, "CELLULARSHEAVES_MPI_PREFIX", ""))
const PROJECT = Base.active_project()
const SCRIPTS = ["grid_boxes.jl"]

# Whether the launcher starts a rank at all, before reporting test failures.
launches = try
    MPI.mpiexec(exe -> success(pipeline(`$exe -n 1 $PREFIX $(Base.julia_cmd()) --startup-file=no -e 0`; stdout, stderr)))
catch err
    skip("the MPI launcher does not run: $(sprint(showerror, err))")
end
launches || skip("the MPI launcher could not start a rank.")

@testset "MPI" begin
    for script in SCRIPTS, n in RANKS
        @testset "$script on $n ranks" begin
            path = joinpath(@__DIR__, script)
            ok = MPI.mpiexec() do exe                       # the launcher's environment applies inside the block
                cmd = `$exe -n $n $PREFIX $(Base.julia_cmd()) --project=$PROJECT --startup-file=no $path`
                success(pipeline(cmd; stdout, stderr))
            end
            @test ok
        end
    end
end
