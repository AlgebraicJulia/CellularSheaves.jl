# Distributed benchmarks for the implicit (matrix-free) HJB policy iteration
# (GridPolicyIteration with an MPI communicator): one box of the n⁴ grid per
# MPI rank, KernelAbstractions kernels on each rank's CPU threads.
#
# For each grid size and preconditioner the solve runs HJB_REPEATS times (the
# nodes are shared with other jobs, so the minimum is the least disturbed
# measurement); rank 0 writes one CSV row with the min, median and max wall
# time (the slowest rank's) and the phase breakdown of the last run.
#
#   red_black_sgs  BiCGStab + the global red–black SGS sweep: the same operator
#                  for every number of ranks (two halo exchanges per sweep)
#   ras            BiCGStab + restricted additive Schwarz on the boxes extended
#                  by one point, one local red–black sweep per box (one halo
#                  exchange per application)
#
# Launch with an MPI launcher, e.g. inside one node:
#   julia --project=<env with MPI> -e 'using MPI; mpiexec(exe -> run(`$exe -n 16 julia -t 1 --project=<env> docs/scripts/hjb_mpi_benchmarks.jl`))'
# Environment: HJB_SIZES (default 33), HJB_PRECONDITIONERS (red_black_sgs,ras),
# HJB_REPEATS (3), HJB_OUT (docs/figures/hjb). Set
# JULIA_THREAD_SLEEP_THRESHOLD=infinite when running several threads per rank.
using MPI
using CellularSheaves
using CellularSheaves.ControlSheaves.DoubleIntegratorHJB
using Printf
using Statistics: median

MPI.Init()
const COMM = MPI.COMM_WORLD
const BOXES = mpi_boxes(COMM)
const RANK = MPI.Comm_rank(COMM)
const NRANKS = MPI.Comm_size(COMM)
const THREADS = Threads.nthreads()
const SIZES = parse.(Int, split(get(ENV, "HJB_SIZES", "33"), ","))
const PRECONDITIONERS = Symbol.(split(get(ENV, "HJB_PRECONDITIONERS", "red_black_sgs,ras"), ","))
const REPEATS = parse(Int, get(ENV, "HJB_REPEATS", "3"))
const OUT = get(ENV, "HJB_OUT", joinpath(@__DIR__, "..", "figures", "hjb"))

problem(n) = HJBProblem(StateGrid(fill(-2.0, 4), fill(2.0, 4), fill(n, 4)); control_bound = 1.0, constraint = :disc)

for pre in PRECONDITIONERS                                     # compile everything
    solve(problem(9), GridPolicyIteration(communicator = BOXES, preconditioner = pre))
end

const HEADER = "n,unknowns,preconditioner,ranks,threads,boxes,seconds_min,seconds_median,seconds_max," *
               "assembly,linear,improvement,policy_iterations,inner_iterations,converged"
io = RANK == 0 ? (mkpath(OUT); open(joinpath(OUT, "mpi_r$(NRANKS)_t$(THREADS).csv"), "w")) : devnull
println(io, HEADER)
for n in SIZES, pre in PRECONDITIONERS
    prob = problem(n)
    alg = GridPolicyIteration(communicator = BOXES, preconditioner = pre)
    times = Float64[]
    local sol
    for _ in 1:REPEATS
        GC.gc()
        MPI.Barrier(COMM)
        t = @elapsed sol = solve(prob, alg)
        push!(times, MPI.Allreduce(t, max, COMM))
    end
    boxes = join(balanced_ranks(NRANKS, (n, n, n, n)), "x")
    s = sol.seconds
    if RANK == 0
        @printf("n=%2d (%8d unknowns) %-13s ranks %3d (%s) threads %2d  %8.2f s [median %.2f, max %.2f]  linear %.2f  PI %2d  inner %5d  conv %s\n",
            n, n^4, pre, NRANKS, boxes, THREADS, minimum(times), median(times), maximum(times), s.linear,
            sol.iterations, sum(sol.linear_iterations), sol.converged)
        flush(stdout)
    end
    println(io, join((n, n^4, pre, NRANKS, THREADS, boxes, minimum(times), median(times), maximum(times),
        s.assembly, s.linear, s.improvement, sol.iterations, sum(sol.linear_iterations), sol.converged), ","))
    flush(io)
end
RANK == 0 && close(io)
MPI.Finalize()
