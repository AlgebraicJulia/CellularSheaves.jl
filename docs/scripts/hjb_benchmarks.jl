# Benchmarks for the double-integrator HJB solver (ControlSheaves.DoubleIntegratorHJB).
#
# Planar problem with a disc thrust bound (‖u‖ ≤ 1), discount 0.5, on n⁴ grids
# over [-2, 2]⁴, solved by policy iteration (a semismooth Newton method). The
# methods differ only in how each Newton step A_u V = b_u is solved:
#
#   direct      sparse LU of the whole grid (partial pivoting)
#   gmres-gs    serial GMRES, Gauss–Seidel preconditioner          (serial baseline)
#   gmres-sgs   serial GMRES, symmetric Gauss–Seidel preconditioner (serial baseline)
#   multicolor  Schwarz on b⁴ overlapping boxes, multicolor sweeps  (threaded)
#   ras-gmres   Schwarz on b⁴ boxes, GMRES with a RAS preconditioner (threaded)
#
# HJB_PHASE=correctness  every method on small grids, compared with the direct LU
# HJB_PHASE=runtime      iterative methods only, at the thread count Julia was
#                        started with; run once per thread count for scaling
#
# BLAS runs on one thread, so the only parallelism is the decomposition's.
#
# Run with, e.g.:
#   HJB_PHASE=correctness julia -t 8 --project=. docs/scripts/hjb_benchmarks.jl
#   HJB_PHASE=runtime HJB_SIZES=21,25,29 julia -t 4 --project=. docs/scripts/hjb_benchmarks.jl
# Writes docs/figures/hjb/correctness.csv and docs/figures/hjb/runtime_t<threads>.csv.

using CellularSheaves
using CellularSheaves.ControlSheaves.DoubleIntegratorHJB
using CellularSheaves.ControlSheaves.DoubleIntegratorHJB: solve
using CellularSheaves.NetworkSheaves.SchwarzMethods: SchwarzGMRES, ParallelSweep
using LinearAlgebra
using Printf

BLAS.set_num_threads(1)
const OUT = joinpath(@__DIR__, "..", "figures", "hjb")
mkpath(OUT)
const PHASE = get(ENV, "HJB_PHASE", "correctness")
const THREADS = Threads.nthreads()
const STATES = ([1.0, 0.0, 0.0, 0.0], [0.5, 0.5, 0.0, 0.0], [0.0, 0.0, 0.5, -0.5])

blocks(n) = fill(max(2, round(Int, n / 7)), 4)        # about 7 points per box per dimension

function method(name, n)
    name == "direct" && return DirectPolicyEvaluation()
    name == "gmres-gs" && return KrylovPolicyEvaluation(preconditioner = :gauss_seidel)
    name == "gmres-sgs" && return KrylovPolicyEvaluation(preconditioner = :symmetric_gauss_seidel)
    name == "multicolor" && return SchwarzPolicyEvaluation(blocks(n))
    name == "ras-gmres" && return SchwarzPolicyEvaluation(blocks(n);
        algorithm = SchwarzGMRES(sweep = ParallelSweep(), tol = 1e-10, maxiter = 500))
    error("unknown method $name")
end

struct Run
    n::Int
    method::String
    threads::Int
    total::Float64
    sol::HJBSolution
end

function run(n, name)
    grid = StateGrid(fill(-2.0, 4), fill(2.0, 4), fill(n, 4))
    prob = HJBProblem(grid; control_bound = 1.0, constraint = :disc)
    GC.gc()
    t = @elapsed sol = solve(prob, PolicyIteration(evaluation = method(name, n)))
    s = sol.seconds
    @printf("n=%2d (%7d unknowns) %-10s t=%2d %8.2f s  [assembly %.1f, setup %.1f, linear %.1f, improvement %.1f]  PI %2d  inner %5d  conv %s\n",
        n, length(grid), name, THREADS, t, s.assembly, s.setup, s.linear, s.improvement,
        sol.iterations, sum(sol.linear_iterations), sol.converged)
    flush(stdout)
    return Run(n, name, THREADS, t, sol)
end

row(r::Run, extra...) = join((r.n, r.n^4, repr(r.method), r.threads, r.total, values(r.sol.seconds)...,
    r.sol.iterations, sum(r.sol.linear_iterations), r.sol.converged, extra...), ",")
const HEADER = "n,unknowns,method,threads,seconds,assembly,setup,linear,improvement,policy_iterations,inner_iterations,converged"

for name in ("direct", "gmres-gs", "gmres-sgs", "multicolor", "ras-gmres")    # compile everything
    run(9, name)
end

if PHASE == "correctness"
    sizes = parse.(Int, split(get(ENV, "HJB_SIZES", "13,17"), ","))
    open(joinpath(OUT, "correctness.csv"), "w") do io
        println(io, HEADER, ",max_rel_diff_vs_direct,V1,V2,V3")
        for n in sizes
            ref = run(n, "direct")
            scale = maximum(abs, ref.sol.values)
            for name in ("direct", "gmres-gs", "gmres-sgs", "multicolor", "ras-gmres")
                r = name == "direct" ? ref : run(n, name)
                diff = maximum(abs, r.sol.values - ref.sol.values) / scale
                @printf("    max |V − V_direct| / max|V_direct| = %.1e\n", diff)
                println(io, row(r, diff, (value_at(r.sol, x) for x in STATES)...))
                flush(io)
            end
        end
    end
else
    sizes = parse.(Int, split(get(ENV, "HJB_SIZES", "21,25,29"), ","))
    names = THREADS == 1 ? ("gmres-gs", "gmres-sgs", "multicolor", "ras-gmres") : ("multicolor", "ras-gmres")
    open(joinpath(OUT, "runtime_t$(THREADS).csv"), "w") do io
        println(io, HEADER)
        for n in sizes, name in names
            println(io, row(run(n, name)))
            flush(io)
        end
    end
end
