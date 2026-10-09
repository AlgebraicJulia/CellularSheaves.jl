# Benchmarks for the double-integrator HJB solver (ControlSheaves.DoubleIntegratorHJB).
#
# Policy iteration on the planar problem with a disc thrust bound (‖u‖ ≤ 1),
# discount 0.5, on n⁴ grids over [-2, 2]⁴. Each policy evaluation is solved by
#   - a sparse LU of the whole grid (partial pivoting), while it fits, and
#   - Schwarz on b⁴ overlapping boxes: multicolor sweeps, and GMRES
#     preconditioned by a parallel sweep,
# reporting wall time, policy iterations, Schwarz iterations, the process's peak
# memory so far (Schwarz runs go first at each size, the direct LU last) and
# the value at a few states (for grid convergence).
#
# Run with:  julia -t 8 --project=. docs/scripts/hjb_benchmarks.jl
# Writes docs/figures/hjb/benchmarks.csv.

using CellularSheaves
using CellularSheaves.ControlSheaves.DoubleIntegratorHJB
using CellularSheaves.ControlSheaves.DoubleIntegratorHJB: solve
using CellularSheaves.NetworkSheaves.SchwarzMethods: SchwarzIteration, SchwarzGMRES, MulticolorSweep, ParallelSweep
using Printf

const OUT = joinpath(@__DIR__, "..", "figures", "hjb")
mkpath(OUT)
const STATES = ([1.0, 0.0, 0.0, 0.0], [0.5, 0.5, 0.0, 0.0], [0.0, 0.0, 0.5, -0.5])
const DIRECT_MAX = parse(Int, get(ENV, "HJB_DIRECT_MAX", "21"))
const SIZES = parse.(Int, split(get(ENV, "HJB_SIZES", "13,17,21,25,29")))

struct HJBRow
    n::Int
    unknowns::Int
    method::String
    seconds::Float64
    policy_iterations::Int
    schwarz_iterations::Int
    converged::Bool
    maxrss_gb::Float64
    values::Vector{Float64}
end

function run(n, method, evaluation)
    grid = StateGrid(fill(-2.0, 4), fill(2.0, 4), fill(n, 4))
    prob = HJBProblem(grid; control_bound = 1.0, constraint = :disc)
    GC.gc()
    t = @elapsed sol = solve(prob, PolicyIteration(evaluation = evaluation))
    row = HJBRow(n, length(grid), method, t, sol.iterations, sum(sol.linear_iterations), sol.converged,
        Sys.maxrss() / 2^30, [value_at(sol, x) for x in STATES])
    @printf("n=%d (%d unknowns) %-28s %8.1f s  PI %2d  Schwarz %5d  conv %s  maxrss %.1f GB  V = %s\n",
        n, row.unknowns, method, t, row.policy_iterations, row.schwarz_iterations, row.converged,
        row.maxrss_gb, join((@sprintf("%.4f", v) for v in row.values), ", "))
    flush(stdout)
    return row
end

blocks(n) = fill(max(2, round(Int, n / 7)), 4)        # about 7 points per box per dimension

run(9, "warm-up", DirectPolicyEvaluation())
run(9, "warm-up", SchwarzPolicyEvaluation(blocks(9)))
run(9, "warm-up", SchwarzPolicyEvaluation(blocks(9); algorithm = SchwarzGMRES(sweep = ParallelSweep(), tol = 1e-10, maxiter = 500)))

rows = HJBRow[]
for n in SIZES
    b = blocks(n)

    push!(rows, run(n, "Schwarz multicolor $(b[1])⁴", SchwarzPolicyEvaluation(b)))
    push!(rows, run(n, "Schwarz GMRES $(b[1])⁴", SchwarzPolicyEvaluation(b;
        algorithm = SchwarzGMRES(sweep = ParallelSweep(), tol = 1e-10, maxiter = 500))))
    n <= DIRECT_MAX && push!(rows, run(n, "direct LU", DirectPolicyEvaluation()))
    open(joinpath(OUT, "benchmarks.csv"), "w") do io
        println(io, "n,unknowns,method,seconds,policy_iterations,schwarz_iterations,converged,peak_rss_so_far_gb,V1,V2,V3")
        for r in rows
            println(io, join((r.n, r.unknowns, repr(r.method), r.seconds, r.policy_iterations,
                r.schwarz_iterations, r.converged, r.maxrss_gb, r.values...), ","))
        end
    end
end
println("wrote ", joinpath(OUT, "benchmarks.csv"))
