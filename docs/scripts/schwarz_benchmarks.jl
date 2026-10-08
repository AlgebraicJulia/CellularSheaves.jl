# Benchmarks for the Schwarz solvers against a direct solve.
#
# Poisson problems from SchwarzModelProblems (unit square at several sizes and
# the notched rectangle), decomposed into boxes of 32 × 32 grid points with
# overlap 2. For each problem we compare
#
#   - direct sparse solves (CliqueTrees ChordalLDLt and CHOLMOD),
#   - the Galerkin tower used alone: one coarse solve on the pushforward
#     (Nicolaides) coarse space, i.e. the hierarchical approximation,
#   - Schwarz iterations: one-level Dirichlet, two-level (tower + sweeps),
#     optimized Robin (multicolor and parallel),
#   - sheaf ADMM (Hanks, Riess et al.) at ρ = p*,
#   - two-level Schwarz-preconditioned CG,
#
# reporting setup and solve wall time, iterations, error against the direct
# solution, neighbour-exchange rounds and the largest local problem.
#
# Run with:  julia -t 8 --project=docs docs/scripts/schwarz_benchmarks.jl
# Results are written to docs/figures/schwarz/benchmarks.csv and timing.svg.

get!(ENV, "GKSwstype", "100")

using CellularSheaves
using CliqueTrees.Multifrontal
using LinearAlgebra
using SparseArrays
using Plots
using Printf

const OUT = joinpath(@__DIR__, "..", "figures", "schwarz")
mkpath(OUT)
const TOL = 1e-8
const BOX = 32

# A benchmark problem: a grid domain cut into px × py boxes.
struct BenchmarkCase
    name::String
    domain::GridDomain
    px::Int
    py::Int
    slow_methods::Bool          # run one-level Schwarz and ADMM (slow on large cases)
end

# One measured method on one case.
struct BenchmarkRow
    case::String
    dofs::Int
    subdomains::Int
    method::String
    setup_s::Float64
    solve_s::Float64
    iterations::Int
    rel_error::Float64
    rounds::Int
    max_local_dofs::Int
end

function ldlt_solve(M, b)
    c = M.P' \ b
    z = M.L \ c
    w = M.D \ z
    y = M.L' \ w
    return M.P \ y
end

# Neighbour-exchange rounds per iteration: a multicolor sweep needs one round
# per color, a parallel sweep or an ADMM step one round, a coarse correction or
# a CG step one more (global) round.
rounds(::MulticolorSweep, dd) = length(dd.colors)
rounds(::ParallelSweep, dd) = 1

# Wall time and result of `f()`. Runs shorter than two seconds are repeated
# three times and the fastest is kept, to suppress timer and GC noise.
function timed(f)
    GC.gc()
    t = @elapsed result = f()
    t < 2 || return t, result
    for _ in 1:2
        GC.gc()
        t = min(t, @elapsed f())
    end
    return t, result
end

function measure(case::BenchmarkCase)
    dom = case.domain
    A = poisson_matrix(dom)
    n = size(A, 1)
    f = ones(n)
    parts = box_partition(dom, case.px, case.py)
    rows = BenchmarkRow[]
    row(method, setup, solve, its, u, nrounds, local_dofs) =
        push!(rows, BenchmarkRow(case.name, n, maximum(parts), method, setup, solve, its,
            norm(u - u_ref) / norm(u_ref), nrounds, local_dofs))

    t_setup, F = timed(() -> ldlt!(ChordalLDLt(A), RowMaximum()))
    t_solve, u_ref = timed(() -> ldlt_solve(F, f))
    row("direct (ChordalLDLt)", t_setup, t_solve, 0, u_ref, 0, n)
    t_setup, C = timed(() -> cholesky(A))
    t_solve, u_chol = timed(() -> C \ f)
    row("direct (CHOLMOD)", t_setup, t_solve, 0, u_chol, 0, n)

    h = dom.h
    pstar = optimized_robin_parameter(5h) / h                 # overlap 2: L = 5h
    doms = overlapping_subdomains(A, parts; overlap=2)
    t_dd, dd = timed(() -> SchwarzDecomposition(A, doms; owner=parts))
    t_robin, ddr = timed(() -> SchwarzDecomposition(A, doms; owner=parts, transmission=RobinTransmission(pstar)))
    t_coarse, coarse = timed(() -> TruncatedPushforwardCoarseSpace(dd))
    local_dofs = maximum(length, dd.cover.stalks)
    prob = SchwarzProblem(dd, f)

    t_solve, u_tower = timed(() -> coarse.basis * ldlt_solve(coarse.factor, coarse.basis' * f))
    row("Galerkin tower alone (one coarse solve)", t_coarse, t_solve, 1, u_tower, 1, coarse_dimension(coarse))

    methods = [
            ("Schwarz, one-level (multicolor)", prob, t_dd, SchwarzIteration(sweep=MulticolorSweep(), tol=TOL, maxiter=20_000)),
            ("Schwarz, two-level: tower + sweeps", prob, t_dd + t_coarse,
                SchwarzIteration(sweep=MulticolorSweep(), coarse=coarse, tol=TOL, maxiter=20_000)),
            ("Robin Schwarz p* (multicolor)", SchwarzProblem(ddr, f), t_robin, SchwarzIteration(sweep=MulticolorSweep(), tol=TOL, maxiter=20_000)),
            ("Robin Schwarz p* (parallel)", SchwarzProblem(ddr, f), t_robin, SchwarzIteration(sweep=ParallelSweep(), tol=TOL, maxiter=20_000)),
            ("Robin Schwarz p* + tower (parallel)", SchwarzProblem(ddr, f), t_robin + t_coarse,
                SchwarzIteration(sweep=ParallelSweep(), coarse=coarse, tol=TOL, maxiter=20_000))]
    case.slow_methods || filter!(m -> !startswith(first(m), "Schwarz, one-level"), methods)
    for (method, problem, setup, alg) in methods
        t, r = timed(() -> solve(problem, alg))
        nrounds = r.iterations * (rounds(alg.sweep, problem.decomposition) + (alg.coarse === nothing ? 0 : 1))
        row(method, setup, t, r.converged ? r.iterations : -1, r.u, nrounds, local_dofs)
    end

    if case.slow_methods
        t, r = timed(() -> solve(prob, SheafADMM(rho=pstar, tol=TOL, maxiter=20_000)))
        row("sheaf ADMM, ρ = p*", t_dd, t, r.converged ? r.iterations : -1, r.u, r.iterations, local_dofs)
    end

    t, r = timed(() -> solve(prob, SchwarzCG(coarse=coarse, tol=TOL, maxiter=5_000)))
    row("two-level Schwarz CG", t_dd + t_coarse, t, r.converged ? r.iterations : -1, r.u, 2r.iterations, local_dofs)
    return rows
end

cases = [
    BenchmarkCase("square 64²", unit_square(63), 2, 2, true),
    BenchmarkCase("square 128²", unit_square(127), 4, 4, true),
    BenchmarkCase("square 256²", unit_square(255), 8, 8, true),
    BenchmarkCase("square 512²", unit_square(511), 16, 16, false),
    BenchmarkCase("notched 128×64", notched_rectangle(63), 4, 2, true),
    BenchmarkCase("notched 256×128", notched_rectangle(127), 8, 4, true),
]

measure(BenchmarkCase("warm-up", unit_square(31), 2, 2, true))      # compile everything first
rows = reduce(vcat, [measure(c) for c in cases])

open(joinpath(OUT, "benchmarks.csv"), "w") do io
    println(io, "case,dofs,subdomains,method,setup_s,solve_s,iterations,rel_error,rounds,max_local_dofs")
    for r in rows
        println(io, join((r.case, r.dofs, r.subdomains, r.method, r.setup_s, r.solve_s, r.iterations,
            r.rel_error, r.rounds, r.max_local_dofs), ","))
    end
end

for c in cases
    selected = filter(r -> r.case == c.name, rows)
    first_row = first(selected)
    direct = first_row.setup_s + first_row.solve_s
    @printf("\n%s: %d dofs, %d subdomains, threads = %d\n", c.name, first_row.dofs, first_row.subdomains, Threads.nthreads())
    @printf("%-40s %9s %9s %9s %6s %9s %7s %8s\n", "method", "setup s", "solve s", "× direct", "iters", "error", "rounds", "local")
    for r in selected
        r.iterations < 0 && (@printf("%-40s %9.4f %9.4f %9s %6s %9s %7d %8d\n", r.method, r.setup_s, r.solve_s, "-", "diverged", "-", r.rounds, r.max_local_dofs); continue)
        @printf("%-40s %9.4f %9.4f %9.1f %6d %9.1e %7d %8d\n", r.method, r.setup_s, r.solve_s,
            (r.setup_s + r.solve_s) / direct, r.iterations, r.rel_error, r.rounds, r.max_local_dofs)
    end
end

squares = filter(c -> startswith(c.name, "square"), cases)
plt = plot(; xscale=:log10, yscale=:log10, xlabel="dofs", ylabel="setup + solve time (s)",
    legend=:outerright, size=(950, 450), title="Poisson on the unit square, 32×32 boxes, overlap 2")
for method in unique(r.method for r in rows)
    selected = [r for r in rows if r.method == method && startswith(r.case, "square") && r.iterations >= 0]
    length(selected) >= 3 || continue
    plot!(plt, [r.dofs for r in selected], [r.setup_s + r.solve_s for r in selected];
        label=method, marker=:circle, lw=2)
end
savefig(plt, joinpath(OUT, "timing.svg"))
println("\nwrote ", joinpath(OUT, "benchmarks.csv"))
