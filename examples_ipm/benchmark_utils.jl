######################################################################
# benchmark_utils.jl
#
# Shared utilities for IPM benchmarks. Provides standardized optimizer
# selection with both commercial (Mosek) and open-source (OSQP/Clarabel)
# options, always running primal and dual forms.
#
# Usage in benchmark files:
#   include("../benchmark_utils.jl")
#   opts = parse_benchmark_args(ARGS)
#   optimizer, dual_optimizer = get_optimizers(; opensource=opts.opensource, cones=...)
#
# Command line:
#   julia bench.jl                    # default: opensource=true (Clarabel/OSQP)
#   julia bench.jl --opensource       # same as above
#   julia bench.jl --mosek            # use Mosek + Mosek-D
######################################################################

using Dualization
using JuMP: optimizer_with_attributes
using Printf
using Statistics: median

# ---- Command-line argument parsing ---------------------------------

"""
    parse_benchmark_args(args = ARGS)

Parse command-line arguments for benchmark scripts.

Returns a NamedTuple with:
- `open::Bool` - use open-source solvers (default: true)
- `mosek::Bool` - use Mosek (overrides open)
- `nruns::Int` - number of benchmark runs (default: 5)
- `nwarmup::Int` - number of warmup runs (default: 2)
- `tol::Float64` - solver tolerance for all solvers (default: 1e-8)

Examples:
    julia bench.jl                    # open=true (Clarabel/OSQP)
    julia bench.jl --open             # same as above
    julia bench.jl --mosek            # use Mosek
    julia bench.jl --nruns=10         # 10 benchmark runs
    julia bench.jl --tol=1e-5         # looser tolerance for speed comparison
"""
function parse_benchmark_args(args = ARGS)
    open = true
    mosek = false
    nruns = 5
    nwarmup = 2
    tol = 1e-8

    for arg in args
        if arg == "--open" || arg == "--open=true"
            open = true
        elseif arg == "--open=false"
            open = false
        elseif arg == "--mosek"
            mosek = true
            open = false
        elseif startswith(arg, "--nruns=")
            nruns = parse(Int, split(arg, "=")[2])
        elseif startswith(arg, "--nwarmup=")
            nwarmup = parse(Int, split(arg, "=")[2])
        elseif startswith(arg, "--tol=")
            tol = parse(Float64, split(arg, "=")[2])
        end
    end

    return (open = open && !mosek, mosek = mosek, nruns = nruns, nwarmup = nwarmup, tol = tol)
end

"""Print benchmark configuration."""
function print_benchmark_config(opts; lp_only::Bool = false)
    solver_name = opts.mosek ? "Mosek" : "Clarabel"
    println("Solver: $solver_name + Dualized")
    println("Runs: $(opts.nruns), Warmup: $(opts.nwarmup), Tolerance: $(opts.tol)")
    println()
end

# ---- Cone detection ------------------------------------------------

"""Detect cone types used in a problem."""
function detect_cones(prob)
    has_soc = false
    has_sdp = false
    has_exp = false
    has_pos = false
    has_cofree = false

    for cone in prob.cones
        name = string(typeof(cone))
        if occursin("SecondOrderCone", name)
            has_soc = true
        elseif occursin("SemidefiniteCone", name)
            has_sdp = true
        elseif occursin("ExponentialCone", name)
            has_exp = true
        elseif occursin("PositiveCone", name)
            has_pos = true
        elseif occursin("CofreeCone", name)
            has_cofree = true
        end
    end

    return (soc = has_soc, sdp = has_sdp, exp = has_exp, pos = has_pos, cofree = has_cofree)
end

"""Check if problem is LP-compatible (only PositiveCone and CofreeCone)."""
is_lp_compatible(cones) = !cones.soc && !cones.sdp && !cones.exp

# ---- Optimizer selection -------------------------------------------

"""
    get_optimizers(opts; lp_only::Bool = false)

Returns (primal_optimizer, dual_optimizer) pair based on parsed options,
with termination tolerances pinned to `opts.tol` so all solvers chase
the same accuracy.

- `opts.mosek=true`: Mosek + Dualization.dual_optimizer(Mosek)
- `opts.open=true`: Clarabel + Dualization.dual_optimizer(Clarabel)

The `lp_only` parameter is kept for API compatibility but ignored —
Clarabel handles both LP/QP and conic problems efficiently.

Requires the appropriate solver package to be loaded in the calling file:
  - MosekTools for Mosek
  - Clarabel for all open-source problems
"""
function get_optimizers(opts; lp_only::Bool = false)
    tol = opts.tol
    if opts.mosek
        primal = optimizer_with_attributes(Mosek.Optimizer,
            "MSK_DPAR_INTPNT_TOL_PFEAS"      => tol,   # LP/QP path
            "MSK_DPAR_INTPNT_TOL_DFEAS"      => tol,
            "MSK_DPAR_INTPNT_TOL_REL_GAP"    => tol,
            "MSK_DPAR_INTPNT_CO_TOL_PFEAS"   => tol,   # conic path
            "MSK_DPAR_INTPNT_CO_TOL_DFEAS"   => tol,
            "MSK_DPAR_INTPNT_CO_TOL_REL_GAP" => tol)
    else
        primal = optimizer_with_attributes(Clarabel.Optimizer,
            "tol_feas"    => tol,
            "tol_gap_abs" => tol,
            "tol_gap_rel" => tol)
    end
    dual = Dualization.dual_optimizer(primal)
    return primal, dual
end

"""Dual wrapper for any optimizer."""
dual_of(opt) = Dualization.dual_optimizer(opt)

# ---- Benchmark timing utilities ------------------------------------

"""Time a function over multiple runs, returning median time in ms."""
function bench_time_ms(f; nwarmup::Int = 2, nruns::Int = 5)
    for _ in 1:nwarmup
        f()
    end
    times = [(@elapsed f()) for _ in 1:nruns]
    return 1000.0 * median(times)
end

"""Time a function over multiple runs, returning minimum time in ms."""
function bench_min_ms(f; nwarmup::Int = 2, nruns::Int = 5)
    for _ in 1:nwarmup
        f()
    end
    times = [(@elapsed f()) for _ in 1:nruns]
    return 1000.0 * minimum(times)
end

"""Time a function over multiple runs, returning mean time in ms."""
function bench_mean_ms(f; nwarmup::Int = 2, nruns::Int = 5)
    for _ in 1:nwarmup
        f()
    end
    total = 0.0
    for _ in 1:nruns
        total += @elapsed f()
    end
    return 1000.0 * total / nruns
end

# ---- Standard benchmark runner -------------------------------------

"""
    run_comparison(label, ipm_solve, jump_build;
                   optimizer, dual_optimizer, nwarmup, nruns, verbose)

Standard benchmark comparison: IPM vs primal optimizer vs dual optimizer.
Returns (t_ipm, t_primal, t_dual) in milliseconds.
"""
function run_comparison(label::String, ipm_solve::Function, jump_build::Function;
                        optimizer, dual_optimizer = nothing,
                        nwarmup::Int = 2, nruns::Int = 5, verbose::Bool = true)
    # IPM timing
    for _ in 1:nwarmup
        ipm_solve()
    end
    t_ipm = 0.0
    for _ in 1:nruns
        t_ipm += @elapsed ipm_solve()
    end
    t_ipm = 1000.0 * t_ipm / nruns

    # Primal timing
    t_primal = Inf
    if optimizer !== nothing
        m = jump_build(optimizer)
        for _ in 1:nwarmup
            optimize!(m)
        end
        t_primal = 0.0
        for _ in 1:nruns
            m = jump_build(optimizer)
            t_primal += @elapsed optimize!(m)
        end
        t_primal = 1000.0 * t_primal / nruns
    end

    # Dual timing
    t_dual = Inf
    if dual_optimizer !== nothing
        m = jump_build(dual_optimizer)
        for _ in 1:nwarmup
            optimize!(m)
        end
        t_dual = 0.0
        for _ in 1:nruns
            m = jump_build(dual_optimizer)
            t_dual += @elapsed optimize!(m)
        end
        t_dual = 1000.0 * t_dual / nruns
    end

    if verbose
        rat_p = t_primal < Inf ? t_primal / t_ipm : NaN
        rat_d = t_dual < Inf ? t_dual / t_ipm : NaN
        @printf("%-20s  %7.2f  %7.2f  %7.2f  %6.2fx  %6.2fx\n",
                label, t_ipm, t_primal, t_dual, rat_p, rat_d)
    end

    return t_ipm, t_primal, t_dual
end
