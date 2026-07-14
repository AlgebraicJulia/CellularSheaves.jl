# =============================================================================
# utils.jl — shared utilities for IPM examples
# =============================================================================

using LinearAlgebra
using Printf
using Statistics: median
using JuMP
using Clarabel
using Dualization
import MathOptInterface as MOI

using CellularSheaves.IPM
using CellularSheaves.IPM: OPTIMAL, NEAR_OPTIMAL, PositiveCone, SemidefiniteCone
using CellularSheaves.BlockSparseArrays: block, nvtxs, blocksparse, colrange

# -----------------------------------------------------------------------------
# Command-line parsing
# -----------------------------------------------------------------------------

function parse_args(args = ARGS)
    tiny = "--tiny" in args
    quick = "--quick" in args || tiny
    mosek = "--mosek" in args
    tol = 1e-8
    for arg in args
        if startswith(arg, "--tol=")
            tol = parse(Float64, split(arg, "=")[2])
        end
    end
    return (; tiny, quick, mosek, tol)
end

# -----------------------------------------------------------------------------
# Optimizer factories
# -----------------------------------------------------------------------------

function clarabel_opt(; tol = 1e-8)
    optimizer_with_attributes(Clarabel.Optimizer,
        "verbose" => false,
        "tol_feas" => tol,
        "tol_gap_abs" => tol,
        "tol_gap_rel" => tol)
end

using MosekTools

function mosek_opt(; tol = 1e-8)
    optimizer_with_attributes(Mosek.Optimizer,
        "QUIET" => true,
        "MSK_IPAR_NUM_THREADS" => 1,
        "MSK_DPAR_INTPNT_CO_TOL_PFEAS" => tol,
        "MSK_DPAR_INTPNT_CO_TOL_DFEAS" => tol,
        "MSK_DPAR_INTPNT_CO_TOL_REL_GAP" => tol,
        "MSK_DPAR_INTPNT_TOL_PFEAS" => tol,
        "MSK_DPAR_INTPNT_TOL_DFEAS" => tol,
        "MSK_DPAR_INTPNT_TOL_REL_GAP" => tol)
end

# -----------------------------------------------------------------------------
# Problem statistics
# -----------------------------------------------------------------------------

function problem_stats(prob)
    N0 = size(prob.B, 2)
    N1 = size(prob.B, 1)
    block_dims = [size(block(prob.Q, v, v, v), 1) for v in 1:nvtxs(prob.B)]
    med_block = length(block_dims) > 0 ? median(block_dims) : 0.0
    return (N0 = N0, N1 = N1, med_block = med_block)
end

# -----------------------------------------------------------------------------
# Timing utilities
# -----------------------------------------------------------------------------

"Min-over-runs timing with adaptive nruns for fast rows."
function timed_min(f; nruns::Int = 5, nwarmup::Int = 2)
    for _ in 1:nwarmup
        f()
    end
    t1 = @elapsed f()
    n = t1 < 0.1 ? max(nruns, 10) : (t1 > 20.0 ? 1 : nruns)
    n == 1 && return t1
    return min(t1, minimum(@elapsed(f()) for _ in 1:(n - 1)))
end

# -----------------------------------------------------------------------------
# IPM measurement
# -----------------------------------------------------------------------------

function ipm_objective(prob, res)
    0.5 * dot(res.p, prob.Q * res.p) + dot(prob.c, res.p)
end

function measure_ipm(prob, settings; nruns = 5, nwarmup = 2)
    res = nothing
    for _ in 1:nwarmup
        res = solve(prob, settings)
    end
    ok = res.status in (OPTIMAL, NEAR_OPTIMAL)
    ok || return (t = NaN, status = string(res.status), obj = NaN)
    t = timed_min(() -> (res = solve(prob, settings)); nruns)
    return (t = t, status = string(res.status), obj = ipm_objective(prob, res))
end

# -----------------------------------------------------------------------------
# JuMP baseline measurement
# -----------------------------------------------------------------------------

function measure_jump(build_m; nruns = 5, nwarmup = 2)
    m = nothing
    try
        for _ in 1:nwarmup
            m = build_m()
            optimize!(m)
        end
    catch e
        return (t = NaN, status = "ERROR", obj = NaN)
    end
    st = termination_status(m)
    m = build_m()
    t = timed_min(() -> optimize!(m); nruns)
    obj = try objective_value(m) catch; NaN end
    return (t = t, status = string(st), obj = obj)
end

# Format time with color: yellow for ALMOST_OPTIMAL, red for anything else non-optimal
const YELLOW = "\e[33m"
const RED = "\e[31m"
const RESET = "\e[0m"

function fmt_time(m; width=8)
    if !isfinite(m.t)
        return rpad("—", width)
    end
    tstr = @sprintf("%.1fms", m.t * 1000)
    if m.status == "OPTIMAL"
        return lpad(tstr, width)
    elseif m.status == "ALMOST_OPTIMAL"
        return YELLOW * lpad(tstr, width) * RESET
    else
        return RED * lpad(tstr, width) * RESET
    end
end

# -----------------------------------------------------------------------------
# Log-log slope fitting
# -----------------------------------------------------------------------------

function loglog_slope(dofs, ts)
    keep = [i for i in eachindex(ts) if isfinite(ts[i]) && ts[i] > 0]
    length(keep) < 3 && return NaN
    x = log.(Float64.(dofs[keep]))
    y = log.(ts[keep])
    x̄, ȳ = sum(x) / length(x), sum(y) / length(y)
    return sum((x .- x̄) .* (y .- ȳ)) / sum((x .- x̄) .^ 2)
end
