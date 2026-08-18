# =============================================================================
# utils.jl — shared utilities for IPM examples
# =============================================================================

using Printf
using JuMP
using Clarabel
using MosekTools
using BenchmarkTools

using CellularSheaves.IPM: solve, OPTIMAL, NEAR_OPTIMAL

const YELLOW = "\e[33m"
const RED = "\e[31m"
const RESET = "\e[0m"

# -----------------------------------------------------------------------------
# Command-line parsing
# -----------------------------------------------------------------------------

function parse_args(args = ARGS)
    clarabel = "--clarabel" in args
    mosek = "--mosek" in args
    hsd = !("--no-hsd" in args)
    tol = 1e-8

    for arg in args
        if startswith(arg, "--tol=")
            tol = parse(Float64, split(arg, "=")[2])
        end
    end

    return (; clarabel, mosek, hsd, tol)
end

# -----------------------------------------------------------------------------
# Optimizer factories
# -----------------------------------------------------------------------------

function clarabel_opt(; tol = 1e-8)
    optimizer = optimizer_with_attributes(Clarabel.Optimizer,
        "verbose" => false,
        "tol_feas" => tol,
        "tol_gap_abs" => tol,
        "tol_gap_rel" => tol)

    return optimizer
end

function mosek_opt(; tol = 1e-8)
    optimizer = optimizer_with_attributes(Mosek.Optimizer,
        "QUIET" => true,
        "MSK_IPAR_NUM_THREADS" => 1,
        "MSK_DPAR_INTPNT_CO_TOL_PFEAS" => tol,
        "MSK_DPAR_INTPNT_CO_TOL_DFEAS" => tol,
        "MSK_DPAR_INTPNT_CO_TOL_REL_GAP" => tol,
        "MSK_DPAR_INTPNT_TOL_PFEAS" => tol,
        "MSK_DPAR_INTPNT_TOL_DFEAS" => tol,
        "MSK_DPAR_INTPNT_TOL_REL_GAP" => tol)

    return optimizer
end

# -----------------------------------------------------------------------------
# Measurement — @belapsed owns warmup, tuning, and adaptive sampling
# -----------------------------------------------------------------------------

function measure_ipm(prob, settings)
    res = solve(prob, settings)

    if res.status != OPTIMAL && res.status != NEAR_OPTIMAL
        t = NaN
        obj = NaN
    else
        t = @belapsed solve($prob, $settings)
        obj = res.pobj
    end

    status = string(res.status)
    return (t = t, status = status, obj = obj)
end

function measure_jump(build_m)
    m = build_m()

    try
        optimize!(m)
    catch
        return (t = NaN, status = "ERROR", obj = NaN)
    end

    t = @belapsed optimize!($m)

    local obj

    try
        obj = objective_value(m)
    catch
        obj = NaN
    end

    return (t = t, status = string(termination_status(m)), obj = obj)
end

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

