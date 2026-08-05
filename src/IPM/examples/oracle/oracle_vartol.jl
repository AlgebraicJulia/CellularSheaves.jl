# Does the ORACLE (optimal-α) solve improve under vartol forcing?
# For each problem/solver, sweep α (fix_alpha) under the current absolute μ-schedule (vartol=false) and
# under Zanetti–Gondzio relative forcing (vartol=true), and compare the summed optimal cost — Σ ncraig
# over the chosen (cheapest-α) rows — plus iteration count and final status.
#
# usage: julia --project=examples examples/oracle/oracle_vartol.jl [key1 key2 ...]

using CellularSheaves
using CellularSheaves.IPM: HSDSettings, IPMSettings, init
import CellularSheaves.IPM as IPM
using LinearAlgebra, SparseArrays

const EX = dirname(@__DIR__)
gp(x) = x isa Tuple ? x[1] : x
for k in ("e01", "e02", "e04", "07", "e08", "e15", "e11", "e13")
    include(joinpath(EX, k == "07" ? "e07.jl" : "$k.jl"))
end
include(joinpath(@__DIR__, "oracle.jl"))   # solve_logged, DEFAULT_ALPHA_GRID

function buildprob(key)
    key == "e01" && return gp(build_merton_chain(merton_instance(; nx=21, olap=5, nsteps=10)))
    key == "e02" && return gp(build_qm_problem(generate_qm_instance(; N=10)))
    key == "e04" && return gp(build_sos_spline(sos_instance(; n=13)))
    key == "07"  && return gp(build_poisson_tv(poisson_instance(; N=128, Tsz=16, m=16, k=-1, K=6, R=12, q=3, seed=3)))
    key == "e11" && (i = equalizer_instance(); return gp(build_equalizer(make_system(i, 8), i; form=:exact, white=true)))
    key == "e13" && return gp(build_landing(landing_instance(), 60.0, 20))
    error("no recipe for $key")
end

function oraclecost(prob, solv, vt, tol0)
    tol = 1e-10
    settings = solv == "ipm" ?
        IPMSettings{Float64}(; feas_tol=tol, gap_tol=tol, itmax=300, fix_alpha=true, vartol=vt, tol0=tol0) :
        HSDSettings{Float64}(; feas_tol=tol, gap_tol=tol, itmax=300, fix_alpha=true, vartol=vt, tol0=tol0)
    s0 = init(prob, settings)
    _, records = solve_logged(s0, collect(DEFAULT_ALPHA_GRID))
    chosen = filter(r -> hasproperty(r, :chosen) && r.chosen === true, records)
    cost = isempty(chosen) ? 0 : sum(r.ncraig for r in chosen)
    iters = isempty(records) ? 0 : maximum(r.iter for r in records)
    status = isempty(chosen) ? "NONE" : string(last(chosen).ipm_status)
    return (cost=cost, iters=iters, status=status)
end

# short status tag; factorizations = IPM iterations (the dominant cost) reported first, solves second.
short(s) = s == "NUMERICAL_FAILURE" ? "FAIL" : s == "NEAR_OPTIMAL" ? "NOPT" : s == "OPTIMAL" ? "OPT" : s == "CONTINUE" ? "conv" : s
cell(r) = "fac=$(rpad(r.iters,2)) slv=$(rpad(r.cost,3)) $(short(r.status))"

keys = length(ARGS) ≥ 1 ? ARGS : ["e02", "e04", "e11", "e13"]
TOL0S = [1e-3, 1e-2, 1e-1]

println("fac = factorizations (= IPM iters, DOMINANT cost); slv = Σ F⁻¹ applications (cleanup)")
println("-"^118)
for key in keys
    prob = buildprob(key)
    for solv in ("ipm", "hsd")
        a = oraclecost(prob, solv, false, 1e-3)
        print(rpad("$key $solv", 10), "abs[", cell(a), "]")
        for t0 in TOL0S
            b = oraclecost(prob, solv, true, t0)
            print("  var@", t0, "[", cell(b), "]")
        end
        println()
    end
end
