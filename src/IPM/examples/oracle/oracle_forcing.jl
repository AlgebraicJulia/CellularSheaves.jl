# Oracle A/B across forcing rules: absolute (V0), vartol (V1), EW Choice 2 γ=1 α=φ (C2a), γ=0.9 α=2 (C2b).
# For each problem/solver, sweep α (fix_alpha) under each rule; report factorizations (= IPM iters, the
# dominant cost) and Σ F⁻¹ over the chosen (cheapest-α) rows. NOTE: the oracle picks a different α-path per
# rule, so cross-rule totals carry the α-jitter confound — read the systematic trends, not single numbers.
#   usage: julia --project=examples examples/oracle/oracle_forcing.jl [key1 key2 ...]

using CellularSheaves
using CellularSheaves.IPM: HSDSettings, IPMSettings, init
import CellularSheaves.IPM as IPM
using LinearAlgebra, SparseArrays

const EX = dirname(@__DIR__)
gp(x) = x isa Tuple ? x[1] : x
for k in ("e01", "e02", "e04", "07", "e08", "e15", "e11", "e13")
    include(joinpath(EX, k == "07" ? "e07.jl" : "$k.jl"))
end
include(joinpath(@__DIR__, "oracle.jl"))

function buildprob(key)
    key == "e01" && return gp(build_merton_chain(merton_instance(; nx=21, olap=5, nsteps=10)))
    key == "e02" && return gp(build_qm_problem(generate_qm_instance(; N=10)))
    key == "e04" && return gp(build_sos_spline(sos_instance(; n=13)))
    key == "07"  && return gp(build_poisson_tv(poisson_instance(; N=128, Tsz=16, m=16, k=-1, K=6, R=12, q=3, seed=3)))
    key == "e11" && (i = equalizer_instance(); return gp(build_equalizer(make_system(i, 8), i; form=:exact, white=true)))
    key == "e13" && return gp(build_landing(landing_instance(), 60.0, 20))
    error("no recipe for $key")
end

function run1(prob, solv, forcing)
    tol = 1e-10
    settings = solv == "ipm" ?
        IPMSettings{Float64}(; feas_tol=tol, gap_tol=tol, itmax=300, fix_alpha=true, forcing=forcing) :
        HSDSettings{Float64}(; feas_tol=tol, gap_tol=tol, itmax=300, fix_alpha=true, forcing=forcing)
    _, records = solve_logged(init(prob, settings), collect(DEFAULT_ALPHA_GRID))
    chosen = filter(r -> hasproperty(r, :chosen) && r.chosen === true, records)
    fac = isempty(records) ? 0 : maximum(r.iter for r in records)
    slv = isempty(chosen) ? 0 : sum(r.ncraig for r in chosen)
    st = isempty(chosen) ? "NONE" : string(last(chosen).ipm_status)
    return (fac=fac, slv=slv, st=st)
end

short(s) = s == "NUMERICAL_FAILURE" ? "FAIL" : s == "NEAR_OPTIMAL" ? "NOPT" : s == "OPTIMAL" ? "OPT" : "conv"
cell(r) = "fac=$(rpad(r.fac,2)) slv=$(rpad(r.slv,3)) $(short(r.st))"
const MODES = [(0, "abs "), (1, "V   "), (2, "gN  "), (4, "fix ")]

keys = length(ARGS) ≥ 1 ? ARGS : ["e01", "e02", "e04", "07", "e11", "e13"]
println("fac = factorizations (= IPM iters, dominant); slv = Σ F⁻¹ (cleanup).")
println("abs=absolute μ-schedule  V=vartol  gN=gnorm δ·‖g‖ (δ=0.05)  fix=fixtol(1e-8)")
println("-"^118)
for key in keys
    prob = buildprob(key)
    for solv in ("ipm", "hsd")
        print(rpad("$key $solv", 10))
        for (f, nm) in MODES
            print(nm, "[", cell(run1(prob, solv, f)), "]  ")
        end
        println()
    end
end
