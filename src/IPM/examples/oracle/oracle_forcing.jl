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
for k in ("e01","e02","e03","e04","e05","07","e08","e09","e10","e11","e12","e13","e14","e15")
    include(joinpath(EX, k == "07" ? "e07.jl" : "$k.jl"))
end
include(joinpath(@__DIR__, "oracle.jl"))

function buildprob(key)
    key == "e01" && return gp(build_merton_chain(merton_instance(; nx=21, olap=5, nsteps=10)))
    key == "e02" && return gp(build_qm_problem(generate_qm_instance(; N=10)))
    key == "e03" && return gp(build_tensor_problem(generate_tensor_instance(; Mx=2, My=2, n=4, k=2, N=200)))
    key == "e04" && return gp(build_sos_spline(sos_instance(; n=13)))
    key == "e05" && return gp(build_galerkin(galerkin_instance(; N=512)))
    key == "07"  && return gp(build_poisson_tv(poisson_instance(; N=128, Tsz=16, m=16, k=-1, K=6, R=12, q=3, seed=3)))
    if key == "e08"
        dt = 4.0/64
        return gp(build_exec(exec_instance(; n=1, T=64); Σ=fill(dt,1,1), X0=[1.5], η=[1.0/sqrt(dt)], γ=1.0))
    end
    if key == "e09"
        inst = smoother_instance(); sys = make_system(inst); _, y, _ = simulate(sys, inst, 100)
        return gp(build_smoother(sys, inst, y))
    end
    key == "e10" && return gp(build_torus(make_system(tl_instance()), 3))
    key == "e11" && (i = equalizer_instance(); return gp(build_equalizer(make_system(i, 8), i; form=:exact, white=true)))
    key == "e12" && (i = soundzone_instance(); return gp(build_soundzone(make_system(i, 6), i)))
    key == "e13" && return gp(build_landing(landing_instance(), 60.0, 20))
    key == "e14" && return gp(build_hodge(make_system(8)))
    if key == "e15"
        train = [[0.0, 60.0*(v-1), 0.0, 0.0, 0.0, 0.0] for v in 1:V_BASE]
        return gp(build_cw(train, [zeros(6) for _ in 1:V_BASE], N_BASE; gam=1e-1, kap=1e-1))
    end
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
const MODES = [(0, "abs "), (1, "V   "), (2, "gN  "), (3, "vfl "), (4, "fix ")]

keys = length(ARGS) ≥ 1 ? ARGS :
    ["e01","e02","e03","e04","e05","07","e08","e09","e10","e11","e12","e13","e14","e15"]
solvs = split(get(ENV, "FORCING_SOLVS", "ipm"), ",")   # IPM-only by default; set FORCING_SOLVS=ipm,hsd for both
println("fac = factorizations (= IPM iters, dominant); slv = Σ F⁻¹ (cleanup).")
println("abs=absolute μ-schedule  V=vartol  gN=gnorm(fl)  vfl=vfloor  fix=fixtol(feas_tol)  [floored = max(·,feas_tol)]")
println("-"^118)
for key in keys
    prob = try buildprob(key) catch e; println(rpad("$key", 10), "BUILD FAILED ($(typeof(e)))"); nothing end
    prob === nothing && continue
    for solv in solvs
        print(rpad("$key $solv", 10))
        for (f, nm) in MODES
            r = try run1(prob, solv, f) catch e; (fac=0, slv=0, st="ERR") end
            print(nm, "[", cell(r), "]  ")
        end
        println()
    end
end
