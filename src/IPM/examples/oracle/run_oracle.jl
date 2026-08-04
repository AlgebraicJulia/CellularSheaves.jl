# Per-iteration α oracle driver. Runs solve_logged on one (problem, tol, solver) and dumps the full
# per-(iter, α) record to oracle/<key>_<tol>_<solver>.csv + a small meta.json sidecar.
#
# Usage:  julia --project=examples examples/oracle/run_oracle.jl <key> <tol> <ipm|hsd>
#   e.g.  julia --project=examples examples/oracle/run_oracle.jl e01 1e-10 ipm
#
# Parsimonious buildprob: only the starter fleet used this campaign. Add a line per new key.
key   = ARGS[1]
tol   = parse(Float64, ARGS[2])
solvs = length(ARGS) ≥ 3 ? (ARGS[3] == "both" ? ["ipm", "hsd"] : [ARGS[3]]) : ["ipm"]

using CellularSheaves
using CellularSheaves.IPM: HSDSettings, IPMSettings, init
import CellularSheaves.IPM as IPM

const EX  = dirname(@__DIR__)                                   # examples/
const OUT = @__DIR__                                            # examples/oracle/
gp(x) = x isa Tuple ? x[1] : x

include(joinpath(@__DIR__, "oracle.jl"))

incfile = key == "07" ? "$EX/e07.jl" : "$EX/$key.jl"
include(incfile)

function buildprob(key)
    key == "e01" && return gp(build_merton_chain(merton_instance(; nx=21, olap=5, nsteps=10)))
    key == "e02" && return gp(build_qm_problem(generate_qm_instance(; N=10)))
    key == "e04" && return gp(build_sos_spline(sos_instance(; n=13)))
    key == "07"  && return gp(build_poisson_tv(poisson_instance(; N=128, Tsz=16, m=16, k=-1, K=6, R=12, q=3, seed=3)))
    if key == "e08"
        dt = 4.0 / 64; inst = exec_instance(; n=1, T=64)
        return gp(build_exec(inst; Σ=fill(dt, 1, 1), X0=[1.5], η=[1.0 / sqrt(dt)], γ=1.0))
    end
    if key == "e15"
        train = [[0.0, 60.0 * (v - 1), 0.0, 0.0, 0.0, 0.0] for v in 1:V_BASE]
        return gp(build_cw(train, [zeros(6) for _ in 1:V_BASE], N_BASE; gam=1e-1, kap=1e-1))
    end
    error("no recipe for $key")
end

using LinearAlgebra, SparseArrays
prob = buildprob(key)                                          # solver-agnostic; built once
grid = collect(DEFAULT_ALPHA_GRID)
gitshort = try readchomp(`git -C $EX rev-parse --short HEAD`) catch; "unknown" end
jval(v::Real) = isfinite(v) ? string(v) : "null"
jval(v) = "\"$(string(v))\""

for solv in solvs
    settings = solv == "ipm" ?
        IPMSettings{Float64}(; feas_tol=tol, gap_tol=tol, itmax=300, fix_alpha=true) :
        HSDSettings{Float64}(; feas_tol=tol, gap_tol=tol, itmax=300, fix_alpha=true)
    s0 = init(prob, settings)
    final, records = solve_logged(s0, grid)

    niter = isempty(records) ? 0 : maximum(r.iter for r in records)
    chosenrows = filter(r -> r.chosen, records)
    final_status = isempty(chosenrows) ? "NONE" : string(last(chosenrows).ipm_status)
    base = "$(key)_$(ARGS[2])_$(solv)"
    write_oracle_csv(joinpath(OUT, "$base.csv"), records)

    setget(k) = hasproperty(settings, k) ? getproperty(settings, k) : nothing
    meta = [
        "problem"=>key, "solver"=>solv, "tol"=>ARGS[2], "git"=>gitshort,
        "m"=>size(s0.B, 1), "n"=>size(s0.B, 2), "normB"=>norm(s0.B),
        "rgmin"=>setget(:rgmin), "rgmax"=>setget(:rgmax),
        "feas_tol"=>setget(:feas_tol), "gap_tol"=>setget(:gap_tol),
        "refine_itmax"=>setget(:refine_itmax), "refine_stall"=>setget(:refine_stall),
        "final_status"=>final_status, "niter"=>niter, "nrows"=>length(records),
    ]
    open(joinpath(OUT, "$base.meta.json"), "w") do io
        println(io, "{" * join(("\"$k\":$(jval(v))" for (k, v) in meta), ",") * "}")
    end
    println("wrote $(length(records)) records ($niter iters, status=$final_status) -> $base")
end
