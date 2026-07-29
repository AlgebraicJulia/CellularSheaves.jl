# Per-iteration α oracle driver. Runs solve_logged on one (problem, tol, solver) and dumps the full
# per-(iter, α) record — score + μ/ρ/force_tol/floor_tol/σ²_min + all three solves' res0/res1 — to
# oracle/<key>_<tol>_<solver>.csv (this directory).
#
# Usage:  julia --project=examples examples/oracle/run_oracle.jl <key> <tol> <hsd|ipm>
#   e.g.  julia --project=examples examples/oracle/run_oracle.jl e08 1e-8 hsd
key = ARGS[1]; tol = parse(Float64, ARGS[2]); solv = length(ARGS) ≥ 3 ? ARGS[3] : "hsd"
using CellularSheaves
using CellularSheaves.IPM: HSDSettings, IPMSettings, init, solve_logged, write_oracle_csv
import CellularSheaves.IPM as IPM
const EX = dirname(@__DIR__)          # examples/
const OUT = @__DIR__                  # examples/oracle/
gp(x) = x isa Tuple ? x[1] : x

incfile = key == "X03" ? "$EX/adversarial/X03.jl" : key == "X04" ? "$EX/adversarial/X04.jl" :
          key == "06" ? "$EX/e06.jl" : key == "07" ? "$EX/e07.jl" : key == "13" ? "$EX/e13.jl" : "$EX/$key.jl"
include(incfile)

function buildprob(key)
    key == "X03" && return gp(build_narrow(; npos=12, nfree=12, spread=8.0))
    key == "X04" && return gp(build_degenerate(; n=16, degen=8, rowscale=8.0))
    key == "06"  && return gp(build_sqrtloss(sqrtloss_instance(; P=12)))
    key == "07"  && return gp(build_poisson_tv(poisson_instance(; N=128, Tsz=16, m=16, k=-1, K=6, R=12, q=3, seed=3)))
    key == "13"  && return gp(build_landing(landing_instance(), 60.0, 20))
    key == "e02" && return gp(build_qm_problem(generate_qm_instance(; N=10)))
    key == "e03" && return gp(build_tensor_problem(generate_tensor_instance(; Mx=2, My=2, n=4, k=2, N=200)))
    key == "e05" && return gp(build_galerkin(galerkin_instance(; N=512)))
    key == "e10" && return gp(build_torus(make_system(tl_instance()), 3))
    key == "e11" && (i = equalizer_instance(); return gp(build_equalizer(make_system(i, 8), i; form=:exact, white=true)))
    key == "e12" && (i = soundzone_instance(); return gp(build_soundzone(make_system(i, 6), i)))
    key == "e14" && return gp(build_hodge(make_system(8)))
    key == "e01" && return gp(build_merton_chain(merton_instance(; nx=21, olap=5, nsteps=10)))
    key == "e04" && return gp(build_sos_spline(sos_instance(; n=13)))
    if key == "e08"
        dt = 4.0 / 64; inst = exec_instance(; n=1, T=64)
        return gp(build_exec(inst; Σ=fill(dt, 1, 1), X0=[1.5], η=[1.0 / sqrt(dt)], γ=1.0))
    end
    if key == "e09"
        inst = smoother_instance(); sys = make_system(inst); _, y, _ = simulate(sys, inst, 100)
        return gp(build_smoother(sys, inst, y))
    end
    if key == "e15"
        train = [[0.0, 60.0 * (v - 1), 0.0, 0.0, 0.0, 0.0] for v in 1:V_BASE]
        return gp(build_cw(train, [zeros(6) for _ in 1:V_BASE], N_BASE; gam=1e-1, kap=1e-1))
    end
    error("no recipe for $key")
end

settings = solv == "ipm" ?
    IPMSettings{Float64}(feas_tol=tol, gap_tol=tol, itmax=300) :
    HSDSettings{Float64}(feas_tol=tol, gap_tol=tol, itmax=300)

s0 = init(buildprob(key), settings)
final, records = solve_logged(s0)
path = joinpath(OUT, "$(key)_$(ARGS[2])_$(solv).csv")
write_oracle_csv(path, records)
println("wrote $(length(records)) records ($(maximum(r.iter for r in records)) iters) -> $path")
