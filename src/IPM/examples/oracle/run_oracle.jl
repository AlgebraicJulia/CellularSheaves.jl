# Per-iteration α oracle driver. Runs solve_logged on one (problem, tol, solver) and dumps the full
# per-(iter, α) record — score + μ/ρ/force_tol/floor_tol/σ²_min + all three solves' res0/res1 — to
# oracle/<key>_<tol>_<solver>.csv (this directory).
#
# Usage:  julia --project=examples examples/oracle/run_oracle.jl <key> <tol> <hsd|ipm> [dial]
#   e.g.  julia --project=examples examples/oracle/run_oracle.jl e08 1e-8 hsd
#   [dial] (adversarials only): X04/X04b → eps, X03/X03b → spread, X06/X06b → corner_r. CSV names
#   pick up the dial (e.g. X04_1e-8_hsd_eps1e-4.csv) so sweeps don't clobber each other.
key = ARGS[1]; tol = parse(Float64, ARGS[2]); solv = length(ARGS) ≥ 3 ? ARGS[3] : "hsd"
dial = length(ARGS) ≥ 4 ? parse(Float64, ARGS[4]) : nothing
using CellularSheaves
using CellularSheaves.IPM: HSDSettings, IPMSettings, init, solve_logged, write_oracle_csv
import CellularSheaves.IPM as IPM
const EX = dirname(@__DIR__)          # examples/
const OUT = get(ENV, "ORACLE_OUT", @__DIR__)   # examples/oracle/ (override with ORACLE_OUT, e.g. oracle2 for v3)
gp(x) = x isa Tuple ? x[1] : x

xbase = startswith(key, "X") ? key[1:min(3, length(key))] : key      # X04b → X04
incfile = startswith(key, "X") ? "$EX/adversarial/$xbase.jl" :
          key == "06" ? "$EX/e06.jl" : key == "07" ? "$EX/e07.jl" : key == "13" ? "$EX/e13.jl" : "$EX/$key.jl"
include(incfile)

function buildprob(key)
    # X-instances return the raw (prob, meta) tuple so run_oracle can record generator params + severity.
    key == "X03"  && return build_narrow(; spread = dial === nothing ? 8.0 : dial)   # col-spread [I1]
    key == "X03b" && return build_narrow_twin()                                      # benign (spread=0)
    key == "X04"  && return build_degenerate(; eps = dial === nothing ? 1e-4 : dial) # near-dep rows [I2] PRIMARY
    key == "X04b" && return build_degenerate_twin()                                  # benign (independent rows)
    key == "X06"  && return build_corner_soc(; corner_r = dial === nothing ? 0.1 : dial)  # SOC corner
    key == "X06b" && return build_corner_soc_twin()                                  # benign (all interior)
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

using LinearAlgebra, SparseArrays
raw = buildprob(key)
prob = gp(raw)
bmeta = (startswith(key, "X") && raw isa Tuple && length(raw) ≥ 2) ? raw[2] : nothing  # X generator/severity meta
s0 = init(prob, settings)
# v3: uniform 1e0–1e18 half-decade grid (37 pts) for ALL problems — de-censors window ceilings.
grid = collect(CellularSheaves.IPM.DEFAULT_ALPHA_GRID)
final, records = solve_logged(s0, grid)
niter = isempty(records) ? 0 : maximum(r.iter for r in records)
chosenrows = filter(r -> r.chosen, records)
final_status = isempty(chosenrows) ? "NONE" : string(last(chosenrows).ipm_status)
dsuf = dial === nothing ? "" : "_dial$(ARGS[4])"
base = "$(key)_$(ARGS[2])_$(solv)$(dsuf)"
path = joinpath(OUT, "$base.csv")
write_oracle_csv(path, records)

# ---- metadata sidecar (once per run): anchoring, normalization, settings, shape, identity ----
jval(v::AbstractString) = "\"$v\""
jval(v::Bool) = string(v)
jval(v::Real) = isfinite(v) ? string(v) : "null"
jval(v::AbstractVector) = "[" * join(jval.(v), ",") * "]"
jval(v::AbstractDict) = "{" * join(("\"$k\":$(jval(x))" for (k, x) in v), ",") * "}"
jval(v) = "\"$(string(v))\""
setget(k) = hasproperty(settings, k) ? getproperty(settings, k) : nothing
nH = sqrt(size(s0.B, 2)); nQ = norm(Symmetric(s0.Q, :L)); nB = norm(s0.B)
cones = Dict{String,Int}()
for k in s0.K; t = string(nameof(typeof(k))); cones[t] = get(cones, t, 0) + 1; end
gitshort = try readchomp(`git -C $EX rev-parse --short HEAD`) catch; "unknown" end
gitdirty = try !isempty(readchomp(`git -C $EX status --porcelain`)) catch; true end
meta = [
    "problem"=>key, "solver"=>solv, "tol"=>ARGS[2], "git"=>gitshort, "git_dirty"=>gitdirty,
    "grid"=>collect(grid),
    "raug"=>setget(:raug), "aaug"=>setget(:aaug), "alpha_anchor"=>s0.α[],
    "nH"=>nH, "nQ"=>nQ, "nB"=>nB,
    "nc"=>final.nc[], "ng"=>final.ng[], "mu_1"=>first(final.hist.μ),
    "feas_tol"=>setget(:feas_tol), "gap_tol"=>setget(:gap_tol),
    "forcing_frac"=>setget(:forcing_frac), "forcing_ceil"=>setget(:forcing_ceil),
    "refine_itmax"=>setget(:refine_itmax), "refine_stall"=>setget(:refine_stall),
    "floor_patience"=>setget(:floor_patience), "rgmin"=>setget(:rgmin), "rgmax"=>setget(:rgmax),
    "itmax"=>setget(:itmax), "step_frac"=>setget(:step_frac),
    "m"=>size(s0.B, 1), "n"=>size(s0.B, 2), "nu"=>s0.ν, "cones"=>cones,
    "final_status"=>final_status, "niter"=>niter,
]
# X-instances: fold in generator params + severity diagnostics (scalars only; skip pstar/Bd/Qd arrays).
if bmeta !== nothing
    for k in keys(bmeta)
        v = getproperty(bmeta, k)
        v isa Union{Real, Bool, AbstractString} && push!(meta, "x_$(k)"=>v)
    end
end
open(joinpath(OUT, "$base.meta.json"), "w") do io
    println(io, "{" * join(("\"$k\":$(jval(v))" for (k, v) in meta), ",") * "}")
end
println("wrote $(length(records)) records ($(maximum(r.iter for r in records)) iters) + meta -> $base")
