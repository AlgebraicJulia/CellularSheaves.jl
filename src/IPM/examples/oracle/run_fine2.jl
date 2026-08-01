# Two-edge fine-grid sweep driver (fine2 spec). One (problem, solver) per process — same builder as
# run_oracle.jl — running solve_logged_fine2 (brackets BOTH window edges, adaptive per-edge widening).
# Writes fine2/<key>_<tol>_<solver>.csv (+ meta.json) in the oracle4 schema + `grid` + `edge` columns,
# and appends a per-run acceptance line to fine2/_runs.tsv for the aggregator (build_fine2_summary.jl).
#
# Usage:  julia --project=examples examples/oracle/run_fine2.jl <key> <tol> <hsd|ipm> [dial]
key = ARGS[1]; tol = parse(Float64, ARGS[2]); solv = length(ARGS) ≥ 3 ? ARGS[3] : "hsd"
dial = length(ARGS) ≥ 4 ? parse(Float64, ARGS[4]) : nothing
using CellularSheaves
using CellularSheaves.IPM: HSDSettings, IPMSettings, init, solve_logged_fine2, write_oracle_csv
import CellularSheaves.IPM as IPM
using LinearAlgebra, SparseArrays, Printf
const EX  = dirname(@__DIR__)
const OUT = joinpath(EX, "fine2"); mkpath(OUT)
gp(x) = x isa Tuple ? x[1] : x

xbase = startswith(key, "X") ? key[1:min(3, length(key))] : key
incfile = startswith(key, "SOS") ? "$EX/e04.jl" :
          startswith(key, "X") ? "$EX/adversarial/$xbase.jl" :
          key == "06" ? "$EX/e06.jl" : key == "07" ? "$EX/e07.jl" : key == "13" ? "$EX/e13.jl" : "$EX/$key.jl"
include(incfile)

function buildprob(key)
    ms = match(r"^SOS_P(\d+)_n(\d+)_s(\d+)$", key)
    ms !== nothing && return gp(build_sos_spline(sos_instance(; P=parse(Int,ms[1]), n=parse(Int,ms[2]), seed=parse(Int,ms[3]))))
    key == "X03"  && return build_narrow(; spread = dial === nothing ? 8.0 : dial)
    key == "X06b" && return build_corner_soc_twin()
    key == "X08b" && return build_nonstrict_twin()
    key == "X09"  && return build_ceiling_base(; S = dial === nothing ? 1e2 : dial)
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
        dt = 4.0/64; inst = exec_instance(; n=1, T=64)
        return gp(build_exec(inst; Σ=fill(dt,1,1), X0=[1.5], η=[1.0/sqrt(dt)], γ=1.0))
    end
    if key == "e09"
        inst = smoother_instance(); sys = make_system(inst); _, y, _ = simulate(sys, inst, 100)
        return gp(build_smoother(sys, inst, y))
    end
    if key == "e15"
        train = [[0.0, 60.0*(v-1), 0.0, 0.0, 0.0, 0.0] for v in 1:V_BASE]
        return gp(build_cw(train, [zeros(6) for _ in 1:V_BASE], N_BASE; gam=1e-1, kap=1e-1))
    end
    error("no recipe for $key")
end

settings = solv == "ipm" ?
    IPMSettings{Float64}(; feas_tol=tol, gap_tol=tol, itmax=300) :
    HSDSettings{Float64}(; feas_tol=tol, gap_tol=tol, itmax=300)

prob = gp(buildprob(key)); s0 = init(prob, settings)
grid = collect(IPM.DEFAULT_ALPHA_GRID)
t0 = time()
final, records, notes = solve_logged_fine2(s0, grid)
wall = round(time() - t0, digits=1)

isparam = startswith(key, "SOS_P")
base = isparam ? "$(key)_$(solv)" : "$(key)_$(ARGS[2])_$(solv)"
write_oracle_csv(joinpath(OUT, "$base.csv"), records)

# ---- per-run acceptance metrics ----
byiter = Dict{Int,Vector}(); for r in records; push!(get!(byiter, r.iter, []), r); end
swnotes = [n for n in notes if n.did_fine]
edgefrac(e) = (v = [n.bracketed for n in swnotes if n.edge == e]; isempty(v) ? 1.0 : count(v)/length(v))
floorfrac = edgefrac("floor"); ceilfrac = edgefrac("ceil")
detmax = isempty(swnotes) ? 0.0 : maximum(n.det_reldiff for n in swnotes)
maxw = isempty(swnotes) ? 0.0 : maximum(n.width for n in swnotes)
# consistency (edge-aware): the floor lattice should be non-increasing in α (flag interior RISES);
# the ceiling lattice sits on the rising side, so it should be non-decreasing (flag interior DROPS).
incons = Tuple{Int,String}[]
for n in swnotes
    fr = sort([r for r in byiter[n.iter] if r.grid == "fine" && r.edge == n.edge], by = r -> r.alpha)
    length(fr) < 3 && continue
    nc = [r.ncraig for r in fr]
    bad = n.edge == "ceil" ? any(nc[k] < nc[k-1] for k in 2:length(nc)-1) :
                             any(nc[k] > nc[k-1] for k in 2:length(nc)-1)
    bad && push!(incons, (n.iter, n.edge))
end
inconsstr = [ "$(it):$(ed)" for (it, ed) in incons ]
niter = isempty(records) ? 0 : maximum(r.iter for r in records)
nsw = length(unique(n.iter for n in swnotes))
chosen = filter(r -> r.chosen, records)
fstatus = isempty(chosen) ? "NONE" : string(last(chosen).ipm_status)

# ---- meta sidecar ----
jval(v::AbstractString) = "\"$v\""; jval(v::Bool) = string(v)
jval(v::Real) = isfinite(v) ? string(v) : "null"; jval(v::AbstractVector) = "[" * join(jval.(v), ",") * "]"
jval(v) = "\"$(string(v))\""
setget(k) = hasproperty(settings, k) ? getproperty(settings, k) : nothing
nB = norm(s0.B); gitshort = try readchomp(`git -C $EX rev-parse --short HEAD`) catch; "unknown" end
meta = ["problem"=>key, "solver"=>solv, "tol"=>ARGS[2], "git"=>gitshort, "grid"=>collect(grid),
    "fine_step"=>0.1, "fine_widths"=>[1.0,1.5,2.5,3.5], "edges"=>["floor","ceil"],
    "raug"=>setget(:raug), "aaug"=>setget(:aaug), "alpha_anchor"=>s0.α[], "nB"=>nB,
    "feas_tol"=>setget(:feas_tol), "gap_tol"=>setget(:gap_tol), "forcing_frac"=>setget(:forcing_frac),
    "m"=>size(s0.B,1), "n"=>size(s0.B,2), "nu"=>s0.ν, "final_status"=>fstatus, "niter"=>niter,
    "nswept"=>nsw, "floor_bracket_frac"=>floorfrac, "ceil_bracket_frac"=>ceilfrac,
    "det_maxdiff"=>detmax, "max_width"=>maxw, "consistency_exceptions"=>inconsstr,
    "grid_note"=>"Two-edge fine sweep: coarse half-decade grid + 0.1-decade lattices around BOTH the floor (smallest α, ncraig≤nmin+1) and ceiling (largest α) plateau edges, per-edge adaptive width. Columns grid=coarse|fine, edge=coarse|floor|ceil. Superset of oracle4 schema.",
]
open(joinpath(OUT, "$base.meta.json"), "w") do io
    println(io, "{" * join(("\"$k\":$(jval(v))" for (k, v) in meta), ",") * "}")
end

open(joinpath(OUT, "_runs.tsv"), "a") do io
    @printf(io, "%s\t%s\t%d\t%d\t%.2f\t%.2f\t%.1e\t%.1f\t%s\t%.1f\t%d\n",
        key, solv, niter, nsw, floorfrac, ceilfrac, detmax, maxw,
        isempty(inconsstr) ? "none" : join(inconsstr, ","), wall, length(records))
end
println("wrote $base: $(length(records)) rows, $niter iters, $nsw swept, floor=$(round(floorfrac,digits=2)) ceil=$(round(ceilfrac,digits=2)) detmax=$detmax maxw=$maxw ($(wall)s)")
