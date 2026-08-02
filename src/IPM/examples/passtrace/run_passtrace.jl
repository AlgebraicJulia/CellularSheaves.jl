# Per-pass refinement instrumentation driver (E-series "targeted rerun" spec). One (problem, tol, solver)
# per process, run TWICE — run1 (production: refine_stall=0.5, refine_itmax=10) and run2 (stall trigger
# inert: refine_stall=Inf so it NEVER fires, refine_itmax=20 cap). Same coarse α oracle sweep as the
# existing oracle CSVs (solve_logged on DEFAULT_ALPHA_GRID); the schema is the oracle superset plus the
# new per-pass columns: p_/c_/w_ × {dres,pres,tres,nkry}1..6, plus tau_stall / refine_cap / normB2.
# Both runs share the SAME builder as run_oracle.jl / run_fine2.jl.
#
# Usage:  julia --project=examples examples/passtrace/run_passtrace.jl <key> <tol> <hsd|ipm> [dial]
key = ARGS[1]; tol = parse(Float64, ARGS[2]); solv = length(ARGS) ≥ 3 ? ARGS[3] : "hsd"
dial = length(ARGS) ≥ 4 ? parse(Float64, ARGS[4]) : nothing
using CellularSheaves
using CellularSheaves.IPM: HSDSettings, IPMSettings, init, solve_logged, write_oracle_csv
import CellularSheaves.IPM as IPM
using LinearAlgebra, SparseArrays, Printf
const EX  = dirname(@__DIR__)
const OUT = joinpath(EX, "passtrace"); mkpath(OUT)
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
    key == "X06"  && return build_corner_soc(; corner_r = dial === nothing ? 0.1 : dial)
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

mkset(; kw...) = solv == "ipm" ?
    IPMSettings{Float64}(; feas_tol=tol, gap_tol=tol, itmax=300, kw...) :
    HSDSettings{Float64}(; feas_tol=tol, gap_tol=tol, itmax=300, kw...)

prob = gp(buildprob(key))
grid = collect(IPM.DEFAULT_ALPHA_GRID)
isparam = startswith(key, "SOS_P")
stub = isparam ? "$(key)_$(solv)" : "$(key)_$(ARGS[2])_$(solv)"

# run1 = production defaults; run2 = stall trigger strictly inert (never fires) + higher pass cap.
RUNS = (("run1", (;)),                                          # defaults: refine_stall=0.5, refine_itmax=10
        ("run2", (; refine_stall = Inf, refine_itmax = 20)))    # inert stall + cap 20

jval(v::AbstractString) = "\"$v\""; jval(v::Bool) = string(v)
jval(v::Real) = isfinite(v) ? string(v) : "null"; jval(v::AbstractVector) = "[" * join(jval.(v), ",") * "]"
jval(v) = "\"$(string(v))\""
gitshort = try readchomp(`git -C $EX rev-parse --short HEAD`) catch; "unknown" end

for (label, kw) in RUNS
    settings = mkset(; kw...)
    s0 = init(prob, settings)
    t0 = time()
    final, records = solve_logged(s0, grid)
    wall = round(time() - t0, digits=1)
    niter = isempty(records) ? 0 : maximum(r.iter for r in records)
    chosen = filter(r -> r.chosen, records)
    fstatus = isempty(chosen) ? "NONE" : string(last(chosen).ipm_status)
    # how far the per-pass trace was actually exercised: max ppass and count of rows that hit the cap.
    ppmax = isempty(records) ? 0 : maximum(max(getproperty(r, :ppass), getproperty(r, :cpass),
                hasproperty(r, :wpass) && r.wpass !== nothing ? r.wpass : 0) for r in records)
    capped = count(r -> string(r.pstat) == "REFINE_ITMAX" || string(r.cstat) == "REFINE_ITMAX", records)

    base = "$(stub)_$(label)"
    write_oracle_csv(joinpath(OUT, "$base.csv"), records)

    setget(k) = hasproperty(settings, k) ? getproperty(settings, k) : nothing
    nB = norm(s0.B)
    meta = ["problem"=>key, "solver"=>solv, "tol"=>ARGS[2], "run"=>label, "git"=>gitshort,
        "refine_stall"=>setget(:refine_stall), "refine_itmax"=>setget(:refine_itmax),
        "passtrace_len"=>Int(IPM.PASSTRACE_LEN), "grid"=>collect(grid),
        "raug"=>setget(:raug), "aaug"=>setget(:aaug), "alpha_anchor"=>s0.α[], "nB"=>nB, "normB2"=>nB^2,
        "feas_tol"=>setget(:feas_tol), "gap_tol"=>setget(:gap_tol), "forcing_frac"=>setget(:forcing_frac),
        "m"=>size(s0.B,1), "n"=>size(s0.B,2), "nu"=>s0.ν,
        "final_status"=>fstatus, "niter"=>niter, "nrows"=>length(records),
        "max_ppass"=>ppmax, "n_itmax_rows"=>capped, "wall_s"=>wall,
        "schema_note"=>"oracle superset + per-pass trace. New: <role>_dres/pres/tres/nkry 1..6 (role=p/c/w). dres/pres/tres = ENTRY residual (dual/primal/τ split) of refinement pass i, at the top of refinement iter i, BEFORE that pass's CRAIG; NaN = loop never reached pass i; tres NaN on IPM 2-row. nkry = CRAIG iters of pass i (=solve return−1); −1 = pass i did not fire. Base invocation (pass 0): entry=r0_{p,c,w}, CRAIG=pbase−1, post-base residual=pass-1 entry (p_dres1/p_pres1). tau_stall/refine_cap = stall threshold & pass cap in force. normB2 = ordinary ‖B‖² (NOT sigma2max, which is the preconditioned B F⁻¹ Bᵀ Ritz σ̂²max). r1_{p,c,w} = post-CRAIG base residual (equals r0 when base CRAIG did nothing).",
    ]
    open(joinpath(OUT, "$base.meta.json"), "w") do io
        println(io, "{" * join(("\"$k\":$(jval(v))" for (k, v) in meta), ",") * "}")
    end
    open(joinpath(OUT, "_runs.tsv"), "a") do io
        @printf(io, "%s\t%s\t%s\t%s\t%d\t%d\t%d\t%d\t%s\t%.1f\n",
            key, solv, label, fstatus, niter, length(records), ppmax, capped,
            string(setget(:refine_stall)), wall)
    end
    println("wrote $base: $(length(records)) rows, $niter iters, maxppass=$ppmax capped=$capped status=$fstatus ($(wall)s)")
end
