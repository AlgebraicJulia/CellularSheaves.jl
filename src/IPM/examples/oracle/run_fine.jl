# Fine-grid floor sweep driver (fine1 spec). Runs solve_logged_fine on a small fixed set of (problem,
# solver) pairs, writing one CSV + meta.json each into examples/fine1/, in the oracle4 schema + a `grid`
# ("coarse"|"fine") column. Computes the acceptance tests (bracketing / determinism / consistency) and
# writes fine1/SUMMARY.md. Auto-widens a problem's fine range to ±1.5 decades and re-runs if bracketing
# fails on > 1/3 of its swept iterations.
#
# Usage:  julia --project=examples examples/oracle/run_fine.jl
using CellularSheaves
using CellularSheaves.IPM: HSDSettings, IPMSettings, init, solve_logged_fine, write_oracle_csv
import CellularSheaves.IPM as IPM
using LinearAlgebra, SparseArrays, Printf

const EX  = dirname(@__DIR__)                                   # examples/
const OUT = joinpath(EX, "fine1")
mkpath(OUT)
gp(x) = x isa Tuple ? x[1] : x

# ---- problem recipes (subset of run_oracle.jl, only the keys this sweep needs) ----
include("$EX/e01.jl"); include("$EX/e15.jl"); include("$EX/e06.jl"); include("$EX/adversarial/X08.jl")
function buildprob(key)
    key == "e01"  && return gp(build_merton_chain(merton_instance(; nx=21, olap=5, nsteps=10)))
    key == "06"   && return gp(build_sqrtloss(sqrtloss_instance(; P=12)))
    key == "X08b" && return build_nonstrict_twin()
    if key == "e15"
        train = [[0.0, 60.0 * (v - 1), 0.0, 0.0, 0.0, 0.0] for v in 1:V_BASE]
        return gp(build_cw(train, [zeros(6) for _ in 1:V_BASE], N_BASE; gam=1e-1, kap=1e-1))
    end
    error("no recipe for $key")
end

# ---- the pair list (≤5; resolution, not fleet size) ----
const PAIRS = [("e01","ipm"), ("e15","ipm"), ("X08b","ipm"), ("e01","hsd"), ("06","ipm")]
const TOL   = 1e-10

mkset(solv) = solv == "ipm" ?
    IPMSettings{Float64}(; feas_tol=TOL, gap_tol=TOL, itmax=300) :
    HSDSettings{Float64}(; feas_tol=TOL, gap_tol=TOL, itmax=300)

# ---- meta sidecar (mirrors run_oracle.jl) ----
jval(v::AbstractString) = "\"$v\""
jval(v::Bool) = string(v)
jval(v::Real) = isfinite(v) ? string(v) : "null"
jval(v::AbstractVector) = "[" * join(jval.(v), ",") * "]"
jval(v) = "\"$(string(v))\""
function writemeta(path, key, solv, s0, final, settings, grid, hw, fstep, final_status, niter)
    setget(k) = hasproperty(settings, k) ? getproperty(settings, k) : nothing
    nH = sqrt(size(s0.B, 2)); nQ = norm(Symmetric(s0.Q, :L)); nB = norm(s0.B)
    cones = Dict{String,Int}()
    for k in s0.K; t = string(nameof(typeof(k))); cones[t] = get(cones, t, 0) + 1; end
    gitshort = try readchomp(`git -C $EX rev-parse --short HEAD`) catch; "unknown" end
    meta = [
        "problem"=>key, "solver"=>solv, "tol"=>string(TOL), "git"=>gitshort,
        "grid"=>collect(grid), "fine_halfwidth"=>hw, "fine_step"=>fstep,
        "raug"=>setget(:raug), "aaug"=>setget(:aaug), "alpha_anchor"=>s0.α[],
        "nH"=>nH, "nQ"=>nQ, "nB"=>nB, "nc"=>final.nc[], "ng"=>final.ng[], "mu_1"=>first(final.hist.μ),
        "feas_tol"=>setget(:feas_tol), "gap_tol"=>setget(:gap_tol),
        "forcing_frac"=>setget(:forcing_frac), "forcing_ceil"=>setget(:forcing_ceil),
        "refine_itmax"=>setget(:refine_itmax), "refine_stall"=>setget(:refine_stall),
        "rgmin"=>setget(:rgmin), "rgmax"=>setget(:rgmax), "itmax"=>setget(:itmax),
        "m"=>size(s0.B, 1), "n"=>size(s0.B, 2), "nu"=>s0.ν, "cones"=>cones,
        "final_status"=>final_status, "niter"=>niter,
        "grid_note"=>"Two-stage sweep: coarse = half-decade grid (as oracle4); fine = $(fstep)-decade lattice over log α_floor ± $(hw), coincident-with-coarse points skipped. Column `grid`=coarse|fine. Spectral harvest columns unchanged from oracle4.",
    ]
    open(path, "w") do io
        println(io, "{" * join(("\"$k\":$(jval(v))" for (k, v) in meta), ",") * "}")
    end
end

# ---- acceptance metrics from a problem's records + notes ----
function acceptance(records, notes)
    byiter = Dict{Int,Vector}()
    for r in records; push!(get!(byiter, r.iter, []), r); end
    swept = [n for n in notes if n.did_fine]
    # 1. bracketing: BOTH landmarks interior — ≥2 fine rows at plateau (nmin), ≥2 at nmin+1
    #    (full min+1 tread), and ≥2 at nmin+2 (second transition captured).
    brack_ok = 0
    for n in swept
        frows = [r for r in byiter[n.iter] if r.grid == "fine"]
        isempty(frows) && continue
        plat = minimum(r.ncraig for r in frows)
        c0 = count(r -> r.ncraig == plat,     frows)
        c1 = count(r -> r.ncraig == plat + 1, frows)
        c2 = count(r -> r.ncraig >= plat + 2, frows)
        (c0 >= 2 && c1 >= 2 && c2 >= 2) && (brack_ok += 1)
    end
    brackfrac = isempty(swept) ? 1.0 : brack_ok / length(swept)
    # per-iteration bracketing shortfalls (which tread level had <2 fine rows) + max in-window pbase
    shortfalls = String[]; pbmax = 0
    for n in swept
        frows = [r for r in byiter[n.iter] if r.grid == "fine"]
        isempty(frows) && continue
        pbmax = max(pbmax, maximum(r.pbase for r in frows))
        plat = minimum(r.ncraig for r in frows)
        miss = String[]
        count(r -> r.ncraig == plat,     frows) < 2 && push!(miss, "nmin")
        count(r -> r.ncraig == plat + 1, frows) < 2 && push!(miss, "nmin+1")
        count(r -> r.ncraig >= plat + 2, frows) < 2 && push!(miss, "nmin+2")
        isempty(miss) || push!(shortfalls, "it$(n.iter):<2@" * join(miss, "&"))
    end
    # 2. determinism: worst det_reldiff over swept iterations
    detmax = isempty(swept) ? 0.0 : maximum(n.det_reldiff for n in swept)
    # 3. consistency: ncraig non-increasing in α across fine lattice (ignoring the top-edge rise)
    incons = Int[]
    for n in swept
        frows = sort([r for r in byiter[n.iter] if r.grid == "fine"], by = r -> r.alpha)
        length(frows) < 2 && continue
        nc = [r.ncraig for r in frows]
        # drop a single trailing (high-α) rise: find last index where it stops decreasing from the top
        viol = any(nc[k] > nc[k-1] for k in 2:length(nc)-1)   # interior increase (edges allowed to rise)
        viol && push!(incons, n.iter)
    end
    return (; nswept = length(swept), brackfrac, detmax, incons, shortfalls, pbmax)
end

# ---- run one (key, solver), optionally widened ----
function runpair(key, solv, hw)
    settings = mkset(solv); prob = gp(buildprob(key)); s0 = init(prob, settings)
    grid = collect(IPM.DEFAULT_ALPHA_GRID)
    final, records, notes = solve_logged_fine(s0, grid; halfwidth = hw)
    return s0, final, settings, grid, records, notes
end

function main()
summ = IOBuffer()
gitshort = try readchomp(`git -C $EX rev-parse --short HEAD`) catch; "unknown" end
println(summ, "# fine1 — fine-grid floor sweep\n")
println(summ, "git: `$gitshort`  ·  tol: $TOL  ·  coarse: half-decade (37 pts)  ·  fine: 0.1-decade lattice\n")
println(summ, "| problem | solver | iters | swept | fine½w | bracket-frac | det-maxΔ | consistency-exceptions |")
println(summ, "|---|---|---|---|---|---|---|---|")

SKIPLOG = String[]
NOTELOG = String[]
t0 = time()
for (key, solv) in PAIRS
    hw = 1.0
    s0 = final = settings = grid = records = notes = acc = nothing
    for attempt in 1:2
        s0, final, settings, grid, records, notes = runpair(key, solv, hw)
        acc = acceptance(records, notes)
        (acc.brackfrac >= 2/3 || attempt == 2) && break
        hw = 1.5     # widen and re-run once
    end
    tag = "$(key)_$(string(TOL))_$(solv)"
    write_oracle_csv(joinpath(OUT, "$tag.csv"), records)
    niter = isempty(records) ? 0 : maximum(r.iter for r in records)
    chosen = filter(r -> r.chosen, records)
    fstatus = isempty(chosen) ? "NONE" : string(last(chosen).ipm_status)
    writemeta(joinpath(OUT, "$tag.meta.json"), key, solv, s0, final, settings, grid, hw, 0.1, fstatus, niter)
    incs = isempty(acc.incons) ? "none" : join(acc.incons, ",")
    @printf(summ, "| %s | %s | %d | %d | ±%.1f | %.2f | %.1e | %s |\n",
            key, solv, niter, acc.nswept, hw, acc.brackfrac, acc.detmax, incs)
    # skipped-fine notes
    skipped = [(n.iter, n.reason) for n in notes if !n.did_fine]
    println("wrote $tag: $(length(records)) rows, $niter iters, $(acc.nswept) swept, brack=$(round(acc.brackfrac,digits=2)), detmax=$(acc.detmax), skipped=$(length(skipped))")
    for (it, rs) in skipped; push!(SKIPLOG, "$tag iter $it: $rs"); end
    sf = isempty(acc.shortfalls) ? "all 3 treads bracketed" : "shortfalls " * join(acc.shortfalls, ", ")
    push!(NOTELOG, "$tag: max in-window pbase = $(acc.pbmax); $sf")
end
wall = round(time() - t0, digits=1)

println(summ, "\n## Skipped fine passes\n")
for line in SKIPLOG; println(summ, "- $line"); end
isempty(SKIPLOG) && println(summ, "(none)")
println(summ, "\n## Acceptance tests\n")
println(summ, "1. **Bracketing** — fraction of swept iterations with ≥2 fine rows at each of plateau (nmin), nmin+1, and nmin+2, so both transitions and the full min+1 tread are interior to the range (see table). Problems below 2/3 were auto-widened to ±1.5 decades and re-run.")
println(summ, "2. **Determinism** — max relative disagreement between α_floor re-solved on the frozen iterate and its coarse row, over (ncraig, pbase, r0_p). Must be ≤ 1e-12 (see det-maxΔ).")
println(summ, "3. **Consistency** — iterations where ncraig rises in the interior of the fine α-lattice (edges allowed). Reported, not fixed — floor-staircase data.")
println(summ, "\n## Structural notes (data, not analysis)\n")
println(summ, "Max in-window predictor pbase (>1 ⇒ the plateau does genuine base-solve work) and per-iteration bracketing shortfalls (which tread level had <2 fine rows):\n")
for line in NOTELOG; println(summ, "- $line"); end
println(summ, "\nwall-clock: $(wall)s")

open(joinpath(OUT, "SUMMARY.md"), "w") do io; write(io, String(take!(summ))); end
println("\n=== wrote fine1/ (", length(PAIRS), " pairs) in $(wall)s ===")
end

main()
