# Shifted-regime snapshot capture (snapshots3 spec). Captures the SAME iterate at three αs each that
# bracket the ρ-shift transition, so the vector explaining the shifted-regime dual residual can be
# identified. Extra meta: rho_engaged, rho_value, norm_dp, norm_dy, r0_*, *_res0_dual, nc, ng.
# No solver changes; capture is a guarded observer; ρ is read, never altered.
#
# Usage:  julia --project=examples examples/snapshots/run_snapshots3.jl
using CellularSheaves
using CellularSheaves.IPM: HSDSettings, IPMSettings, init, snapshot_capture_at
import CellularSheaves.IPM as IPM
using LinearAlgebra, SparseArrays, Printf

const EX  = dirname(@__DIR__)
const OUT = joinpath(EX, "snapshots3"); mkpath(OUT)
gp(x) = x isa Tuple ? x[1] : x
const GIT = try readchomp(`git -C $EX rev-parse HEAD`) catch; "unknown" end

include("$EX/e07.jl"); include("$EX/e15.jl")
function buildprob(key)
    key == "07"  && return gp(build_poisson_tv(poisson_instance(; N=128, Tsz=16, m=16, k=-1, K=6, R=12, q=3, seed=3)))
    key == "e15" && (train = [[0.0, 60.0*(v-1), 0.0, 0.0, 0.0, 0.0] for v in 1:V_BASE];
                     return gp(build_cw(train, [zeros(6) for _ in 1:V_BASE], N_BASE; gam=1e-1, kap=1e-1)))
    error("no recipe for $key")
end
mkset(solv) = solv == "ipm" ?
    IPMSettings{Float64}(; feas_tol=1e-10, gap_tol=1e-10, itmax=300) :
    HSDSettings{Float64}(; feas_tol=1e-10, gap_tol=1e-10, itmax=300)

function write_mtx(path, S::SparseMatrixCSC)
    open(path, "w") do io
        println(io, "%%MatrixMarket matrix coordinate real general")
        println(io, "$(size(S,1)) $(size(S,2)) $(nnz(S))")
        rows = rowvals(S); vals = nonzeros(S)
        for j in 1:size(S,2), k in nzrange(S, j); @printf(io, "%d %d %.17g\n", rows[k], j, vals[k]); end
    end
end
read_mtx(path) = begin
    lines = readlines(path); h = findfirst(l -> !startswith(strip(l), "%") && !isempty(strip(l)), lines)
    m, n, _ = parse.(Int, split(strip(lines[h]))); I = Int[]; J = Int[]; V = Float64[]
    for l in lines[h+1:end]; isempty(strip(l)) && continue; p = split(strip(l))
        push!(I, parse(Int, p[1])); push!(J, parse(Int, p[2])); push!(V, parse(Float64, p[3])); end
    sparse(I, J, V, m, n)
end
function write_vectors(path, cols)
    L = maximum(length(v) for (_, v) in cols)
    open(path, "w") do io
        println(io, join((k for (k, _) in cols), ","))
        for i in 1:L; println(io, join(((i <= length(v) ? @sprintf("%.17g", v[i]) : "") for (_, v) in cols), ",")); end
    end
end
# ρ-aware fidelity: rebuild F = A/α + B'B (+ ρI if rho≠0), backsolve, compare ‖g-Bx0‖ to r0.
function fidelity(A, B, f, g, y0, alpha, r0; rho=0.0)
    β = 1 / alpha
    F = β .* A .+ B' * B
    rho != 0 && (F = F + rho * I)
    rhs = β .* f .+ B' * (g .+ (y0 === nothing ? zero(g) : β .* y0))
    x0 = F \ rhs
    return abs(norm(g .- B * x0) - r0) / max(r0, eps())
end
jstr(v::AbstractString) = "\"$v\""; jstr(v::Bool) = string(v)
jstr(v::Real) = isfinite(v) ? string(v) : "null"; jstr(v) = "\"$(string(v))\""
write_meta(path, d) = open(path, "w") do io
    println(io, "{" * join(("\"$k\":$(jstr(v))" for (k, v) in d), ",") * "}")
end
function rolevecs(snaps, solv)
    solv == "ipm" && return (p = snaps[findfirst(x -> x.role == :p, snaps)],
                             c = snaps[findfirst(x -> x.role == :c, snaps)], w = nothing)
    nts = [x for x in snaps if x.role == :newton]
    return (p = nts[1], c = nts[2], w = snaps[findfirst(x -> x.role == :w, snaps)])
end

# capture one (problem, iter) at one α; determine rho_engaged by whether the ρ-free reconstruction
# reproduces r0_p (≤1%); write the directory. Returns a summary NamedTuple.
function capture_one!(rows, key, solv, K, logα)
    prob = gp(buildprob(key)); s0 = init(prob, mkset(solv))
    α = exp10(logα)
    cap = snapshot_capture_at(s0, K, α)
    tagname = @sprintf("%s_%s_it%d_a%.1f", key, solv, K, logα)
    if cap === nothing
        push!(rows, (; tagname, key, solv, K, logα, ok=false, note="terminated before it$K"))
        println("SKIP $tagname"); return
    end
    rv = rolevecs(cap.snaps, solv); A = rv.p.A; B = rv.p.Bmat; n, m = size(A,1), size(B,1)
    row = cap.row; ρ = row.ρ                      # applied shift (0 when unshifted) — now authoritative
    engaged = ρ > 0
    ρuse = ρ
    # fidelity reconstructs with the SAME F the solver factored: ρ=0 unshifted, ρ>0 shifted.
    fp = fidelity(A, B, rv.p.f, rv.p.g, rv.p.y0, cap.alpha, row.r0_p; rho=ρ)
    fc = fidelity(A, B, rv.c.f, rv.c.g, rv.c.y0, cap.alpha, row.r0_c; rho=ρ)
    r0_w = hasproperty(row, :r0_w) ? row.r0_w : nothing
    fw = solv == "ipm" ? nothing : fidelity(A, B, rv.w.f, rv.w.g, rv.w.y0, cap.alpha, r0_w; rho=ρ)
    dir = joinpath(OUT, tagname); mkpath(dir)
    write_mtx(joinpath(dir, "A.mtx"), A); write_mtx(joinpath(dir, "B.mtx"), B)
    zerom = zeros(m)
    cols = ["fp"=>rv.p.f, "gp"=>rv.p.g, "y0p"=>(rv.p.y0===nothing ? zerom : rv.p.y0),
            "fc"=>rv.c.f, "gc"=>rv.c.g, "y0c"=>(rv.c.y0===nothing ? zerom : rv.c.y0)]
    solv != "ipm" && append!(cols, ["fw"=>rv.w.f, "gw"=>rv.w.g, "y0w"=>(rv.w.y0===nothing ? zerom : rv.w.y0)])
    write_vectors(joinpath(dir, "vectors.csv"), cols)
    s0d = hasproperty(row, :pres0_d) ? row.pres0_d : NaN
    c0d = hasproperty(row, :cres0_d) ? row.cres0_d : NaN
    w0d = hasproperty(row, :wres0_d) ? row.wres0_d : nothing
    meta = Any["problem"=>key, "solver"=>solv, "iter"=>K, "git_commit"=>GIT,
        "alpha_used"=>cap.alpha, "log10_alpha"=>logα, "rho_engaged"=>engaged, "rho_value"=>ρuse,
        "force_tol"=>cap.force_tol, "floor_tol"=>cap.floor_tol, "mu"=>row.μ,
        "nc"=>cap.entry.nc, "ng"=>cap.entry.ng, "norm_dp"=>cap.norm_dp, "norm_dy"=>cap.norm_dy]
    solv != "ipm" && append!(meta, ["tau"=>cap.entry.tau, "kappa"=>cap.entry.kappa])
    append!(meta, ["r0_p"=>row.r0_p, "r0_c"=>row.r0_c])
    solv != "ipm" && push!(meta, "r0_w"=>r0_w)
    append!(meta, ["p_res0_dual"=>s0d, "c_res0_dual"=>c0d])
    solv != "ipm" && push!(meta, "w_res0_dual"=>w0d)
    append!(meta, ["n_rows_A"=>n, "n_rows_B"=>m, "nnz_A"=>nnz(A), "nnz_B"=>nnz(B)])
    write_meta(joinpath(dir, "meta.json"), meta)
    push!(rows, (; tagname, key, solv, K, logα, ok=true, engaged, ρuse, n, m, nnzA=nnz(A), nnzB=nnz(B), fp, fc, fw, note=""))
    println("captured $tagname: ρ_engaged=$engaged ρ=$(round(ρuse,sigdigits=3)) fid p=$(round(fp,sigdigits=2)) c=$(round(fc,sigdigits=2))" * (fw===nothing ? "" : " w=$(round(fw,sigdigits=2))"))
    return
end

# probe a (problem, iteration) for its ρ-transition: scan log α, return (below, above) straddling the flip.
function probe_transition(key, solv, K; lo=8.0, hi=13.0, step=0.1)
    prob = gp(buildprob(key)); s0 = init(prob, mkset(solv))
    below = nothing; above = nothing
    for logα in lo:step:hi
        cap = snapshot_capture_at(s0, K, exp10(logα))
        cap === nothing && continue
        if cap.row.ρ > 0            # shift engaged ⇒ we are at/above the transition
            above = round(logα, digits=2); break
        else
            below = round(logα, digits=2)
        end
    end
    return below, above
end

rows = NamedTuple[]
trans = Dict{String,Tuple}()
t0 = time()
# For each (problem, iteration): probe the actual ρ-transition, then capture the spec's alphas UNION the
# probed transition bracket UNION the deep row — so every solve is bracketed tightly regardless of guess.
for (key, solv, K, spec, deep, lo, hi) in (
        ("07",  "ipm", 4, [10.7, 10.8], 12.0, 10.4, 12.5),
        ("07",  "hsd", 3, Float64[],    12.0,  8.0, 11.0),
        ("e15", "ipm", 3, [11.8, 11.9], 13.0, 11.4, 12.5),
    )
    bl, ab = probe_transition(key, solv, K; lo, hi)
    trans[key * "_" * solv] = (bl, ab)
    println("$key $solv it$K transition: below=$bl above=$ab")
    alphas = sort(unique(vcat(spec, filter(!isnothing, [bl, ab]), deep)))
    for la in alphas; capture_one!(rows, key, solv, K, la); end
end
wall = round(time() - t0, digits=1)

# ---- SUMMARY ----
io = IOBuffer()
println(io, "# snapshots3 — shifted-regime captures\n")
println(io, "git: `$(GIT[1:min(10,end)])`  ·  captured $(count(r->r.ok, rows))/$(length(rows))  ·  wall $(wall)s")
println(io, "Probed ρ-transitions (last-unshifted → first-shifted log α):")
for (k, (bl, ab)) in sort(collect(trans))
    println(io, "- $k : $bl → $ab")
end
println(io, "")
println(io, "One directory per (problem, iteration, α). A/B are α-independent (identical across a problem's three captures); vectors/meta are per-α.\n")
println(io, "| problem | solver | it | log α | ρ_engaged | ρ_value | A (n×n) | B (m×n) | fid p | fid c | fid w |")
println(io, "|---|---|---|---|---|---|---|---|---|---|---|")
for r in rows
    if r.ok
        fw = r.fw === nothing ? "—" : @sprintf("%.1e", r.fw)
        @printf(io, "| %s | %s | %d | %.1f | %s | %.2e | %d×%d | %d×%d | %.1e | %.1e | %s |\n",
            r.key, r.solv, r.K, r.logα, r.engaged, r.ρuse, r.n, r.n, r.m, r.n, r.fp, r.fc, fw)
    else
        println(io, "| $(r.key) | $(r.solv) | $(r.K) | $(r.logα) | | | | | | $(r.note) |")
    end
end
allf = [f for r in rows if r.ok for f in (r.fp, r.fc, r.fw) if f !== nothing]
println(io, "\n**Fidelity: max relative Δ = $(isempty(allf) ? "n/a" : @sprintf("%.2e", maximum(allf)))** (must ≤ 1e-2; ρ-shift applied where ρ_engaged).")
# transition bracketing check
println(io, "\n## Transition bracketing (test 3)\n")
for (k,s,K) in (("07","ipm",4),("07","hsd",3),("e15","ipm",3))
    grp = sort([r for r in rows if r.ok && r.key==k && r.solv==s], by=r->r.logα)
    isempty(grp) && continue
    # the adjacent unshifted→shifted pair actually captured
    li = findlast(r -> !r.engaged, grp); ui = findfirst(r -> r.engaged, grp)
    ok = li !== nothing && ui !== nothing && grp[li].logα < grp[ui].logα
    lo = first(grp); hi = last(grp)
    br = "ρ_engaged $(lo.engaged)@$(lo.logα) → $(hi.engaged)@$(hi.logα)"
    tight = ok ? "tight bracket $(grp[li].logα)(false)→$(grp[ui].logα)(true) ✓" : "⚠ no false→true pair captured"
    println(io, "- $k $s it$K: $br; $tight")
end
println(io, "\n## Notes\n- rho_engaged = (applied shift row.ρ > 0). init_uzw! now reports the ACTUAL shift (0 when the unshifted factorization succeeds; the ladder rung otherwise) instead of rgmin. Fidelity rebuilds F with that same ρ (acceptance test 1). ρ is report-only in the solver, so the change is behavior-preserving.\n- Extra meta fields per spec: rho_engaged, rho_value, norm_dp, norm_dy, r0_{p,c,w}, {p,c,w}_res0_dual, nc, ng.\n- Capture point / A-B-vector conventions identical to snapshots/ and snapshots2/.\n- Existing sweeps (fine2/, oracle*/, snapshots/, snapshots2/) predate the ρ fix and still carry the old rgmin ρ-column convention.")
write(joinpath(OUT, "SUMMARY.md"), String(take!(io)))
println("\n=== wrote snapshots3/ ($(count(r->r.ok, rows))/$(length(rows))) in $(wall)s ===")
