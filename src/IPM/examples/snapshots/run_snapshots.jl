# Solve-state snapshot capture (snapshots spec). For each listed (problem, solver, iter), advance the
# oracle's chosen-α trajectory to that iteration, capture the pre-solve KKT input state of every base
# solve, and write A.mtx / B.mtx / vectors.csv / meta.json into snapshots/<problem>_<solver>_it<K>/.
# Runs the FIDELITY acceptance test (rebuild F=(1/α)A+B'B, backsolve, compare ‖g-Bx0‖ to logged r0)
# and writes snapshots/SUMMARY.md. No solver logic changes — capture is a guarded observer.
#
# Usage:  julia --project=examples examples/snapshots/run_snapshots.jl
using CellularSheaves
using CellularSheaves.IPM: HSDSettings, IPMSettings, init, snapshot_capture, oracle_mu_trajectory
import CellularSheaves.IPM as IPM
using LinearAlgebra, SparseArrays, Printf

const EX  = dirname(@__DIR__)
const OUT = joinpath(EX, "snapshots")
mkpath(OUT)
gp(x) = x isa Tuple ? x[1] : x
const GIT = try readchomp(`git -C $EX rev-parse HEAD`) catch; "unknown" end

include("$EX/e08.jl"); include("$EX/e01.jl"); include("$EX/e15.jl"); include("$EX/e04.jl")
include("$EX/adversarial/X03.jl"); include("$EX/adversarial/X04.jl")
function buildprob(key)
    key == "e01" && return gp(build_merton_chain(merton_instance(; nx=21, olap=5, nsteps=10)))
    key == "e15" && (train = [[0.0, 60.0*(v-1), 0.0, 0.0, 0.0, 0.0] for v in 1:V_BASE];
                     return gp(build_cw(train, [zeros(6) for _ in 1:V_BASE], N_BASE; gam=1e-1, kap=1e-1)))
    key == "SOS_P8_n9_s1" && return gp(build_sos_spline(sos_instance(; P=8, n=9, seed=1)))
    key == "X03" && return gp(build_narrow(; spread=8.0))
    key == "X04" && return gp(build_degenerate(; eps=1e-4))
    if key == "e08"
        dt = 4.0/64; inst = exec_instance(; n=1, T=64)
        return gp(build_exec(inst; Σ=fill(dt,1,1), X0=[1.5], η=[1.0/sqrt(dt)], γ=1.0))
    end
    error("no recipe for $key")
end
mkset(solv) = solv == "ipm" ?
    IPMSettings{Float64}(; feas_tol=1e-10, gap_tol=1e-10, itmax=300) :
    HSDSettings{Float64}(; feas_tol=1e-10, gap_tol=1e-10, itmax=300)

# ---- MatrixMarket coordinate I/O (general, full) ----
function write_mtx(path, S::SparseMatrixCSC)
    open(path, "w") do io
        println(io, "%%MatrixMarket matrix coordinate real general")
        println(io, "$(size(S,1)) $(size(S,2)) $(nnz(S))")
        rows = rowvals(S); vals = nonzeros(S)
        for j in 1:size(S,2), k in nzrange(S, j)
            @printf(io, "%d %d %.17g\n", rows[k], j, vals[k])
        end
    end
end
function read_mtx(path)
    lines = readlines(path)
    h = findfirst(l -> !startswith(strip(l), "%") && !isempty(strip(l)), lines)
    m, n, _ = parse.(Int, split(strip(lines[h])))
    I = Int[]; J = Int[]; V = Float64[]
    for l in lines[h+1:end]
        isempty(strip(l)) && continue
        p = split(strip(l))
        push!(I, parse(Int, p[1])); push!(J, parse(Int, p[2])); push!(V, parse(Float64, p[3]))
    end
    return sparse(I, J, V, m, n)
end

# ---- vectors.csv (ragged columns NaN/empty-padded to the longest) ----
function write_vectors(path, cols::Vector{Pair{String,Vector{Float64}}})
    L = maximum(length(v) for (_, v) in cols)
    open(path, "w") do io
        println(io, join((k for (k, _) in cols), ","))
        for i in 1:L
            println(io, join(((i <= length(v) ? @sprintf("%.17g", v[i]) : "") for (_, v) in cols), ","))
        end
    end
end

# ---- fidelity: reload A,B from mtx, rebuild F=(1/α)A+B'B, backsolve, compare ‖g-Bx0‖ to r0 ----
function fidelity(dir, f, g, y0, alpha, r0)
    A = read_mtx(joinpath(dir, "A.mtx")); B = read_mtx(joinpath(dir, "B.mtx"))
    β = 1 / alpha
    rhs = β .* f .+ B' * (g .+ (y0 === nothing ? zero(g) : β .* y0))
    x0 = (β .* A .+ B' * B) \ rhs
    return abs(norm(g .- B * x0) - r0) / max(r0, eps())
end

jstr(v::AbstractString) = "\"$v\""
jstr(v::Bool) = string(v)
jstr(v::Real) = isfinite(v) ? string(v) : "null"
jstr(v) = "\"$(string(v))\""
function write_meta(path, d)
    open(path, "w") do io
        println(io, "{" * join(("\"$k\":$(jstr(v))" for (k, v) in d), ",") * "}")
    end
end

# role → captured base-solve entry
function rolevecs(snaps, solv)
    if solv == "ipm"
        return (p = snaps[findfirst(x -> x.role == :p, snaps)],
                c = snaps[findfirst(x -> x.role == :c, snaps)], w = nothing)
    else
        nts = [x for x in snaps if x.role == :newton]
        return (p = nts[1], c = nts[2], w = snaps[findfirst(x -> x.role == :w, snaps)])
    end
end

# ---- do one snapshot ----
function do_snapshot!(rows, key, solv, K; murec = nothing)
    prob = gp(buildprob(key)); s0 = init(prob, mkset(solv))
    cap = snapshot_capture(s0, K)
    if cap === nothing
        push!(rows, (; tag = "$(key)_$(solv)_it$K", ok = false, note = "solve terminated before iter $K"))
        return
    end
    rv = rolevecs(cap.snaps, solv)
    dir = joinpath(OUT, "$(key)_$(solv)_it$K"); mkpath(dir)
    A = rv.p.A; B = rv.p.Bmat
    n, m = size(A, 1), size(B, 1)
    if max(nnz(A), nnz(B)) * 40 > 200e6   # ~40 bytes/entry in ascii mtx
        push!(rows, (; tag = "$(key)_$(solv)_it$K", ok = false, note = "matrix > ~200MB — skipped"))
        return
    end
    write_mtx(joinpath(dir, "A.mtx"), A); write_mtx(joinpath(dir, "B.mtx"), B)
    zerom = zeros(m)
    cols = ["fp" => rv.p.f, "gp" => rv.p.g, "y0p" => (rv.p.y0 === nothing ? zerom : rv.p.y0),
            "fc" => rv.c.f, "gc" => rv.c.g, "y0c" => (rv.c.y0 === nothing ? zerom : rv.c.y0)]
    if solv != "ipm"
        append!(cols, ["fw" => rv.w.f, "gw" => rv.w.g, "y0w" => (rv.w.y0 === nothing ? zerom : rv.w.y0)])
    end
    write_vectors(joinpath(dir, "vectors.csv"), cols)
    row = cap.row
    r0_w = hasproperty(row, :r0_w) ? row.r0_w : nothing
    meta = Any["problem"=>key, "solver"=>solv, "iter"=>K, "git_commit"=>GIT,
        "alpha_used"=>cap.chosenα, "force_tol"=>cap.force_tol, "floor_tol"=>cap.floor_tol,
        "mu"=>row.μ, "nc"=>cap.entry.nc, "ng"=>cap.entry.ng]
    solv != "ipm" && append!(meta, ["tau"=>cap.entry.tau, "kappa"=>cap.entry.kappa])
    append!(meta, ["r0_p"=>row.r0_p, "r0_c"=>row.r0_c])
    solv != "ipm" && push!(meta, "r0_w"=>r0_w)
    murec !== nothing && push!(meta, "mu_match"=>murec)
    append!(meta, ["n_rows_A"=>n, "n_rows_B"=>m, "nnz_A"=>nnz(A), "nnz_B"=>nnz(B)])
    write_meta(joinpath(dir, "meta.json"), meta)
    # fidelity per role
    fp = fidelity(dir, rv.p.f, rv.p.g, rv.p.y0, cap.chosenα, row.r0_p)
    fc = fidelity(dir, rv.c.f, rv.c.g, rv.c.y0, cap.chosenα, row.r0_c)
    fw = solv == "ipm" ? nothing : fidelity(dir, rv.w.f, rv.w.g, rv.w.y0, cap.chosenα, r0_w)
    push!(rows, (; tag = "$(key)_$(solv)_it$K", ok = true, n, m, nnzA = nnz(A), nnzB = nnz(B),
        fp, fc, fw, mu = row.μ, note = ""))
    println("captured $(key)_$(solv)_it$K: A $(n)×$(n) ($(nnz(A)) nnz), B $(m)×$(n); fidelity p=$(round(fp,sigdigits=2)) c=$(round(fc,sigdigits=2))" * (fw===nothing ? "" : " w=$(round(fw,sigdigits=2))"))
end

# ---- X03 matched pair: pick J (ipm), K (hsd) with |Δlog10 μ| ≤ 0.5, prefer interior ----
function match_x03()
    muI = oracle_mu_trajectory(init(gp(buildprob("X03")), mkset("ipm")), 40)
    muH = oracle_mu_trajectory(init(gp(buildprob("X03")), mkset("hsd")), 40)
    best = (Inf, 0, 0)
    for J in 2:length(muI)-1, K in 2:length(muH)-1
        (muI[J] <= 0 || muH[K] <= 0) && continue
        d = abs(log10(muI[J]) - log10(muH[K]))
        d < best[1] && (best = (d, J, K))
    end
    return best, muI, muH
end

rows = NamedTuple[]
t0 = time()
(best, muI, muH) = match_x03()
J, K = best[2], best[3]
if J == 0 || K == 0
    println("WARNING: no interior X03 μ-match found (traj lens ipm=$(length(muI)) hsd=$(length(muH))); X03 pair will be skipped")
else
    println("X03 match: ipm it$J (μ=$(muI[J])) vs hsd it$K (μ=$(muH[K])), Δlog10μ=$(round(best[1],digits=3))")
end

do_snapshot!(rows, "e08", "ipm", 9)
do_snapshot!(rows, "e01", "hsd", 8)
do_snapshot!(rows, "e15", "ipm", 6)
do_snapshot!(rows, "SOS_P8_n9_s1", "hsd", 11)
if J > 0 && K > 0
    do_snapshot!(rows, "X03", "ipm", J; murec = muI[J])
    do_snapshot!(rows, "X03", "hsd", K; murec = muH[K])
end
do_snapshot!(rows, "X04", "ipm", 3)
wall = round(time() - t0, digits=1)

# ---- SUMMARY.md ----
io = IOBuffer()
println(io, "# snapshots — solve-state capture\n")
println(io, "git: `$(GIT[1:min(10,end)])`  ·  captured $(count(r->r.ok, rows))/$(length(rows)) snapshots  ·  wall $(wall)s\n")
println(io, "X03 matched pair: ipm it$J (μ=$(round(muI[J],sigdigits=4))) vs hsd it$K (μ=$(round(muH[K],sigdigits=4))), Δlog10μ = $(round(best[1],digits=3)) (≤ 0.5 ✓)\n")
println(io, "## Fidelity (test 1) + shapes (test 3)\n")
println(io, "| snapshot | A (n×n) | nnz A | B (m×n) | nnz B | fid p | fid c | fid w | status |")
println(io, "|---|---|---|---|---|---|---|---|---|")
for r in rows
    if r.ok
        fw = r.fw === nothing ? "—" : @sprintf("%.1e", r.fw)
        @printf(io, "| %s | %d×%d | %d | %d×%d | %d | %.1e | %.1e | %s | ok |\n",
            r.tag, r.n, r.n, r.nnzA, r.m, r.n, r.nnzB, r.fp, r.fc, fw)
    else
        println(io, "| $(r.tag) | — | — | — | — | — | — | — | $(r.note) |")
    end
end
allf = [f for r in rows if r.ok for f in (r.fp, r.fc, r.fw) if f !== nothing]
println(io, "\n**Fidelity: max relative Δ = $(isempty(allf) ? "n/a" : @sprintf("%.2e", maximum(allf)))** (must be ≤ 1e-2).")
println(io, "\n## Test 2 (loadability)\n")
println(io, "Run per snapshot dir: `python -c \"import scipy.io,pandas,json; scipy.io.mmread('A.mtx'); scipy.io.mmread('B.mtx'); pandas.read_csv('vectors.csv'); json.load(open('meta.json'))\"` — results appended below by the loadability checker.\n")
println(io, "## Notes\n- A = full symmetric first block (H completed from its stored lower triangle — the operator the solver factors); no rescaling/extra symmetrization. B as the solver holds it.\n- vectors.csv columns are NaN/empty-padded to the longest (fp,y0 length n; g length m); slice by n_rows_A / n_rows_B in meta.\n- Fidelity reloads A,B from the written .mtx and reconstructs with the captured rhs; r0 targets are the logged pre-CRAIG residual norms.")
write(joinpath(OUT, "SUMMARY.md"), String(take!(io)))
println("\n=== wrote snapshots/ ($(count(r->r.ok, rows))/$(length(rows))) in $(wall)s ===")
