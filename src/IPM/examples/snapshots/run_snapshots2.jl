# Targeted solve-state snapshot capture (snapshots2 spec). Reads the Python-produced _selection.tsv
# (worst/best floor/ceiling-error iterations + late iterations), and for each capture advances the
# oracle chosen-α trajectory to the selected iteration, capturing every role's pre-solve KKT input
# state. Same machinery/fidelity as run_snapshots.jl. Writes snapshots2/<key>_<solv>_it<K>/ + SUMMARY.md.
#
# Usage:  julia --project=examples examples/snapshots/run_snapshots2.jl
using CellularSheaves
using CellularSheaves.IPM: HSDSettings, IPMSettings, init, snapshot_capture
import CellularSheaves.IPM as IPM
using LinearAlgebra, SparseArrays, Printf

const EX  = dirname(@__DIR__)
const OUT = joinpath(EX, "snapshots2"); mkpath(OUT)
gp(x) = x isa Tuple ? x[1] : x
const GIT = try readchomp(`git -C $EX rev-parse HEAD`) catch; "unknown" end

include("$EX/e01.jl"); include("$EX/e02.jl"); include("$EX/e08.jl"); include("$EX/e09.jl")
include("$EX/e15.jl"); include("$EX/e04.jl"); include("$EX/adversarial/X06.jl"); include("$EX/adversarial/X09.jl")
function buildprob(key)
    key == "e01" && return gp(build_merton_chain(merton_instance(; nx=21, olap=5, nsteps=10)))
    key == "e02" && return gp(build_qm_problem(generate_qm_instance(; N=10)))
    key == "X06" && return build_corner_soc(; corner_r=0.1)
    key == "X09" && return build_ceiling_base(; S=1e2)
    key == "SOS_P8_n9_s1" && return gp(build_sos_spline(sos_instance(; P=8, n=9, seed=1)))
    key == "e15" && (train = [[0.0, 60.0*(v-1), 0.0, 0.0, 0.0, 0.0] for v in 1:V_BASE];
                     return gp(build_cw(train, [zeros(6) for _ in 1:V_BASE], N_BASE; gam=1e-1, kap=1e-1)))
    if key == "e08"
        dt = 4.0/64; inst = exec_instance(; n=1, T=64)
        return gp(build_exec(inst; Σ=fill(dt,1,1), X0=[1.5], η=[1.0/sqrt(dt)], γ=1.0))
    end
    if key == "e09"
        inst = smoother_instance(); sys = make_system(inst); _, y, _ = simulate(sys, inst, 100)
        return gp(build_smoother(sys, inst, y))
    end
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
function read_mtx(path)
    lines = readlines(path)
    h = findfirst(l -> !startswith(strip(l), "%") && !isempty(strip(l)), lines)
    m, n, _ = parse.(Int, split(strip(lines[h])))
    I = Int[]; J = Int[]; V = Float64[]
    for l in lines[h+1:end]
        isempty(strip(l)) && continue
        p = split(strip(l)); push!(I, parse(Int, p[1])); push!(J, parse(Int, p[2])); push!(V, parse(Float64, p[3]))
    end
    return sparse(I, J, V, m, n)
end
function write_vectors(path, cols)
    L = maximum(length(v) for (_, v) in cols)
    open(path, "w") do io
        println(io, join((k for (k, _) in cols), ","))
        for i in 1:L
            println(io, join(((i <= length(v) ? @sprintf("%.17g", v[i]) : "") for (_, v) in cols), ","))
        end
    end
end
function fidelity(dir, f, g, y0, alpha, r0)
    A = read_mtx(joinpath(dir, "A.mtx")); B = read_mtx(joinpath(dir, "B.mtx"))
    β = 1 / alpha
    rhs = β .* f .+ B' * (g .+ (y0 === nothing ? zero(g) : β .* y0))
    x0 = (β .* A .+ B' * B) \ rhs
    return abs(norm(g .- B * x0) - r0) / max(r0, eps())
end
jstr(v::AbstractString) = "\"$v\""; jstr(v::Bool) = string(v)
jstr(v::Real) = isfinite(v) ? string(v) : "null"; jstr(v) = "\"$(string(v))\""
write_meta(path, d) = open(path, "w") do io
    println(io, "{" * join(("\"$k\":$(jstr(v))" for (k, v) in d), ",") * "}")
end
function rolevecs(snaps, solv)
    if solv == "ipm"
        return (p = snaps[findfirst(x -> x.role == :p, snaps)],
                c = snaps[findfirst(x -> x.role == :c, snaps)], w = nothing)
    else
        nts = [x for x in snaps if x.role == :newton]
        return (p = nts[1], c = nts[2], w = snaps[findfirst(x -> x.role == :w, snaps)])
    end
end

function do_capture!(rows, group, tag, key, solv, K, errtype, errval)
    prob = gp(buildprob(key)); s0 = init(prob, mkset(solv))
    cap = snapshot_capture(s0, K)
    tagname = "$(key)_$(solv)_it$K"
    if cap === nothing
        push!(rows, (; group, tag, tagname, ok=false, note="solve terminated before iter $K"))
        println("SKIP $tagname: terminated before iter $K"); return
    end
    rv = rolevecs(cap.snaps, solv)
    dir = joinpath(OUT, tagname); mkpath(dir)
    A = rv.p.A; B = rv.p.Bmat; n, m = size(A, 1), size(B, 1)
    write_mtx(joinpath(dir, "A.mtx"), A); write_mtx(joinpath(dir, "B.mtx"), B)
    zerom = zeros(m)
    cols = ["fp"=>rv.p.f, "gp"=>rv.p.g, "y0p"=>(rv.p.y0===nothing ? zerom : rv.p.y0),
            "fc"=>rv.c.f, "gc"=>rv.c.g, "y0c"=>(rv.c.y0===nothing ? zerom : rv.c.y0)]
    solv != "ipm" && append!(cols, ["fw"=>rv.w.f, "gw"=>rv.w.g, "y0w"=>(rv.w.y0===nothing ? zerom : rv.w.y0)])
    write_vectors(joinpath(dir, "vectors.csv"), cols)
    row = cap.row; r0_w = hasproperty(row, :r0_w) ? row.r0_w : nothing
    meta = Any["problem"=>key, "solver"=>solv, "iter"=>K, "git_commit"=>GIT, "group"=>group, "select_tag"=>tag,
        "select_err_type"=>errtype, "select_err"=>errval,
        "alpha_used"=>cap.chosenα, "force_tol"=>cap.force_tol, "floor_tol"=>cap.floor_tol,
        "mu"=>row.μ, "nc"=>cap.entry.nc, "ng"=>cap.entry.ng]
    solv != "ipm" && append!(meta, ["tau"=>cap.entry.tau, "kappa"=>cap.entry.kappa])
    append!(meta, ["r0_p"=>row.r0_p, "r0_c"=>row.r0_c])
    solv != "ipm" && push!(meta, "r0_w"=>r0_w)
    append!(meta, ["n_rows_A"=>n, "n_rows_B"=>m, "nnz_A"=>nnz(A), "nnz_B"=>nnz(B)])
    write_meta(joinpath(dir, "meta.json"), meta)
    fp = fidelity(dir, rv.p.f, rv.p.g, rv.p.y0, cap.chosenα, row.r0_p)
    fc = fidelity(dir, rv.c.f, rv.c.g, rv.c.y0, cap.chosenα, row.r0_c)
    fw = solv == "ipm" ? nothing : fidelity(dir, rv.w.f, rv.w.g, rv.w.y0, cap.chosenα, r0_w)
    push!(rows, (; group, tag, tagname, ok=true, n, m, nnzA=nnz(A), nnzB=nnz(B), fp, fc, fw,
        errtype, errval, note=""))
    println("captured $tagname (G$group $tag $errtype=$errval): A $(n)×$(n), B $(m)×$(n); fid p=$(round(fp,sigdigits=2)) c=$(round(fc,sigdigits=2))" * (fw===nothing ? "" : " w=$(round(fw,sigdigits=2))"))
end

# ---- read selection and capture ----
sel = readlines(joinpath(OUT, "_selection.tsv"))[2:end]
rows = NamedTuple[]
t0 = time()
for line in sel
    isempty(strip(line)) && continue
    p = split(line, '\t')
    group = parse(Int, p[1]); tag = p[2]; key = p[3]; solv = p[4]; K = parse(Int, p[5])
    errtype = p[6]; errval = tryparse(Float64, p[7])
    do_capture!(rows, group, tag, key, solv, K, errtype, errval === nothing ? NaN : errval)
end
wall = round(time() - t0, digits=1)

# ---- SUMMARY.md ----
io = IOBuffer()
println(io, "# snapshots2 — targeted assumption/law correspondence captures\n")
println(io, "git: `$(GIT[1:min(10,end)])`  ·  captured $(count(r->r.ok, rows))/$(length(rows))  ·  wall $(wall)s\n")
println(io, "## Captures (selection reproducibility + fidelity + shapes)\n")
println(io, "| group | tag | problem | solver | iter | sel-err | A (n×n) | nnz A | B (m×n) | nnz B | fid p | fid c | fid w |")
println(io, "|---|---|---|---|---|---|---|---|---|---|---|---|---|")
for r in rows
    if r.ok
        parts = split(r.tagname, "_"); K = parts[end]
        prob = join(parts[1:end-2], "_"); solv = parts[end-1]
        fw = r.fw === nothing ? "—" : @sprintf("%.1e", r.fw)
        selerr = isnan(r.errval) ? "—" : @sprintf("%s %+.3f", r.errtype, r.errval)
        @printf(io, "| %d | %s | %s | %s | %s | %s | %d×%d | %d | %d×%d | %d | %.1e | %.1e | %s |\n",
            r.group, r.tag, prob, solv, replace(K,"it"=>""), selerr, r.n, r.n, r.nnzA, r.m, r.n, r.nnzB, r.fp, r.fc, fw)
    else
        println(io, "| $(r.group) | $(r.tag) | $(r.tagname) | | | | | | | | | | $(r.note) |")
    end
end
allf = [f for r in rows if r.ok for f in (r.fp, r.fc, r.fw) if f !== nothing]
println(io, "\n**Fidelity: max relative Δ = $(isempty(allf) ? "n/a" : @sprintf("%.2e", maximum(allf)))** (must ≤ 1e-2).")
println(io, "\n## Notes\n- Selection per snapshots2 spec: W0 = α with every role base==1 & refn==0; probe = middle of W0; pred_floor = maxᵣ α‖r0ᵣ‖/ε, pred_ceil = minᵣ αε/‖s0ᵣ‖ (‖s0‖ = role *_res0_dual, ε = max(force_tol,floor_tol)); err = log10(pred) − log10(obs). Sweeps: examples/fine2/.\n- G1 worst-error, G2 best-error controls, G3 late iterations (last-W0 / last-overall; e08 & e01 had them coincide → one capture each).\n- Capture point, A/B/vector conventions, and fidelity identical to snapshots/ (run 1).")
write(joinpath(OUT, "SUMMARY.md"), String(take!(io)))
println("\n=== wrote snapshots2/ ($(count(r->r.ok, rows))/$(length(rows))) in $(wall)s ===")
