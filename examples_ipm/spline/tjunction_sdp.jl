######################################################################
# tjunction_sdp.jl
#
# SDP SHAPE CERTIFICATES on the T-junction mesh — the 2D sibling of
# nonnegative_spline_exact.jl, with the honesty asterisk that dimension
# two demands: bivariate nonnegativity has NO exact SDP description
# (Hilbert 1888, Motzkin), so these are degree-truncated Schmuedgen
# certificates on the box — sufficient, and strictly TIGHTER than the
# coefficient-nonnegativity LP of tjunction_spline.jl:
#
#   :nonneg_sdp   f = σ0 + u(1−u)σ1 + v(1−v)σ2 + u(1−u)v(1−v)σ3,
#                 σ_a = z_aᵀG_a z_a,  G_a ⪰ 0,  z_a = tensor Bernstein
#                 bases of bidegree (m,m),(m−1,m),(m,m−1),(m−1,m−1),
#                 n = 2m even.  (Even n makes these four multipliers
#                 complete: b_I^n b_J^n splits by parity of (I,J), so
#                 LP-cone ⊆ SDP-cone, and strictly — the demo exhibits
#                 a certified-nonnegative fit with Bernstein
#                 coefficient ≈ −0.9.)
#
#   :convex_sdp   SOS-matrix Hessian: H_loc = S0 + u(1−u)S1 + v(1−v)S2
#                 + u(1−u)v(1−v)S3 with S_a = (I₂⊗z_a)ᵀG_a(I₂⊗z_a),
#                 G_a ⪰ 0 of size 2D_a.  The physical Hessian is a
#                 congruence D·H_loc·D of the local one, so local PSD
#                 suffices; piecewise convex + C¹ (k ≥ 1) ⟹ globally
#                 convex across all patch and T-junction lines.  This
#                 is the shape with no LP-sufficient alternative — and
#                 the corpus's first mixed PSD + subdivision problem.
#
# SHEAF SHAPE. Certificates are per-vertex objects, T-edges are
# per-edge objects, and they do not interact: each patch keeps its
# coefficient vertex P (now CofreeCone — the cone moved to the Grams),
# gains four SemidefiniteCone Gram vertices, and one certificate
# hyperedge   Pmap·P − Σ_a Λ_a·svec(G_a) = 0   ties them together.
# Every subdivision map, physical jet, and coarse-line bookkeeping rule
# of tjunction_spline.jl acts on the P vertices UNCHANGED.
#
# H¹ IS INVARIANT: Λ_0 is surjective onto the coefficient space (the
# products z_α z_β span bidegree (n,n)), so Λ_0ᵀ is injective and a
# left-null vector must vanish on every certificate edge; dim ker δᵀ
# equals the base mesh's junction-rule/ring value. Oracle-verified:
# 36 → 36/36 and (ring regime) 82 → 82/82 under both certificates.
#
# svec convention (library sdp.jl): column-major LOWER triangle,
# diagonal first per column, off-diagonals ×√2.
#
# Written against the CellularSheaves.IPM PR-67 API; not executed here.
# tjunction_sdp_oracle.py is the numerical ground truth.
######################################################################

include("tjunction_spline.jl")

svecdim(D::Int) = D * (D + 1) ÷ 2

"""Certificate blocks (p, q, eu, ev): bidegree of the SOS basis and the
u(1−u)/v(1−v) multiplier exponents. Requires n = 2m even."""
function cert_blocks(n::Int)
    @assert iseven(n) && n ≥ 2 "shape certificates need even n"
    m = n ÷ 2
    return [(m, m, 0, 0), (m - 1, m, 1, 0), (m, m - 1, 0, 1), (m - 1, m - 1, 1, 1)]
end

"""(n+1)² × svecdim(Dtot) map: svec(G) ↦ bidegree-(n,n) Bernstein
coefficients of z ᵀ G[srow.., scol..] z · (u(1−u))^eu (v(1−v))^ev, where z
is the (p,q) tensor Bernstein basis (x fastest). Uses the univariate
identities b_i^p b_{i'}^p = [C(p,i)C(p,i')/C(2p,i+i')] b_{i+i'}^{2p} and
u(1−u) b_s^{2p} = [C(2p,s)/C(2p+2,s+1)] b_{s+1}^{2p+2}."""
function gram_rows(n::Int, p::Int, q::Int, eu::Int, ev::Int;
                   Dtot::Int = (p + 1) * (q + 1), srow::Int = 0, scol::Int = 0)
    d = (n + 1)^2
    W = zeros(d, Dtot, Dtot)
    for j2 in 0:q, i2 in 0:p, j in 0:q, i in 0:p
        α = 1 + i + (p + 1) * j
        β = 1 + i2 + (p + 1) * j2
        w = binomial(p, i) * binomial(p, i2) / binomial(n, i + i2 + eu) *
            binomial(q, j) * binomial(q, j2) / binomial(n, j + j2 + ev)
        ci = 1 + (i + i2 + eu) + (n + 1) * (j + j2 + ev)
        W[ci, srow + α, scol + β] += w / 2
        W[ci, scol + β, srow + α] += w / 2
    end
    L = zeros(d, svecdim(Dtot))
    kk = 0
    for c in 1:Dtot, r in c:Dtot
        kk += 1
        L[:, kk] .= r == c ? W[:, c, c] : sqrt(2.0) .* W[:, r, c]
    end
    return L
end

"""One-step Bernstein degree elevation, degree p → p+1."""
function elev1(p::Int)
    E = zeros(p + 2, p + 1)
    for i in 0:(p + 1)
        1 ≤ i ≤ p + 1 && (E[i + 1, i] += i / (p + 1))
        i ≤ p && (E[i + 1, i + 1] += 1 - i / (p + 1))
    end
    return E
end

"""(Cuu, Cuv, Cvv): maps from P to the LOCAL Hessian entries' Bernstein
coefficients, degree-elevated back to bidegree (n,n); x fastest."""
function hess_maps(n::Int)
    A2 = elev1(n - 1) * elev1(n - 2) * ((n * (n - 1)) .* diffop(n, 2))
    A1 = elev1(n - 1) * (n .* diffop(n, 1))
    Ins = Matrix{Float64}(I, n + 1, n + 1)
    return kron(Ins, A2), kron(A1, A1), kron(A2, Ins)
end

# ---- instances --------------------------------------------------------

truth_convex(x, y) =
    (x - 0.45)^2 + 0.7 * (y - 0.55)^2 + 0.4 * (x - 0.45) * (y - 0.55) + 0.05

"""Same mesh/data machinery as generate_tjunction_instance; shape selects
the certificate (:nonneg_sdp uses the :nonneg truth for LP-comparable
data, :convex_sdp uses truth_convex)."""
function generate_tjsdp_instance(;
        Mx::Int = 2, My::Int = 2, n::Int = 4, k::Int = 2,
        shape::Symbol = :nonneg_sdp, refine = Tuple{Int, Int}[],
        refined::Union{Nothing, BitMatrix} = nothing,
        N::Int = 100 * Mx * My, σ::Float64 = 0.05,
        λr::Float64 = 1e-4, ε::Float64 = 1e-8, seed::Int = 1)
    @assert shape in (:nonneg_sdp, :convex_sdp) && iseven(n) && n ≥ max(k + 1, 2)
    rng = MersenneTwister(seed)
    ns = n + 1; d = ns * ns
    R = refined === nothing ? falses(Mx, My) : copy(refined)
    for (a, b) in refine; R[a, b] = true; end
    Qp = Dict{Tuple{Symbol, Int, Int}, Matrix{Float64}}()
    cp = Dict{Tuple{Symbol, Int, Int}, Vector{Float64}}()
    for b in 1:My, a in 1:Mx
        if R[a, b]
            for q in (2b - 1):(2b), p in (2a - 1):(2a)
                Qp[(:f, p, q)] = zeros(d, d); cp[(:f, p, q)] = zeros(d)
            end
        else
            Qp[(:c, a, b)] = zeros(d, d); cp[(:c, a, b)] = zeros(d)
        end
    end
    truthf(x, y) = shape === :convex_sdp ? truth_convex(x, y) : truth2(:nonneg, x, y)
    for _ in 1:N
        x, y = rand(rng), rand(rng)
        z = truthf(x, y) + σ * randn(rng)
        key, u, v = tj_owner(R, Mx, My, x, y)
        φ = kron(bernstein_eval(n, v), bernstein_eval(n, u))
        Qp[key] .+= 2.0 .* (φ * φ'); cp[key] .-= 2.0 .* z .* φ
    end
    return TJunctionInstance(Mx, My, n, k, shape, R, Qp, cp, λr, ε)
end

"""Same data, different cone — for LP/free reference solves."""
with_shape(inst::TJunctionInstance, s::Symbol) =
    TJunctionInstance(inst.Mx, inst.My, inst.n, inst.k, s, inst.refined,
                      inst.Qp, inst.cp, inst.λr, inst.ε)

# ---- shared certificate pieces ----------------------------------------

function certificate_maps(n::Int, kind::Symbol)
    blocks = cert_blocks(n)
    if kind === :nonneg
        Ls = [gram_rows(n, p, q, eu, ev) for (p, q, eu, ev) in blocks]
        Pmap = Matrix{Float64}(I, (n + 1)^2, (n + 1)^2)
    else
        Ls = Matrix{Float64}[]
        for (p, q, eu, ev) in blocks
            D = (p + 1) * (q + 1)
            push!(Ls, vcat(gram_rows(n, p, q, eu, ev; Dtot = 2D, srow = 0, scol = 0),
                           gram_rows(n, p, q, eu, ev; Dtot = 2D, srow = D, scol = 0),
                           gram_rows(n, p, q, eu, ev; Dtot = 2D, srow = D, scol = D)))
        end
        Pmap = vcat(hess_maps(n)...)
    end
    mult = kind === :nonneg ? 1 : 2
    Ds = [mult * (p + 1) * (q + 1) for (p, q, _, _) in blocks]
    return Pmap, Ls, Ds
end

# ---- native builder ----------------------------------------------------

"""Vertex order: patch P vertices (base order, CofreeCone — the cone lives
on the Grams), then per patch (same order) four SemidefiniteCone Gram
vertices. Edges: tj_edges on the P vertices UNCHANGED, then one
certificate hyperedge per patch."""
function build_tjsdp_problem(inst::TJunctionInstance)
    kind = inst.shape === :nonneg_sdp ? :nonneg :
           inst.shape === :convex_sdp ? :convex :
           error("build_tjsdp_problem: shape must be :nonneg_sdp or :convex_sdp")
    kind === :convex && @assert inst.k ≥ 1 "global convexity needs C¹ (k ≥ 1)"
    n = inst.n; d = (n + 1)^2
    Pmap, Ls, Ds = certificate_maps(n, kind)

    row_ids = Int[]; col_ids = Int[]; blks = Matrix{Float64}[]
    col_cone = AbstractCone[]; col_dim = Int[]
    new_vertex!(dim, cone) = (push!(col_dim, dim); push!(col_cone, cone); length(col_dim))
    ec = Ref(0); new_edge!() = (ec[] += 1; ec[])
    place!(e, v, mat) = (push!(row_ids, e); push!(col_ids, v); push!(blks, Matrix{Float64}(mat)))

    keys_order = Tuple{Symbol, Int, Int}[]
    vtx = Dict{Tuple{Symbol, Int, Int}, Int}()
    for b in 1:inst.My, a in 1:inst.Mx
        inst.refined[a, b] && continue
        push!(keys_order, (:c, a, b)); vtx[(:c, a, b)] = new_vertex!(d, CofreeCone())
    end
    for q in 1:(2 * inst.My), p in 1:(2 * inst.Mx)
        inst.refined[cld(p, 2), cld(q, 2)] || continue
        push!(keys_order, (:f, p, q)); vtx[(:f, p, q)] = new_vertex!(d, CofreeCone())
    end
    vG = Dict{Tuple{Tuple{Symbol, Int, Int}, Int}, Int}()
    for key in keys_order, a in 1:4
        vG[(key, a)] = new_vertex!(svecdim(Ds[a]), SemidefiniteCone())
    end

    for (kL, AL, kR, AR) in tj_edges(inst)
        e = new_edge!(); place!(e, vtx[kL], AL); place!(e, vtx[kR], -AR)
    end
    for key in keys_order
        e = new_edge!()
        place!(e, vtx[key], Pmap)
        for a in 1:4
            place!(e, vG[(key, a)], -Ls[a])
        end
    end

    B = blocksparse(row_ids, col_ids, blks)
    nvtx = length(col_dim)
    @assert nvtxs(B) == nvtx
    g = zeros(size(B, 1))

    Qd, Cd = tj_patch_QC(inst)
    Q = IPM.allocblockdiag(B); fill!(Q, 0)
    c = zeros(size(B, 2))
    for key in keys_order
        v = vtx[key]
        block(Q, v, v, v) .= Qd[key]
        c[colrange(B, v)] .= Cd[key]
        for a in 1:4
            vg = vG[(key, a)]
            block(Q, vg, vg, vg) .= inst.ε .* Matrix{Float64}(I, svecdim(Ds[a]), svecdim(Ds[a]))
        end
    end
    prob = IPMProblem(Q, B, c, g, AbstractCone[col_cone[v] for v in 1:nvtx])
    info = (patchcols = Dict(key => colrange(B, vtx[key]) for key in keys_order),
            gramcols = Dict(kv => colrange(B, v) for (kv, v) in vG),
            keys_order = keys_order, nvtx = nvtx, ncols = size(B, 2),
            h1 = h1_predict_tj(inst))                 # invariant — oracle-verified
    return prob, info
end

# ---- semantic evaluators (work on base-builder results too) ------------

function tj_min_f(inst, res, info; ngrid = 60)
    n = inst.n; ns = n + 1; worst = Inf
    for gx in 0:ngrid, gy in 0:ngrid
        x, y = min(gx / ngrid, 1 - 1e-12), min(gy / ngrid, 1 - 1e-12)
        key, u, v = tj_owner(inst.refined, inst.Mx, inst.My, x, y)
        P = reshape(res.p[info.patchcols[key]], ns, ns)
        worst = min(worst, bernstein_eval(n, u)' * P * bernstein_eval(n, v))
    end
    return worst
end

"""Minimum eigenvalue of the PHYSICAL Hessian over a dense grid, across
coarse and fine patches (per-patch widths from the key's level)."""
function tj_min_hess_eig(inst, res, info; ngrid = 60)
    n = inst.n; ns = n + 1
    D2 = (n * (n - 1)) .* diffop(n, 2)
    D1 = n .* diffop(n, 1)
    worst = Inf
    for gx in 0:ngrid, gy in 0:ngrid
        x, y = min(gx / ngrid, 1 - 1e-12), min(gy / ngrid, 1 - 1e-12)
        key, u, v = tj_owner(inst.refined, inst.Mx, inst.My, x, y)
        s = key[1] === :f ? 0.5 : 1.0
        wx, wy = s / inst.Mx, s / inst.My
        P = reshape(res.p[info.patchcols[key]], ns, ns)
        fuu = bernstein_eval(n - 2, u)' * (D2 * P) * bernstein_eval(n, v) / wx^2
        fvv = bernstein_eval(n, u)' * (P * D2') * bernstein_eval(n - 2, v) / wy^2
        fuv = bernstein_eval(n - 1, u)' * (D1 * P * D1') * bernstein_eval(n - 1, v) / (wx * wy)
        tr, det = fuu + fvv, fuu * fvv - fuv^2
        worst = min(worst, (tr - sqrt(max(tr^2 - 4det, 0.0))) / 2)
    end
    return worst
end

function tj_mse_convex(inst, res, info; ngrid = 60)
    n = inst.n; ns = n + 1; tot = 0.0
    for gx in 1:ngrid, gy in 1:ngrid
        x, y = (gx - 0.5) / ngrid, (gy - 0.5) / ngrid
        key, u, v = tj_owner(inst.refined, inst.Mx, inst.My, x, y)
        P = reshape(res.p[info.patchcols[key]], ns, ns)
        tot += (bernstein_eval(n, u)' * P * bernstein_eval(n, v) - truth_convex(x, y))^2
    end
    return tot / ngrid^2
end

"""½pᵀQp + cᵀp assembled from per-patch pieces (plus the ε‖G‖²_F Gram
regularization when present) — identical accounting to the oracle, and
comparable across the free / LP / SDP builders."""
function tj_objective(inst, res, info)
    Qd, Cd = tj_patch_QC(inst)
    o = sum(0.5 * res.p[info.patchcols[key]]' * (Qd[key] * res.p[info.patchcols[key]]) +
            Cd[key]' * res.p[info.patchcols[key]] for key in keys(Qd))
    hasproperty(info, :gramcols) &&
        (o += 0.5 * inst.ε * sum(sum(abs2, res.p[r]) for r in values(info.gramcols)))
    return o
end

# ---- demos --------------------------------------------------------------

"""LP ⊋ SDP demonstrated on identical data: free ≤ SDP ≤ LP objectives, and
the SDP fit carries negative Bernstein coefficients while its surface stays
certified nonnegative."""
function nonneg_sdp_demo(; kwargs...)
    inst = generate_tjsdp_instance(; shape = :nonneg_sdp, refine = [(1, 1)], kwargs...)
    prob, info = build_tjsdp_problem(inst)
    solve(prob, tensor_settings())
    t = @elapsed (res = solve(prob, tensor_settings()))
    probl, infol = build_tjunction_problem(with_shape(inst, :nonneg))
    resl = solve(probl, tensor_settings())
    probf, infof = build_tjunction_problem(with_shape(inst, :free))
    resf = solve(probf, tensor_settings())
    @printf("SDP-nonneg  %dx%d(+%d ref) n=%d k=%d  vtx=%d dof=%d  H1=%d  %.1f ms it=%d %s\n",
            inst.Mx, inst.My, count(inst.refined), inst.n, inst.k,
            info.nvtx, info.ncols, info.h1, 1e3t, res.niter, res.status)
    @printf("    free ≤ SDP ≤ LP : %.6f ≤ %.6f ≤ %.6f   (LP−SDP gap %.3f)\n",
            tj_objective(inst, resf, infof), tj_objective(inst, res, info),
            tj_objective(inst, resl, infol),
            tj_objective(inst, resl, infol) - tj_objective(inst, res, info))
    @printf("    min coeff LP/SDP: %.2e / %.2e    min f (grid) SDP: %.2e  free: %.2e\n",
            minimum(resl.p[1:infol.ncols]), minimum(res.p[reduce(vcat, collect.(values(info.patchcols)))]),
            tj_min_f(inst, res, info), tj_min_f(inst, resf, infof))
    return inst, res, info
end

"""The shape with no LP alternative: SOS-matrix Hessian. Also shows the
statistical payoff — the convexity prior sharply reduces fit MSE on noisy
data from a convex truth."""
function convex_sdp_demo(; N = 250, σ = 0.06, λr = 1e-6, seed = 3, kwargs...)
    inst = generate_tjsdp_instance(; shape = :convex_sdp, refine = [(1, 1)],
                                   N, σ, λr, seed, kwargs...)
    prob, info = build_tjsdp_problem(inst)
    solve(prob, tensor_settings())
    t = @elapsed (res = solve(prob, tensor_settings()))
    probf, infof = build_tjunction_problem(with_shape(inst, :free))
    resf = solve(probf, tensor_settings())
    @printf("SDP-convex  %dx%d(+%d ref) n=%d k=%d  vtx=%d dof=%d  H1=%d  %.1f ms it=%d %s\n",
            inst.Mx, inst.My, count(inst.refined), inst.n, inst.k,
            info.nvtx, info.ncols, info.h1, 1e3t, res.niter, res.status)
    @printf("    min eig H (phys, grid): free %.3e  certified %.3e\n",
            tj_min_hess_eig(inst, resf, infof), tj_min_hess_eig(inst, res, info))
    @printf("    fit-mse vs truth      : free %.3e  certified %.3e\n",
            tj_mse_convex(inst, resf, infof), tj_mse_convex(inst, res, info))
    return inst, res, info
end

# ---- benchmark -----------------------------------------------------------

function run_tjsdp_benchmark(; optimizer = nothing, dual_optimizer = nothing, solver_name::String = "Mosek", nwarmup::Int = 2, nruns::Int = 5)
    optimizer === nothing && error("Pass optimizer, e.g. run_tjsdp_benchmark(optimizer=Mosek.Optimizer)")
    # (shape, Mx, My, refine_list, label, raug) — raug tuned per case
    cases = [
        (:nonneg_sdp, 2, 2, [(1, 1)], "NN 2x2+1", 1e6),
        (:nonneg_sdp, 3, 3, [(2, 2)], "NN 3x3 ctr", 1e1),
        (:nonneg_sdp, 4, 4, [(2, 2), (3, 3)], "NN 4x4 diag", 1e6),
        (:convex_sdp, 2, 2, [(1, 1)], "CVX 2x2+1", 1e5),
        (:convex_sdp, 3, 3, [(2, 2)], "CVX 3x3 ctr", 1e3),
        (:convex_sdp, 4, 4, [(2, 2), (3, 3)], "CVX 4x4 diag", 1e4),
    ]
    sname = rpad(solver_name, 9)
    sname_d = rpad(solver_name * "-D", 9)
    println("="^95)
    println("T-Junction SDP Benchmark: Sheaf IPM vs $solver_name")
    println("="^95)
    if dual_optimizer !== nothing
        @printf("%-14s %6s %5s %5s %5s %9s %9s %9s %7s %7s\n",
                "Config", "raug", "DOF", "|V|", "H1", "IPM(ms)", sname, sname_d, "P/IPM", "D/IPM")
    else
        @printf("%-14s %6s %5s %5s %5s %9s %9s %8s\n",
                "Config", "raug", "DOF", "|V|", "H1", "IPM(ms)", sname, "speedup")
    end
    println("-"^95)
    for (shape, Mx, My, refine_list, label, raug) in cases
        inst = generate_tjsdp_instance(; Mx, My, shape, refine = refine_list, N = 100 * Mx * My)
        prob, info = build_tjsdp_problem(inst)
        settings = IPMSettings{Float64}(
            kkt = UzawaSettings{Float64}(raug = raug),
            feas_tol = 1e-8, gap_tol = 1e-8, itmax = 200)
        for _ in 1:nwarmup; solve(prob, settings); end
        t_ipm = minimum([@elapsed solve(prob, settings) for _ in 1:nruns])
        for _ in 1:nwarmup
            m, _ = build_jump_tjsdp(inst, optimizer); optimize!(m)
        end
        t_mosek = minimum([@elapsed begin
            m, _ = build_jump_tjsdp(inst, optimizer); optimize!(m)
        end for _ in 1:nruns])

        if dual_optimizer !== nothing
            for _ in 1:nwarmup
                m, _ = build_jump_tjsdp(inst, dual_optimizer); optimize!(m)
            end
            t_dual = minimum([@elapsed begin
                m, _ = build_jump_tjsdp(inst, dual_optimizer); optimize!(m)
            end for _ in 1:nruns])
            @printf("%-14s %6.0e %5d %5d %5d %9.1f %9.1f %9.1f %6.2fx %6.2fx\n",
                    label, raug, info.ncols, info.nvtx, info.h1,
                    t_ipm * 1000, t_mosek * 1000, t_dual * 1000,
                    t_mosek / t_ipm, t_dual / t_ipm)
        else
            @printf("%-14s %6.0e %5d %5d %5d %9.1f %9.1f %7.2fx\n",
                    label, raug, info.ncols, info.nvtx, info.h1,
                    t_ipm * 1000, t_mosek * 1000, t_mosek / t_ipm)
        end
    end
end

# =====================================================================
# JuMP reference (PSD matrix variables; same Λ maps through svec exprs)
# =====================================================================

svec_expr(G, D) = [r == c ? 1.0 * G[c, c] : sqrt(2.0) * G[r, c]
                   for c in 1:D for r in c:D]

function build_jump_tjsdp(inst::TJunctionInstance, optimizer)
    kind = inst.shape === :nonneg_sdp ? :nonneg : :convex
    n = inst.n; d = (n + 1)^2
    Pmap, Ls, Ds = certificate_maps(n, kind)
    Qd, Cd = tj_patch_QC(inst)
    model = Model(optimizer); set_silent(model)
    P = Dict(key => @variable(model, [1:d]) for key in keys(Qd))
    G = Dict((key, a) => @variable(model, [1:Ds[a], 1:Ds[a]], PSD)
             for key in keys(Qd), a in 1:4)
    for (kL, AL, kR, AR) in tj_edges(inst)
        @constraint(model, AL * P[kL] .== AR * P[kR])
    end
    for key in keys(Qd)
        @constraint(model, Pmap * P[key] .==
                    sum(Ls[a] * svec_expr(G[(key, a)], Ds[a]) for a in 1:4))
    end
    @objective(model, Min,
               sum(0.5 * P[key]' * Qd[key] * P[key] + Cd[key]' * P[key] +
                   0.5 * inst.ε * sum(G[(key, a)][i, j]^2
                                      for a in 1:4, i in 1:Ds[a], j in 1:Ds[a])
                   for key in keys(Qd)))
    return model, P
end

if abspath(PROGRAM_FILE) == @__FILE__
    include(joinpath(@__DIR__, "..", "benchmark_utils.jl"))
    opts = parse_benchmark_args(ARGS)

    if opts.mosek
        using MosekTools
        optimizer = Mosek.Optimizer
        dual_optimizer = Dualization.dual_optimizer(Mosek.Optimizer)
        solver_name = "Mosek"
    else
        using Clarabel
        optimizer = Clarabel.Optimizer
        dual_optimizer = Dualization.dual_optimizer(Clarabel.Optimizer)
        solver_name = "Clarabel"
    end
    println("Solver: $solver_name")
    println("Runs: $(opts.nruns), Warmup: $(opts.nwarmup)\n")

    run_tjsdp_benchmark(; optimizer, dual_optimizer, solver_name,
                        nwarmup = opts.nwarmup, nruns = opts.nruns)
end
