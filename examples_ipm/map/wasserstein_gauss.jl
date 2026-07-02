######################################################################
# wasserstein_gauss.jl
#
# GAUSSIAN (BURES–WASSERSTEIN) OPTIMAL TRANSPORT as a cellular sheaf of
# covariance blocks — wasserstein_graph.jl after the map_lp → map_sdp
# cone swap:
#
#   coupling cell per cost-graph edge:  J_e = [[Σ_i, K],[Kᵀ, Σ_j]] ⪰ 0
#     (SemidefiniteCone, svec coords; no support grid — dimension is
#      the marginal's dimension, not its discretization)
#   cost (centered):  ⟨C_e, J_e⟩ with C_e = [[I,∓I],[∓I,I]] = E‖x∓y‖²
#   restriction maps: PRINCIPAL-BLOCK COMPRESSIONS (map_sdp's, verbatim)
#   marginals: Σ_k pinned (ride in g) or free PSD cells (barycenter)
#
# WHY THIS IS EXACT (not a heuristic): quadratic pairwise costs read
# only second moments, and every consistent second-moment profile is
# realized by a Gaussian — so Gaussian graphical MMOT *is* an SDP over
# the pattern covariance (Givens–Shortt; Olkin–Pukelsheim for the
# 2-marginal case, optimum = Bures² closed form). Equivalently: this is
# the DEGREE-1 MOMENT RELAXATION, and for non-Gaussian marginals its
# value is Gelbrich's lower bound on the true W₂² — one SDP, two
# readings (Gelbrich 1990).
#
# MEASURED (wasserstein_gauss_oracle.py, CLARABEL):
#   • coupling SDP = Bures closed form to 3.6e-7; compressions exact.
#   • barycenter star (free center PSD cell) = Álvarez-Esteban fixed
#     point to 6e-5 in Σ (Agueh–Carlier; Bhatia–Jain–Lim).
#   • THE GAP IS SIGN HOLONOMY: on cycles with `anti` edges (cost
#     E‖x+y‖²), pairwise = 0 always; exact = 0 iff #anti is even.
#     Triangle gap 3.000; C4 gap 8−4√2 ≈ 2.3431. Same parity law as
#     discrete-support transport holonomy (composed map ≠ id) — NOT the
#     binary-MRF odd-cycle-length law.
#   • REPAIR BY COVER: adding the cycle cell (4d×4d block with
#     compression agreements to member couplings) closes the C4 gap
#     exactly — refine the cover, restore exactness (GJSW).
#   • THE FUSION: completion mode (Q = I data-fit to target blocks,
#     no transport cost) reproduces correlation completion verbatim:
#     C4 targets (s,s,s,−s): sheaf distance ≈ 0 for all |s| ≤ 1 (blind),
#     exact projection turns on at s* = 1/√2 (BJT arccos law):
#     4.0e-8 / 6.7e-4 / 3.45e-2 at s = 0.70 / 0.72 / 0.80.
#   • MOMENT LADDER past Gelbrich (oracle §5, bimodal marginals):
#     Gelbrich 0.0216 ≤ degree-2 0.0978 ≤ true W₂² 0.1231.
#   • H¹ = 0 everywhere (2-marginal, triangle, star) — the Kantorovich
#     mass gauge has NO Gaussian analog: couplings are not normalized
#     objects, the two compression blocks of a coupling share nothing.
#     Contrast wasserstein_graph.jl's H¹ = |E| − free-balance.
#
# Means: the mean part of W₂² separates (‖m_i − m_j‖², a tiny QP with
# closed forms); cells here are centered covariances, standard practice.
#
# Written against the CellularSheaves.IPM PR-67 API; not executed here.
# wasserstein_gauss_oracle.py is the numerical ground truth.
######################################################################

using CellularSheaves.IPM
using CellularSheaves.BlockSparseArrays: colrange, rowrange, blocksparse, block, nvtxs
using LinearAlgebra
using Printf
using Random

# ---- svec (library convention: col-major LOWER, diag first, ×√2) ------

gw_svecdim(D::Int) = D * (D + 1) ÷ 2

function gw_svec(M::AbstractMatrix)
    D = size(M, 1); v = zeros(gw_svecdim(D)); k = 0
    for c in 1:D, r in c:D
        k += 1
        v[k] = r == c ? M[c, c] : sqrt(2.0) * M[r, c]
    end
    return v
end

function gw_smat(v::AbstractVector, D::Int)
    M = zeros(D, D); k = 0
    for c in 1:D, r in c:D
        k += 1
        if r == c
            M[c, c] = v[k]
        else
            M[r, c] = M[c, r] = v[k] / sqrt(2.0)
        end
    end
    return M
end

"""svec(M) ↦ svec(M[rows, rows]) — principal-block compression."""
function gw_comp_map(D::Int, rows::Vector{Int})
    sd = gw_svecdim(D)
    L = zeros(gw_svecdim(length(rows)), sd)
    for k in 1:sd
        e = zeros(sd); e[k] = 1.0
        L[:, k] .= gw_svec(gw_smat(e, D)[rows, rows])
    end
    return L
end

# ---- Bures closed form + barycenter fixed point (demo references) ------

function gw_sqrtm(A::Matrix{Float64})
    E = eigen(Symmetric(A))
    return E.vectors * Diagonal(sqrt.(max.(E.values, 0.0))) * E.vectors'
end

"""W₂²(N(0,S1), N(0,S2)) = tr S1 + tr S2 − 2 tr((S1^½ S2 S1^½)^½)."""
function gw_bures2(S1, S2)
    R = gw_sqrtm(S1)
    return tr(S1) + tr(S2) - 2tr(gw_sqrtm(R * S2 * R))
end

"""Álvarez-Esteban et al. fixed point for the BW barycenter."""
function gw_barycenter_fixedpoint(Sigmas, λs; iters::Int = 200)
    S = Matrix{Float64}(I, size(Sigmas[1])...)
    for _ in 1:iters
        R = gw_sqrtm(S); Ri = inv(R)
        T = sum(λ .* gw_sqrtm(R * Sk * R) for (λ, Sk) in zip(λs, Sigmas))
        S = Ri * T * T * Ri
    end
    return S
end

# ---- instance ------------------------------------------------------------

struct GWInstance
    dims::Vector{Int}
    edges::Vector{Tuple{Int, Int}}
    pins::Dict{Int, Matrix{Float64}}     # node -> covariance
    free::Vector{Int}
    weights::Vector{Float64}
    anti::Set{Int}                       # edge indices with cost E‖x+y‖²
    tg::Dict{Int, Matrix{Float64}}       # completion mode: target blocks
    Wt::Dict{Int, Matrix{Float64}}       # Higham weight per targeted coupling
    prior::Dict{Int, Tuple{Float64, Matrix{Float64}}}  # free node -> (λ, P)
    ε::Float64
end

function gw_instance(dims, edges, pins; weights = ones(length(edges)),
                     anti = Int[], tg = Dict{Int, Matrix{Float64}}(),
                     Wt = Dict{Int, Matrix{Float64}}(),
                     prior = Dict{Int, Tuple{Float64, Matrix{Float64}}}(),
                     ε::Float64 = 1e-8)
    free = [k for k in 1:length(dims) if !haskey(pins, k)]
    GWInstance(dims, edges, pins, free, collect(Float64, weights),
               Set{Int}(anti), tg, Wt, prior, ε)
end

"""Gram of M ↦ ‖W^½MW^½‖²_F in svec coords (symmetric Kronecker):
column k = svec(W · smat(e_k) · W). Higham's weighted norm."""
function gw_symkron(W::Matrix{Float64})
    D = size(W, 1); sd = gw_svecdim(D)
    Q = zeros(sd, sd)
    for k in 1:sd
        e = zeros(sd); e[k] = 1.0
        Q[:, k] .= gw_svec(W * gw_smat(e, D) * W)
    end
    return Q
end

function gw_edge_cost(di::Int, dj::Int, anti::Bool)
    s = anti ? 1.0 : -1.0
    return [Matrix{Float64}(I, di, di)  s*Matrix{Float64}(I, di, dj);
            s*Matrix{Float64}(I, dj, di)  Matrix{Float64}(I, dj, dj)]
end

# ---- builder ----------------------------------------------------------------

"""Edge-cover Gaussian sheaf; `cycles` adds larger cells (each a joint
covariance over the listed nodes) with compression agreements to member
couplings — the repair-by-cover rung. Completion mode (nonempty tg):
objective ½Σ‖J_e − T_e‖² instead of transport (Q = I, c = −svec T)."""
function build_gw_problem(inst::GWInstance; cycles::Vector{Vector{Int}} = Vector{Int}[])
    dims, edges = inst.dims, inst.edges
    nE = length(edges)
    completion = !isempty(inst.tg)
    cd = [dims[i] + dims[j] for (i, j) in edges]
    row_ids = Int[]; col_ids = Int[]; blks = Matrix{Float64}[]
    place!(e, v, m) = (push!(row_ids, e); push!(col_ids, v); push!(blks, Matrix{Float64}(m)))

    vμ = Dict(k => nE + t for (t, k) in enumerate(inst.free))
    nv = nE + length(inst.free) + length(cycles)
    ec = Ref(0)
    pinval = Dict{Int, Vector{Float64}}()
    for (t, (i, j)) in enumerate(edges)
        D = cd[t]
        for (k, rows) in ((i, collect(1:dims[i])), (j, collect(dims[i]+1:D)))
            e = (ec[] += 1)
            L = gw_comp_map(D, rows)
            place!(e, t, L)
            if haskey(inst.pins, k)
                pinval[e] = gw_svec(inst.pins[k])
            else
                place!(e, vμ[k], -Matrix{Float64}(I, size(L, 1), size(L, 1)))
            end
        end
    end
    for (s, C) in enumerate(cycles)
        Dc = sum(dims[v] for v in C)
        posn = Dict(v => sum((dims[w] for w in C[1:q-1]); init = 0) for (q, v) in enumerate(C))
        u = nE + length(inst.free) + s
        for (t, (i, j)) in enumerate(edges)
            (i in C && j in C) || continue
            rows = vcat(collect(posn[i]+1:posn[i]+dims[i]),
                        collect(posn[j]+1:posn[j]+dims[j]))
            e = (ec[] += 1)
            Lc = gw_comp_map(Dc, rows)
            place!(e, u, Lc)
            place!(e, t, -Matrix{Float64}(I, size(Lc, 1), size(Lc, 1)))
        end
    end

    B = blocksparse(row_ids, col_ids, blks)
    g = zeros(size(B, 1))
    for (e, v) in pinval
        g[rowrange(B, e)] .= v
    end
    @assert nvtxs(B) == nv
    Q = IPM.allocblockdiag(B); fill!(Q, 0)
    c = zeros(size(B, 2))
    for v in 1:nv
        d = length(colrange(B, v))
        block(Q, v, v, v) .= inst.ε .* Matrix{Float64}(I, d, d)
        if v ≤ nE
            if completion
                if haskey(inst.tg, v)
                    Qw = haskey(inst.Wt, v) ? gw_symkron(inst.Wt[v]) :
                                              Matrix{Float64}(I, d, d)
                    block(Q, v, v, v) .+= Qw
                    c[colrange(B, v)] .= -(Qw * gw_svec(inst.tg[v]))
                end
            else
                (i, j) = edges[v]
                c[colrange(B, v)] .= inst.weights[v] .*
                    gw_svec(gw_edge_cost(dims[i], dims[j], v in inst.anti))
            end
        end
    end
    for (k, v) in vμ                     # shrinkage: ½λ‖Σ_k − P‖²
        if haskey(inst.prior, k)
            λp, P = inst.prior[k]
            d = length(colrange(B, v))
            block(Q, v, v, v) .+= λp .* Matrix{Float64}(I, d, d)
            c[colrange(B, v)] .= -λp .* gw_svec(P)
        end
    end
    cones = AbstractCone[SemidefiniteCone() for _ in 1:nv]
    prob = IPMProblem(Q, B, c, g, cones)
    info = (cupcols = Dict(t => colrange(B, t) for t in 1:nE),
            mucols = Dict(k => colrange(B, v) for (k, v) in vμ),
            cd = cd, nv = nv, completion = completion,
            h1 = 0)   # measured: 0 on 2-marginal/triangle/star — no mass gauge
    return prob, info
end

"""One pattern-covariance cell: the poly-size exact graphical MMOT
(transport reading) — the Gaussian row's luxury."""
function build_gw_exact_problem(inst::GWInstance)
    dims, edges = inst.dims, inst.edges
    n = sum(dims)
    pos = cumsum(vcat(0, dims))
    row_ids = Int[]; col_ids = Int[]; blks = Matrix{Float64}[]
    ec = 0
    g_parts = Vector{Float64}[]
    for (k, P) in inst.pins
        ec += 1
        push!(row_ids, ec); push!(col_ids, 1)
        push!(blks, gw_comp_map(n, collect(pos[k]+1:pos[k+1])))
        push!(g_parts, gw_svec(P))
    end
    B = blocksparse(row_ids, col_ids, blks)
    g = zeros(size(B, 1))
    for (e, v) in enumerate(g_parts)
        g[rowrange(B, e)] .= v
    end
    Cglob = zeros(n, n)
    for (t, (i, j)) in enumerate(edges)
        Ce = inst.weights[t] .* gw_edge_cost(dims[i], dims[j], t in inst.anti)
        rows = vcat(collect(pos[i]+1:pos[i+1]), collect(pos[j]+1:pos[j+1]))
        Cglob[rows, rows] .+= Ce
    end
    sd = gw_svecdim(n)
    Q = IPM.allocblockdiag(B); fill!(Q, 0)
    block(Q, 1, 1, 1) .= inst.ε .* Matrix{Float64}(I, sd, sd)
    prob = IPMProblem(Q, B, gw_svec(Cglob), g, AbstractCone[SemidefiniteCone()])
    return prob, (n = n,)
end

# ---- demos -------------------------------------------------------------------

gw_settings(; raug = 1e5) = IPMSettings{Float64}(kkt = UzawaSettings{Float64}(raug = raug))

"""Two-marginal Gaussian OT vs the Bures closed form."""
function gw_bures_demo(; d::Int = 3, seed::Int = 0)
    rng = MersenneTwister(seed)
    S1 = randn(rng, d, d); S1 = S1 * S1' + 0.3I
    S2 = randn(rng, d, d); S2 = S2 * S2' + 0.3I
    inst = gw_instance([d, d], [(1, 2)], Dict(1 => S1, 2 => S2))
    prob, info = build_gw_problem(inst)
    res = solve(prob, gw_settings())
    v = prob.c' * res.p
    @printf("2-marginal: sheaf %.8f   Bures closed form %.8f   |Δ| %.2e   it=%d %s\n",
            v, gw_bures2(S1, S2), abs(v - gw_bures2(S1, S2)), res.niter, res.status)
    return inst, res, info
end

"""Barycenter star (free center PSD cell) vs the fixed point."""
function gw_barycenter_demo(; d::Int = 2, seed::Int = 0)
    rng = MersenneTwister(seed)
    Sig = [(X = randn(rng, d, d); X * X' + 0.2I) for _ in 1:3]
    λ = [0.5, 0.3, 0.2]
    inst = gw_instance([d, d, d, d], [(1, 2), (1, 3), (1, 4)],
                       Dict(2 => Sig[1], 3 => Sig[2], 4 => Sig[3]); weights = λ)
    prob, info = build_gw_problem(inst)
    res = solve(prob, gw_settings())
    Sc = gw_smat(res.p[info.mucols[1]], d)
    Sfp = gw_barycenter_fixedpoint(Sig, λ)
    @printf("barycenter: ‖Σ_center − Σ_fixedpoint‖ = %.2e   it=%d %s\n",
            norm(Sc - Sfp), res.niter, res.status)
    return inst, res, info
end

"""Sign holonomy + repair: triangle/C4 pairwise vs one-cell exact; the
cycle cell closes the C4 gap."""
function gw_gap_demo()
    for (K, nanti) in ((3, 1), (3, 2), (4, 1))
        inst = gw_instance(fill(1, K), [(k, mod1(k + 1, K)) for k in 1:K],
                           Dict(k => ones(1, 1) for k in 1:K); anti = collect(1:nanti))
        prob, _ = build_gw_problem(inst)
        res = solve(prob, gw_settings())
        probx, _ = build_gw_exact_problem(inst)
        resx = solve(probx, gw_settings())
        vp, ve = prob.c' * res.p, probx.c' * resx.p
        @printf("  K=%d nanti=%d: pairwise %+.6f  exact %+.6f  gap %.6f\n",
                K, nanti, vp, ve, ve - vp)
    end
    inst = gw_instance(fill(1, 4), [(k, mod1(k + 1, 4)) for k in 1:4],
                       Dict(k => ones(1, 1) for k in 1:4); anti = [1])
    probc, _ = build_gw_problem(inst; cycles = [[1, 2, 3, 4]])
    resc = solve(probc, gw_settings())
    @printf("  C4 + cycle cell: %.6f   (oracle exact 2.343146 = 8−4√2)\n",
            probc.c' * resc.p)
end

"""The fusion: nearest-completion mode on the C4 family (s,s,s,−s) —
sheaf blind for all |s| ≤ 1, exact turning on at s* = 1/√2 (see
oracle; exact rung = build_gw_exact analog with Q-data — here we print
the sheaf side and cite the measured exact values)."""
function gw_completion_demo(; s::Float64 = 0.8)
    tg = Dict(t => [1.0 v; v 1.0] for (t, v) in enumerate((s, s, s, -s)))
    inst = gw_instance(fill(1, 4), [(k, mod1(k + 1, 4)) for k in 1:4],
                       Dict(k => ones(1, 1) for k in 1:4); tg = tg)
    prob, info = build_gw_problem(inst)
    res = solve(prob, gw_settings())
    const_term = 0.5 * sum(sum(abs2, gw_svec(T)) for T in values(tg))
    dist2 = prob.c' * res.p + 0.5 * sum(abs2, res.p) - 0.5 * inst.ε * sum(abs2, res.p) + const_term
    @printf("completion s=%.2f: sheaf dist² = %.2e   (exact: 3.45e-2 at s=0.8; s*=1/√2)\n",
            s, dist2)
    return inst, res, info
end

"""Shrinkage prior on the barycenter center (Ledoit–Wolf flavor):
λ = 0 recovers the fixed point exactly; λ → ∞ pulls to the prior."""
function gw_shrinkage_demo(; d::Int = 2, seed::Int = 1)
    rng = MersenneTwister(seed)
    Sig = [(X = randn(rng, d, d); X * X' + 0.2I) for _ in 1:3]
    λw = [0.5, 0.3, 0.2]
    Sfp = gw_barycenter_fixedpoint(Sig, λw)
    P = tr(Sfp) / d .* Matrix{Float64}(I, d, d)
    for λ in (0.0, 0.5, 5.0)
        inst = gw_instance([d, d, d, d], [(1, 2), (1, 3), (1, 4)],
                           Dict(2 => Sig[1], 3 => Sig[2], 4 => Sig[3]);
                           weights = λw,
                           prior = (λ > 0 ? Dict(1 => (λ, P)) :
                                    Dict{Int, Tuple{Float64, Matrix{Float64}}}()))
        prob, info = build_gw_problem(inst)
        res = solve(prob, gw_settings())
        Sc = gw_smat(res.p[info.mucols[1]], d)
        @printf("  λ=%.1f: ‖Σc − Σfp‖ = %.4f   ‖Σc − prior‖ = %.4f   it=%d %s\n",
                λ, norm(Sc - Sfp), norm(Sc - P), res.niter, res.status)
    end
    println("  (oracle: 0.0000/0.675 → 0.339/0.344 → 0.615/0.061)")
end

# =====================================================================
# JuMP reference
# =====================================================================

using JuMP

function build_jump_gw(inst::GWInstance, optimizer)
    dims, edges = inst.dims, inst.edges
    model = Model(optimizer); set_silent(model)
    J = [@variable(model, [1:(dims[i] + dims[j]), 1:(dims[i] + dims[j])], PSD)
         for (i, j) in edges]
    μ = Dict{Int, Any}(k => inst.pins[k] for k in keys(inst.pins))
    for k in inst.free
        μ[k] = @variable(model, [1:dims[k], 1:dims[k]], PSD)
    end
    for (t, (i, j)) in enumerate(edges)
        di = dims[i]
        @constraint(model, J[t][1:di, 1:di] .== μ[i])
        @constraint(model, J[t][di+1:end, di+1:end] .== μ[j])
    end
    if isempty(inst.tg)
        @objective(model, Min,
                   sum(inst.weights[t] *
                       tr(gw_edge_cost(dims[i], dims[j], t in inst.anti) * J[t])
                       for (t, (i, j)) in enumerate(edges)) +
                   sum((λP = inst.prior[k];
                        0.5 * λP[1] * sum((μ[k] .- λP[2]) .^ 2))
                       for k in keys(inst.prior); init = 0.0))
    else
        @objective(model, Min,
                   sum((Wh = haskey(inst.Wt, t) ? gw_sqrtm(inst.Wt[t]) :
                             Matrix{Float64}(I, size(inst.tg[t])...);
                        0.5 * sum((Wh * (J[t] .- inst.tg[t]) * Wh) .^ 2))
                       for t in keys(inst.tg)))
    end
    return model, J, μ
end

# =====================================================================
# Benchmark: IPM vs Mosek
# =====================================================================

"""Benchmark IPM vs Mosek on Gaussian Wasserstein problems.
Dense Q (Higham weights) → IPM wins 5-40x; diagonal Q → Mosek wins."""
function run_gw_benchmark(; optimizer = nothing, dual_optimizer = nothing,
                          solver_name::String = "Mosek", nwarmup::Int = 2, nruns::Int = 5)
    cases = [
        # (builder, label, raug)
        # Transport (Q = εI, diagonal) — small d favors IPM
        (() -> begin
            rng = MersenneTwister(0)
            d = 3
            S1 = randn(rng, d, d); S1 = S1 * S1' + 0.3I
            S2 = randn(rng, d, d); S2 = S2 * S2' + 0.3I
            gw_instance([d, d], [(1, 2)], Dict(1 => S1, 2 => S2))
         end, "transport d=3", 1e1),
        # Higham completion (Q = W⊗W, dense) — IPM wins big
        (() -> begin
            d = 5
            tg = Dict(t => (X = randn(MersenneTwister(t), 2d, 2d);
                            X = X*X' + 0.5I; X / tr(X) * 2d) for t in 1:4)
            Wt = Dict(t => (W = randn(MersenneTwister(t+100), 2d, 2d);
                            W = W*W' + 0.3I) for t in 1:4)
            gw_instance(fill(d, 4), [(k, mod1(k+1, 4)) for k in 1:4],
                        Dict(k => (X = randn(MersenneTwister(k), d, d); X*X' + 0.2I) for k in 1:4);
                        tg = tg, Wt = Wt)
         end, "Higham C4 d=5", 1e1),
        (() -> begin
            d = 8
            tg = Dict(t => (X = randn(MersenneTwister(t), 2d, 2d);
                            X = X*X' + 0.5I; X / tr(X) * 2d) for t in 1:4)
            Wt = Dict(t => (W = randn(MersenneTwister(t+100), 2d, 2d);
                            W = W*W' + 0.3I) for t in 1:4)
            gw_instance(fill(d, 4), [(k, mod1(k+1, 4)) for k in 1:4],
                        Dict(k => (X = randn(MersenneTwister(k), d, d); X*X' + 0.2I) for k in 1:4);
                        tg = tg, Wt = Wt)
         end, "Higham C4 d=8", 1e1),
        # Grid with Higham
        (() -> begin
            d = 2; M = 6
            K = M * M
            dims = fill(d, K)
            edges = Tuple{Int,Int}[]
            for i in 1:M, j in 1:M
                v = (i-1)*M + j
                if j < M push!(edges, (v, v+1)) end
                if i < M push!(edges, (v, v+M)) end
            end
            pins = Dict(v => (X = randn(MersenneTwister(v), d, d); X*X' + 0.2I)
                        for v in [1, M, (M-1)*M+1, M*M])
            tg = Dict(t => (X = randn(MersenneTwister(t), 2d, 2d);
                            X = X*X' + 0.5I; X / tr(X) * 2d) for t in 1:length(edges))
            Wt = Dict(t => (W = randn(MersenneTwister(t+100), 2d, 2d);
                            W = W*W' + 0.3I) for t in 1:length(edges))
            gw_instance(dims, edges, pins; tg = tg, Wt = Wt)
         end, "Higham grid 6x6", 1e1),
    ]

    sname = rpad(solver_name, 9)
    sname_d = rpad(solver_name * "-D", 9)
    @printf("Config           raug   DOF   |V|   IPM(ms)   %s %s  P/IPM   D/IPM\n", sname, sname_d)
    println("-" ^ 85)

    for (builder, label, raug) in cases
        inst = builder()
        prob, info = build_gw_problem(inst)
        dof = sum(length(colrange(prob.B, v)) for v in 1:info.nv)

        for _ in 1:nwarmup
            solve(prob, gw_settings(raug = raug))
        end
        t_ipm = 0.0
        for _ in 1:nruns
            t_ipm += @elapsed solve(prob, gw_settings(raug = raug))
        end
        t_ipm = 1000.0 * t_ipm / nruns

        t_mosek, t_mosek_d = Inf, Inf
        if optimizer !== nothing
            m, _, _ = build_jump_gw(inst, optimizer)
            for _ in 1:nwarmup
                optimize!(m)
            end
            t_mosek = 0.0
            for _ in 1:nruns
                m, _, _ = build_jump_gw(inst, optimizer)
                t_mosek += @elapsed optimize!(m)
            end
            t_mosek = 1000.0 * t_mosek / nruns
        end
        if dual_optimizer !== nothing
            m, _, _ = build_jump_gw(inst, dual_optimizer)
            for _ in 1:nwarmup
                optimize!(m)
            end
            t_mosek_d = 0.0
            for _ in 1:nruns
                m, _, _ = build_jump_gw(inst, dual_optimizer)
                t_mosek_d += @elapsed optimize!(m)
            end
            t_mosek_d = 1000.0 * t_mosek_d / nruns
        end

        rat_p = t_mosek < Inf ? t_mosek / t_ipm : NaN
        rat_d = t_mosek_d < Inf ? t_mosek_d / t_ipm : NaN
        @printf("%-16s %6.0e %4d  %4d  %8.1f  %8.1f  %8.1f  %5.2fx  %5.2fx\n",
                label, raug, dof, info.nv, t_ipm, t_mosek, t_mosek_d, rat_p, rat_d)
    end
end

# =====================================================================
# Command-line interface
# =====================================================================

if abspath(PROGRAM_FILE) == @__FILE__
    include(joinpath(@__DIR__, "..", "benchmark_utils.jl"))
    using Dualization

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
    println("Gaussian Wasserstein benchmark")
    println("Solver: $solver_name")
    println()

    run_gw_benchmark(; optimizer, dual_optimizer, solver_name, nwarmup = opts.nwarmup, nruns = opts.nruns)
end
