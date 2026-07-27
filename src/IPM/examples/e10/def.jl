# =============================================================================
# e10/def.jl — Compositional Lyapunov certification on a K×K torus of subsystems
#          (Agler-split dense SDP with Gramian-weighted dense Q)
#
# Usage:
#   julia --project e10/run.jl              # Clarabel baseline
#   julia --project e10/run.jl --mosek      # + Mosek
#   julia --project e10/run.jl --quick      # quick (smaller) sweep
#
# Problem. A K×K torus of coupled subsystems, each with DENSE internal
# dynamics A_loc ∈ R^{d×d} (skew + dense SPD damping) and dense
# nearest-neighbor couplings ε·{Cp, Cm, Dp, Dm}, certified stable
# compositionally:
#
#   min  Σ_v [ tr(P_v) + (γ/2)‖W^{1/2} P_v W^{1/2}‖_F² ]
#   s.t. A′P + PA + I + S = 0,   P = diag(P_v),  P_v ⪰ 0,
#        S = Σ_edges embed(S_e),  S_e ⪰ 0  (2d×2d, one per torus edge),
#
# with W the local controllability Gramian (A_loc W + W A_loc′ + I = 0) —
# a canonical, dynamics-derived, fully dense weight. The quadratic term is
# the flagship (γ = 0.1): it makes the program strictly convex in P, the
# certificate UNIQUE (Tikhonov selection on the linear SDP's optimal
# face), and cross-solver gates comparable ON STATES — a first for the
# suite's SDP examples. γ → 0 recovers the pure linear certificate at
# measured rate γ^1.00.
#
# WHY THIS EXAMPLE EXISTS — the conjunction of the suite's two measured
# win-conditions, deliberately placed on the 2D topology where the
# quantum-marginal torus failed. (1) DENSE Q: dense 10-dim svec Gram
# blocks on every P_v stalk — dense Q is free for the sheaf solver (b³
# per block regardless) and a pure tax on baselines (dense cliques in
# their KKT fill, or dense factor rows after conic reformulation).
# (2) SMALL, BALANCED BLOCKS: vertex stalks svec 10 (P_v) and 36 (S_e)
# against row blocks of 10–16, so b_v ≈ b_e — neither side is thin, and
# no solver inherits a cubed-constant advantage (the quantum torus's
# 136-dim stalks against 10-dim rows handed the constraint-side
# factorizer a ~2500× smaller cubic constant; this example is the
# controlled correction, and the falsifiable prediction of the thin-side
# model is that the IPM is competitive HERE). Additionally, coker B = 0
# BY THE PRIVATE-SLOT THEOREM — every row group pins a fresh S_e block,
# so the sum-split never manufactures redundancy — verified by numerical
# rank below; this example is a redundancy CONTROL (vs a coboundary-constrained companion).
#
# Two new species. Sum-decomposition edges (vs E09's agreements,
# E02's partial traces, E03's jets): the Agler split of the banded slack S
# into clique-supported PSD pieces, which on the torus's non-chordal
# cyclic band pattern is a strict RESTRICTION — that is what a
# compositional certificate is, and the conservatism is measured, not
# hidden. A spectral (2D-DFT / group-averaging) oracle: the monolithic
# comparator (P block-diag, one big LMI) reduces exactly to per-wavenumber
# d-dim Hermitian LMIs at ANY K, so the conservatism gap curve is
# computable everywhere (measured per-cell at γ = 0: 8.11% at K = 3 →
# 3.36% → 2.90% → 3.10% at K = 8, converging to its continuum limit with
# mild wavenumber-aliasing wobble). E02's non-chordal argument, relocated
# to the constraint side.
#
# Sheaf structure. Stalks: K² certificate blocks P_v (SemidefiniteCone,
# svec dim 10) + 2K² coupling slacks Sh_v, Sv_v (SemidefiniteCone, svec
# dim 36; edge indexed by its left/bottom endpoint). Edges: per vertex,
# one vertex equation (10 rows; touches P_v and FOUR incident S stalks —
# the Lyapunov map A_loc′(·) + (·)A_loc is a dense 10×10 block) and two
# off-diagonal coupling equations (16 rows each; dense d×d maps through
# ε·C's, plus the S_e off-diagonal extraction). rhs: svec(−I) on vertex
# equations. Q: γ·(Gramian Gram) + ridge on P stalks (dense SPD), ridge
# on S stalks. c: svec(I) on P stalks.
#
# Validation (tools/torus_lyapunov_oracle.py, self-contained numpy +
# scipy + cvxpy, fully executed; d = 4, ε = 0.15, γ = 0.1; note the
# Julia instance below is a different random draw — all gates re-derive
# their certificates in-process):
#   * setup: Gramian residual 8.9e-16; W dense SPD; global A Hurwitz.
#   * COKER = 0 verified by rank at K = 3, 4 (B is 378×738 / 672×1312,
#     aspect 0.512, full row rank) — the private-slot theorem.
#   * Agler-split FEASIBILITY at flagship coupling (conservatism can mean
#     "no certificate exists"; here one does, v̄ = 1.50529852).
#   * SYMMETRY ORACLE: |v(K) − K²·v̄| = 1.7e-12 / 3.2e-12 at K = 3 / 4,
#     and STATE UNIQUENESS max|P_v − P̄| = 3.3e-9 / 6.2e-9. (K = 2 is
#     degenerate — wraparound doubles couplings — documented, not gated.)
#   * SPECTRAL ORACLE: DFT vs direct monolithic |dv| = 5.1e-8 at K = 3;
#     conservatism curve as quoted above.
#   * DYNAMICAL CERTIFICATE: V = x′Px monotone decreasing along
#     ẋ = Ax and V(t) ≤ V(0)·exp(−t/λmax(P)) with max ratio 1.000.
#   * INSTABILITY LOCALIZATION: destabilize one subsystem → global A
#     unstable, program infeasible, and the minimal per-vertex relaxation
#     concentrates EXACTLY there (t_broken = 1.036, elsewhere 0.0e0;
#     healthy control all zeros).
#   * DEFORMATION: rate γ^1.00 (predicted 1).
#   * FORMULATION CERTIFICATE: svec assembly (this file's layout) vs
#     matrix-native |dv| = 3.5e-11, states 1.8e-8.
#   * CROSS-SOLVER ON STATES: Clarabel vs SCS |dP| = 1.1e-8.
#
# Dials: torus size K (primary; DOF = 82K², N1 = 42K²) and subsystem
# dim d (secondary — fattens BOTH sides together, preserving b_v ≈ b_e);
# γ and ε are physics knobs.
#
# STATUS: authored against the E02 build API; the oracle is fully
# executed but this file has NOT yet been run against the IPM. First-run
# checklist: (1) this is the falsifiable test of the thin-side account —
# dense Q + small balanced blocks on a 2D torus is predicted
# IPM-competitive; if it loses like the quantum torus did, the model
# needs revision; (2) B has FULL row rank by design (the redundancy
# control for the coboundary case); (3) confirm raug (1e2); (4) vertex equations touch
# five stalks on one edge id — first five-block rows in the suite;
# (5) fill the sample-run table.
# =============================================================================

using AppleAccelerate
using LinearAlgebra, SparseArrays, Printf, Random

include("../utils.jl")

using CellularSheaves.IPM: SemidefiniteCone
using CellularSheaves.BlockSparseArrays: rowrange

const OPTS = parse_args(ARGS)
const TOL = OPTS.tol
const NRUNS = 5
const RIDGE = 1e-9

# -----------------------------------------------------------------------------
# Problem definition
# -----------------------------------------------------------------------------

struct TLInstance
    d::Int; eps::Float64; γ::Float64; seed::Int
end

tl_instance(; d = 4, eps = 0.15, γ = 0.1, seed = 11) =
    TLInstance(d, eps, γ, seed)

svecdim(D::Int) = D * (D + 1) ÷ 2

function tosvec(M::AbstractMatrix)
    D = size(M, 1); v = zeros(svecdim(D)); k = 0
    for c in 1:D, r in c:D
        k += 1
        v[k] = r == c ? M[c, c] : sqrt(2.0) * M[r, c]
    end
    v
end

function tosmat(v::AbstractVector, D::Int)
    M = zeros(D, D); k = 0
    for c in 1:D, r in c:D
        k += 1
        if r == c
            M[c, c] = v[k]
        else
            M[r, c] = M[c, r] = v[k] / sqrt(2.0)
        end
    end
    M
end

function make_system(inst::TLInstance)
    d = inst.d
    rng = MersenneTwister(inst.seed)
    G = randn(rng, d, d)
    S0 = randn(rng, d, d); S0 .-= S0'
    Aloc = S0 .- (Matrix{Float64}(I, d, d) .+ 0.3 .* (G * G'))
    Cp = randn(rng, d, d) ./ sqrt(d)
    Cm = randn(rng, d, d) ./ sqrt(d)
    Dp = randn(rng, d, d) ./ sqrt(d)
    Dm = randn(rng, d, d) ./ sqrt(d)
    W = lyap(Aloc, Matrix{Float64}(I, d, d))     # A W + W A' + I = 0
    Wh = Matrix(sqrt(Symmetric(W)))
    sp = svecdim(d); ss = svecdim(2 * d)
    return (; d, eps = inst.eps, γ = inst.γ, Aloc, Cp, Cm, Dp, Dm, W, Wh,
            sp, ss)
end

"Matrix of a linear map given on svec coordinates."
function lin_map(f, din::Int)
    cols = Vector{Vector{Float64}}(undef, din)
    for k in 1:din
        e = zeros(din); e[k] = 1.0
        cols[k] = vec(f(e))
    end
    reduce(hcat, cols)
end

function build_maps(sys)
    d, sp, ss = sys.d, sys.sp, sys.ss
    LY = lin_map(e -> tosvec(sys.Aloc' * tosmat(e, d) .+ tosmat(e, d) * sys.Aloc), sp)
    GG = lin_map(e -> tosvec(sys.Wh * tosmat(e, d) * sys.Wh), sp)
    E11 = lin_map(e -> tosvec(tosmat(e, 2 * d)[1:d, 1:d]), ss)
    E22 = lin_map(e -> tosvec(tosmat(e, 2 * d)[d + 1:2 * d, d + 1:2 * d]), ss)
    EOD = lin_map(e -> tosmat(e, 2 * d)[1:d, d + 1:2 * d], ss)
    odL(C) = lin_map(e -> sys.eps .* (C' * tosmat(e, d)), sp)
    odR(C) = lin_map(e -> sys.eps .* (tosmat(e, d) * C), sp)
    return (; LY, GG, E11, E22, EOD,
            Hv = odR(sys.Cp), Hw = odL(sys.Cm),
            Vv = odR(sys.Dp), Vw = odL(sys.Dm))
end

function global_A(sys, K::Int)
    d = sys.d
    widx(x, y) = (mod1(x, K) - 1) * K + mod1(y, K)
    A = zeros(K * K * d, K * K * d)
    bl(v) = ((v - 1) * d + 1):(v * d)
    for x in 1:K, y in 1:K
        v = widx(x, y)
        A[bl(v), bl(v)] .= sys.Aloc
        A[bl(v), bl(widx(x + 1, y))] .+= sys.eps .* sys.Cp
        A[bl(v), bl(widx(x - 1, y))] .+= sys.eps .* sys.Cm
        A[bl(v), bl(widx(x, y + 1))] .+= sys.eps .* sys.Dp
        A[bl(v), bl(widx(x, y - 1))] .+= sys.eps .* sys.Dm
    end
    A
end

# -----------------------------------------------------------------------------
# Build (K² P stalks + 2K² S stalks; vertex + coupling equations)
# -----------------------------------------------------------------------------
#
# Stalk order: P_1..P_{K²} | Sh_1..Sh_{K²} | Sv_1..Sv_{K²}, edges indexed
# by their left/bottom endpoint. K = 0 means the single-cell reduced
# program is requested (build_cell): same machinery, one of everything.

function build_torus(sys, K::Int; γ = sys.γ)
    W = K * K
    d, sp, ss = sys.d, sys.sp, sys.ss
    mp = build_maps(sys)
    widx(x, y) = (mod1(x, K) - 1) * K + mod1(y, K)
    pstalk(v) = v
    hstalk(v) = W + v
    vstalk(v) = 2 * W + v

    row_ids, col_ids, blks = Int[], Int[], Matrix{Float64}[]
    rhs_val = Dict{Int, Vector{Float64}}()
    place!(e, v, mat) = (push!(row_ids, e); push!(col_ids, v);
                         push!(blks, Matrix{Float64}(mat)))
    e = 0
    for x in 1:K, y in 1:K
        v = widx(x, y)
        vr = widx(x + 1, y); vu = widx(x, y + 1)
        vl = widx(x - 1, y); vd = widx(x, y - 1)
        # vertex equation (five stalks on one edge id)
        e += 1
        place!(e, pstalk(v), mp.LY)
        place!(e, hstalk(v), mp.E11)
        place!(e, vstalk(v), mp.E11)
        place!(e, hstalk(vl), mp.E22)
        place!(e, vstalk(vd), mp.E22)
        rhs_val[e] = tosvec(-Matrix{Float64}(I, d, d))
        # horizontal coupling equation
        e += 1
        place!(e, pstalk(v), mp.Hv)
        place!(e, pstalk(vr), mp.Hw)
        place!(e, hstalk(v), mp.EOD)
        # vertical coupling equation
        e += 1
        place!(e, pstalk(v), mp.Vv)
        place!(e, pstalk(vu), mp.Vw)
        place!(e, vstalk(v), mp.EOD)
    end

    B = blocksparse(row_ids, col_ids, blks)
    g = zeros(size(B, 1))
    for (edge, val) in rhs_val
        g[rowrange(B, edge)] .= val
    end
    Q = IPM.allocblockdiag(B); fill!(Q, 0)
    Qp = γ .* (mp.GG' * mp.GG) .+ RIDGE .* Matrix{Float64}(I, sp, sp)
    Isr = RIDGE .* Matrix{Float64}(I, ss, ss)
    c = zeros(size(B, 2))
    ci = tosvec(Matrix{Float64}(I, d, d))
    for v in 1:W
        block(Q, pstalk(v), pstalk(v), pstalk(v)) .= Qp
        block(Q, hstalk(v), hstalk(v), hstalk(v)) .= Isr
        block(Q, vstalk(v), vstalk(v), vstalk(v)) .= Isr
        c[colrange(B, pstalk(v))] .= .-ci
    end
    cones = vcat(IPM.AbstractCone[SemidefiniteCone() for _ in 1:W],
                 IPM.AbstractCone[SemidefiniteCone() for _ in 1:(2 * W)])
    prob = IPMProblem(Q, B, c, g, cones)
    prob, (; K, W, γ, sys)
end

"Single-cell reduced program (the translation-invariant limit)."
function build_cell(sys; γ = sys.γ)
    d, sp, ss = sys.d, sys.sp, sys.ss
    mp = build_maps(sys)
    row_ids, col_ids, blks = Int[], Int[], Matrix{Float64}[]
    place!(e, v, mat) = (push!(row_ids, e); push!(col_ids, v);
                         push!(blks, Matrix{Float64}(mat)))
    # stalks: 1 = P, 2 = Sh, 3 = Sv
    place!(1, 1, mp.LY)
    place!(1, 2, mp.E11 .+ mp.E22)
    place!(1, 3, mp.E11 .+ mp.E22)
    place!(2, 1, mp.Hv .+ mp.Hw)
    place!(2, 2, mp.EOD)
    place!(3, 1, mp.Vv .+ mp.Vw)
    place!(3, 3, mp.EOD)
    B = blocksparse(row_ids, col_ids, blks)
    g = zeros(size(B, 1))
    g[rowrange(B, 1)] .= tosvec(-Matrix{Float64}(I, d, d))
    Q = IPM.allocblockdiag(B); fill!(Q, 0)
    block(Q, 1, 1, 1) .= γ .* (mp.GG' * mp.GG) .+ RIDGE .* Matrix{Float64}(I, sp, sp)
    block(Q, 2, 2, 2) .= RIDGE .* Matrix{Float64}(I, ss, ss)
    block(Q, 3, 3, 3) .= RIDGE .* Matrix{Float64}(I, ss, ss)
    c = zeros(size(B, 2))
    c[colrange(B, 1)] .= .-tosvec(Matrix{Float64}(I, d, d))
    prob = IPMProblem(Q, B, c, g,
                      IPM.AbstractCone[SemidefiniteCone() for _ in 1:3])
    prob
end

extract_Ps(prob, ctx, pv) =
    [tosmat(pv[colrange(prob.B, v)], ctx.sys.d) for v in 1:ctx.W]

function tl_objective(ctx, pv, prob)
    sys = ctx.sys
    F = 0.0
    for v in 1:ctx.W
        P = tosmat(pv[colrange(prob.B, v)], sys.d)
        F += tr(P) + 0.5 * ctx.γ * sum(abs2, sys.Wh * P * sys.Wh)
        F += 0.5 * RIDGE * sum(abs2, pv[colrange(prob.B, v)])
    end
    for j in (ctx.W + 1):(3 * ctx.W)
        F += 0.5 * RIDGE * sum(abs2, pv[colrange(prob.B, j)])
    end
    F
end

# -----------------------------------------------------------------------------
# Baseline (same conic program via JuMP, matrix-native)
# -----------------------------------------------------------------------------

function build_jump_torus(sys, K::Int, optimizer; γ = sys.γ)
    model = Model(optimizer); set_silent(model)
    d = sys.d
    widx(x, y) = (mod1(x, K) - 1) * K + mod1(y, K)
    P = Dict(v => @variable(model, [1:d, 1:d], PSD) for v in 1:K * K)
    Sh = Dict(v => @variable(model, [1:2 * d, 1:2 * d], PSD) for v in 1:K * K)
    Sv = Dict(v => @variable(model, [1:2 * d, 1:2 * d], PSD) for v in 1:K * K)
    Id = Matrix{Float64}(I, d, d)
    for x in 1:K, y in 1:K
        v = widx(x, y)
        vr = widx(x + 1, y); vu = widx(x, y + 1)
        vl = widx(x - 1, y); vd = widx(x, y - 1)
        @constraint(model,
            sys.Aloc' * P[v] .+ P[v] * sys.Aloc .+ Id
            .+ Sh[v][1:d, 1:d] .+ Sv[v][1:d, 1:d]
            .+ Sh[vl][d + 1:2 * d, d + 1:2 * d]
            .+ Sv[vd][d + 1:2 * d, d + 1:2 * d] .== 0)
        @constraint(model,
            sys.eps .* (sys.Cm' * P[vr] .+ P[v] * sys.Cp)
            .+ Sh[v][1:d, d + 1:2 * d] .== 0)
        @constraint(model,
            sys.eps .* (sys.Dm' * P[vu] .+ P[v] * sys.Dp)
            .+ Sv[v][1:d, d + 1:2 * d] .== 0)
    end
    ridge_term(M, n) = sum(M[i, j]^2 for i in 1:n, j in 1:n)
    @objective(model, Min,
        sum(tr(P[v])
            + 0.5 * γ * sum((sys.Wh * P[v] * sys.Wh)[i, j]^2
                            for i in 1:d, j in 1:d)
            + 0.5 * RIDGE * ridge_term(P[v], d)
            + 0.5 * RIDGE * ridge_term(Sh[v], 2 * d)
            + 0.5 * RIDGE * ridge_term(Sv[v], 2 * d)
            for v in 1:K * K))
    model, P
end

# -----------------------------------------------------------------------------
# Gate tests
# -----------------------------------------------------------------------------

function test_maps(sys)
    rng = MersenneTwister(0)
    mp = build_maps(sys)
    X = randn(rng, sys.d, sys.d); X = X + X'
    e1 = maximum(abs.(tosmat(mp.LY * tosvec(X), sys.d)
                      .- (sys.Aloc' * X .+ X * sys.Aloc)))
    e2 = maximum(abs.(tosmat(mp.GG * tosvec(X), sys.d)
                      .- sys.Wh * X * sys.Wh))
    gram = maximum(abs.(sys.Aloc * sys.W .+ sys.W * sys.Aloc'
                        .+ Matrix{Float64}(I, sys.d, sys.d)))
    @assert e1 < 1e-12 && e2 < 1e-12 && gram < 1e-12
    println("  [PASS] map identities (Lyapunov $(round(e1, sigdigits = 2)), ",
        "Gramian-Gram $(round(e2, sigdigits = 2))); Gramian residual ",
        "$(round(gram, sigdigits = 2)); W dense SPD ",
        "(min eig $(round(eigmin(Symmetric(sys.W)), digits = 3)))")
end

function test_coker(sys)
    prob, _ = build_torus(sys, 3)
    Bd = Matrix(sparse(prob.B))
    N1, N0 = size(Bd)
    r = rank(Bd)
    @assert r == N1 "coker ≠ 0: rank $r of $N1 rows (private-slot theorem!)"
    println("  [PASS] coker B = 0 (private-slot theorem): B is $(N1)×$(N0) ",
        "(aspect $(round(N1 / N0, digits = 3)), wide), FULL row rank — ",
        "a redundancy control")
end

function test_stability(sys)
    A = global_A(sys, 3)
    s = maximum(real.(eigvals(A)))
    @assert s < 0 "global A not Hurwitz: $s"
    println("  [PASS] global A Hurwitz (max Re eig = $(round(s, digits = 3)))")
end

function test_symmetry(sys, settings)
    probc = build_cell(sys)
    resc = solve(probc, settings)
    @assert resc.status in (OPTIMAL, NEAR_OPTIMAL) "cell: $(resc.status)"
    vbar = ipm_objective(probc, resc)
    Pbar = tosmat(resc.p[colrange(probc.B, 1)], sys.d)
    prob, ctx = build_torus(sys, 3)
    res = solve(prob, settings)
    @assert res.status in (OPTIMAL, NEAR_OPTIMAL) "K=3: $(res.status)"
    v3 = ipm_objective(prob, res)
    d = abs(v3 - 9 * vbar)
    Ps = extract_Ps(prob, ctx, res.p)
    dev = maximum(maximum(abs.(P .- Pbar)) for P in Ps)
    @assert d < 1e-4 "symmetry oracle: $d (oracle: 1.7e-12)"
    @assert dev < 1e-4 "state uniqueness: $dev (oracle: 3.3e-9)"
    println("  [PASS] symmetry oracle: |v(3) − 9v̄| = $(round(d, sigdigits = 2)); ",
        "state uniqueness max|P_v − P̄| = $(round(dev, sigdigits = 2)) ",
        "(v̄ = $(round(vbar, digits = 6)); Agler split feasible)")
    prob, ctx, res, Ps
end

function test_dynamical(sys, Ps)
    K = 3; d = sys.d
    A = global_A(sys, K)
    Pg = zeros(K * K * d, K * K * d)
    for v in 1:K * K
        Pg[((v - 1) * d + 1):(v * d), ((v - 1) * d + 1):(v * d)] .= Ps[v]
    end
    λ = eigmax(Symmetric(Pg))
    rng = MersenneTwister(0)
    worst = 0.0
    ts = range(0.0, 6.0; length = 25)
    for _ in 1:4
        x0 = randn(rng, K * K * d)
        V = [begin xt = exp(A .* t) * x0; xt' * Pg * xt end for t in ts]
        @assert all(diff(V) .< 1e-10) "V not monotone"
        worst = max(worst, maximum(V ./ (V[1] .* exp.(-collect(ts) ./ λ))))
    end
    @assert worst <= 1.0 + 1e-9
    println("  [PASS] dynamical certificate: V = x′Px monotone decreasing; ",
        "V(t) ≤ V(0)exp(−t/λmax(P)) (max ratio $(round(worst, digits = 3)))")
end

function test_deformation(sys, settings)
    vals = Float64[]
    for γ in (0.4, 0.1)
        probc = build_cell(sys; γ = γ)
        resc = solve(probc, settings)
        @assert resc.status in (OPTIMAL, NEAR_OPTIMAL)
        push!(vals, ipm_objective(probc, resc))
    end
    probc0 = build_cell(sys; γ = 0.0)
    resc0 = solve(probc0, settings)
    @assert resc0.status in (OPTIMAL, NEAR_OPTIMAL)
    v0 = ipm_objective(probc0, resc0)
    rate = log(abs(vals[2] - v0) / abs(vals[1] - v0)) / log(0.1 / 0.4)
    @assert 0.7 < rate < 1.3 "deformation rate γ^$(round(rate, digits = 2)) (oracle: 1.00)"
    println("  [PASS] deformation γ → 0: rate γ^$(round(rate, digits = 2)) ",
        "(predicted 1); pure linear certificate v₀ = $(round(v0, digits = 6))")
end

function test_objective_identity(prob, ctx, res)
    o_direct = tl_objective(ctx, res.p, prob)
    o_ipm = ipm_objective(prob, res)
    rel = abs(o_direct - o_ipm) / max(abs(o_direct), 1.0)
    @assert rel < 1e-6 "assembly mismatch: $o_direct vs $o_ipm"
    println("  [PASS] objective identity (rel $(round(rel, sigdigits = 2)))")
end

function test_ipm_vs_clarabel(sys, ctx, Ps)
    model, Pj = build_jump_torus(sys, ctx.K, clarabel_opt(; tol = TOL))
    optimize!(model)
    dx = maximum(maximum(abs.(value.(Pj[v]) .- Ps[v])) for v in 1:ctx.W)
    @assert dx < 1e-4 "IPM vs Clarabel states: $dx (oracle cross-state: 1.1e-8)"
    println("  [PASS] IPM vs Clarabel ON STATES (unique optimum): ",
        "max|ΔP| = $(round(dx, sigdigits = 2))")
end

# -----------------------------------------------------------------------------
# Sweep
# -----------------------------------------------------------------------------

sizes() = OPTS.tiny ? [4] :
    OPTS.quick ? [4, 6, 8] :
    [4, 6, 8, 12]

function run()
    println("\n", "=" ^ 78)
    println("  E10: compositional Lyapunov on a K×K torus ",
        "(dials: K; d fattens both sides)")
    println("=" ^ 78)

    ipm_settings = IPMSettings{Float64}(
        feas_tol = TOL, gap_tol = TOL, itmax = 200)
    hsd_settings = HSDSettings{Float64}(
        feas_tol = TOL, gap_tol = TOL, itmax = 200)

    inst = tl_instance()
    sys = make_system(inst)

    println("\n  Gate tests (d = 4, ε = 0.15, γ = 0.1):")
    test_maps(sys)
    test_coker(sys)
    test_stability(sys)
    prob3, ctx3, res3, Ps3 = test_symmetry(sys, ipm_settings)
    test_dynamical(sys, Ps3)
    test_deformation(sys, ipm_settings)
    test_objective_identity(prob3, ctx3, res3)
    test_ipm_vs_clarabel(sys, ctx3, Ps3)
    println()

    cla_opt = clarabel_opt(; tol = TOL)
    msk_opt = OPTS.mosek ? mosek_opt(; tol = TOL) : nothing

    println("  Dial: torus size K (DOF = 82K², N1 = 42K², coker = 0, ",
        "b_v ≈ b_e)")
    rows = []
    for K in sizes()
        prob, ctx = build_torus(sys, K)
        stats = problem_stats(prob)
        @printf("  K=%-3d dof=%-6d n1=%-6d blk=%-4.0f  ",
            K, stats.N0, stats.N1, stats.med_block)
        m_ipm = measure_ipm(prob, ipm_settings; nruns = NRUNS)
        m_hsd = measure_ipm(prob, hsd_settings; nruns = NRUNS)
        m_cla = measure_jump(() -> first(build_jump_torus(sys, K, cla_opt));
                             nruns = NRUNS)
        # Mosek benchmarked in DUAL form (Dualization.jl): ~1.4x faster than primal here.
        m_msk = msk_opt !== nothing ?
            measure_jump(() -> first(build_jump_torus(sys, K, Dualization.dual_optimizer(msk_opt)));
                         nruns = NRUNS) :
            (t = NaN, status = "", obj = NaN)
        ratio(b, base) = isfinite(b.t) && isfinite(base.t) ? b.t / base.t : NaN
        fmt_ratio(b, base) = isnan(ratio(b, base)) ? "—" : @sprintf("%.2fx", ratio(b, base))
        println("IPM ", fmt_time(m_ipm), "  HSD ", fmt_time(m_hsd),
            " (", fmt_ratio(m_hsd, m_ipm), ")  Cla ", fmt_time(m_cla),
            " (", fmt_ratio(m_cla, m_ipm), ")  Msk ", fmt_time(m_msk),
            " (", fmt_ratio(m_msk, m_ipm), ")")
        push!(rows, (dof = stats.N0, ipm = m_ipm, hsd = m_hsd, cla = m_cla, msk = m_msk))
    end
    dofs = [r.dof for r in rows]
    println("\nFitted log-log slopes (time vs DOF):")
    for (name, get_t) in [("IPM", r -> r.ipm.t), ("HSD", r -> r.hsd.t), ("Clarabel", r -> r.cla.t),
                          ("Mosek", r -> r.msk.t)]
        sl = loglog_slope(dofs, [get_t(r) for r in rows])
        print("  $name: ", isnan(sl) ? "n/a" : @sprintf("DOF^%.2f", sl))
    end
    println()

    # ---- secondary dial: subsystem dim d at fixed K (fattens BOTH sides;
    # baselines pay per nonzero ~d², the IPM per block — prediction: the
    # Mosek ratio improves with d)
    Kd = OPTS.quick || OPTS.tiny ? 4 : 6
    println("\n  Fatness dial: subsystem dim d at fixed K = $Kd")
    for dd in (OPTS.tiny ? (4,) : (4, 6, 8))
        instd = tl_instance(; d = dd)
        sysd = make_system(instd)
        probd, ctxd = build_torus(sysd, Kd)
        statsd = problem_stats(probd)
        @printf("  d=%-3d dof=%-6d n1=%-6d blk=%-4.0f  ",
            dd, statsd.N0, statsd.N1, statsd.med_block)
        m_ipm = measure_ipm(probd, ipm_settings; nruns = NRUNS)
        m_hsd = measure_ipm(probd, hsd_settings; nruns = NRUNS)
        m_cla = measure_jump(() -> first(build_jump_torus(sysd, Kd, cla_opt));
                             nruns = NRUNS)
        # Mosek benchmarked in DUAL form (Dualization.jl): ~1.4x faster than primal here.
        m_msk = msk_opt !== nothing ?
            measure_jump(() -> first(build_jump_torus(sysd, Kd, Dualization.dual_optimizer(msk_opt)));
                         nruns = NRUNS) :
            (t = NaN, status = "", obj = NaN)
        ratio(b, base) = isfinite(b.t) && isfinite(base.t) ? b.t / base.t : NaN
        fmt_ratio(b, base) = isnan(ratio(b, base)) ? "—" : @sprintf("%.2fx", ratio(b, base))
        println("IPM ", fmt_time(m_ipm), "  HSD ", fmt_time(m_hsd),
            " (", fmt_ratio(m_hsd, m_ipm), ")  Cla ", fmt_time(m_cla),
            " (", fmt_ratio(m_cla, m_ipm), ")  Msk ", fmt_time(m_msk),
            " (", fmt_ratio(m_msk, m_ipm), ")")
    end
end

