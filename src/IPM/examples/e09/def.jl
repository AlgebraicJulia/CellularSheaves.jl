# =============================================================================
# e09/def.jl — Robust fixed-interval smoothing under burst-corrupted measurements
#          (dense-Q + SOC habitat: the Kalman smoother meets E06's sqrt loss)
#
# Usage:
#   julia --project e09/run.jl              # Clarabel baseline
#   julia --project e09/run.jl --mosek      # + Mosek
#   julia --project e09/run.jl --quick      # quick (smaller) sweep
#
# Problem. Linear-Gaussian state-space smoothing over a horizon T,
#
#     x_{t+1} = A x_t + w_t,  w_t ~ N(0, W),   x_0 ~ N(μ0, P0)
#     y_t     = C x_t + v_t,  v_t ~ N(0, V),   t = 0..T,
#
# with the SQUARE-ROOT measurement loss (E06's loss, dynamic edition):
#
#     min_x  ½‖L0⁻¹(x_0 − μ0)‖² + ½ Σ_t ‖Lw⁻¹(x_{t+1} − A x_t)‖²
#            + λ Σ_t ‖Lv⁻¹(y_t − C x_t)‖         (norms, NOT squares).
#
# The sqrt loss bounds each timestep's pull on the trajectory at λ, so
# measurement BURSTS — contiguous windows where every sensor is corrupted
# (telemetry dropouts, EMI) — cannot drag the estimate: the chain coasts
# across the outage on dynamics alone. With the quadratic loss the same
# program IS the Kalman/RTS smoother, which is what makes this family an
# analytic oracle. Lineage: Rauch–Tung–Striebel (1965); Van Loan (1978);
# Aravkin–Burke–Ljung–Lozano–Pillonetto, "Generalized Kalman smoothing"
# (Automatica 2017).
#
# Dense by nature. A = exp(Ac·Δt) for a damped spring-mass chain — the
# matrix exponential is dense even though Ac is tridiagonal; W comes from
# the Van Loan block exponential (dense for the same reason); C is dense
# sensing rows; V = σ²(FF' + I) is a dense factor model (common-mode
# sensor drift, cf. E08's factor Σ), so M = Lv⁻¹C is dense by nature.
# Measured densities at 1e-8 of the max entry: A 0.986, W 1.0, V 1.0,
# M 1.0.
#
# Sheaf structure (dense-Q + SOC — a habitat cell no prior example
# occupies: E08 is dense-Q + pow, E06 is SOC with Q ≡ 0). Vertices t = 1..T
# are OVERLAPPING PAIR STALKS z_t = (x_{t-1}, x_t) ∈ R^{2n} (CofreeCone)
# carrying the dense process Gram
#     G = [A'ΩA  −A'Ω; −ΩA  Ω],   Ω = W⁻¹
# (+ prior P0⁻¹ on the first block of z_1) — E01's overlapping-patch
# discipline for coupled quadratics — with n-row agreement edges between
# consecutive pairs. Each timestep hangs one SOC stalk (u_τ, w_τ, s_τ) of
# dim m+2 off the chain: m dense tie rows w_τ = M x_τ − b_τ (nonzero rhs,
# E06's pattern) and a one-row pin s_τ = δ; the objective adds λ·u_τ.
# δ = 0 is the flagship sqrt loss; δ > 0 with λ = δ is the pseudo-Huber
# deformation ψ_δ(r) = δ(‖(r, δ)‖ − δ) → ‖r‖²/2, whose limit is EXACTLY
# the RTS smoother. No Gram whitening is needed (unlike E06): the residual
# dim m is small and the tie is dense by nature already.
#
# Validation (tools/robust_smoother_oracle.py, self-contained
# numpy + scipy + cvxpy, fully executed; figures at n = 12, m = 4,
# T = 200, λ = 1.5, 3 bursts of 8 steps):
#   * ANALYTIC ORACLE (the suite's third, after E01's MC and E08's ODE):
#     three independent algorithms agree — RTS recursion vs banded
#     information-form solve 2.0e-14; conic quadratic MAP vs RTS 4.6e-15.
#   * FORMULATION CERTIFICATE: the pair-stalk + agreement + SOC-tie
#     assembly (exactly what this file builds) equals the native norm-sum
#     formulation to ‖Δx‖∞ = 7.4e-7 (rel Δv 1.7e-12); the pair-Gram split
#     identity holds at random points to 4.7e-16.
#   * DEFORMATION: pseudo-Huber x(δ) → RTS at measured rate δ^-1.98
#     (predicted −2), ‖Δx‖∞ 1.4e-1 → 2.3e-3 over δ = 2 → 16.
#   * ESCALATION (corruption ×1 → ×128): robust clean-region state RMSE
#     saturates 0.268 → 0.339 while RTS blows up 0.261 → 2.158; overall
#     RMSE at ×128: robust 0.554 vs RTS 5.644 (10.2×). Clean-data price
#     +2.7%.
#   * SUPPORT RECOVERY + INFLUENCE CAP: the top-24 whitened residuals are
#     exactly the 24 corrupted steps (margin ×2.8); the quadratic loss
#     lets corrupted steps pull with ‖r‖ up to 95.9 vs the robust cap
#     λ = 1.5 (64×).
#   * KKT: subgradient-selection stationarity residual 1.35e-4
#     (70 interpolation-active steps ‖w_τ‖ = 0, E05's discipline).
#   * CLASSICAL (DARE): boundary-perturbation influence on interior states
#     decays at e^-0.2881/step vs the closed-loop ρ((I−KC)A) = 0.7519
#     predicted by the Riccati equation (log −0.2851) — ratio 1.010.
#   * CROSS-SOLVER: Clarabel vs SCS ‖Δx‖∞ = 2.1e-6.
# The gate tests below re-derive the fast subset in-process, including the
# analytic RTS gate via native filter + information-form references.
#
# Dial: horizon T at fixed n = 12, m = 4 (vertex count grows, per-block
# work constant, cf. E06/E08). DOF = 30T + 6; N1 = 17T − 7.
#
# STATUS: authored against the E06/E08 build API; the oracle is fully
# executed but this file has NOT yet been run against the IPM. First-run
# checklist: (1) confirm raug (1e2 chosen to match E08's dense-Q setting);
# (2) watch cone boundary behavior at interpolation-active steps where
# ‖w_τ‖ → 0 (the E06 regime); (3) fill the sample-run table below.
# =============================================================================

using AppleAccelerate
using LinearAlgebra, SparseArrays, Printf, Random

include("../utils.jl")

using CellularSheaves.IPM: SecondOrderCone, CofreeCone
using CellularSheaves.BlockSparseArrays: rowrange

const OPTS = parse_args(ARGS)
const TOL = OPTS.tol
const NRUNS = 5
const RIDGE = 1e-8

# -----------------------------------------------------------------------------
# Problem definition
# -----------------------------------------------------------------------------

struct SmootherInstance
    K::Int; m::Int                          # masses (n = 2K states), sensors
    dt::Float64; k::Float64; kg::Float64    # step, coupling, ground stiffness
    damp::Float64; q::Float64               # damping, disturbance intensity
    sig::Float64; sigf::Float64             # sensor noise scale, factor scale
    λ::Float64                              # sqrt-loss weight
    burst_len::Int; corrupt_scale::Float64
    seed_sys::Int; seed_sim::Int
end

smoother_instance(; K = 6, m = 4, dt = 0.5, k = 4.0, kg = 1.0, damp = 0.4,
                  q = 0.5, sig = 0.12, sigf = 0.7, λ = 1.5, burst_len = 8,
                  corrupt_scale = 32.0, seed_sys = 7, seed_sim = 3) =
    SmootherInstance(K, m, dt, k, kg, damp, q, sig, sigf, λ, burst_len,
                     corrupt_scale, seed_sys, seed_sim)

"Spring-mass chain + Van Loan discretization + dense factor V. All dense by nature."
function make_system(inst::SmootherInstance)
    K, m, n = inst.K, inst.m, 2 * inst.K
    rng = MersenneTwister(inst.seed_sys)
    Kmat = SymTridiagonal(fill(2 * inst.k + inst.kg, K), fill(-inst.k, K - 1))
    IK = Matrix{Float64}(I, K, K)
    Ac = [zeros(K, K) IK; (-Matrix(Kmat)) (-inst.damp .* IK)]
    Wc = zeros(n, n); Wc[K + 1:n, K + 1:n] .= inst.q .* IK
    # Van Loan (1978): F = exp([-Ac Wc; 0 Ac']·dt); A = F22', W = A·F12
    F = exp([-Ac Wc; zeros(n, n) Ac'] .* inst.dt)
    A = Matrix(F[n + 1:2n, n + 1:2n]')
    W = A * F[1:n, n + 1:2n]; W = Symmetric(0.5 .* (W .+ W'))
    C = randn(rng, m, n) ./ sqrt(n)
    Fn = inst.sigf .* randn(rng, m, 2)
    V = Symmetric(inst.sig^2 .* (Fn * Fn' .+ Matrix{Float64}(I, m, m)))
    P0 = Symmetric(0.5 .* Matrix{Float64}(I, n, n))
    μ0 = zeros(n)
    Lw = Matrix(cholesky(W).L); Lv = Matrix(cholesky(V).L)
    L0 = Matrix(cholesky(P0).L)
    M = Lv \ C                              # dense by nature
    return (; n, m, A, W = Matrix(W), C, V = Matrix(V), P0 = Matrix(P0), μ0,
            Lw, Lv, L0, M, Ω = inv(Matrix(W)), Ω0 = inv(Matrix(P0)))
end

"Truth + measurements y_0..y_T; contiguous bursts of scaled sensor noise."
function simulate(sys, inst::SmootherInstance, T::Int; corrupt_scale = inst.corrupt_scale)
    rng = MersenneTwister(inst.seed_sim)
    n, m = sys.n, sys.m
    x = zeros(T + 1, n); y = zeros(T + 1, m)
    x[1, :] .= sys.μ0 .+ sys.L0 * randn(rng, n)
    for t in 1:T
        x[t + 1, :] .= sys.A * x[t, :] .+ sys.Lw * randn(rng, n)
    end
    nburst = max(3, T ÷ 64)
    lo, hi = 5, T - inst.burst_len - 5
    starts = sort(randperm(rng, hi - lo + 1)[1:nburst] .+ (lo - 1))
    for i in 2:nburst
        starts[i] = max(starts[i], starts[i - 1] + inst.burst_len + 4)
    end
    corrupted = falses(T + 1)
    for s in starts
        s = min(s, T - inst.burst_len)
        corrupted[s:s + inst.burst_len - 1] .= true
    end
    for t in 0:T
        sc = corrupted[t + 1] ? corrupt_scale : 1.0
        y[t + 1, :] .= sys.C * x[t + 1, :] .+ sc .* (sys.Lv * randn(rng, m))
    end
    return x, y, corrupted
end

# -----------------------------------------------------------------------------
# Native references: RTS recursion and banded information form (no IPM)
# -----------------------------------------------------------------------------

"Kalman filter + Rauch–Tung–Striebel backward pass (closed form)."
function rts_smoother(sys, y)
    n = sys.n; T = size(y, 1) - 1
    xf = zeros(T + 1, n); Pf = [zeros(n, n) for _ in 1:T + 1]
    xp = zeros(T + 1, n); Pp = [zeros(n, n) for _ in 1:T + 1]
    xpred, Ppred = copy(sys.μ0), copy(sys.P0)
    for t in 1:T + 1
        xp[t, :] .= xpred; Pp[t] .= Ppred
        S = sys.C * Ppred * sys.C' .+ sys.V
        K = (S \ (sys.C * Ppred))'
        xf[t, :] .= xpred .+ K * (y[t, :] .- sys.C * xpred)
        Pf[t] .= (I - K * sys.C) * Ppred; Pf[t] .= 0.5 .* (Pf[t] .+ Pf[t]')
        xpred = sys.A * xf[t, :]; Ppred = sys.A * Pf[t] * sys.A' .+ sys.W
    end
    xs = copy(xf)
    for t in T:-1:1
        G = (Pp[t + 1] \ (sys.A * Pf[t]))'
        xs[t, :] .= xf[t, :] .+ G * (xs[t + 1, :] .- xp[t + 1, :])
    end
    return xs
end

"Block-tridiagonal information-form MAP J x = h (independent algorithm)."
function info_smoother(sys, y)
    n = sys.n; T = size(y, 1) - 1; N = (T + 1) * n
    Ωv = inv(sys.V)
    J = zeros(N, N); h = zeros(N)
    sl(t) = (t * n + 1):((t + 1) * n)       # t = 0..T
    J[sl(0), sl(0)] .+= sys.Ω0; h[sl(0)] .+= sys.Ω0 * sys.μ0
    AtΩ = sys.A' * sys.Ω
    for t in 0:T - 1
        J[sl(t), sl(t)] .+= AtΩ * sys.A
        J[sl(t + 1), sl(t + 1)] .+= sys.Ω
        J[sl(t), sl(t + 1)] .+= -AtΩ
        J[sl(t + 1), sl(t)] .+= -AtΩ'
    end
    CtΩv = sys.C' * Ωv
    for t in 0:T
        J[sl(t), sl(t)] .+= CtΩv * sys.C
        h[sl(t)] .+= CtΩv * y[t + 1, :]
    end
    return reshape(J \ h, n, T + 1)'
end

# -----------------------------------------------------------------------------
# Build (pair stalks + agreement edges + SOC ties)
# -----------------------------------------------------------------------------
#
# Stalk order: pair stalks 1..T (dim 2n, z_t = (x_{t-1}, x_t)) | SOC stalks
# vsoc(τ) = T + 1 + τ for τ = 0..T (dim m+2, (u, w, s)). x_τ is hosted by
# pair vertex τ+1 (first half) for τ < T and by pair vertex T (second half)
# for τ = T.

function build_smoother(sys, inst::SmootherInstance, y;
                        weight = inst.λ, δ = 0.0)
    n, m = sys.n, sys.m
    T = size(y, 1) - 1
    d = m + 2
    vpair(t) = t                            # t = 1..T
    vsoc(τ) = T + 1 + τ                     # τ = 0..T
    hostcols(τ) = τ < T ? (1:n) : (n + 1:2n)
    hostvtx(τ) = τ < T ? τ + 1 : T

    row_ids, col_ids, blocks = Int[], Int[], Matrix{Float64}[]
    rhs_val = Dict{Int, Vector{Float64}}()
    In = Matrix{Float64}(I, n, n)
    e = 0

    # ---- agreement: second half of z_t == first half of z_{t+1}
    for t in 1:T - 1
        e += 1
        A1 = zeros(n, 2n); A1[:, n + 1:2n] .= In
        A2 = zeros(n, 2n); A2[:, 1:n] .= -In
        push!(row_ids, e); push!(col_ids, vpair(t)); push!(blocks, A1)
        push!(row_ids, e); push!(col_ids, vpair(t + 1)); push!(blocks, A2)
    end

    # ---- ties: w_τ − M x_τ = −b_τ (dense by nature, nonzero rhs)
    for τ in 0:T
        e += 1
        Au = zeros(m, d); Au[:, 2:m + 1] .= Matrix{Float64}(I, m, m)
        Ax = zeros(m, 2n); Ax[:, hostcols(τ)] .= -sys.M
        push!(row_ids, e); push!(col_ids, vsoc(τ)); push!(blocks, Au)
        push!(row_ids, e); push!(col_ids, hostvtx(τ)); push!(blocks, Ax)
        rhs_val[e] = -(sys.Lv \ y[τ + 1, :])
    end

    # ---- pins: s_τ = δ (δ = 0: flagship sqrt loss; δ > 0: pseudo-Huber)
    for τ in 0:T
        e += 1
        Ap = zeros(1, d); Ap[1, d] = 1.0
        push!(row_ids, e); push!(col_ids, vsoc(τ)); push!(blocks, Ap)
        rhs_val[e] = [δ]
    end

    B = blocksparse(row_ids, col_ids, blocks)

    g_rhs = zeros(size(B, 1))
    for (edge, val) in rhs_val
        g_rhs[rowrange(B, edge)] .= val
    end

    # ---- dense pair Grams (+ prior on z_1) — the dense-Q habitat
    AtΩ = sys.A' * sys.Ω
    G = [(AtΩ * sys.A) (-AtΩ); (-AtΩ') sys.Ω]
    Q = IPM.allocblockdiag(B); fill!(Q, 0)
    for t in 1:T
        Qt = copy(G)
        t == 1 && (Qt[1:n, 1:n] .+= sys.Ω0)
        for j in 1:2n
            Qt[j, j] += RIDGE
        end
        block(Q, t, t, t) .= Qt
    end
    c = zeros(size(B, 2))
    c[colrange(B, 1)[1:n]] .= (sys.Ω0 * sys.μ0)
    for τ in 0:T
        c[first(colrange(B, vsoc(τ)))] = -weight
    end

    K_cones = vcat(IPM.AbstractCone[CofreeCone() for _ in 1:T],
                   IPM.AbstractCone[SecondOrderCone() for _ in 1:T + 1])
    prob = IPMProblem(Q, B, c, g_rhs, K_cones)

    ctx = (; n, m, T, d, weight, δ, y, sys)
    return prob, ctx
end

"State trajectory (T+1) × n from the solution vector."
function extract_x(prob, ctx, p::AbstractVector)
    X = zeros(ctx.T + 1, ctx.n)
    for τ in 0:ctx.T - 1
        X[τ + 1, :] .= p[colrange(prob.B, τ + 1)][1:ctx.n]
    end
    X[ctx.T + 1, :] .= p[colrange(prob.B, ctx.T)][ctx.n + 1:2 * ctx.n]
    return X
end

"Whitened measurement residuals r_τ = M x_τ − b_τ, one per timestep."
residuals(ctx, X) =
    [ctx.sys.M * X[τ + 1, :] .- (ctx.sys.Lv \ ctx.y[τ + 1, :]) for τ in 0:ctx.T]

"Direct evaluation of the flagship objective (assembly check)."
function smoother_objective(ctx, X)
    sys = ctx.sys
    F = 0.5 * sum(abs2, sys.L0 \ (X[1, :] .- sys.μ0))
    for t in 1:ctx.T
        F += 0.5 * sum(abs2, sys.Lw \ (X[t + 1, :] .- sys.A * X[t, :]))
    end
    if ctx.δ == 0.0
        F += ctx.weight * sum(norm, residuals(ctx, X))
    else
        F += sum(ctx.weight * norm(vcat(r, ctx.δ)) for r in residuals(ctx, X))
    end
    return F
end

# -----------------------------------------------------------------------------
# Baseline (same conic program via JuMP)
# -----------------------------------------------------------------------------

function build_jump_smoother(prob, ctx, optimizer)
    model = Model(optimizer)
    set_silent(model)
    T, n, d = ctx.T, ctx.n, ctx.d
    xs = Vector{Vector{VariableRef}}(undef, T + (T + 1))
    for t in 1:T
        xs[t] = @variable(model, [1:2n])
    end
    for j in 1:T + 1                        # SOC stalks: same convention
        v = @variable(model, [1:d])
        @constraint(model, v in JuMP.SecondOrderCone())
        xs[T + j] = v
    end
    x = reduce(vcat, xs)
    Qs = sparse(prob.Q); Bs = sparse(prob.B)
    @objective(model, Min, 0.5 * x' * Qs * x - prob.c' * x)
    @constraint(model, Bs * x .== prob.g)
    return model, x
end

# -----------------------------------------------------------------------------
# Gate tests
# -----------------------------------------------------------------------------

density(Mx; tol = 1e-8) = count(abs.(Mx) .> tol * maximum(abs.(Mx))) / length(Mx)

function test_dense_by_nature(sys)
    ds = (density(sys.A), density(sys.W), density(sys.V), density(sys.M))
    @assert minimum(ds) > 0.95 "dense-by-nature broken: $ds"
    println("  [PASS] dense by nature: density@1e-8 A $(round(ds[1], digits = 3)) ",
        "W $(round(ds[2], digits = 3)) V $(round(ds[3], digits = 3)) ",
        "M $(round(ds[4], digits = 3))")
end

function test_references(sys, inst)
    _, y, _ = simulate(sys, inst, 100; corrupt_scale = 1.0)
    e = maximum(abs.(rts_smoother(sys, y) .- info_smoother(sys, y)))
    @assert e < 1e-9 "RTS vs information form: $e (oracle: 2.0e-14)"
    println("  [PASS] analytic references: RTS recursion vs information form ",
        "$(round(e, sigdigits = 2))")
end

function test_objective_identity(prob, ctx, res)
    X = extract_x(prob, ctx, res.p)
    o_direct = smoother_objective(ctx, X)
    o_ipm = ipm_objective(prob, res)
    rel = abs(o_direct - o_ipm) / abs(o_direct)
    @assert rel < 1e-5 "assembly mismatch: $o_direct vs $o_ipm"
    println("  [PASS] objective identity: F(x) = $(round(o_direct, digits = 4)) ",
        "(rel $(round(rel, sigdigits = 2)))")
    return X
end

function test_deformation(sys, inst, settings)
    _, y, _ = simulate(sys, inst, 100; corrupt_scale = 1.0)
    xs = rts_smoother(sys, y)
    errs = Float64[]
    for δ in (4.0, 16.0)
        prob, ctx = build_smoother(sys, inst, y; weight = δ, δ = δ)
        res = solve(prob, settings)
        @assert res.status in (OPTIMAL, NEAR_OPTIMAL) "δ=$δ: $(res.status)"
        push!(errs, maximum(abs.(extract_x(prob, ctx, res.p) .- xs)))
    end
    ratio = errs[2] / errs[1]
    rate = log(ratio) / log(4.0)
    @assert errs[2] < 2e-2 "pseudo-Huber δ=16 too far from RTS: $(errs[2])"
    @assert ratio < 0.15 "deformation rate broken: δ^$(round(rate, digits = 2)) (oracle: δ^-1.98)"
    println("  [PASS] deformation to RTS: ‖x(δ) − RTS‖∞ $(round(errs[1], sigdigits = 2)) ",
        "→ $(round(errs[2], sigdigits = 2)) over δ = 4 → 16 (rate δ^$(round(rate, digits = 2)))")
end

function test_robustness(sys, inst, settings)
    vals = Dict{Float64, NTuple{2, Float64}}()
    for sc in (1.0, 32.0)
        x, y, corr = simulate(sys, inst, 100; corrupt_scale = sc)
        clean = .!corr
        rmse(X) = sqrt(sum(abs2, (X .- x)[clean, :]) / (count(clean) * sys.n))
        prob, ctx = build_smoother(sys, inst, y)
        res = solve(prob, settings)
        @assert res.status in (OPTIMAL, NEAR_OPTIMAL) "×$sc: $(res.status)"
        vals[sc] = (rmse(extract_x(prob, ctx, res.p)), rmse(rts_smoother(sys, y)))
    end
    (r1, k1), (r32, k32) = vals[1.0], vals[32.0]
    @assert r32 < 0.6 * k32 "robustness gate: robust $r32 vs RTS $k32"
    @assert r32 < 1.6 * r1 "sqrt loss did not saturate: $r32 vs $r1"
    @assert r1 < 1.3 * k1 "clean-data price too high: $r1 vs $k1"
    println("  [PASS] robustness: clean-region RMSE ×1 robust $(round(r1, digits = 4)) ",
        "RTS $(round(k1, digits = 4)); ×32 robust $(round(r32, digits = 4)) ",
        "RTS $(round(k32, digits = 4))")
end

function test_support_recovery(ctx, X, corrupted)
    rn = norm.(residuals(ctx, X))
    k = count(corrupted)
    top = sortperm(rn; rev = true)[1:k]
    tp = count(τ -> corrupted[τ], top)
    margin = minimum(rn[corrupted]) / maximum(rn[.!corrupted])
    @assert tp == k "support recovery: $tp/$k (oracle: 24/24)"
    @assert margin > 1.5 "support margin: $margin (oracle: 2.8)"
    println("  [PASS] support recovery: top-$k residuals are exactly the corrupted ",
        "steps (margin ×$(round(margin, digits = 1)); influence cap ",
        "$(round(maximum(rn) / ctx.weight, digits = 0))× vs quadratic)")
end

function test_ipm_vs_clarabel(prob, ctx, X_ipm)
    model, xv = build_jump_smoother(prob, ctx, clarabel_opt(; tol = TOL))
    optimize!(model)
    X_cla = extract_x(prob, ctx, value.(xv))
    dx = maximum(abs.(X_ipm .- X_cla))
    @assert dx < 1e-3 "IPM vs Clarabel mismatch: $dx"
    println("  [PASS] IPM vs Clarabel (same conic program): ‖Δx‖∞ = ",
        "$(round(dx, sigdigits = 2))")
end

# -----------------------------------------------------------------------------
# Sweep
# -----------------------------------------------------------------------------

sizes() = OPTS.tiny ? [64] :
    OPTS.quick ? [64, 128, 256] :
    [64, 128, 256, 512, 1024]

function run()
    println("\n", "=" ^ 78)
    println("  E09: robust Kalman smoothing, burst corruption ",
        "(dial: horizon T; n = 12, m = 4)")
    println("=" ^ 78)

    ipm_settings = IPMSettings{Float64}(
        kkt = UzawaSettings{Float64}(raug = 1e4),
        feas_tol = TOL, gap_tol = TOL, itmax = 200)
    hsd_settings = HSDSettings{Float64}(
        kkt = UzawaSettings{Float64}(raug = 1e4),
        feas_tol = TOL, gap_tol = TOL, itmax = 200)

    inst = smoother_instance()
    sys = make_system(inst)

    println("\n  Gate tests (n = 12, m = 4, T = 100):")
    test_dense_by_nature(sys)
    test_references(sys, inst)
    x_true, y, corrupted = simulate(sys, inst, 100)
    prob, ctx = build_smoother(sys, inst, y)
    res = solve(prob, ipm_settings)
    @assert res.status in (OPTIMAL, NEAR_OPTIMAL) "IPM status $(res.status)"
    X = test_objective_identity(prob, ctx, res)
    test_deformation(sys, inst, ipm_settings)
    test_robustness(sys, inst, ipm_settings)
    test_support_recovery(ctx, X, corrupted)
    test_ipm_vs_clarabel(prob, ctx, X)
    println()

    cla_opt = clarabel_opt(; tol = TOL)
    msk_opt = OPTS.mosek ? mosek_opt(; tol = TOL) : nothing

    rows = []
    for T in sizes()
        _, yT, _ = simulate(sys, inst, T)
        prob, ctx = build_smoother(sys, inst, yT)
        stats = problem_stats(prob)

        @printf("  T=%-5d dof=%-6d n1=%-5d blk=%-4.0f  ",
            T, stats.N0, stats.N1, stats.med_block)

        m_ipm = measure_ipm(prob, ipm_settings; nruns = NRUNS)
        m_hsd = measure_ipm(prob, hsd_settings; nruns = NRUNS)
        m_cla = measure_jump(() -> first(build_jump_smoother(prob, ctx, cla_opt)); nruns = NRUNS)
        m_msk = msk_opt !== nothing ?
            measure_jump(() -> first(build_jump_smoother(prob, ctx, msk_opt)); nruns = NRUNS) :
            (t = NaN, status = "", obj = NaN)

        ratio(b, base) = isfinite(b.t) && isfinite(base.t) ? b.t / base.t : NaN
        fmt_ratio(b, base) = isnan(ratio(b, base)) ? "—" : @sprintf("%.2fx", ratio(b, base))

        println("IPM ", fmt_time(m_ipm), "  HSD ", fmt_time(m_hsd), " (", fmt_ratio(m_hsd, m_ipm),
            ")  Cla ", fmt_time(m_cla), " (", fmt_ratio(m_cla, m_ipm),
            ")  Msk ", fmt_time(m_msk), " (", fmt_ratio(m_msk, m_ipm), ")")

        push!(rows, (size = T, dof = stats.N0, ipm = m_ipm, hsd = m_hsd, cla = m_cla, msk = m_msk))
    end

    dofs = [r.dof for r in rows]
    println("\nFitted log-log slopes (time vs DOF):")
    for (name, get_t) in [("IPM", r -> r.ipm.t), ("HSD", r -> r.hsd.t), ("Clarabel", r -> r.cla.t), ("Mosek", r -> r.msk.t)]
        sl = loglog_slope(dofs, [get_t(r) for r in rows])
        print("  $name: ", isnan(sl) ? "n/a" : @sprintf("DOF^%.2f", sl))
    end
    println()
end

