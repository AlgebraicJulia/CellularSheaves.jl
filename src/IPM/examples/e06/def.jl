# =============================================================================
# e06/def.jl — Square-root-loss spectral regression (whitened SOC chain)
#
# Usage:
#   julia --project e06/run.jl              # Clarabel only
#   julia --project e06/run.jl --mosek      # + Mosek
#   julia --project e06/run.jl --quick      # small sweep
#
# Problem. Fit a smooth function with P Chebyshev patches (degree n, C^k jets)
# under BLOCK-CORRUPTED noise — a few patches heavily contaminated:
#
#     min_θ  Σ_p ||E_p θ_p − y_p||₂      (norms, NOT squares)
#     s.t.   C^k jets at patch interfaces.
#
# The square-root loss bounds each patch's influence: the gradient of ||r_p||
# is a unit vector, so a corrupted patch's pull on the global fit saturates.
# Oracle escalation curve (clean-region RMSE, corruption ×1 → ×128):
#     sqrt-loss 0.0217 → 0.0233 (saturates);  LS 0.0216 → 0.1893 (blows up);
# and the price of robustness on clean data is ~0.5%.
#
# Whitened formulation (the recipe for keeping epigraph SOCPs in the dense-
# block habitat). Naively, t_p ≥ ||E_p θ_p − y_p|| needs a residual variable
# r_p tied to θ by identity blocks — internally-sparse maps, E5's disease.
# Instead, with E'E = LL', b = L⁻¹E'y, c = y'y − ||b||² ≥ 0:
#     ||E θ − y||² = ||L'θ − b||² + c,
# so t ≥ ||Eθ − y||  ⟺  (t, L'θ − b, √c) ∈ SOC(n+3). Eliminating
# θ = L⁻ᵀ(w + b) leaves ONE SOC stalk (t, w, s) per patch: jets become dense
# maps through L⁻ᵀ with nonzero right-hand side, s is pinned to √c by a
# one-row block, and there are no free stalks and no identity blocks. Every
# map is dense by nature; the bloat factor is ~1.
#
# Sheaf structure: a path of P SecondOrderCone stalks (dim n+3), (k+1)-row
# dense jet edges (row-normalized, g ≠ 0), one pin row per stalk. Q = tiny
# ridge on the w slots. Dial: chain length P — vertex count grows, per-block
# work constant.
#
# Validation (tools/sqrtloss_oracle.py):
#   * whitening identity exact to 5.9e-16 (random-θ algebra check)
#   * whitened conic == direct norm-sum solve: ||Δf||∞ = 4.9e-5,
#     rel objective diff 6.5e-10
#   * robustness escalation curve as above; saturation and LS blowup gated
#   * stationarity certificate (gradient + jet-multiplier span): 1.2e-5
#   * scaling P = 8..256 all solved at tol 1e-8
# The gate tests below re-derive the fast subset in-process.
# =============================================================================

using AppleAccelerate
using LinearAlgebra, SparseArrays, Printf, Random

include("../utils.jl")

const OPTS = parse_args(ARGS)
const TOL = OPTS.tol
const NRUNS = 5

# -----------------------------------------------------------------------------
# Problem definition
# -----------------------------------------------------------------------------

struct SqrtLossInstance
    P::Int; n::Int; k::Int; npts::Int
    noise::Float64; corrupt_frac::Float64; corrupt_scale::Float64
    seed::Int
end

sqrtloss_instance(; P = 64, n = 24, k = 2, npts = 200, noise = 0.05,
                  corrupt_frac = 1 / 6, corrupt_scale = 8.0, seed = 1) =
    SqrtLossInstance(P, n, k, npts, noise, corrupt_frac, corrupt_scale, seed)

f_true(x) = sin(2π * (3x + 4x^2)) * (1 + 0.5x) + 1.5x

function make_data(inst::SqrtLossInstance)
    rng = MersenneTwister(inst.seed)
    knots = collect(range(0.0, 1.0; length = inst.P + 1))
    ncor = max(1, round(Int, inst.corrupt_frac * inst.P))
    corrupted = sort(randperm(rng, inst.P)[1:ncor])
    xs = Vector{Vector{Float64}}(undef, inst.P)
    ys = Vector{Vector{Float64}}(undef, inst.P)
    for p in 1:inst.P
        x = sort(rand(rng, inst.npts)) .* (knots[p + 1] - knots[p]) .+ knots[p]
        s = inst.noise * (p in corrupted ? inst.corrupt_scale : 1.0)
        xs[p] = x
        ys[p] = f_true.(x) .+ s .* randn(rng, inst.npts)
    end
    return knots, xs, ys, corrupted
end

# ---- Chebyshev machinery -----------------------------------------------------

function chebval(u::Real, c::AbstractVector{Float64})
    b1 = 0.0; b2 = 0.0
    for j in length(c):-1:2
        b1, b2 = 2u * b1 - b2 + c[j], b1
    end
    return u * b1 - b2 + c[1]
end

function chebder(c::AbstractVector{Float64})
    nn = length(c) - 1
    nn == 0 && return [0.0]
    d = zeros(nn + 2)
    for kk in nn:-1:1
        d[kk] = d[kk + 2] + 2kk * c[kk + 1]
    end
    d[1] /= 2
    return d[1:nn]
end

function chebder(c::AbstractVector{Float64}, r::Int)
    for _ in 1:r; c = chebder(c); end
    return c
end

function chebvander(u::AbstractVector, n::Int)
    V = zeros(length(u), n + 1)
    V[:, 1] .= 1.0
    n ≥ 1 && (V[:, 2] .= u)
    for j in 3:n + 1
        V[:, j] .= 2.0 .* u .* V[:, j - 1] .- V[:, j - 2]
    end
    V
end

function jet_rows(n::Int, k::Int, h::Float64)
    Rp = zeros(k + 1, n + 1); Rm = zeros(k + 1, n + 1)
    for r in 0:k, j in 1:n + 1
        e = zeros(n + 1); e[j] = 1.0
        d = r == 0 ? e : chebder(e, r)
        Rp[r + 1, j] = chebval(+1.0, d) * (2.0 / h)^r
        Rm[r + 1, j] = chebval(-1.0, d) * (2.0 / h)^r
    end
    return Rp, Rm
end

eval_mat(x, knots, p, n) =
    chebvander(2.0 .* (x .- knots[p]) ./ (knots[p + 1] - knots[p]) .- 1.0, n)

# -----------------------------------------------------------------------------
# Build (whitened)
# -----------------------------------------------------------------------------

const RIDGE = 1e-8

function build_sqrtloss(inst::SqrtLossInstance)
    knots, xs, ys, corrupted = make_data(inst)
    P, n, k = inst.P, inst.n, inst.k
    d = n + 3
    Ls = Vector{Matrix{Float64}}(undef, P)
    bs = Vector{Vector{Float64}}(undef, P)
    cs = zeros(P)
    Es = Vector{Matrix{Float64}}(undef, P)
    for p in 1:P
        E = eval_mat(xs[p], knots, p, n)
        L = cholesky(Symmetric(E' * E)).L
        b = L \ (E' * ys[p])
        Ls[p] = Matrix(L); bs[p] = b
        cs[p] = max(ys[p]' * ys[p] - b' * b, 0.0)
        Es[p] = E
    end
    LinvT = [Matrix(inv(L)') for L in Ls]

    row_ids, col_ids, blocks = Int[], Int[], Matrix{Float64}[]
    g = Float64[]
    e = 0
    Rp, Rm = jet_rows(n, k, 1.0 / P)
    for p in 1:P - 1                        # jets: dense through L^{-T}, g ≠ 0
        A1 = Rp * LinvT[p]; A2 = Rm * LinvT[p + 1]
        rhs = -A1 * bs[p] .+ A2 * bs[p + 1]
        B1 = zeros(k + 1, d); B2 = zeros(k + 1, d)
        for r in 1:k + 1
            s = norm(vcat(A1[r, :], A2[r, :]))
            B1[r, 2:n + 2] .= A1[r, :] ./ s
            B2[r, 2:n + 2] .= -A2[r, :] ./ s
            push!(g, rhs[r] / s)
        end
        e += 1
        push!(row_ids, e); push!(col_ids, p); push!(blocks, B1)
        push!(row_ids, e); push!(col_ids, p + 1); push!(blocks, B2)
    end
    for p in 1:P                            # s_p = sqrt(c_p)
        Bp = zeros(1, d); Bp[1, d] = 1.0
        e += 1
        push!(row_ids, e); push!(col_ids, p); push!(blocks, Bp)
        push!(g, sqrt(cs[p]))
    end
    B = blocksparse(row_ids, col_ids, blocks)

    Q = IPM.allocblockdiag(B); fill!(Q, 0)
    c = zeros(size(B, 2))
    for p in 1:P
        Qp = zeros(d, d)
        for j in 2:n + 2
            Qp[j, j] = RIDGE
        end
        block(Q, p, p, p) .= Qp
        c[first(colrange(B, p))] = -1.0     # minimize Σ t_p
    end
    K = [IPM.SecondOrderCone() for _ in 1:P]
    prob = IPMProblem(Q, B, c, Vector{Float64}(g), K)

    ctx = (; P, n, k, d, knots, xs, ys, corrupted, Ls, bs, cs, Es)
    return prob, ctx
end

"Recover per-patch Chebyshev coefficients: θ_p = L_p^{-T}(w_p + b_p)."
function extract_thetas(prob, ctx, p::AbstractVector)
    [ctx.Ls[q]' \ (p[colrange(prob.B, q)][2:ctx.n + 2] .+ ctx.bs[q]) for q in 1:ctx.P]
end

function eval_fit(ctx, thetas, xq::AbstractVector)
    yq = similar(xq)
    for (t, xx) in enumerate(xq)
        p = clamp(1 + floor(Int, xx * ctx.P), 1, ctx.P)
        u = 2.0 * (xx - ctx.knots[p]) / (ctx.knots[p + 1] - ctx.knots[p]) - 1.0
        yq[t] = chebval(u, thetas[p])
    end
    yq
end

"Direct objective Σ ||E_p θ_p − y_p|| (assembly check)."
sqrt_objective(ctx, thetas) =
    sum(norm(ctx.Es[p] * thetas[p] .- ctx.ys[p]) for p in 1:ctx.P)

"Equality-constrained least squares (dense KKT): the non-robust baseline."
function solve_ls(ctx)
    P, n, k = ctx.P, ctx.n, ctx.k
    nc = n + 1; N = P * nc
    G = zeros(N, N); b = zeros(N)
    for p in 1:P
        sl = (p - 1) * nc + 1:p * nc
        G[sl, sl] .= 2.0 .* (ctx.Es[p]' * ctx.Es[p])
        b[sl] .= 2.0 .* (ctx.Es[p]' * ctx.ys[p])
    end
    Rp, Rm = jet_rows(n, k, 1.0 / P)
    nA = (P - 1) * (k + 1)
    A = zeros(nA, N); rr = 0
    for p in 1:P - 1, r in 1:k + 1
        rr += 1
        A[rr, (p - 1) * nc + 1:p * nc] .= Rp[r, :]
        A[rr, p * nc + 1:(p + 1) * nc] .= -Rm[r, :]
    end
    KKT = [G A'; A zeros(nA, nA)]
    sol = KKT \ vcat(b, zeros(nA))
    return [sol[(p - 1) * nc + 1:p * nc] for p in 1:P]
end

# -----------------------------------------------------------------------------
# Baselines (same whitened conic program)
# -----------------------------------------------------------------------------

function build_jump_sqrtloss(prob, ctx, optimizer)
    model = Model(optimizer)
    set_silent(model)
    xs = Vector{Vector{VariableRef}}(undef, ctx.P)
    for p in 1:ctx.P
        v = @variable(model, [1:ctx.d])
        @constraint(model, v in JuMP.SecondOrderCone())
        xs[p] = v
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

function test_whitening(ctx)
    rng = MersenneTwister(7)
    worst = 0.0
    for _ in 1:20
        p = rand(rng, 1:ctx.P)
        θ = randn(rng, ctx.n + 1)
        lhs = sum(abs2, ctx.Es[p] * θ .- ctx.ys[p])
        rhs = sum(abs2, ctx.Ls[p]' * θ .- ctx.bs[p]) + ctx.cs[p]
        worst = max(worst, abs(lhs - rhs) / max(lhs, 1.0))
    end
    @assert worst < 1e-10 "whitening identity broken: $worst"
    println("  [PASS] whitening identity: worst rel err $(round(worst, sigdigits = 2))")
end

function test_objective_identity(prob, ctx, res)
    thetas = extract_thetas(prob, ctx, res.p)
    o_direct = sqrt_objective(ctx, thetas)
    o_ipm = ipm_objective(prob, res)
    rel = abs(o_direct - o_ipm) / abs(o_direct)
    @assert rel < 1e-4 "assembly mismatch: $o_direct vs $o_ipm"
    println("  [PASS] objective identity: Σ‖r_p‖ = $(round(o_direct, digits = 5)) ",
        "(rel $(round(rel, sigdigits = 2)))")
    return thetas
end

function test_robustness(settings)
    xg = collect(range(0.0, 1.0; length = 4001))
    rms_clean(ctx, f) = begin
        mask = trues(length(xg))
        for p in ctx.corrupted
            mask .&= .!((xg .≥ ctx.knots[p]) .& (xg .≤ ctx.knots[p + 1]))
        end
        sqrt(sum(abs2, (f .- f_true.(xg))[mask]) / count(mask))
    end
    vals = Dict{Float64, Tuple{Float64, Float64}}()
    for cs_ in (1.0, 32.0)
        inst = sqrtloss_instance(; P = 12, corrupt_scale = cs_)
        prob, ctx = build_sqrtloss(inst)
        res = solve(prob, settings)
        @assert res.status in (OPTIMAL, NEAR_OPTIMAL)
        f_sq = eval_fit(ctx, extract_thetas(prob, ctx, res.p), xg)
        f_ls = eval_fit(ctx, solve_ls(ctx), xg)
        vals[cs_] = (rms_clean(ctx, f_sq), rms_clean(ctx, f_ls))
    end
    (s1, l1), (s32, l32) = vals[1.0], vals[32.0]
    @assert s32 < 0.6 * l32 "robustness gate: sqrt $s32 vs LS $l32"
    @assert s32 < 1.3 * s1 "sqrt-loss did not saturate: $s32 vs $s1"
    @assert s1 < 1.3 * l1 "clean-data price too high: $s1 vs $l1"
    println("  [PASS] robustness: clean-region RMSE ×1 sqrt $(round(s1, digits = 4)) ",
        "LS $(round(l1, digits = 4)); ×32 sqrt $(round(s32, digits = 4)) ",
        "LS $(round(l32, digits = 4))")
end

function test_ipm_vs_clarabel(prob, ctx, thetas)
    model, xv = build_jump_sqrtloss(prob, ctx, clarabel_opt(; tol = TOL))
    optimize!(model)
    th_c = extract_thetas(prob, ctx, value.(xv))
    xg = collect(range(0.0, 1.0; length = 4001))
    df = maximum(abs.(eval_fit(ctx, thetas, xg) .- eval_fit(ctx, th_c, xg)))
    o_i = sqrt_objective(ctx, thetas); o_c = sqrt_objective(ctx, th_c)
    @assert df < 1e-3 "IPM vs Clarabel: $df"
    @assert abs(o_i - o_c) < 1e-4 * abs(o_c)
    println("  [PASS] IPM vs Clarabel: ‖Δf‖∞ = $(round(df, sigdigits = 2)); ",
        "obj $(round(o_i, digits = 5)) vs $(round(o_c, digits = 5))")
end

# -----------------------------------------------------------------------------
# Sweep
# -----------------------------------------------------------------------------

sizes() = OPTS.tiny ? [32] :
    OPTS.quick ? [32, 64, 128] :
    [32, 64, 128, 256, 512, 1024]

function run()
    println("\n", "=" ^ 78)
    println("  E06: sqrt-loss spectral regression (dial: chain length P; n = 24)")
    println("=" ^ 78)

    ipm_settings = IPMSettings{Float64}(
        kkt = UzawaSettings{Float64}(raug = 1e4),
        feas_tol = TOL, gap_tol = TOL, itmax = 200)
    hsd_settings = HSDSettings{Float64}(
        kkt = UzawaSettings{Float64}(raug = 1e4),
        feas_tol = TOL, gap_tol = TOL, itmax = 200)

    println("\n  Gate tests:")
    inst = sqrtloss_instance(; P = 12)
    prob, ctx = build_sqrtloss(inst)
    test_whitening(ctx)
    res = solve(prob, ipm_settings)
    @assert res.status in (OPTIMAL, NEAR_OPTIMAL) "IPM status $(res.status)"
    thetas = test_objective_identity(prob, ctx, res)
    test_robustness(ipm_settings)
    test_ipm_vs_clarabel(prob, ctx, thetas)
    println()

    cla_opt = clarabel_opt(; tol = TOL)
    msk_opt = OPTS.mosek ? mosek_opt(; tol = TOL) : nothing

    rows = []
    for P in sizes()
        inst = sqrtloss_instance(; P = P)
        prob, ctx = build_sqrtloss(inst)
        stats = problem_stats(prob)

        @printf("  P=%-5d dof=%-6d n1=%-5d blk=%-4.0f  ",
            P, stats.N0, stats.N1, stats.med_block)

        m_ipm = measure_ipm(prob, ipm_settings; nruns = NRUNS)
        m_hsd = measure_ipm(prob, hsd_settings; nruns = NRUNS)
        m_cla = measure_jump(() -> first(build_jump_sqrtloss(prob, ctx, cla_opt)); nruns = NRUNS)
        m_msk = msk_opt !== nothing ?
            measure_jump(() -> first(build_jump_sqrtloss(prob, ctx, msk_opt)); nruns = NRUNS) :
            (t = NaN, status = "", obj = NaN)

        ratio(b, base) = isfinite(b.t) && isfinite(base.t) ? b.t / base.t : NaN
        fmt_ratio(b, base) = isnan(ratio(b, base)) ? "—" : @sprintf("%.2fx", ratio(b, base))

        println("IPM ", fmt_time(m_ipm), "  HSD ", fmt_time(m_hsd), " (", fmt_ratio(m_hsd, m_ipm),
            ")  Cla ", fmt_time(m_cla), " (", fmt_ratio(m_cla, m_ipm),
            ")  Msk ", fmt_time(m_msk), " (", fmt_ratio(m_msk, m_ipm), ")")

        push!(rows, (size = P, dof = stats.N0, ipm = m_ipm, hsd = m_hsd, cla = m_cla, msk = m_msk))
    end

    dofs = [r.dof for r in rows]
    println("\nFitted log-log slopes (time vs DOF):")
    for (name, get_t) in [("IPM", r -> r.ipm.t), ("HSD", r -> r.hsd.t), ("Clarabel", r -> r.cla.t), ("Mosek", r -> r.msk.t)]
        sl = loglog_slope(dofs, [get_t(r) for r in rows])
        print("  $name: ", isnan(sl) ? "n/a" : @sprintf("DOF^%.2f", sl))
    end
    println()
end

