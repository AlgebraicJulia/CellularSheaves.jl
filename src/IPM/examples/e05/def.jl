# =============================================================================
# e05/def.jl — Galerkin-compressed nonlocal TV (dense-restriction-map SOC habitat)
#
# Usage:
#   julia --project e05/run.jl              # Clarabel only
#   julia --project e05/run.jl --mosek      # + Mosek
#   julia --project e05/run.jl --quick      # small sweep
#
# The same variational problem as E5 (Gilboa–Osher isotropic nonlocal TV),
# re-discretized: u is a piecewise polynomial — P non-overlapping tiles, m
# local Chebyshev coefficients each, C^k jets at tile interfaces — instead of
# a nodal vector. E5 sits at the boundary of this solver's habitat (its
# restriction maps are internally sparse; block-dense storage inflates B by
# ~85x and scalar-sparse solvers rightly win). Here every map is dense BY
# NATURE: point evaluation of a spectral representation touches all m
# coefficients, so z-rows, jets, and the per-tile data Gram are honestly
# dense, and the bloat factor is ~1 by construction, not by decree. The
# compression is bought with a genuinely smaller problem (u-DOF ÷4 at m=32),
# and every solver receives the same compressed conic program.
#
#     min_θ  1/2 Σ_t ||E_t θ_t − y_t||² + λ Σ_i ||G_i(θ)||₂
#     s.t.   C^k jets match at tile interfaces
#
# Sheaf structure. P tile stalks (CofreeCone, dim m) with dense data Grams
# E'E on the diagonal of Q; N node stalks (SecondOrderCone, dim K+1) holding
# (t_i, z_i) with λ on the t slot; edges: (k+1)×m dense jet blocks between
# adjacent tiles (row-normalized) and K-row z-blocks z_i = √w (eval_i −
# eval_j) θ touching one or two tile stalks (contributions to a shared tile
# are merged). No halos, no copies, no ownership weights — the coefficient
# representation removes all of E5's duplication machinery.
#
# Statistical bonus (measured, not assumed): the truncated basis filters
# fine-scale noise the TV term leaves behind — on the reference instance the
# compressed estimator BEATS the nodal one (MSE 0.0275 at m=32 vs 0.0344
# nodal vs 0.0501 best-λ local TV).
#
# Validation (tools/galerkin_nltv_oracle.py):
#   * graph gate: clean-signal matches 36% vs 5% chance (6.8x enrichment)
#   * nodal reference certified by subgradient-selection KKT check (7.9e-6)
#   * compression family: F monotone in m; resolution threshold visible at
#     m ≈ π·(Tsz/period) as predicted; MSE optimal at m = 32–40
#   * completeness certificate: full basis + no jets reproduces the nodal
#     optimum exactly (excess +1.3e-6 at Tsz = m = 16)
#   * sheaf assembly == independent monolithic solve: ||Δu||∞ = 5.0e-6
#   * conditioning bound: equispaced Chebyshev evaluation is well-conditioned
#     for m ≤ 48 (cond 5.4 at m=32, 1.9e2 at m=48) and blows up past m ≈ 64 —
#     the m dial must respect this
# The gate tests below re-derive the fast subset in-process, including a
# small-scale completeness certificate against a nodal JuMP solve.
# =============================================================================

using AppleAccelerate
using LinearAlgebra, SparseArrays, Printf, Random
using Statistics: median

include("../utils.jl")

const OPTS = parse_args(ARGS)
const TOL = OPTS.tol
const NRUNS = 5

# -----------------------------------------------------------------------------
# Problem definition
# -----------------------------------------------------------------------------

struct GalerkinNLTVInstance
    N::Int; Tsz::Int; m::Int; k::Int      # nodes, tile size, coefs/tile, C^k
    K::Int; R::Int; q::Int                # nonlocal neighbors, radius, patch
    λ::Float64; noise::Float64; seed::Int
end

galerkin_instance(; N = 1024, Tsz = 128, m = 32, k = 2, K = 12, R = 64,
                  q = 4, λ = 0.12, noise = 0.35, seed = 1) =
    GalerkinNLTVInstance(N, Tsz, m, k, K, R, q, λ, noise, seed)

"Period-16 texture under a smooth (raised-cosine) envelope: C^k-representable."
function make_signal(N::Int, noise::Float64, seed::Int)
    rng = MersenneTwister(seed)
    ramp(s) = 0.5 - 0.5 * cospi(clamp(s, 0.0, 1.0))
    u_true = zeros(N); x = 0:N - 1
    for (t, xx) in enumerate(x)
        env = 0.5 + 0.5 * ramp((128 - abs(xx % 256 - 128)) / 48.0)
        base = 0.8 * ramp((256 - abs(xx % 512 - 256)) / 64.0)
        u_true[t] = env * sin(2π * xx / 16) + base
    end
    return u_true, u_true .+ noise .* randn(rng, N)
end

# ---- nonlocal graph: identical to E5's (see e5.jl) --------------------------

function patch_matrix(y::Vector{Float64}, q::Int)
    N = length(y)
    yp = vcat(reverse(y[2:q + 1]), y, reverse(y[N - q:N - 1]))
    P = zeros(N, 2q + 1)
    for i in 1:N
        P[i, :] .= yp[i:i + 2q]
    end
    return P
end

function nonlocal_graph(y::Vector{Float64}, K::Int, q::Int, R::Int)
    N = length(y)
    P = patch_matrix(y, q)
    h = median(abs.(diff(y)))
    nbrs = zeros(Int, N, K); wts = zeros(N, K)
    for i in 1:N
        js = [j for j in max(1, i - R):min(N, i + R) if abs(j - i) > 2]
        d2 = [sum(abs2, view(P, j, :) .- view(P, i, :)) / (2q + 1) for j in js]
        order = partialsortperm(d2, 1:K)
        nbrs[i, :] .= js[order]
        w = exp.(-d2[order] ./ (2 * h^2))
        wts[i, :] .= w .* (K / sum(w))
    end
    return nbrs, wts
end

# ---- Chebyshev machinery (local basis per tile) ------------------------------

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

"Evaluation functional at global node g in tile (lo, hi), m coefficients."
function eval_row(g::Int, lo::Int, hi::Int, m::Int)
    u = 2.0 * (g - lo) / (hi - lo) - 1.0
    row = zeros(m)
    for j in 1:m
        e = zeros(m); e[j] = 1.0
        row[j] = chebval(u, e)
    end
    return row
end

"C^k jet functionals at local ±1, tile width h."
function jet_rows(m::Int, k::Int, h::Float64)
    Rp = zeros(k + 1, m); Rm = zeros(k + 1, m)
    for r in 0:k, j in 1:m
        e = zeros(m); e[j] = 1.0
        d = r == 0 ? e : chebder(e, r)
        Rp[r + 1, j] = chebval(+1.0, d) * (2.0 / h)^r
        Rm[r + 1, j] = chebval(-1.0, d) * (2.0 / h)^r
    end
    return Rp, Rm
end

# -----------------------------------------------------------------------------
# Build
# -----------------------------------------------------------------------------

function build_galerkin(inst::GalerkinNLTVInstance)
    u_true, y = make_signal(inst.N, inst.noise, inst.seed)
    nbrs, wts = nonlocal_graph(y, inst.K, inst.q, inst.R)
    N, Tsz, m, k, K, λ = inst.N, inst.Tsz, inst.m, inst.k, inst.K, inst.λ
    P = N ÷ Tsz
    @assert P * Tsz == N && P ≥ 2
    tiles = [((t - 1) * Tsz + 1, t * Tsz) for t in 1:P]        # inclusive
    tile_of(g) = min((g - 1) ÷ Tsz + 1, P)
    h = Float64(Tsz)

    Emats = [vcat([eval_row(g, tiles[t][1], tiles[t][2], m)' for g in
                   tiles[t][1]:tiles[t][2]]...) for t in 1:P]

    row_ids, col_ids, blocks = Int[], Int[], Matrix{Float64}[]
    e = 0
    if k ≥ 0
        Rp, Rm = jet_rows(m, k, h)
        for t in 1:P - 1
            R1 = copy(Rp); R2 = copy(Rm)
            for r in 1:k + 1
                s = norm(vcat(R1[r, :], R2[r, :]))
                R1[r, :] ./= s; R2[r, :] ./= s
            end
            e += 1
            push!(row_ids, e); push!(col_ids, t); push!(blocks, R1)
            push!(row_ids, e); push!(col_ids, t + 1); push!(blocks, -R2)
        end
    end
    for i in 1:N                          # z-rows; merge same-tile contributions
        ti = tile_of(i)
        ei = eval_row(i, tiles[ti][1], tiles[ti][2], m)
        acc = Dict{Int, Matrix{Float64}}()
        acc[ti] = zeros(K, m)
        An = zeros(K, K + 1)
        for kk in 1:K
            j = nbrs[i, kk]; tj = tile_of(j)
            ej = eval_row(j, tiles[tj][1], tiles[tj][2], m)
            sw = sqrt(wts[i, kk])
            acc[ti][kk, :] .+= sw .* ei
            haskey(acc, tj) || (acc[tj] = zeros(K, m))
            acc[tj][kk, :] .-= sw .* ej
            An[kk, 1 + kk] = -1.0
        end
        e += 1
        for (t, A) in sort(collect(acc); by = first)
            push!(row_ids, e); push!(col_ids, t); push!(blocks, A)
        end
        push!(row_ids, e); push!(col_ids, P + i); push!(blocks, An)
    end
    B = blocksparse(row_ids, col_ids, blocks)

    Q = IPM.allocblockdiag(B); fill!(Q, 0)
    c = zeros(size(B, 2))
    for t in 1:P
        G = Emats[t]' * Emats[t]
        ridge = 1e-10 * tr(G) / m
        block(Q, t, t, t) .= Symmetric(G .+ ridge .* Matrix(I, m, m))
        c[colrange(B, t)] .= Emats[t]' * y[tiles[t][1]:tiles[t][2]]
    end
    for i in 1:N
        c[first(colrange(B, P + i))] = -λ
    end
    g = zeros(size(B, 1))
    K_cones = vcat(IPM.AbstractCone[IPM.CofreeCone() for _ in 1:P],
                   IPM.AbstractCone[IPM.SecondOrderCone() for _ in 1:N])
    prob = IPMProblem(Q, B, c, g, K_cones)

    ctx = (; N, P, Tsz, m, k, K, λ, tiles, nbrs, wts, y, u_true, Emats,
             constterm = 0.5 * (y' * y))
    return prob, ctx
end

function extract_u(prob, ctx, p::AbstractVector)
    u = zeros(ctx.N)
    for g in 1:ctx.N
        t = min((g - 1) ÷ ctx.Tsz + 1, ctx.P)
        θ = p[colrange(prob.B, t)]
        u[g] = eval_row(g, ctx.tiles[t][1], ctx.tiles[t][2], ctx.m)' * θ
    end
    return u
end

"Direct evaluation of F(u(θ)) — assembly check."
function galerkin_objective(ctx, u::AbstractVector)
    F = 0.5 * sum(abs2, u .- ctx.y)
    for i in 1:ctx.N
        F += ctx.λ * norm(sqrt.(ctx.wts[i, :]) .* (u[i] .- u[ctx.nbrs[i, :]]))
    end
    return F
end

# -----------------------------------------------------------------------------
# Baselines (same conic program) and reference solves for gates
# -----------------------------------------------------------------------------

function build_jump_galerkin(prob, ctx, optimizer)
    model = Model(optimizer)
    set_silent(model)
    xs = Vector{Vector{VariableRef}}(undef, ctx.P + ctx.N)
    for t in 1:ctx.P
        xs[t] = @variable(model, [1:ctx.m])
    end
    for i in 1:ctx.N
        v = @variable(model, [1:ctx.K + 1])
        @constraint(model, v in JuMP.SecondOrderCone())
        xs[ctx.P + i] = v
    end
    x = reduce(vcat, xs)
    Qs = sparse(prob.Q); Bs = sparse(prob.B)
    @objective(model, Min, 0.5 * x' * Qs * x - prob.c' * x)
    @constraint(model, Bs * x .== prob.g)
    return model, x
end

"Nodal NLTV via JuMP (E5's formulation, monolithic): completeness reference."
function solve_nodal_nltv(y, nbrs, wts, λ)
    N, K = size(nbrs)
    model = Model(clarabel_opt(; tol = 1e-9))
    set_silent(model)
    @variable(model, u[1:N])
    @variable(model, t[1:N])
    for i in 1:N
        z = [sqrt(wts[i, kk]) * (u[i] - u[nbrs[i, kk]]) for kk in 1:K]
        @constraint(model, vcat(t[i], z) in JuMP.SecondOrderCone())
    end
    @objective(model, Min, 0.5 * sum((u .- y).^2) + λ * sum(t))
    optimize!(model)
    return value.(u), objective_value(model)
end

function solve_local_tv(y::Vector{Float64}, λ::Float64)
    N = length(y)
    model = Model(clarabel_opt(; tol = 1e-8))
    set_silent(model)
    @variable(model, u[1:N])
    @variable(model, s[1:N - 1])
    @constraint(model, [i = 1:N - 1], s[i] ≥ u[i + 1] - u[i])
    @constraint(model, [i = 1:N - 1], s[i] ≥ u[i] - u[i + 1])
    @objective(model, Min, 0.5 * sum((u .- y).^2) + λ * sum(s))
    optimize!(model)
    return value.(u)
end

# -----------------------------------------------------------------------------
# Gate tests
# -----------------------------------------------------------------------------

function test_graph(ctx)
    q = 4
    Pt = patch_matrix(ctx.u_true, q)
    dch = Float64[]
    for i in 1:ctx.N, kk in 1:ctx.K
        push!(dch, sum(abs2, view(Pt, ctx.nbrs[i, kk], :) .- view(Pt, i, :)) / (2q + 1))
    end
    dall = Float64[]
    for i in 1:7:ctx.N, j in max(1, i - 64):min(ctx.N, i + 64)
        abs(j - i) > 2 && push!(dall, sum(abs2, view(Pt, j, :) .- view(Pt, i, :)) / (2q + 1))
    end
    fc = count(<(0.05), dch) / length(dch)
    fa = count(<(0.05), dall) / length(dall)
    @assert fc > 3 * fa "graph no better than chance"
    println("  [PASS] graph: clean matches $(round(100fc))% vs $(round(100fa))% ",
        "chance ($(round(fc / fa, digits = 1))x)")
end

function test_objective_identity(prob, ctx, res)
    u = extract_u(prob, ctx, res.p)
    o_direct = galerkin_objective(ctx, u)
    o_ipm = ipm_objective(prob, res) + ctx.constterm
    rel = abs(o_direct - o_ipm) / abs(o_direct)
    @assert rel < 1e-5 "assembly mismatch: $o_direct vs $o_ipm"
    println("  [PASS] objective identity: F(u) = $(round(o_direct, digits = 5)) ",
        "(rel $(round(rel, sigdigits = 2)))")
    return u
end

function test_model_value(ctx, u_g)
    mse(u) = sum(abs2, u .- ctx.u_true) / ctx.N
    m_y = mse(ctx.y); m_g = mse(u_g)
    best_tv = minimum(mse(solve_local_tv(ctx.y, λ)) for λ in (0.05, 0.1, 0.2, 0.4, 0.8))
    @assert m_g < best_tv < m_y "MSE ordering violated"
    println("  [PASS] model value: MSE noisy $(round(m_y, digits = 4)) > ",
        "best local TV $(round(best_tv, digits = 4)) > Galerkin $(round(m_g, digits = 4))")
end

function test_completeness(settings)
    # full basis + no jets on a small well-conditioned instance == nodal NLTV
    inst = galerkin_instance(; N = 128, Tsz = 16, m = 16, k = -1, K = 6, R = 12, q = 3)
    prob, ctx = build_galerkin(inst)
    res = solve(prob, settings)
    @assert res.status in (OPTIMAL, NEAR_OPTIMAL)
    u_g = extract_u(prob, ctx, res.p)
    F_g = galerkin_objective(ctx, u_g)
    u_n, F_n = solve_nodal_nltv(ctx.y, ctx.nbrs, ctx.wts, ctx.λ)
    @assert abs(F_g - F_n) < 1e-4 * abs(F_n) "completeness violated: $F_g vs $F_n"
    @assert maximum(abs.(u_g .- u_n)) < 1e-3 "completeness: solutions differ"
    println("  [PASS] completeness (Tsz=m=16, no jets): F = $(round(F_g, digits = 6)) ",
        "== nodal $(round(F_n, digits = 6))")
end

function test_ipm_vs_clarabel(prob, ctx, u_ipm)
    model, xv = build_jump_galerkin(prob, ctx, clarabel_opt(; tol = TOL))
    optimize!(model)
    u_cla = extract_u(prob, ctx, value.(xv))
    du = maximum(abs.(u_ipm .- u_cla))
    o_i = galerkin_objective(ctx, u_ipm); o_c = galerkin_objective(ctx, u_cla)
    @assert du < 1e-4 "IPM vs Clarabel mismatch: $du"
    @assert abs(o_i - o_c) < 1e-5 * abs(o_c)
    println("  [PASS] IPM vs Clarabel: ‖Δu‖∞ = $(round(du, sigdigits = 2)); ",
        "obj $(round(o_i, digits = 5)) vs $(round(o_c, digits = 5))")
end

# -----------------------------------------------------------------------------
# Sweep
# -----------------------------------------------------------------------------

sizes() = OPTS.tiny ? [512] :
    OPTS.quick ? [512, 1024, 2048] :
    [512, 1024, 2048, 4096, 8192]

function run()
    println("\n", "=" ^ 78)
    println("  E05: Galerkin-compressed nonlocal TV (dial: N; Tsz = 128, m = 32, K = 12)")
    println("=" ^ 78)

    ipm_settings = IPMSettings{Float64}(
        feas_tol = TOL, gap_tol = TOL, itmax = 200)
    hsd_settings = HSDSettings{Float64}(
        feas_tol = TOL, gap_tol = TOL, itmax = 200)

    println("\n  Gate tests (N = 512):")
    inst = galerkin_instance(; N = 512)
    prob, ctx = build_galerkin(inst)
    test_graph(ctx)
    res = solve(prob, ipm_settings)
    @assert res.status in (OPTIMAL, NEAR_OPTIMAL) "IPM status $(res.status)"
    u_ipm = test_objective_identity(prob, ctx, res)
    test_model_value(ctx, u_ipm)
    test_completeness(ipm_settings)
    test_ipm_vs_clarabel(prob, ctx, u_ipm)
    println()

    cla_opt = clarabel_opt(; tol = TOL)
    msk_opt = OPTS.mosek ? mosek_opt(; tol = TOL) : nothing

    rows = []
    for N in sizes()
        inst = galerkin_instance(; N = N)
        prob, ctx = build_galerkin(inst)
        stats = problem_stats(prob)

        @printf("  N=%-5d dof=%-6d n1=%-5d blk=%-4.0f  ",
            N, stats.N0, stats.N1, stats.med_block)

        m_ipm = measure_ipm(prob, ipm_settings; nruns = NRUNS)
        m_hsd = measure_ipm(prob, hsd_settings; nruns = NRUNS)
        m_cla = measure_jump(() -> first(build_jump_galerkin(prob, ctx, cla_opt)); nruns = NRUNS)
        m_msk = msk_opt !== nothing ?
            measure_jump(() -> first(build_jump_galerkin(prob, ctx, msk_opt)); nruns = NRUNS) :
            (t = NaN, status = "", obj = NaN)

        ratio(b, base) = isfinite(b.t) && isfinite(base.t) ? b.t / base.t : NaN
        fmt_ratio(b, base) = isnan(ratio(b, base)) ? "—" : @sprintf("%.2fx", ratio(b, base))

        println("IPM ", fmt_time(m_ipm), "  HSD ", fmt_time(m_hsd), " (", fmt_ratio(m_hsd, m_ipm),
            ")  Cla ", fmt_time(m_cla), " (", fmt_ratio(m_cla, m_ipm),
            ")  Msk ", fmt_time(m_msk), " (", fmt_ratio(m_msk, m_ipm), ")")

        push!(rows, (size = N, dof = stats.N0, ipm = m_ipm, hsd = m_hsd, cla = m_cla, msk = m_msk))
    end

    dofs = [r.dof for r in rows]
    println("\nFitted log-log slopes (time vs DOF):")
    for (name, get_t) in [("IPM", r -> r.ipm.t), ("HSD", r -> r.hsd.t), ("Clarabel", r -> r.cla.t), ("Mosek", r -> r.msk.t)]
        sl = loglog_slope(dofs, [get_t(r) for r in rows])
        print("  $name: ", isnan(sl) ? "n/a" : @sprintf("DOF^%.2f", sl))
    end
    println()
end

