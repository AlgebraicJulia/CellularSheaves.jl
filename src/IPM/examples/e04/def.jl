# =============================================================================
# e04/def.jl — SOS-certified shape-constrained regression (dense-SDP habitat)
#
# Usage:
#   julia --project e04/run.jl              # Clarabel only
#   julia --project e04/run.jl --mosek      # + Mosek
#   julia --project e04/run.jl --quick      # small sweep
#
# Problem. Nonparametric regression with a CERTIFIED-nonnegative estimator:
# fit noisy samples of f ≥ 0 with a piecewise polynomial of odd degree n on P
# intervals, C^k at the knots, subject to p(x) ≥ 0 everywhere — exactly, not
# on a sample grid. By Markov–Lukács, p ≥ 0 on an interval iff
#     p(u) = (1+u)·v(u)ᵀS v(u) + (1−u)·v(u)ᵀT v(u),   S, T ⪰ 0,
# (local coordinate u ∈ [−1,1], v = Chebyshev basis of degree (n−1)/2), so the
# constraint is two PSD Gram matrices per interval and the fit is a QP over
# SDP cones. This is E03's reference formulation (Papp–Alizadeh 2014) done
# exactly: E03's nonnegative control net is a sufficient-condition surrogate;
# here nonnegativity is certified and the sampled-grid relaxation provably
# converges to this optimum from below (oracle check [4]).
#
# Sheaf structure. One vertex per interval; the stalk is svec(M) of the merged
# Gram M = [[S, X], [Xᵀ, T]] ⪰ 0 (size n+1). The cross block X is free: the
# projection of {M ⪰ 0} onto (S, T) is exactly {S ⪰ 0} × {T ⪰ 0} (take X = 0),
# so the merged single-cone stalk represents the Gram pair without loss — and
# with the ridge the optimizer sets X = 0 itself (oracle: ‖X‖ = 0 at optimum).
# The per-vertex Q is the data Gram + roughness pulled back through the
# Gram→coefficient map: genuinely dense, rank n+1 plus an ε ridge (the SOS
# parameterization is non-unique; ε = 1e-8·tr/dim biases the estimator by
# ~1e-4 per decade — negligible against statistical error). Edges are the
# knots: C^k jet-matching rows (row-normalized; jets scale like n^{2r}),
# non-trivial restriction maps in the tensor_spline style. Chordal
# decomposition sees nothing here: the coupling is through jet functionals of
# Gram entries, not shared matrix indices.
#
# Local Chebyshev basis, not monomial: at n = 25 the effective conditioning of
# the pulled-back Gram is 4.0e2 (Chebyshev) vs 2.9e15 (monomial) — Papp–Yildiz
# 2019 make the same point for interpolant bases at high degree.
#
# Validation (spline_sos_oracle.py, self-contained numpy + cvxpy/Clarabel):
#   * ML representation spans degree-n exactly; positive polys representable
#   * constraint live: unconstrained fit dips −0.05; SOS fit certified ≥ −5e-8
#   * exactness: sampled-grid relaxation objective ↗ SOS objective
#     (gap 5.9e-2 → 3.5e-6 as the grid refines 11 → 641)
#   * merged-stalk formulation ≡ explicit-pair reference: ‖Δp‖∞ = 2e-6
#   * sandwich at benchmark tol 1e-8, n = 9..41: all solved, certified min ≥ 0,
#     obj within [lb, lb + 1.4e-2] of the sampled lower bound
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

struct SOSSplineInstance
    P::Int           # intervals (vertices)
    n::Int           # polynomial degree per interval (odd)
    k::Int           # continuity order at knots (C^k)
    npts::Int        # samples per interval
    λ::Float64       # roughness penalty on ∫ (p'')²
    noise::Float64
    seed::Int
end

sos_instance(; P = 8, n = 13, k = 2, npts = max(60, 8 * (n + 1)),
             λ = 1e-5, noise = 0.08, seed = 1) =
    SOSSplineInstance(P, n, k, npts, λ, noise, seed)

f_true(x) = sin(3π * x)^2 * (1.1 - x)     # ≥ 0, touches 0 at 0, 1/3, 2/3, 1

const RIDGE_EPS = 1e-8

# -----------------------------------------------------------------------------
# Chebyshev machinery (local basis on [−1, 1])
# -----------------------------------------------------------------------------

"Clenshaw evaluation of Σ c[j+1] T_j(u)."
function chebval(u::Real, c::AbstractVector{Float64})
    b1 = 0.0; b2 = 0.0
    for j in length(c):-1:2
        b1, b2 = 2u * b1 - b2 + c[j], b1
    end
    return u * b1 - b2 + c[1]
end

chebvander(u::AbstractVector, n::Int) = begin
    V = zeros(length(u), n + 1)
    V[:, 1] .= 1.0
    n ≥ 1 && (V[:, 2] .= u)
    for j in 3:n + 1
        V[:, j] .= 2.0 .* u .* V[:, j - 1] .- V[:, j - 2]
    end
    V
end

"Coefficients of d/du of a Chebyshev series: d[k-1] = d[k+1] + 2k·c[k], d[0] halved."
function chebder(c::AbstractVector{Float64})
    nn = length(c) - 1
    nn == 0 && return [0.0]
    d = zeros(nn + 2)                       # d[j+1] holds coefficient of T_j
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

"Product of Chebyshev series via T_i·T_j = ½(T_{i+j} + T_{|i−j|})."
function chebmul(a::AbstractVector{Float64}, b::AbstractVector{Float64})
    out = zeros(length(a) + length(b) - 1)
    for (i, ai) in enumerate(a), (j, bj) in enumerate(b)
        ai * bj == 0.0 && continue
        p, q = i - 1, j - 1
        out[p + q + 1] += 0.5 * ai * bj
        out[abs(p - q) + 1] += 0.5 * ai * bj
    end
    return out
end

"Gauss–Legendre nodes/weights on [−1,1] via Golub–Welsch."
function gauss_legendre(ngl::Int)
    β = [kk / sqrt(4kk^2 - 1) for kk in 1:ngl - 1]
    J = SymTridiagonal(zeros(ngl), β)
    F = eigen(J)
    return F.values, 2.0 .* (F.vectors[1, :] .^ 2)
end

# -----------------------------------------------------------------------------
# svec (solver convention: lower triangle, column-major, off-diagonals ×√2)
# -----------------------------------------------------------------------------

svecdim(m::Int) = m * (m + 1) ÷ 2

function ssvec(M::AbstractMatrix)
    m = size(M, 1); v = zeros(svecdim(m)); kk = 0
    for c in 1:m, r in c:m
        kk += 1
        v[kk] = r == c ? M[c, c] : sqrt(2.0) * M[r, c]
    end
    v
end

# -----------------------------------------------------------------------------
# Markov–Lukács maps
# -----------------------------------------------------------------------------

"""
Coefficient map Cmap : svec(M) → Chebyshev coefficients (degree ≤ n),
M = [[S, X], [Xᵀ, T]] of size 2m, m = (n+1)/2, with
p(u) = (1+u)·v(u)ᵀS v(u) + (1−u)·v(u)ᵀT v(u).  X columns are zero.
"""
function merged_coef_map(n::Int)
    @assert isodd(n)
    m = (n + 1) ÷ 2
    Msz = 2m
    nsv = svecdim(Msz)
    Cmap = zeros(n + 1, nsv)
    kk = 0
    for c in 1:Msz, r in c:Msz
        kk += 1
        scale = r == c ? 1.0 : sqrt(2.0)
        if r ≤ m && c ≤ m                              # S block: weight (1+u)
            er = zeros(m); er[r] = 1.0
            ec = zeros(m); ec[c] = 1.0
            pc = chebmul(chebmul(er, ec), [1.0, 1.0]) .* scale
            Cmap[1:length(pc), kk] .+= pc
        elseif r > m && c > m                          # T block: weight (1−u)
            er = zeros(m); er[r - m] = 1.0
            ec = zeros(m); ec[c - m] = 1.0
            pc = chebmul(chebmul(er, ec), [1.0, -1.0]) .* scale
            Cmap[1:length(pc), kk] .+= pc
        end                                            # X block: zero column
    end
    return Cmap, Msz
end

"∫ (p'')² dx over an interval of width h, as a form on Chebyshev coefficients."
function roughness_gram(n::Int, h::Float64; ngl::Int = 64)
    t, wq = gauss_legendre(ngl)
    Vpp = zeros(ngl, n + 1)
    for j in 1:n + 1
        e = zeros(n + 1); e[j] = 1.0
        dd = chebder(e, 2)
        for a in 1:ngl
            Vpp[a, j] = chebval(t[a], dd)
        end
    end
    return (2.0 / h)^4 * (h / 2.0) .* (Vpp' * (wq .* Vpp))
end

"Row functional: r-th x-derivative at local coordinate u (interval width h)."
function jet_row(n::Int, u::Float64, r::Int, h::Float64)
    row = zeros(n + 1)
    for j in 1:n + 1
        e = zeros(n + 1); e[j] = 1.0
        row[j] = chebval(u, chebder(e, r))
    end
    return (2.0 / h)^r .* row
end

# -----------------------------------------------------------------------------
# Build
# -----------------------------------------------------------------------------

function build_sos_spline(inst::SOSSplineInstance)
    @assert isodd(inst.n)
    @assert inst.P ≥ 2
    P, n = inst.P, inst.n
    Cmap, Msz = merged_coef_map(n)
    nsv = svecdim(Msz)
    knots = collect(range(0.0, 1.0; length = P + 1))
    h = 1.0 / P
    rng = MersenneTwister(inst.seed)

    xdat = [sort(rand(rng, inst.npts)) .* h .+ knots[i] for i in 1:P]
    ydat = [f_true.(xdat[i]) .+ inst.noise .* randn(rng, inst.npts) for i in 1:P]

    Qblocks = Vector{Matrix{Float64}}(undef, P)
    cvecs = Vector{Vector{Float64}}(undef, P)
    constterm = 0.0
    R = roughness_gram(n, h)
    for i in 1:P
        u = 2.0 .* (xdat[i] .- knots[i]) ./ h .- 1.0
        E = chebvander(u, n)
        G = E' * E .+ inst.λ .* R
        Qv = 2.0 .* (Cmap' * G * Cmap)
        ridge = RIDGE_EPS * tr(Qv) / nsv
        Qblocks[i] = Qv + ridge * I
        cvecs[i] = -2.0 .* (Cmap' * (E' * ydat[i]))
        constterm += ydat[i]' * ydat[i]
    end
    scale = sum(tr(Q) / nsv for Q in Qblocks) / P     # one COMMON factor: safe
    for i in 1:P
        Qblocks[i] ./= scale
        cvecs[i] ./= scale
    end

    # agreement edges: C^k jets at each interior knot, row-normalized
    row_ids, col_ids, blocks = Int[], Int[], Matrix{Float64}[]
    for e in 1:P - 1
        R1 = vcat([jet_row(n, +1.0, r, h)' for r in 0:inst.k]...) * Cmap
        R2 = vcat([jet_row(n, -1.0, r, h)' for r in 0:inst.k]...) * Cmap
        for r in 1:inst.k + 1
            s = norm(vcat(R1[r, :], R2[r, :]))
            R1[r, :] ./= s; R2[r, :] ./= s
        end
        push!(row_ids, e); push!(col_ids, e); push!(blocks, R1)
        push!(row_ids, e); push!(col_ids, e + 1); push!(blocks, -R2)
    end
    B = blocksparse(row_ids, col_ids, blocks)
    Q = IPM.allocblockdiag(B); fill!(Q, 0)
    for v in 1:P
        block(Q, v, v, v) .= Symmetric(0.5 .* (Qblocks[v] .+ Qblocks[v]'))
    end
    c = vcat(cvecs...)
    g = zeros(size(B, 1))
    K = [SemidefiniteCone() for _ in 1:P]
    prob = IPMProblem(Q, B, c, g, K)

    ctx = (; P, n, k = inst.k, nsv, Msz, Cmap, knots, h, xdat, ydat,
             λ = inst.λ, scale, constterm)
    return prob, ctx
end

"Per-vertex Chebyshev coefficients from a solution vector."
function extract_coefs(prob, ctx, p::AbstractVector)
    [ctx.Cmap * p[colrange(prob.B, v)] for v in 1:ctx.P]
end

function eval_spline(ctx, coefs, xq::AbstractVector)
    yq = similar(xq)
    for (t, xx) in enumerate(xq)
        i = clamp(1 + floor(Int, xx * ctx.P), 1, ctx.P)
        u = 2.0 * (xx - ctx.knots[i]) / ctx.h - 1.0
        yq[t] = chebval(u, coefs[i])
    end
    yq
end

"Full-scale objective from a solution vector (undoes normalization)."
function full_objective(prob, ctx, p::AbstractVector)
    return ctx.scale * (0.5 * dot(p, prob.Q * p) + dot(prob.c, p)) + ctx.constterm
end

# -----------------------------------------------------------------------------
# Baselines (same conic program: merged PSD stalks, jets, quadratic objective)
# -----------------------------------------------------------------------------

function build_jump_sos(prob, ctx, optimizer)
    model = Model(optimizer)
    set_silent(model)
    Msz = ctx.Msz
    Ms = [@variable(model, [1:Msz, 1:Msz] in PSDCone()) for _ in 1:ctx.P]
    svs = [[r == c ? Ms[v][c, c] : sqrt(2.0) * Ms[v][r, c]
            for c in 1:Msz for r in c:Msz] for v in 1:ctx.P]
    x = reduce(vcat, svs)
    Qs = sparse(prob.Q); Bs = sparse(prob.B)
    @objective(model, Min, 0.5 * x' * Qs * x + prob.c' * x)
    @constraint(model, Bs * x .== prob.g)
    return model, x
end

"Sampled-nonnegativity relaxation over raw coefficients (LOWER bound on the
 SOS optimum): plain QP with p(grid) ≥ 0 — solved with Clarabel via JuMP."
function sampled_lower_bound(inst::SOSSplineInstance, ctx; ns = 2001)
    model = Model(clarabel_opt(; tol = TOL))
    set_silent(model)
    P, n = ctx.P, ctx.n
    θ = [@variable(model, [1:n + 1]) for _ in 1:P]
    obj = QuadExpr()
    R = roughness_gram(n, ctx.h)
    for i in 1:P
        u = 2.0 .* (ctx.xdat[i] .- ctx.knots[i]) ./ ctx.h .- 1.0
        E = chebvander(u, n)
        add_to_expression!(obj, θ[i]' * (E' * E .+ ctx.λ .* R) * θ[i])
        add_to_expression!(obj, -2.0 * (E' * ctx.ydat[i])' * θ[i])
    end
    for i in 1:P - 1
        for r in 0:ctx.k
            @constraint(model, jet_row(n, +1.0, r, ctx.h)' * θ[i] ==
                               jet_row(n, -1.0, r, ctx.h)' * θ[i + 1])
        end
    end
    xs = collect(range(0.0, 1.0; length = ns))
    for i in 1:P
        sel = [xx for xx in xs if ctx.knots[i] ≤ xx ≤ ctx.knots[i + 1]]
        isempty(sel) && continue
        u = 2.0 .* (sel .- ctx.knots[i]) ./ ctx.h .- 1.0
        @constraint(model, chebvander(u, n) * θ[i] .≥ 0)
    end
    @objective(model, Min, obj)
    optimize!(model)
    return objective_value(model) + ctx.constterm
end

"Unconstrained fit (equality-constrained least squares via dense KKT)."
function unconstrained_fit(ctx)
    P, n = ctx.P, ctx.n
    nc = n + 1; N = P * nc
    G = zeros(N, N); b = zeros(N)
    R = roughness_gram(n, ctx.h)
    for i in 1:P
        u = 2.0 .* (ctx.xdat[i] .- ctx.knots[i]) ./ ctx.h .- 1.0
        E = chebvander(u, n)
        rows = (i - 1) * nc + 1:i * nc
        G[rows, rows] .= 2.0 .* (E' * E .+ ctx.λ .* R)
        b[rows] .= 2.0 .* (E' * ctx.ydat[i])
    end
    nA = (P - 1) * (ctx.k + 1)
    A = zeros(nA, N); rr = 0
    for i in 1:P - 1, r in 0:ctx.k
        rr += 1
        A[rr, (i - 1) * nc + 1:i * nc] .= jet_row(n, +1.0, r, ctx.h)
        A[rr, i * nc + 1:(i + 1) * nc] .= -jet_row(n, -1.0, r, ctx.h)
    end
    KKT = [G A'; A zeros(nA, nA)]
    sol = KKT \ vcat(b, zeros(nA))
    return [sol[(i - 1) * nc + 1:i * nc] for i in 1:P]
end

# -----------------------------------------------------------------------------
# Gate tests
# -----------------------------------------------------------------------------

function test_representation()
    n = 13
    Cmap, Msz = merged_coef_map(n)
    @assert rank(Cmap) == n + 1 "ML map does not span degree-n polynomials"
    # Slater point: S = T = I/2 ⟹ p(u) = Σ_j T_j(u)² > 0, M = diag ≻ 0,
    # identical on every (uniform) interval ⟹ all jets match.
    m = Msz ÷ 2
    M0 = Matrix(0.5I, Msz, Msz)
    p0 = Cmap * ssvec(M0)
    umin = minimum(chebval(u, p0) for u in range(-1.0, 1.0; length = 2001))
    @assert umin > 0 "Slater point not strictly positive"
    println("  [PASS] ML representation: rank $(rank(Cmap)) == $(n + 1); ",
        "Slater point p_min = $(round(umin, digits = 3)) > 0")
end

function test_constraint_live_and_certified(settings)
    inst = sos_instance(; n = 13)
    prob, ctx = build_sos_spline(inst)
    xg = collect(range(0.0, 1.0; length = 8001))
    fu = eval_spline(ctx, unconstrained_fit(ctx), xg)
    @assert minimum(fu) < -1e-3 "unconstrained fit does not dip negative"
    res = solve(prob, settings)
    @assert res.status in (OPTIMAL, NEAR_OPTIMAL) "IPM status $(res.status)"
    fc = eval_spline(ctx, extract_coefs(prob, ctx, res.p), xg)
    @assert minimum(fc) > -1e-5 "SOS fit not certified nonnegative: $(minimum(fc))"
    println("  [PASS] constraint live: unconstrained min $(round(minimum(fu), digits = 4)) < 0; ",
        "SOS fit min $(round(minimum(fc), sigdigits = 2)) ≥ 0")
    return prob, ctx, res
end

function test_ipm_vs_clarabel(prob, ctx, res)
    model, xv = build_jump_sos(prob, ctx, clarabel_opt(; tol = TOL))
    optimize!(model)
    xg = collect(range(0.0, 1.0; length = 8001))
    f_ipm = eval_spline(ctx, extract_coefs(prob, ctx, res.p), xg)
    coefs_cla = [ctx.Cmap * value.(xv[(v - 1) * ctx.nsv + 1:v * ctx.nsv]) for v in 1:ctx.P]
    f_cla = eval_spline(ctx, coefs_cla, xg)
    derr = maximum(abs.(f_ipm .- f_cla))
    o_ipm = full_objective(prob, ctx, res.p)
    o_cla = full_objective(prob, ctx, value.(xv))
    @assert derr < 1e-3 "IPM vs Clarabel fit mismatch: $derr"
    @assert abs(o_ipm - o_cla) < 1e-3 * abs(o_cla) "objective mismatch"
    println("  [PASS] IPM vs Clarabel: ‖Δp‖∞ = $(round(derr, sigdigits = 2)); ",
        "obj $(round(o_ipm, digits = 5)) vs $(round(o_cla, digits = 5))")
    return o_ipm
end

function test_sandwich(o_ipm)
    inst = sos_instance(; n = 13)
    _, ctx = build_sos_spline(inst)
    lb = sampled_lower_bound(inst, ctx; ns = 2001)
    gap = o_ipm - lb
    @assert gap > -1e-4 * abs(lb) "objective below sampled lower bound: $gap"
    @assert gap < 5e-3 * abs(lb) "sandwich gap suspiciously large: $gap"
    println("  [PASS] sandwich: obj $(round(o_ipm, digits = 5)) ∈ ",
        "[lb, lb + $(round(gap, sigdigits = 2))], lb = $(round(lb, digits = 5))")
end

# -----------------------------------------------------------------------------
# Sweep
# -----------------------------------------------------------------------------

sizes() = OPTS.tiny ? [8] :
    OPTS.quick ? [8, 16, 32] :
    [8, 16, 32, 64, 128]

function run()
    println("\n", "=" ^ 78)
    println("  E04: SOS-certified shape-constrained regression (dial: intervals P; n = 13)")
    println("=" ^ 78)

    settings = IPMSettings{Float64}(
        kkt = UzawaSettings{Float64}(raug = 1e2),
        feas_tol = TOL, gap_tol = TOL, itmax = 200)

    println("\n  Gate tests:")
    test_representation()
    prob13, ctx13, res13 = test_constraint_live_and_certified(settings)
    o_ipm = test_ipm_vs_clarabel(prob13, ctx13, res13)
    test_sandwich(o_ipm)
    println()

    cla_opt = clarabel_opt(; tol = TOL)
    msk_opt = OPTS.mosek ? mosek_opt(; tol = TOL) : nothing

    rows = []
    for P in sizes()
        inst = sos_instance(; P = P, n = 13, λ = 1e-5 / P)
        prob, ctx = build_sos_spline(inst)
        stats = problem_stats(prob)

        @printf("  P=%-3d dof=%-6d n1=%-4d blk=%-4.0f  ",
            P, stats.N0, stats.N1, stats.med_block)

        m_ipm = measure_ipm(prob, settings; nruns = NRUNS)
        m_cla = measure_jump(() -> first(build_jump_sos(prob, ctx, cla_opt)); nruns = NRUNS)
        m_msk = msk_opt !== nothing ?
            measure_jump(() -> first(build_jump_sos(prob, ctx, msk_opt)); nruns = NRUNS) :
            (t = NaN, status = "", obj = NaN)

        ratio(b) = isfinite(b.t) && isfinite(m_ipm.t) ? b.t / m_ipm.t : NaN
        fmt_ratio(b) = isnan(ratio(b)) ? "—" : @sprintf("%.2fx", ratio(b))

        println("IPM ", fmt_time(m_ipm), "  Cla ", fmt_time(m_cla), " (", fmt_ratio(m_cla),
            ")  Msk ", fmt_time(m_msk), " (", fmt_ratio(m_msk), ")")

        push!(rows, (size = P, dof = stats.N0, ipm = m_ipm, cla = m_cla, msk = m_msk))
    end

    dofs = [r.dof for r in rows]
    println("\nFitted log-log slopes (time vs DOF):")
    for (name, get_t) in [("IPM", r -> r.ipm.t), ("Clarabel", r -> r.cla.t), ("Mosek", r -> r.msk.t)]
        sl = loglog_slope(dofs, [get_t(r) for r in rows])
        print("  $name: ", isnan(sl) ? "n/a" : @sprintf("DOF^%.2f", sl))
    end
    println()
end

