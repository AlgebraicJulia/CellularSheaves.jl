# adversarial/X09.jl — "ceiling hunter": a PLANTED, altitude-controlled, μ-descending ceiling.
#
# The ceiling law (adaptive-raug-spectrum T2/T3):  α_max ≈ ρ_stall · λ_min(H|ker B) / (ε_mach ‖B‖²).
# X09 controls the ceiling by planting the minimal restricted curvature λ_min(H|ker B) at a chosen,
# μ-descending altitude.
#
# CONSTRUCTION NOTE (numpy pre-validation, 2026-07-30): the spec's original "one exactly-zero column"
# idea was REFUTED by the proto — a zero column puts eⱼ∈ker B but leaves pⱼ UNPINNED (B doesn't see it),
# so pⱼ blows up (~1e14) and its curvature is uncontrolled, S-independent (knob dead). A descending
# barrier curvature μ/pⱼ² *requires* pⱼ constraint-pinned, which a zero column forbids. The fix that
# preserves every goal — measured λ_min(H|ker B)=μ/S² to the digit, both knobs live — is a PARALLEL
# (identical) COLUMN PAIR: two coords j₁,j₂ with B_{·j₁}=B_{·j₂}. The constraint pins their SUM (→2S, so
# each→S by symmetry ⇒ hⱼⱼ=μ/S², descending, altitude ∝1/S²); their DIFFERENCE eⱼ₁−eⱼ₂∈ker B exactly.
# Equilibration-robust in the way that matters: per-coordinate D-scaling leaves the columns parallel, so
# a low-curvature ker B direction on the two active coords persists — no diagonal scaling removes it.
#
# Manufactured optimum (PositiveCone, well-conditioned B): planted pair active at p*=S, d*=0; remaining
# coords strictly complementary (mixed primal/dual active); Q=0 (base) ⇒ c=−Bᵀy*−d*. Identical columns
# ⇒ the split is non-unique (degenerate) ⇒ the gate checks objective+feasibility, not ‖p−p*‖.
#
# Variants: X09 base = one pair (clean descent ∝μ/S²). X09c = one pair + ridge Qⱼⱼ (descent FLATTENS when
# μ/S²<Qⱼⱼ — tests whether the tracker stops descending when the physics does; hdiag_min is post-Q so the
# tracked h_min should flatten with it). X09b = two pairs, one ridged, so the binding (minimal) curvature
# SWITCHES mid-solve (kinked descent — tests every-stall recalibration). Benign twin: same builder, p*=O(1)
# (via S≈1), so no anomalously-small planted curvature.

using CellularSheaves
using CellularSheaves.IPM: IPMProblem, IPMSettings, HSDSettings, OPTIMAL, NEAR_OPTIMAL, solve!, init,
    PositiveCone
import CellularSheaves.IPM as IPM
using CellularSheaves.BlockSparseArrays: blocksparse
using LinearAlgebra, SparseArrays, Random, Printf

function _scalar_blocksparse(Bd::AbstractMatrix)
    m, nv = size(Bd)
    rows = Int[]; cols = Int[]; blks = Matrix{Float64}[]
    for j in 1:nv, i in 1:m
        v = Bd[i, j]
        iszero(v) && continue
        push!(rows, i); push!(cols, j); push!(blks, fill(float(v), 1, 1))
    end
    return blocksparse(rows, cols, blks)
end

"""
    build_ceiling(; n=50, pairs=[(S=1e2, ridge=0.0)], benign=false, seed=1) -> (prob, meta)

PositiveCone LP with `length(pairs)` planted parallel-column pairs. Pair `k` plants identical columns at
altitude `pairs[k].S` (barrier curvature μ/S², the ceiling driver) plus optional `pairs[k].ridge` on Q
(floors that curvature). `benign=true` sets every planted p*=O(1) (no small curvature). `meta` carries
`pstar`, `S`/`ridge` per pair, and the planted-column indices.
"""
function build_ceiling(; n::Int = 50, pairs = [(S = 1e2, ridge = 0.0)], benign::Bool = false, seed::Int = 1)
    rng = MersenneTwister(seed)
    npair = length(pairs)
    npl = 2 * npair
    @assert n ≥ npl + 2
    m = max(n ÷ 2, 1)
    Bd = randn(rng, m, n); Bd ./= opnorm(Bd)          # well-conditioned, ‖B‖₂ ≈ 1
    pstar = zeros(n); dstar = zeros(n); Qdiag = zeros(n)
    planted = Int[]
    for k in 1:npair
        j1 = 2k - 1; j2 = 2k
        Bd[:, j2] .= @view Bd[:, j1]                  # identical columns ⇒ eⱼ₁−eⱼ₂ ∈ ker B
        val = benign ? 1.0 + rand(rng) : pairs[k].S
        pstar[j1] = pstar[j2] = val                   # active, p=S each (sum pinned by the constraint)
        Qdiag[j1] = Qdiag[j2] = pairs[k].ridge
        push!(planted, j1, j2)
    end
    for j in (npl + 1):n                              # generic strictly-complementary coords
        isodd(j) ? (pstar[j] = 1.0 + rand(rng)) : (dstar[j] = 1.0 + rand(rng))
    end
    ystar = randn(rng, m)
    g = Bd * pstar
    c = Qdiag .* pstar .- Bd' * ystar .- dstar        # c = Qp* − Bᵀy* − d*

    B = _scalar_blocksparse(Bd)
    Q = IPM.allocblockdiag(B); fill!(Q, 0)
    for j in 1:n; iszero(Qdiag[j]) || (IPM.block(Q, j, j, j) .= Qdiag[j]); end
    K = IPM.AbstractCone[PositiveCone() for _ in 1:n]
    meta = (n = n, m = m, npair = npair, S = Float64[p.S for p in pairs], ridge = Float64[p.ridge for p in pairs],
            benign = benign, planted = planted, sigmin = minimum(svdvals(Bd)),
            pstar = pstar, dstar = dstar, Qdiag = Qdiag, Bd = Bd, g = g)
    return IPMProblem(Q, B, c, g, K), meta
end

"X09 base — single clean pair. `S` is the altitude dial (ceiling ∝ μ/S²)."
build_ceiling_base(; S::Float64 = 1e2, kw...) = build_ceiling(; pairs = [(S = S, ridge = 0.0)], kw...)

"X09c — single pair + ridge; the μ/S² descent flattens at `ridge`."
build_ceiling_ridge(; S::Float64 = 1e2, ridge::Float64 = 1e-10, kw...) =
    build_ceiling(; pairs = [(S = S, ridge = ridge)], kw...)

"X09b — two pairs, the high-altitude one ridged, so the binding curvature switches mid-solve (kink)."
build_ceiling_kink(; S1::Float64 = 1e4, S2::Float64 = 1e2, ridge1::Float64 = 1e-9, kw...) =
    build_ceiling(; pairs = [(S = S1, ridge = ridge1), (S = S2, ridge = 0.0)], kw...)

# λ_min(H|ker B) at a surrogate central-path point for parameter μ: planted/primal-active coords sit at
# p*, dual-active coords approach the boundary p(μ)=μ/d*; H = diag(μ/p² ) + Q (post-Q, as the solver sees).
function _lmin_HkerB(meta, μ)
    # planted/primal-active coords sit at p*; dual-active coords are CLAMPED to p=1 (moderate curvature μ)
    # rather than the true path value p=μ/d* (curvature d²/μ, ~1e6). They are never the minimum anyway
    # (large curvature), and clamping keeps the tiny planted eigenvalue μ/S² resolvable in float64.
    p = similar(meta.pstar)
    @inbounds for j in 1:meta.n
        p[j] = meta.pstar[j] > 0 ? meta.pstar[j] : 1.0
    end
    h = μ ./ p .^ 2 .+ meta.Qdiag
    N = nullspace(Matrix(meta.Bd))
    return eigmin(Symmetric(N' * (h .* N)))
end

"Requirement (3): assert the planted ceiling law λ_min(H|ker B) = μ/S²·(1±tol) across μ and both S knobs,
and that the benign twin plants no anomalously-small restricted curvature (stays ≥ O(μ))."
function check_ceiling_law(; seed = 1)
    println("== ceiling law: λ_min(H|ker B) vs μ/S²  (single un-ridged pair) ==")
    for S in (1e2, 1e4)
        _, meta = build_ceiling_base(; S = S, seed = seed)
        for μ in (1e-2, 1e-4, 1e-6)
            lm = _lmin_HkerB(meta, μ); pred = μ / S^2
            @printf("  S=%.0e μ=%.0e | λmin=%.3e  μ/S²=%.3e  ratio=%.3f  %s\n",
                S, μ, lm, pred, lm / pred, abs(lm / pred - 1) ≤ 0.05 ? "OK" : "OFF")
        end
    end
    _, mb = build_ceiling_base(; benign = true, seed = seed)
    for μ in (1e-2, 1e-6)
        lm = _lmin_HkerB(mb, μ)
        @printf("  benign  μ=%.0e | λmin=%.3e  (≥O(μ)? %s — no planted low direction)\n",
            μ, lm, lm ≥ 0.1μ ? "OK" : "OFF")
    end
end

function gate(builder; tol = 1e-8, kw...)
    prob, meta = builder(; kw...)
    obj(p) = 0.5 * dot(p, meta.Qdiag .* p) - dot(prob.c, p)
    condBBt = cond(Matrix(meta.Bd * meta.Bd'))                 # requirement (2): row-space conditioning
    Sstr = join(("$(round(Int,log10(s)))" for s in meta.S), ",")
    @printf("X09 ceiling  n=%d pairs=%d  S=10^{%s}  ridge=%s  σmin(B)=%.2e  cond(BBᵀ)=%.1e %s\n",
        meta.n, meta.npair, Sstr, string(meta.ridge), meta.sigmin, condBBt,
        condBBt ≤ 1e4 ? "(modest ✓ — orthogonal to X04)" : "(!! near-dependent rows leaked in)")
    for (tag, St) in (("HSD", HSDSettings), ("IPM", IPMSettings))
        r = solve!(init(prob, St{Float64}(feas_tol = tol, gap_tol = tol, itmax = 300)))
        objgap = abs(obj(r.p) - obj(meta.pstar)) / (1 + abs(obj(meta.pstar)))
        feas = norm(meta.Bd * r.p .- meta.g, Inf) / (1 + norm(meta.g, Inf))
        # requirement (1): the pair splits evenly (symmetry) — assert on the SOLVED point, not assumed
        split = maximum(abs(r.p[2k-1] - r.p[2k]) / (1 + abs(r.p[2k-1])) for k in 1:meta.npair)
        @printf("  %-3s status=%-15s niter=%-3d  objgap=%.2e  feas=%.2e  split=%.2e  %s\n",
            tag, string(r.status), r.niter, objgap, feas, split,
            (objgap ≤ 1e-6 && feas ≤ 1e-6 && split ≤ 1e-4) ? "PASS (obj+feas+split)" : "CHECK")
    end
    return nothing
end

if abspath(PROGRAM_FILE) == @__FILE__
    check_ceiling_law()
    println("-- X09 base --");         gate(build_ceiling_base; S = 1e2); gate(build_ceiling_base; S = 1e4)
    println("-- X09c ridge --");       gate(build_ceiling_ridge)
    println("-- X09b kink --");        gate(build_ceiling_kink)
    println("-- benign twin --");      gate(build_ceiling_base; benign = true)
end
