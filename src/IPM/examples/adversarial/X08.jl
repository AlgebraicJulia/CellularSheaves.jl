# adversarial/X08.jl — NON-STRICT complementarity [solution-structure pathology]. Per adversarial_spec.md
# queued item 4.
#
# Unlike X03/X04/X06, the pathology here is NOT in the matrix (B is well-conditioned, equilibration-
# invariance is irrelevant) — it is in the SOLUTION. At a strictly-complementary optimum every
# PositiveCone coordinate has exactly one of (pⱼ, dⱼ) zero and the other strictly positive; the central
# path approaches it non-degenerately (pⱼdⱼ ≈ μ with bounded spread). X08 makes a `degen_frac` fraction
# of coordinates NON-strict: BOTH pⱼ* = 0 AND dⱼ* = 0. There the ratio dⱼ/pⱼ (the barrier-Hessian
# diagonal) is 0/0 — the pⱼdⱼ ≈ μ cluster splits and the optimal face acquires positive dimension
# (p* is non-unique). This is the instance the `bar_hdiag_*` columns were added to watch: the mid-cluster
# fraction should collapse as degen_frac rises.
#
# Construction (PositiveCone only, n coords, well-conditioned random B, m = n/2 rows):
#   • degen set (first ⌊degen_frac·n⌋ coords): pⱼ* = dⱼ* = 0        (non-strict complementarity)
#   • remaining coords alternate strictly complementary: primal-active (pⱼ*>0, dⱼ*=0) / dual-active (=0,>0)
# Manufactured optimum: g = B p*, y* random, c = -Bᵀy* - d*  (Q = 0, pure LP). (p*,y*,d*) is a valid KKT
# point; the degen set makes it a DEGENERATE (non-unique) optimum, so the gate checks objective +
# feasibility, not ‖p-p*‖.
#
# Benign twin (`benign=true`): disjoint zero patterns — the degen coords become strictly complementary
# too, so every coordinate has exactly one of (p*,d*) zero. Same B/seed; isolates the non-strictness.

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
    build_nonstrict(; n=50, degen_frac=0.5, benign=false, seed=1) -> (prob, meta)

PositiveCone LP with a manufactured degenerate optimum: a `degen_frac` fraction of coordinates have
pⱼ*=dⱼ*=0 (non-strict complementarity). **`degen_frac` is the dial.** `benign=true` makes those
coordinates strictly complementary instead (disjoint zero patterns) — the matched twin. `meta`: `pstar`,
`dstar`, `both_zero_frac` (realized non-strict fraction), dials.
"""
function build_nonstrict(; n::Int = 50, degen_frac::Float64 = 0.5, benign::Bool = false, seed::Int = 1)
    rng = MersenneTwister(seed)
    m = max(n ÷ 2, 1)
    Bd = randn(rng, m, n)                       # well-conditioned constraint block
    ndeg = round(Int, degen_frac * n)
    pstar = zeros(n); dstar = zeros(n)
    for j in 1:n
        if j ≤ ndeg && !benign
            # non-strict: leave both pstar[j] = dstar[j] = 0
        elseif isodd(j)
            pstar[j] = 1.0 + rand(rng)          # primal-active (dⱼ* = 0)
        else
            dstar[j] = 1.0 + rand(rng)          # dual-active   (pⱼ* = 0)
        end
    end
    ystar = randn(rng, m)
    g = Bd * pstar
    c = .-(Bd' * ystar) .- dstar                 # Q = 0 ⇒ stationarity: c = -Bᵀy* - d*

    B = _scalar_blocksparse(Bd)
    Q = IPM.allocblockdiag(B); fill!(Q, 0)
    K = IPM.AbstractCone[PositiveCone() for _ in 1:n]
    both_zero = count(j -> pstar[j] == 0 && dstar[j] == 0, 1:n)
    meta = (n = n, m = m, degen_frac = degen_frac, ndegen = benign ? 0 : ndeg, benign = benign,
            both_zero_frac = both_zero / n, sigmin = minimum(svdvals(Bd)),
            pstar = pstar, dstar = dstar, Bd = Bd, g = g)
    return IPMProblem(Q, B, c, g, K), meta
end

"Matched benign twin (disjoint zero patterns — strict complementarity everywhere)."
build_nonstrict_twin(; kw...) = build_nonstrict(; benign = true, kw...)

function gate(; tol = 1e-8, kw...)
    prob, meta = build_nonstrict(; kw...)
    obj(p) = -dot(prob.c, p)                     # Q = 0
    @printf("X08 non-strict complementarity  n=%d degen_frac=%.2f  both-zero=%.2f  σmin(B)=%.2e\n",
        meta.n, meta.degen_frac, meta.both_zero_frac, meta.sigmin)
    for (tag, S) in (("HSD", HSDSettings), ("IPM", IPMSettings))
        r = solve!(init(prob, S{Float64}(feas_tol = tol, gap_tol = tol, itmax = 300)))
        objgap = abs(obj(r.p) - obj(meta.pstar)) / (1 + abs(obj(meta.pstar)))
        feas = norm(meta.Bd * r.p .- meta.g, Inf) / (1 + norm(meta.g, Inf))
        @printf("  %-3s status=%-15s niter=%-3d  objgap=%.2e  feas=%.2e  %s\n",
            tag, string(r.status), r.niter, objgap, feas,
            (objgap ≤ 1e-6 && feas ≤ 1e-6) ? "PASS (obj+feas)" : "CHECK")
    end
    return nothing
end

if abspath(PROGRAM_FILE) == @__FILE__
    for df in (0.25, 0.5); gate(; degen_frac = df); end
    println("  -- benign twin --"); gate(; benign = true)
end
