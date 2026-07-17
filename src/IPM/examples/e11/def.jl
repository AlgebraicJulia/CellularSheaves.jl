# =============================================================================
# e11/def.jl — Short multichannel FIR equalizer-bank design under colored noise
#          (regularized MINT: the D = 1 synthesis filter bank; dense-Q + SOC
#          at E01's block sizes, with the gluing maps a Cholesky factor)
#
# Usage:
#   julia --project e11/run.jl              # Clarabel baseline (exact form)
#   julia --project e11/run.jl --mosek      # + Mosek
#   julia --project e11/run.jl --quick      # quick (smaller) sweep
#   julia --project e11/run.jl --budget     # budget-split variant (O(M) coupling, star)
#   julia --project e11/run.jl --tree       # budget-tree variant (O(M) coupling, tree)
#   julia --project e11/run.jl --white      # diagonal-Q ablation (pre-registered)
#   julia --project e11/run.jl --ldial      # block-size dial Lf at fixed M
#
# Problem. M sensors observe a source through reverberant FIR channels h_m
# (length Lh, exponentially decaying, unit norm); design one equalizer f_m
# (length Lf) per channel so the composite response c = Σ_m h_m ∗ f_m = T f
# approximates a pure delay g = e_d, while minimizing the amplified ambient-
# noise power at the output:
#
#     min_f  ½ Σ_m f_m' Q_m f_m,   Q_m = ς_m² Θ + λ I
#     s.t.   ‖T f − g‖ ≤ ε,        T = [Toep(h_1) … Toep(h_M)]
#
# Θ is the autocorrelation Gram of the AR(2) ambient noise field (a shared
# interference peak — dense BY NATURE because the noise is colored), ς_m are
# heterogeneous per-sensor noise scales (mics of different quality/distance —
# what makes the reconstruction budget allocation load-bearing), and λ is the
# white sensor self-noise floor (strict convexity ⟹ unique optimum, E10's
# uniqueness discipline — physics, not decoration). Lineage: MINT, Miyoshi &
# Kaneda (IEEE TASSP 1988); regularized/short multichannel equalization
# (Kodrasi–Doclo lineage); channel shortening (Melsa–Younce–Rohrs 1996);
# Bezout/coprimality (Kailath); eigenfilters (Vaidyanathan–Nguyen 1987).
#
# SHORT-EQUALIZER (TALL) REGIME. The benchmark lives BELOW the MINT length
# threshold Lf* = ⌈(Lh−1)/(M−1)⌉: T is tall with full column rank, exact
# inversion is impossible (ε_LS > 0), and the design trades residual ISI
# against noise amplification — the practically relevant regime (exact MINT
# is notoriously noise-fragile). ε_LS = ‖(I − T G⁻¹T')g‖ is the short-
# equalizer floor, G = T'T. The oracle certifies the threshold classically:
# ε_LS(Lf) drops from 2.2e-2 at Lf* − 1 to 2.6e-14 at Lf* (Miyoshi–Kaneda),
# and a shared spectral null keeps ε_LS ≥ √(2/P) analytically above the
# threshold (Bezout falsification, bound EXACTLY tight at 0.1085).
#
# WHITENING (E06's trick, relocated from objective to CONSTRAINT). With
# G = LL' and b = L⁻¹T'g = L'f̂ (f̂ = G⁻¹T'g the pinv bank):
#
#     ‖T f − g‖² = ‖L'f − b‖² + ε_LS²           (joint identity, exact)
#
# so with ε = hypot(ε_LS, εx) the ball becomes ‖L'f − b‖ ≤ εx: the excess
# budget εx IS the whitened radius. Endpoints: εx → 0 pins f = f̂ (measured
# deformation rate εx^1.00); ε ≥ ‖g‖ = 1 admits f = 0.
#
# SHEAF STRUCTURE (exact slice split — the flagship). The whitened residual
# partitions EXACTLY by Cholesky row blocks, ‖L'f − b‖² = Σ_i ‖(L'f − b)_i‖²,
# giving zero-conservatism decomposition:
#   * M filter stalks f_m ∈ R^Lf (CofreeCone) carrying dense Q_m — density
#     in Q itself, the E03/E08/E09 mechanism;
#   * M slice SOC stalks (s_i, w_i) ∈ R^{Lf+1} with dense triangular ties
#       w_i − Σ_{j≥i} (L')_{ij} f_j = −b_i        (nonzero rhs, E06's pattern)
#     — the gluing maps ARE the Cholesky blocks of the joint convolution
#     Gram, dense by nature (measured densities 1.000 at 1e-8);
#   * one (M+1)-dim hub SOC (h0, ŝ) with one-row s-links ŝ_i = s_i and the
#     pin h0 = εx, enforcing √(Σ s_i²) ≤ εx, hence ‖L'f − b‖ ≤ εx exactly.
# Coupling mass ~ M²Lf²/2 (slice i touches filters i..M): a pre-registered
# falsifiable prediction for the M dial — the per-iteration B cost grows one
# power of M faster than the vertex factorizations, so the M-dial slope
# should show it if Uzawa applications dominate.
#
# BUDGET SPLIT (--budget: the O(M)-coupling variant, E10's Agler move
# relocated to the reconstruction budget). Center the target split at the
# pinv bank, g_m = T_m f̂_m + r̂/M (r̂ ⊥ range(T_m) makes the split identity
# exact); one SOC stalk (s_m, w_m, σ_m) per channel with dense ties
# w_m = L_m'f_m − b_m (G_m = T_m'T_m = L_m L_m', the per-channel
# autocorrelation Gram), pin σ_m = ε_LS/M, and ONE budget row Σ s_m = ε.
# Feasibility CERTIFIES ‖Tf − g‖ ≤ ε by the triangle inequality — a
# compositional certificate whose conservatism is measured, not hidden:
# v_budget/v_joint = 12.09 → 2.04 → 1.18 → 1.04 over εx = 0.2 → 0.003
# (→ 1 at the pinv endpoint), and at the operating point the budget form
# uses only 0.33 of the whitened ball vs the exact form's 1.000. The
# certified budget buys O(M) coupling mass; the benchmark measures both
# sides of the tradeoff.
#
# MEASURED EPILOGUE (post-first-run; the star-vs-tree lesson). The
# --budget split as originally designed aggregates the M per-filter slack
# scalars in ONE arity-M hub SOC of dim M+2. Coupling mass is O(M), but
# the hub is a star: a single vertex of unbounded dimension, (M+2)³
# per-iteration cost, measured slope DOF^1.97 and a 6× LOSS to Clarabel
# at M = 32 (whose presolve sees a nearly separable program plus one
# coupling row: slope 1.04). --tree replaces the star by a binary hypot
# tree of dim-3 combiners s_parent = ‖(s_left, s_right)‖ — exact by
# nested norms, O(M) extra tiny SOCs, bounded degree — and the slope
# drops to DOF^0.99: 36.5 ms at M = 32 (was 585 ms), 1.77× WIN vs
# Clarabel. Diagnosis confirmed: the bottleneck was hub TOPOLOGY, not the
# budget formulation. Habitat rule refined accordingly: bounded coupling
# mass is necessary but not sufficient — vertex arity and max block
# dimension must be bounded too. Recommendation flips with scale:
# flagship (zero conservatism) for M ≤ 16; --tree beyond.
#
# WHY D = 1 (the polyphase exclusion). For a decimated filter bank (D > 1)
# the per-filter constraint Gram is EXACTLY polyphase-block-diagonal:
# (T_m'T_m)[j,j'] = 0 whenever j ≢ j' (mod D) — structural zeros no
# whitening can remove (measured: off-class entries exactly 0.0 at D = 4).
# Structural 1/D block density is sheaf shape without sheaf density, in B
# rather than Q: the decimated bank is out of habitat BY STRUCTURE and is
# pre-registered as a boundary-suite control (predicted bloat ~ D). E14 is
# the D = 1 member of the filter-bank family, where every Gram — Θ, G_m,
# L_m', joint G — is honestly dense (all measured 1.000 at 1e-8).
#
# Habitat claim: dense-Q + SOC (E09's cell) at E01's block sizes — med block
# Lf + 1 ≈ 65 vs E09's 6 — with a new edge species (Cholesky/convolution
# ties) and non-path coupling topology (triangular slice hyperedges + hub).
# Q ≡ 0 on cone stalks, no free auxiliary stalks, bloat ~ 1.
#
# Validation (tools/mint_equalizer_oracle.py, self-contained numpy + scipy
# + cvxpy, fully executed; figures at M = 8, Lh = 1024, Lf = 64, trev = 256,
# AR(2) r = 0.97 φ = 0.55π, λ = 1e-3, εx = 0.05, d = 544; ε_LS = 0.7513):
#   * DENSE BY NATURE + EXCLUSION: densities@1e-8 Θ 1.000, G_m ≥ 1.000,
#     L_m' ≥ 1.000, joint G 1.000; D = 4 control off-class Gram entries
#     EXACTLY 0.
#   * MINT / BEZOUT (classical certificates): threshold ε_LS 2.17e-2 →
#     2.56e-14 at Lf* = 40; shared spectral null keeps ε_LS = 0.1085 ≥ the
#     analytic projection bound 0.1085 (~√(2/P)) — the bound is TIGHT.
#   * FORMULATION CERTIFICATE: whitening identities at random points 0.0 /
#     4.1e-16; the budget stalk assembly == native to ‖Δf‖∞ = 8.6e-7 (rel
#     Δv 1.4e-8); the exact stalk assembly (exactly what this file builds)
#     == native to ‖Δf‖∞ = 4.7e-7 (rel Δv 1.0e-9).
#   * TRS (classical certificate; the suite's species after the Kiefer–Wolfowitz gate of tools/doptimal_oracle.py): the
#     joint problem is a trust-region subproblem — the secular-equation
#     solve (Gay 1981; Moré–Sorensen 1983) matches the conic solver to rel
#     Δv 4.2e-9, with boundary residual 2.1e-17 and stationarity 1.7e-13.
#   * EXACTNESS + CONSERVATISM: v_exact == v_joint to 2.8e-10 with ball use
#     1.000; budget certificate holds with ball use 0.329; conservatism
#     curve as above.
#   * ALLOCATION (budget form): optimized budget beats uniform by 1.9%.
#   * DEFORMATION: ‖f(εx) − f̂‖∞ 1.55e-3 → 3.90e-4 over εx 0.02 → 0.005
#     (rate εx^1.00).
#   * NOISE SHAPING (diagonal-Q control): the white-noise design moves the
#     solution by 12.1% and pays ×1.77 amplified colored-noise power — the
#     --white ablation row, pre-registered (D6/smooth-bary mechanism).
#   * EQUALIZATION: the designed bank cuts output noise power to 0.450× the
#     pinv bank (predicted 0.453× by the quadratic forms) at a certified
#     ISI budget.
#   * CROSS-SOLVER: Clarabel vs SCS ‖Δf‖∞ = 1.1e-7.
# The gate tests below re-derive the fast subset in-process, including the
# TRS comparator via a native eigendecomposition + bisection (no cvxpy).
#
# Dial: channel count M at fixed Lf = 64 with Lh = 2 M Lf (regime-preserving:
# the MINT threshold stays ≈ 2Lf·M/(M−1) ≥ 2Lf, comfortably tall, and the
# aspect of T is constant across the dial). DOF = M(2Lf + 2) + 1;
# N1 = M(Lf + 1) + 1; med block ≈ Lf. --ldial sweeps Lf at fixed M = 8
# (Lh = 16 Lf) — the block-size dial across the med_blk regimes.
#
# STATUS: CONFIRMED. IPM wins vs both Clarabel (3.4-6.8x) and Mosek (1.5-3.2x)
# at M = 4-16. The M² coupling-mass prediction held: IPM slope DOF^2.29 vs
# Mosek's 1.72 — the triangular slice hyperedges are visible. raug = 1e4 OK.
# =============================================================================

using AppleAccelerate
using LinearAlgebra, SparseArrays, Printf, Random
import Metis, CliqueTrees

include("../utils.jl")

using CellularSheaves.IPM: SecondOrderCone, CofreeCone
using CellularSheaves.BlockSparseArrays: rowrange

const OPTS = parse_args(ARGS)
const TOL = OPTS.tol
const NRUNS = 5
const BUDGET = "--budget" in ARGS
const TREE = "--tree" in ARGS
const WHITE = "--white" in ARGS
const LDIAL = "--ldial" in ARGS

# -----------------------------------------------------------------------------
# Problem definition
# -----------------------------------------------------------------------------

struct EqualizerInstance
    Lf::Int                    # equalizer length (block size)
    lhfac::Int                 # Lh = lhfac * M * Lf (regime-preserving)
    trevfrac::Float64          # reverb decay = trevfrac * Lh
    r::Float64; phi::Float64   # AR(2) ambient-noise poles r e^{±iφ}
    lam::Float64               # white sensor self-noise floor
    epsx::Float64              # excess reconstruction budget (whitened radius)
    sglo::Float64; sghi::Float64  # per-sensor noise scale range (geometric)
    seed::Int
end

equalizer_instance(; Lf = 64, lhfac = 2, trevfrac = 0.25, r = 0.97,
                   phi = 0.55π, lam = 1e-3, epsx = 0.05, sglo = 0.5,
                   sghi = 2.0, seed = 7) =
    EqualizerInstance(Lf, lhfac, trevfrac, r, phi, lam, epsx, sglo, sghi, seed)

"Autocorrelation Toeplitz of the AR(2) ambient noise, normalized diag = 1."
function ar2_gram(Lf, r, phi; nimp = 8192)
    a1, a2 = 2r * cos(phi), -r^2
    ψ = zeros(nimp)
    ψ[1] = 1.0
    ψ[2] = a1
    for t in 3:nimp
        ψ[t] = a1 * ψ[t - 1] + a2 * ψ[t - 2]
    end
    ρ = [dot(view(ψ, 1:nimp - k), view(ψ, k + 1:nimp)) for k in 0:Lf - 1]
    ρ ./= ρ[1]
    return [ρ[abs(i - j) + 1] for i in 1:Lf, j in 1:Lf]
end

"Tall convolution block Toep(h) ∈ R^{(Lh+Lf-1) × Lf}."
function toepblock(h, Lf)
    Lh = length(h)
    T = zeros(Lh + Lf - 1, Lf)
    for j in 1:Lf
        T[j:j + Lh - 1, j] .= h
    end
    return T
end

"""Channels + noise model + whitening data. Every Gram is dense by nature:
Θ (colored noise), G_m and L_m' (reverberant autocorrelations, Lf ≤ Lh),
joint G (cross-channel correlations)."""
function make_system(inst::EqualizerInstance, M::Int)
    Lf = inst.Lf
    Lh = inst.lhfac * M * Lf
    rng = MersenneTwister(inst.seed)
    env = exp.(-(0:Lh - 1) ./ (inst.trevfrac * Lh))
    chans = [begin h = randn(rng, Lh) .* env; h ./ norm(h) end for _ in 1:M]
    sigmas = exp.(range(log(inst.sglo), log(inst.sghi); length = M))
    Θ = ar2_gram(Lf, inst.r, inst.phi)
    Qms = [sigmas[m]^2 .* Θ .+ inst.lam .* Matrix{Float64}(I, Lf, Lf)
           for m in 1:M]
    P = Lh + Lf - 1
    d = (Lh + Lf) ÷ 2                       # delay (0-based lag)
    Ts = [toepblock(h, Lf) for h in chans]
    T = hcat(Ts...)
    g = zeros(P); g[d + 1] = 1.0
    G = Symmetric(T' * T)
    chol = cholesky(G)                       # per-block SPD gate, E01 style
    L = Matrix(chol.L)
    fhat = chol \ (T' * g)
    rhat = g - T * fhat
    eLS = norm(rhat)
    b = L \ (T' * g)                         # = L' fhat
    Gms = [Symmetric(Tm' * Tm) for Tm in Ts]
    LTms = [Matrix(cholesky(Gm).U) for Gm in Gms]     # L_m'
    fhm = [fhat[(m - 1) * Lf + 1:m * Lf] for m in 1:M]
    bms = [LTms[m] * fhm[m] for m in 1:M]
    eps = hypot(eLS, inst.epsx)
    return (; M, Lh, Lf, P, d, chans, sigmas, Θ, Qms, Ts, T, g, L, LT = Matrix(L'),
            fhat, fhm, rhat, eLS, b, Gms, bms, LTms, eps)
end

# -----------------------------------------------------------------------------
# Build
# -----------------------------------------------------------------------------
#
# Exact slice split (default). Stalk order: filter stalks 1..M (dim Lf,
# CofreeCone) | slice SOC stalks M+i (dim Lf+1, (s_i, w_i)) | hub SOC 2M+1
# (dim M+1, (h0, ŝ)). Edges: slice ties (Lf rows, w_i − Σ_{j≥i}(L')_{ij}f_j
# = −b_i), s-links (ŝ_i − s_i = 0), hub pin (h0 = εx).
#
# Budget split (--budget). Stalk order: filter stalks 1..M | SOC stalks M+m
# (dim Lf+2, (s_m, w_m, σ_m)). Edges: ties (w_m − L_m'f_m = −b_m), pins
# (σ_m = ε_LS/M), one budget row (Σ_m s_m = ε).

function build_equalizer(sys, inst::EqualizerInstance;
                         form = TREE ? :budget_tree : (BUDGET ? :budget : :exact), white = WHITE)
    M, Lf = sys.M, sys.Lf
    vfil(m) = m
    vsoc(m) = M + m
    row_ids, col_ids, blocks = Int[], Int[], Matrix{Float64}[]
    rhs_val = Dict{Int, Vector{Float64}}()
    ILf = Matrix{Float64}(I, Lf, Lf)
    e = 0

    if form == :exact
        dsoc = Lf + 1
        vhub = 2M + 1
        for i in 1:M                          # ---- slice ties (dense L' blocks)
            e += 1
            Aw = zeros(Lf, dsoc); Aw[:, 2:Lf + 1] .= ILf
            push!(row_ids, e); push!(col_ids, vsoc(i)); push!(blocks, Aw)
            ri = (i - 1) * Lf + 1:i * Lf
            for j in i:M                      # upper-triangular coupling
                cj = (j - 1) * Lf + 1:j * Lf
                push!(row_ids, e); push!(col_ids, vfil(j))
                push!(blocks, -sys.LT[ri, cj])
            end
            rhs_val[e] = -sys.b[ri]
        end
        for i in 1:M                          # ---- s-links ŝ_i = s_i
            e += 1
            Ah = zeros(1, M + 1); Ah[1, 1 + i] = 1.0
            As = zeros(1, dsoc); As[1, 1] = -1.0
            push!(row_ids, e); push!(col_ids, vhub); push!(blocks, Ah)
            push!(row_ids, e); push!(col_ids, vsoc(i)); push!(blocks, As)
        end
        e += 1                                # ---- hub pin h0 = εx
        Ap = zeros(1, M + 1); Ap[1, 1] = 1.0
        push!(row_ids, e); push!(col_ids, vhub); push!(blocks, Ap)
        rhs_val[e] = [inst.epsx]
        K_cones = vcat(IPM.AbstractCone[CofreeCone() for _ in 1:M],
                       IPM.AbstractCone[SecondOrderCone() for _ in 1:M + 1])
        stalks = vcat([(Lf, :free) for _ in 1:M],
                      [(dsoc, :soc) for _ in 1:M], [(M + 1, :soc)])
    elseif form == :budget
        dsoc = Lf + 2
        for m in 1:M                          # ---- ties w_m − L_m'f_m = −b_m
            e += 1
            Aw = zeros(Lf, dsoc); Aw[:, 2:Lf + 1] .= ILf
            push!(row_ids, e); push!(col_ids, vsoc(m)); push!(blocks, Aw)
            push!(row_ids, e); push!(col_ids, vfil(m)); push!(blocks, -sys.LTms[m])
            rhs_val[e] = -sys.bms[m]
        end
        for m in 1:M                          # ---- pins σ_m = ε_LS/M
            e += 1
            Ap = zeros(1, dsoc); Ap[1, dsoc] = 1.0
            push!(row_ids, e); push!(col_ids, vsoc(m)); push!(blocks, Ap)
            rhs_val[e] = [sys.eLS / M]
        end
        e += 1                                # ---- budget row Σ s_m = ε
        for m in 1:M
            As = zeros(1, dsoc); As[1, 1] = 1.0
            push!(row_ids, e); push!(col_ids, vsoc(m)); push!(blocks, As)
        end
        rhs_val[e] = [sys.eps]
        K_cones = vcat(IPM.AbstractCone[CofreeCone() for _ in 1:M],
                       IPM.AbstractCone[SecondOrderCone() for _ in 1:M])
        stalks = vcat([(Lf, :free) for _ in 1:M], [(dsoc, :soc) for _ in 1:M])
    elseif form == :budget_tree
        dsoc = Lf + 2
        for m in 1:M                          # ---- ties w_m − L_m'f_m = −b_m
            e += 1
            Aw = zeros(Lf, dsoc); Aw[:, 2:Lf + 1] .= ILf
            push!(row_ids, e); push!(col_ids, vsoc(m)); push!(blocks, Aw)
            push!(row_ids, e); push!(col_ids, vfil(m)); push!(blocks, -sys.LTms[m])
            rhs_val[e] = -sys.bms[m]
        end
        for m in 1:M                          # ---- pins σ_m = ε_LS/M
            e += 1
            Ap = zeros(1, dsoc); Ap[1, dsoc] = 1.0
            push!(row_ids, e); push!(col_ids, vsoc(m)); push!(blocks, Ap)
            rhs_val[e] = [sys.eLS / M]
        end
        ncomb = M - 1
        vcomb(k) = 2M + k
        nodes = collect(1:M)
        node_is_leaf = trues(M)
        combiner_idx = 0
        combiner_children = Vector{Tuple{Int,Bool,Int,Bool}}()
        while length(nodes) > 1
            new_nodes = Int[]
            new_is_leaf = Bool[]
            i = 1
            while i <= length(nodes)
                if i + 1 <= length(nodes)
                    combiner_idx += 1
                    push!(combiner_children, (nodes[i], node_is_leaf[i],
                                              nodes[i+1], node_is_leaf[i+1]))
                    push!(new_nodes, combiner_idx)
                    push!(new_is_leaf, false)
                    i += 2
                else
                    push!(new_nodes, nodes[i])
                    push!(new_is_leaf, node_is_leaf[i])
                    i += 1
                end
            end
            nodes = new_nodes
            node_is_leaf = new_is_leaf
        end
        root_is_leaf = node_is_leaf[1]
        root_idx = nodes[1]
        for (k, (left, left_leaf, right, right_leaf)) in enumerate(combiner_children)
            e += 1
            Al = zeros(1, 3); Al[1, 2] = 1.0
            push!(row_ids, e); push!(col_ids, vcomb(k)); push!(blocks, Al)
            if left_leaf
                As = zeros(1, dsoc); As[1, 1] = -1.0
                push!(row_ids, e); push!(col_ids, vsoc(left)); push!(blocks, As)
            else
                Ac = zeros(1, 3); Ac[1, 1] = -1.0
                push!(row_ids, e); push!(col_ids, vcomb(left)); push!(blocks, Ac)
            end
            e += 1
            Ar = zeros(1, 3); Ar[1, 3] = 1.0
            push!(row_ids, e); push!(col_ids, vcomb(k)); push!(blocks, Ar)
            if right_leaf
                As = zeros(1, dsoc); As[1, 1] = -1.0
                push!(row_ids, e); push!(col_ids, vsoc(right)); push!(blocks, As)
            else
                Ac = zeros(1, 3); Ac[1, 1] = -1.0
                push!(row_ids, e); push!(col_ids, vcomb(right)); push!(blocks, Ac)
            end
        end
        e += 1
        if root_is_leaf
            Ap = zeros(1, dsoc); Ap[1, 1] = 1.0
            push!(row_ids, e); push!(col_ids, vsoc(root_idx)); push!(blocks, Ap)
        else
            Ap = zeros(1, 3); Ap[1, 1] = 1.0
            push!(row_ids, e); push!(col_ids, vcomb(root_idx)); push!(blocks, Ap)
        end
        rhs_val[e] = [sys.eps]
        K_cones = vcat(IPM.AbstractCone[CofreeCone() for _ in 1:M],
                       IPM.AbstractCone[SecondOrderCone() for _ in 1:M],
                       IPM.AbstractCone[SecondOrderCone() for _ in 1:ncomb])
        stalks = vcat([(Lf, :free) for _ in 1:M],
                      [(dsoc, :soc) for _ in 1:M],
                      [(3, :soc) for _ in 1:ncomb])
    else
        error("unknown form $form")
    end

    B = blocksparse(row_ids, col_ids, blocks)
    g_rhs = zeros(size(B, 1))
    for (edge, val) in rhs_val
        g_rhs[rowrange(B, edge)] .= val
    end

    # ---- dense noise Grams on the filter stalks — the dense-Q habitat
    Q = IPM.allocblockdiag(B); fill!(Q, 0)
    trΘ = tr(sys.Θ) / Lf
    for m in 1:M
        Qm = white ?
            (sys.sigmas[m]^2 * trΘ + inst.lam) .* Matrix{Float64}(I, Lf, Lf) :
            sys.Qms[m]
        block(Q, m, m, m) .= Qm
    end
    c = zeros(size(B, 2))

    prob = IPMProblem(Q, B, c, g_rhs, K_cones)
    ctx = (; M, Lf, form, white, stalks, sys, inst)
    return prob, ctx
end

"Equalizer bank (M × Lf rows) from the solution vector."
function extract_f(prob, ctx, p::AbstractVector)
    F = zeros(ctx.M, ctx.Lf)
    for m in 1:ctx.M
        F[m, :] .= p[colrange(prob.B, m)]
    end
    return F
end

"Direct evaluation of the objective (assembly check)."
equalizer_objective(ctx, F) =
    0.5 * sum(dot(F[m, :], ctx.sys.Qms[m], F[m, :]) for m in 1:ctx.M)

"Composite-response residual ‖T f − g‖ (the quantity ε certifies)."
isi_residual(ctx, F) = norm(ctx.sys.T * vec(F') - ctx.sys.g)

"Whitened-ball utilization ‖L'f − b‖ / εx."
ball_use(ctx, F) = norm(ctx.sys.LT * vec(F') - ctx.sys.b) / ctx.inst.epsx

# -----------------------------------------------------------------------------
# Native TRS comparator (joint ball; no conic solver) — Gay 1981, MS 1983
# -----------------------------------------------------------------------------

"""Joint problem min ½f'Q̄f s.t. ‖Tf − g‖ ≤ hypot(ε_LS, εx) by the secular
equation in whitened coordinates y = L'f − b: y(ν) = −(H + νI)⁻¹Hb with
H = L⁻¹Q̄L⁻ᵀ, ν solving ‖y(ν)‖ = εx. The suite's TRS oracle, in-process."""
function trs_joint(sys, inst::EqualizerInstance)
    M, Lf = sys.M, sys.Lf
    Qb = zeros(M * Lf, M * Lf)
    for m in 1:M
        r = (m - 1) * Lf + 1:m * Lf
        Qb[r, r] .= sys.Qms[m]
    end
    H = Symmetric((sys.L \ Qb) / sys.LT)
    F = eigen(H)
    β = F.vectors' * (H * sys.b)
    @assert norm(sys.b) > inst.epsx "constraint inactive: f = 0 optimal"
    ynorm(ν) = sqrt(sum(abs2, β ./ (F.values .+ ν)))
    lo, hi = 0.0, 1.0
    while ynorm(hi) > inst.epsx
        hi *= 4
    end
    for _ in 1:200
        mid = 0.5 * (lo + hi)
        if ynorm(mid) > inst.epsx
            lo = mid
        else
            hi = mid
        end
    end
    ν = 0.5 * (lo + hi)
    y = F.vectors * (-β ./ (F.values .+ ν))
    f = sys.LT \ (y + sys.b)
    ms = (; ν, boundary = abs(norm(y) - inst.epsx),
          stat = norm(H * (y + sys.b) + ν * y))
    return f, 0.5 * dot(f, Qb, f), ms
end

# -----------------------------------------------------------------------------
# Baseline (same conic program via JuMP)
# -----------------------------------------------------------------------------

function build_jump_equalizer(prob, ctx, optimizer)
    model = Model(optimizer)
    set_silent(model)
    xs = Vector{Vector{VariableRef}}(undef, length(ctx.stalks))
    for (v, (dim, kind)) in enumerate(ctx.stalks)
        x = @variable(model, [1:dim])
        kind == :soc && @constraint(model, x in JuMP.SecondOrderCone())
        xs[v] = x
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

function tridensity(U; tol = 1e-8)
    n = size(U, 1)
    v = [abs(U[i, j]) for j in 1:n for i in 1:j]
    count(v .> tol * maximum(v)) / length(v)
end

function test_dense_by_nature(sys)
    dΘ = density(sys.Θ)
    dG = minimum(density(Matrix(Gm)) for Gm in sys.Gms)
    dL = minimum(tridensity(LTm) for LTm in sys.LTms)
    dJ = density(Matrix(sys.T' * sys.T))
    @assert min(dΘ, dG, dL, dJ) > 0.95 "dense-by-nature broken: $((dΘ, dG, dL, dJ))"
    println("  [PASS] dense by nature: density@1e-8 Θ $(round(dΘ, digits = 3)) ",
        "G_m ≥ $(round(dG, digits = 3)) L_m' ≥ $(round(dL, digits = 3)) ",
        "joint G $(round(dJ, digits = 3))")
end

"The documented boundary: decimation D > 1 makes the Gram polyphase-block-
diagonal with EXACT structural zeros — out of habitat by structure."
function test_polyphase_exclusion(; D = 4, Lh = 64, Lf = 32, seed = 5)
    rng = MersenneTwister(seed)
    h = randn(rng, Lh)
    P = Lh + Lf - 1
    Td = zeros(D * P, Lf)
    for r in 0:D - 1, τ in 0:P - 1, j in 0:Lf - 1
        u = τ - j
        if 0 <= u < Lh && mod(u + r, D) == 0
            Td[r * P + τ + 1, j + 1] = h[u + 1]
        end
    end
    Gd = Td' * Td
    offmax = 0.0
    for j in 1:Lf, jp in 1:Lf
        mod(j - jp, D) != 0 && (offmax = max(offmax, abs(Gd[j, jp])))
    end
    @assert offmax == 0.0 "polyphase exclusion broken: $offmax"
    println("  [PASS] polyphase exclusion (D = $D control): off-class Gram ",
        "entries EXACTLY 0 — structural 1/D density, whitening cannot ",
        "remove it (the decimated bank is a boundary-suite control)")
end

function test_mint_threshold(; M = 4, Lh = 121, trev = 40.0, seed = 3)
    rng = MersenneTwister(seed)
    env = exp.(-(0:Lh - 1) ./ trev)
    ch = [randn(rng, Lh) .* env for _ in 1:M]
    Lstar = ceil(Int, (Lh - 1) / (M - 1))
    els(Lf) = begin
        Tl = hcat([toepblock(h, Lf) for h in ch]...)
        gl = zeros(Lh + Lf - 1); gl[(Lh + Lf) ÷ 2 + 1] = 1.0
        norm(gl - Tl * (Tl \ gl))
    end
    e0, e1 = els(Lstar - 1), els(Lstar)
    @assert e0 > 1e-3 && e1 < 1e-9 "MINT threshold broken: $e0, $e1 (oracle: 2.2e-2, 2.6e-14)"
    println("  [PASS] MINT length threshold: ε_LS $(round(e0, sigdigits = 2)) ",
        "at Lf* − 1 → $(round(e1, sigdigits = 2)) at Lf* = $Lstar ",
        "(Miyoshi–Kaneda)")
end

function test_whitening_identity(sys)
    rng = MersenneTwister(11)
    f = randn(rng, sys.M * sys.Lf)
    lhs = norm(sys.T * f - sys.g)^2
    ej = abs(lhs - (norm(sys.LT * f - sys.b)^2 + sys.eLS^2)) / lhs
    m = max(1, sys.M ÷ 2)
    fm = randn(rng, sys.Lf)
    gm = sys.Ts[m] * sys.fhm[m] + sys.rhat ./ sys.M
    lhs2 = norm(sys.Ts[m] * fm - gm)^2
    es = abs(lhs2 - (norm(sys.LTms[m] * fm - sys.bms[m])^2 +
                     (sys.eLS / sys.M)^2)) / lhs2
    @assert max(ej, es) < 1e-10 "whitening identity broken: $ej, $es"
    println("  [PASS] whitening identities at random points: joint ",
        "$(round(ej, sigdigits = 2)), split $(round(es, sigdigits = 2)) ",
        "(ε_LS = $(round(sys.eLS, digits = 4)), ε = $(round(sys.eps, digits = 4)))")
end

function test_objective_identity(prob, ctx, res)
    F = extract_f(prob, ctx, res.p)
    o_direct = equalizer_objective(ctx, F)
    o_ipm = ipm_objective(prob, res)
    rel = abs(o_direct - o_ipm) / abs(o_direct)
    @assert rel < 1e-5 "assembly mismatch: $o_direct vs $o_ipm"
    println("  [PASS] objective identity: F(f) = $(round(o_direct, digits = 5)) ",
        "(rel $(round(rel, sigdigits = 2)))")
    return F
end

"The analytic gate: IPM (exact form) == native TRS secular solve."
function test_trs(sys, inst, F_ipm, v_ipm)
    f_trs, v_trs, ms = trs_joint(sys, inst)
    @assert ms.boundary < 1e-10 && ms.stat < 1e-6 "MS conditions: $ms"
    dv = abs(v_ipm - v_trs) / abs(v_trs)
    df = maximum(abs.(vec(F_ipm') - f_trs))
    @assert dv < 1e-5 "IPM vs TRS: rel Δv $dv (oracle: 4.2e-9 conic-vs-TRS)"
    println("  [PASS] TRS classical certificate (Gay 1981; Moré–Sorensen 1983): ",
        "IPM == secular solve, rel Δv $(round(dv, sigdigits = 2)), ",
        "‖Δf‖∞ $(round(df, sigdigits = 2)), ν = $(round(ms.ν, digits = 2))")
end

function test_certificate_and_ball(ctx, F)
    isi = isi_residual(ctx, F)
    use = ball_use(ctx, F)
    @assert isi <= ctx.sys.eps * (1 + 1e-5) "certificate broken: $isi > $(ctx.sys.eps)"
    lo = ctx.form == :exact ? 0.99 : 0.15
    hi = ctx.form == :exact ? 1.01 : 0.6
    @assert lo < use < hi "ball use $use outside $(ctx.form) range (oracle: exact 1.000, budget 0.329)"
    println("  [PASS] certificate: ‖Tf − g‖ = $(round(isi, digits = 4)) ≤ ε = ",
        "$(round(ctx.sys.eps, digits = 4)); ball use $(round(use, digits = 3)) ",
        "($(ctx.form) form; oracle: exact 1.000 vs budget 0.329)")
end

function test_conservatism(sys, inst, settings)
    ratios = Float64[]
    for εx in (inst.epsx, inst.epsx / 4)
        inst2 = EqualizerInstance(inst.Lf, inst.lhfac, inst.trevfrac, inst.r,
            inst.phi, inst.lam, εx, inst.sglo, inst.sghi, inst.seed)
        sys2 = (; sys..., eps = hypot(sys.eLS, εx))
        prob, ctx = build_equalizer(sys2, inst2; form = :budget, white = false)
        res = solve(prob, settings)
        @assert res.status in (OPTIMAL, NEAR_OPTIMAL) "εx=$εx: $(res.status)"
        _, v_j, _ = trs_joint(sys2, inst2)
        push!(ratios, ipm_objective(prob, res) / v_j)
    end
    @assert ratios[2] < ratios[1] && ratios[2] < 1.30 "conservatism curve broken: $ratios (oracle: 2.04 → 1.18)"
    println("  [PASS] budget-split conservatism: v_budget/v_joint ",
        "$(round(ratios[1], digits = 3)) → $(round(ratios[2], digits = 3)) ",
        "over εx = $(inst.epsx) → $(inst.epsx / 4) (oracle: 2.04 → 1.18, ",
        "→ 1 at the pinv endpoint)")
end

function test_noise_shaping(sys, inst, settings, F_col)
    prob, ctx = build_equalizer(sys, inst; form = :exact, white = true)
    res = solve(prob, settings)
    @assert res.status in (OPTIMAL, NEAR_OPTIMAL) "white control: $(res.status)"
    F_w = extract_f(prob, ctx, res.p)
    move = norm(F_w - F_col) / norm(F_col)
    excess = equalizer_objective(ctx, F_w) / equalizer_objective(ctx, F_col)
    @assert move > 0.05 && excess > 1.2 "noise shaping not load-bearing: $move, $excess (oracle: 0.121, 1.77)"
    println("  [PASS] noise shaping load-bearing: diagonal-Q (white) control ",
        "moves the solution by $(round(100move, digits = 1))% and pays ",
        "×$(round(excess, digits = 2)) colored-noise power (oracle: 12.1%, ",
        "×1.77) — the D6/smooth-bary ablation, run it with --white")
end

function test_ipm_vs_clarabel(prob, ctx, F_ipm)
    model, xv = build_jump_equalizer(prob, ctx, clarabel_opt(; tol = TOL))
    optimize!(model)
    F_cla = extract_f(prob, ctx, value.(xv))
    df = maximum(abs.(F_ipm .- F_cla))
    @assert df < 1e-3 "IPM vs Clarabel mismatch: $df"
    println("  [PASS] IPM vs Clarabel (same conic program): ‖Δf‖∞ = ",
        "$(round(df, sigdigits = 2))")
end

# -----------------------------------------------------------------------------
# Sweep
# -----------------------------------------------------------------------------

msizes() = OPTS.tiny ? [4] : OPTS.quick ? [4, 8, 16] : [4, 8, 16, 32]
lsizes() = OPTS.tiny ? [32] : OPTS.quick ? [32, 64, 128] : [32, 64, 128, 256]

function run()
    form = TREE ? :budget_tree : (BUDGET ? :budget : :exact)
    println("\n", "=" ^ 78)
    println("  E11: multichannel equalizer bank, regularized MINT ",
        LDIAL ? "(dial: block size Lf; M = 8)" : "(dial: channels M; Lf = 64)",
        " [$(form) form$(WHITE ? ", WHITE ablation" : "")]")
    println("=" ^ 78)

    ipm_settings = IPMSettings{Float64}(
        kkt = UzawaSettings{Float64}(raug = 1e4, elim = CliqueTrees.METIS()),
        feas_tol = TOL, gap_tol = TOL, itmax = 200)
    hsd_settings = HSDSettings{Float64}(
        kkt = UzawaSettings{Float64}(raug = 1e4, elim = CliqueTrees.METIS()),
        feas_tol = TOL, gap_tol = TOL, itmax = 200)

    inst = equalizer_instance()
    sys = make_system(inst, 8)

    println("\n  Gate tests (M = 8, Lh = 1024, Lf = 64, ε_LS = ",
        "$(round(sys.eLS, digits = 4))):")
    test_dense_by_nature(sys)
    test_polyphase_exclusion()
    test_mint_threshold()
    test_whitening_identity(sys)
    prob, ctx = build_equalizer(sys, inst; form = :exact, white = false)
    res = solve(prob, ipm_settings)
    @assert res.status in (OPTIMAL, NEAR_OPTIMAL) "IPM status $(res.status)"
    F = test_objective_identity(prob, ctx, res)
    test_trs(sys, inst, F, ipm_objective(prob, res))
    test_certificate_and_ball(ctx, F)
    test_conservatism(sys, inst, ipm_settings)
    test_noise_shaping(sys, inst, ipm_settings, F)
    test_ipm_vs_clarabel(prob, ctx, F)
    println()

    cla_opt = clarabel_opt(; tol = TOL)
    msk_opt = OPTS.mosek ? mosek_opt(; tol = TOL) : nothing

    rows = []
    for sz in (LDIAL ? lsizes() : msizes())
        M = LDIAL ? 8 : sz
        instk = LDIAL ? equalizer_instance(; Lf = sz) : inst
        sysk = make_system(instk, M)
        probk, ctxk = build_equalizer(sysk, instk; form, white = WHITE)
        stats = problem_stats(probk)

        @printf("  %s=%-4d dof=%-6d n1=%-6d blk=%-4.0f  ",
            LDIAL ? "Lf" : "M", sz, stats.N0, stats.N1, stats.med_block)

        m_ipm = measure_ipm(probk, ipm_settings; nruns = NRUNS)
        m_hsd = measure_ipm(probk, hsd_settings; nruns = NRUNS)
        m_cla = measure_jump(() -> first(build_jump_equalizer(probk, ctxk, cla_opt)); nruns = NRUNS)
        m_msk = msk_opt !== nothing ?
            measure_jump(() -> first(build_jump_equalizer(probk, ctxk, msk_opt)); nruns = NRUNS) :
            (t = NaN, status = "", obj = NaN)

        ratio(b, base) = isfinite(b.t) && isfinite(base.t) ? b.t / base.t : NaN
        fmt_ratio(b, base) = isnan(ratio(b, base)) ? "—" : @sprintf("%.2fx", ratio(b, base))

        println("IPM ", fmt_time(m_ipm), "  HSD ", fmt_time(m_hsd), " (", fmt_ratio(m_hsd, m_ipm),
            ")  Cla ", fmt_time(m_cla), " (", fmt_ratio(m_cla, m_ipm),
            ")  Msk ", fmt_time(m_msk), " (", fmt_ratio(m_msk, m_ipm), ")")

        push!(rows, (size = sz, dof = stats.N0, ipm = m_ipm, hsd = m_hsd, cla = m_cla, msk = m_msk))
    end

    dofs = [r.dof for r in rows]
    println("\nFitted log-log slopes (time vs DOF):")
    for (name, get_t) in [("IPM", r -> r.ipm.t), ("HSD", r -> r.hsd.t), ("Clarabel", r -> r.cla.t), ("Mosek", r -> r.msk.t)]
        sl = loglog_slope(dofs, [get_t(r) for r in rows])
        print("  $name: ", isnan(sl) ? "n/a" : @sprintf("DOF^%.2f", sl))
    end
    println()
end

