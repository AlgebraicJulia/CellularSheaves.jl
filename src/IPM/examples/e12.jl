# =============================================================================
# e12.jl — Distributed sound-zone control: cooperative multizone
#              rendering with local coupling (E11 tiled over a zone graph —
#              dense-Q + SOC at bounded degree, coupling mass O(Z))
#
# Usage:
#   julia --project e12.jl              # Clarabel baseline
#   julia --project e12.jl --mosek      # + Mosek
#   julia --project e12.jl --quick      # quick (smaller) sweep
#   julia --project e12.jl --joint      # per-ball joint SOCs (ablation)
#   julia --project e12.jl --nobudget   # drop amplifier caps (ablation)
#   julia --project e12.jl --ldial      # block-size dial Lf at fixed Z
#
# Problem. Z zones on a path, each with an array of S loudspeakers and C
# control points. Acoustic paths are reverberant FIRs (length Lh, decay τ)
# whose inter-zone amplitude falls off as γ^{|Δz|} with propagation onset
# d0·|Δz|. Every array within radius ρ of a program COOPERATES in rendering
# it (pressure matching drives all nearby speakers per program — Betlehem–
# Zhang–Poulsen–Abhayapala, IEEE SPM 2015; Coleman et al.; car-cabin
# personal audio: Cheer & Elliott), so the stalks are f_{z',p} ∈ R^{S·Lf},
# one per (array z', program p) with |z' − p| ≤ ρ. Per (zone z, program p),
# |p − z| ≤ ρ, the composite response SUMS across arrays,
#
#     r_{z,p} = Σ_{z' ∈ W(p)} A_{z,z'} f_{z',p} − g_{z,p},
#
# with REFERENCE-FIELD bright targets g_{z,z} = A_{z,z} f_ref (virtual-
# source pressure matching; delta targets would demand full dereverberation,
# push the floors to ~0.9 of the signal, and admit f = 0) and dark targets
# 0, constrained per ball ‖r_{z,p}‖ ≤ ε_{z,p} with SIGNAL-SCALE radii
#
#     ε_{z,p} = hypot(floor_{z,p}, η·s_z),   s_z = ‖g_{z,z}‖,
#
# floor_{z,p} = ‖r_{z,p}(f̂)‖ at the global stacked LS point f̂. Objective:
# ½ Σ f'Q_p f with Q_p the RADIATED-POWER Gram — mutual radiation
# resistance sinc(ω·δ·|Δs|) (Morse–Ingard) integrated against program p's
# AR(2) PSD (heterogeneous peak frequencies, the analog of E11's sensor
# scales) — plus the driver floor λI. PER-ARRAY ACOUSTIC-OUTPUT BUDGETS
# (SPL / noise-ordinance caps) couple the programs at each array:
#
#     Σ_{p ∈ P_a} f_{a,p}' Q_p f_{a,p} ≤ β_a²,  β_a = θ·√(power_a(f̂)),
#
# feasible because the effort-minimal solution runs at 0.70–0.79 of f̂'s
# per-array powers (measured; utilization 0.83–0.94 at θ = 0.85).
#
# THREE RECORDED STRUCTURAL OBSERVATIONS (design-process negatives, in the
# polyphase tradition; all measured by the oracle):
#  (i) ONE ARRAY PER PROGRAM makes the per-zone gluing EXACTLY block-
#      diagonal across neighborhood stalks — Cholesky tri-density
#      1/(2ρ+1) = 0.338 measured, structural, whitening-proof. Cooperative
#      rendering restores density 1.000.
# (ii) WITHOUT the budgets the problem decomposes EXACTLY by program (balls
#      never mix programs; the objective is stalk-separable): Z connected
#      components — fatal for a benchmark, and machine-checked by the
#      connectivity certificate (1 component with budgets). Long-range
#      influence is NOT the cure and should not be engineered: the physics
#      is SCREENING — influence decays geometrically at a γ-tracked rate
#      (hop ratio 0.21 at γ = 0.25 vs 0.31 at γ = 0.4, the E09-DARE
#      species). Non-decomposability and reach are different properties.
# (iii) A SIGNAL-POWER budget (I_S ⊗ Toeplitz) would have exactly block-
#      diagonal Cholesky ties, tri-density 1/S — avoided by capping
#      acoustic rather than electrical power.
#
# LOCAL WHITENING + EXACT SLICE SPLIT (E11's machinery as the LOCAL
# structure at every ball). Per reproduction ball, G = M'M = LL' (tall,
# PD), ‖r‖² = ‖L'φ − b‖² + c², rows partition by window-stalk slices:
# k = |W(p)| slice SOCs (s_j, w_j) of dim S·Lf + 1 with dense triangular
# ties, a (k+2)-dim hub (h, ŝ, σ) with pins h = ε, σ = c, and s-links.
# Budget balls are program-partitioned by construction: |P_a| slices
# w_{a,p} = R̂_p'f_{a,p} (dense Cholesky-of-radiation-Gram ties) + a
# (|P_a|+1)-dim hub pinned at β_a. Zero conservatism; bounded arity —
# coupling mass O(Z·ρ²·(S·Lf)²), the pre-registered answer to E11's M-dial
# reading (slope DOF^2.29 from O(M²) global coupling).
#
# ENDPOINT. Radii at the global floors: any point of the intersection is a
# global LS minimizer and the stacked map has full column rank, so η → 0
# pins f̂ EXACTLY (budgets excluded: f̂ sits above the θ < 1 caps) —
# measured rate η^0.98.
#
# Habitat claim: dense-Q + SOC (E09/E11's cell) at med block ≈ S·Lf + 1 =
# 65, bounded-degree coupling on a path of zones. New edge species: window
# Cholesky ties (arity ≤ 2ρ+2) and radiation-Gram budget ties; hub SOCs of
# dim ≤ 2ρ+3. Q ≡ 0 on cone stalks, no free auxiliary stalks, bloat ~ 1.
#
# Validation (tools/soundzone_oracle.py, self-contained numpy + scipy +
# cvxpy, fully executed; figures at Z = 6, S = 4, C = 3, Lf = 16, Lh = 96,
# ρ = 1, γ = 0.25, τ = 32, λ = 1e-3, η = 0.1, θ = 0.85; bright floors
# 0.058–0.108 of signal scale, dark 0.141–0.189):
#   * DENSE BY NATURE + EXCLUSION (i): densities@1e-8 Q_p 1.000, per-ball
#     G 1.000, reproduction ties 1.000, budget ties 1.000; the one-array
#     control has cross-stalk Gram blocks EXACTLY 0, tri-density 0.338.
#   * ENDPOINT: ‖f(η) − f̂‖∞ 8.75e-3 → 2.25e-3 over η 0.02 → 0.005 (rate
#     η^0.98; stacked min eig 7.13).
#   * FORMULATION CERTIFICATE: whitening identities 1.0e-15 (reproduction)
#     / 0.0 (budget partition); the stalk assembly (exactly what this file
#     builds) == native to ‖Δf‖∞ = 1.8e-7, rel Δv = 7.3e-11.
#   * TRS (classical certificate): Z = 1, ρ = 0, no-budget reduces to a
#     single-ball trust-region subproblem — secular == conic to rel Δv
#     4.0e-10; Moré–Sorensen boundary 5.6e-17, stationarity 2.9e-17.
#   * EXACTNESS: slice-split == joint per-ball balls, budgets included
#     (zero conservatism), same figures as the assembly gate.
#   * STRUCTURE: connectivity 6 components without budgets vs 1 with;
#     γ = 0 separates per ARRAY, stacked == joint to 1.2e-9; the non-
#     cooperative baseline violates the dark balls by ×1.8.
#   * SCREENING: hop-1 influence ratio 0.21 (γ = 0.25) vs 0.31 (γ = 0.4).
#   * TRUNCATION REACH: uncontrolled leakage γ^1.79 at distance ρ+1 and
#     γ^2.80 at ρ+2 — increment 1.01 (one power of γ per zone); the base
#     exceeds the naive γ^1 because the dark balls SHADOW farther leakage.
#   * DELIVERABLE: per-speaker signal power matches the spectral-Gram
#     prediction (rel 0.002); response-domain contrast improves by +2.3 to
#     +4.3 dB per zone (mean +3.3) over non-cooperative; utilization
#     0.83–0.94.
#   * CROSS-SOLVER: Clarabel vs SCS ‖Δf‖∞ = 4.2e-8.
# The gate tests below re-derive the fast subset in-process (including the
# TRS comparator natively); screening, truncation, and the simulation
# deliverables are oracle-side.
#
# Dial: zone count Z at fixed S = 4, C = 3, Lf = 16, Lh = 96 (vertex count
# grows, per-block work constant — the chain discipline of E06/E08/E09).
# --ldial sweeps Lf ∈ {16, 32, 48} at Z = 8 with Lh = 6·Lf (block dial).
# γ, ρ, η, θ are physics knobs.
#
# STATUS: authored against the e09/e11 build API; the oracle is fully
# executed but this file has NOT yet been run against the IPM. First-run
# checklist: (1) confirm raug (1e4 per the dense-Q + SOC precedent);
# (2) THE PRE-REGISTERED CLAIM: the Z-dial slope should land near the
# chain examples (~DOF^1.1–1.3), against E11's DOF^2.29 — bounded-degree
# coupling is the whole point; predicted win vs Clarabel and Mosek at all
# Z (E09/E11's cell at blk 65); (3) --joint ablation: per-ball SOCs of dim
# k·S·Lf + 2 ≈ 194 vs slices at 65 — the block-size ablation, E11's
# monolithic-comparator move one level down; (4) --nobudget runs Z
# independent programs (the decomposability control — expect baselines to
# gain relatively if their presolves exploit it, and the IPM to be
# indifferent); (5) note the build performs one dense global LS solve
# (nstk·S·Lf ≈ 6k at Z = 32, a few seconds, untimed); (6) fill the
# sample-run table in run.jl.
# =============================================================================

using AppleAccelerate
using LinearAlgebra, SparseArrays, Printf, Random

include("utils.jl")

using CellularSheaves.IPM: SecondOrderCone, CofreeCone
using CellularSheaves.BlockSparseArrays: rowrange

const OPTS = parse_args(ARGS)
const TOL = OPTS.tol
const NRUNS = 5
const JOINT = "--joint" in ARGS
const NOBUDGET = "--nobudget" in ARGS
const LDIAL = "--ldial" in ARGS

# -----------------------------------------------------------------------------
# Problem definition
# -----------------------------------------------------------------------------

struct SoundZoneInstance
    S::Int; C::Int; Lf::Int
    lhfac::Int                 # Lh = lhfac * Lf
    rho::Int                   # cooperation / constraint radius
    gam::Float64               # inter-zone gain (coupling-strength dial)
    taufrac::Float64           # reverb decay = taufrac * Lh
    d0::Int                    # propagation delay per zone (samples)
    delta::Float64             # inter-driver travel time (samples)
    lam::Float64               # driver-effort floor
    eta::Float64               # reproduction slack (signal-scale)
    theta::Float64             # amplifier-cap fraction of f̂ power
    dr::Int                    # reference-driving delay
    r_ar::Float64              # program AR(2) pole radius
    seed::Int
end

soundzone_instance(; S = 4, C = 3, Lf = 16, lhfac = 6, rho = 1, gam = 0.25,
                   taufrac = 1 / 3, d0 = 4, delta = 2.0, lam = 1e-3,
                   eta = 0.1, theta = 0.85, dr = 4, r_ar = 0.9, seed = 7) =
    SoundZoneInstance(S, C, Lf, lhfac, rho, gam, taufrac, d0, delta, lam,
                      eta, theta, dr, r_ar, seed)

"AR(2) coefficients of program p (heterogeneous peak frequencies)."
ar_coeffs(inst, phis, p) =
    [1.0, -2 * inst.r_ar * cos(phis[p]), inst.r_ar^2]

"Radiated-power Gram: sinc radiation coupling × program-p PSD (quadrature)."
function effort_gram(inst::SoundZoneInstance, phis, p)
    S, Lf = inst.S, inst.Lf
    nf = S * Lf
    nw = 2048
    w = range(1e-4, π; length = nw)
    a = ar_coeffs(inst, phis, p)
    Sw = [1.0 / abs2(a[1] + a[2] * cis(-ω) + a[3] * cis(-2ω)) for ω in w]
    Sw ./= maximum(Sw)
    dw = step(w)
    trap = [i == 1 || i == nw ? 0.5 : 1.0 for i in 1:nw]
    Q = zeros(nf, nf)
    for s in 0:S - 1, sp in s:S - 1
        kern = sinc.(w .* (inst.delta * (sp - s)) ./ π)
        for lag in 0:Lf - 1
            v = sum(trap .* Sw .* kern .* cos.(w .* lag)) * dw / π
            for i in 0:Lf - 1 - lag
                Q[s * Lf + i + 1, sp * Lf + i + lag + 1] = v
                Q[sp * Lf + i + lag + 1, s * Lf + i + 1] = v
                Q[s * Lf + i + lag + 1, sp * Lf + i + 1] = v
                Q[sp * Lf + i + 1, s * Lf + i + lag + 1] = v
            end
        end
    end
    return Q + inst.lam * Matrix{Float64}(I, nf, nf)
end

"Tall convolution block Toep(h) ∈ R^{P × Lf}."
function toepblock(h, Lf)
    Lh = length(h)
    T = zeros(Lh + Lf - 1, Lf)
    for j in 1:Lf
        T[j:j + Lh - 1, j] .= h
    end
    return T
end

"""Channels, targets, whitening, budgets. Every Gram is dense by nature:
Q_p (radiation × colored spectrum), per-ball G (reverberant window cross-
correlations), and both tie families."""
function make_system(inst::SoundZoneInstance, Z::Int; gam = inst.gam)
    S, C, Lf, rho = inst.S, inst.C, inst.Lf, inst.rho
    nf = S * Lf
    Lh = inst.lhfac * Lf
    P = Lh + Lf - 1
    tau = inst.taufrac * Lh
    rng = MersenneTwister(inst.seed)
    raw = Dict{Tuple{Int, Int}, Array{Float64, 3}}()   # (z,zp) -> C×S×Lh
    for z in 0:Z - 1, zp in 0:Z - 1
        H = zeros(C, S, Lh)
        for c in 1:C, s in 1:S
            on = min(inst.d0 * abs(z - zp), Lh - 8)
            n = Lh - on
            H[c, s, on + 1:Lh] .= randn(rng, n) .* exp.(-(0:n - 1) ./ tau)
        end
        raw[z, zp] = H
    end
    rirf(z, zp, c, s) = gam^abs(z - zp) .* raw[z, zp][c, s, :]
    Ablock(z, zp) = vcat([hcat([toepblock(rirf(z, zp, c, s), Lf)
                                for s in 1:S]...) for c in 1:C]...)
    phis = Z == 1 ? [0.5π] : collect(range(0.3π, 0.7π; length = Z))
    Qs = [effort_gram(inst, phis, p) for p in 1:Z]
    RhatT = [Matrix(cholesky(Symmetric(Q)).U) for Q in Qs]   # R̂'
    window(p) = [z for z in 0:Z - 1 if abs(z - p) <= rho]
    programs_at(a) = [p for p in 0:Z - 1 if a in window(p)]
    stalks = [(zp, p) for p in 0:Z - 1 for zp in window(p)]
    sidx = Dict(sk => i for (i, sk) in enumerate(stalks))
    balls = [(z, p) for p in 0:Z - 1 for z in 0:Z - 1 if abs(p - z) <= rho]
    f_ref = zeros(nf)
    for s in 0:S - 1
        f_ref[s * Lf + inst.dr + 1] = 1.0 / S
    end
    Mb = Dict{Tuple{Int, Int}, Matrix{Float64}}()
    gb = Dict{Tuple{Int, Int}, Vector{Float64}}()
    colsb = Dict{Tuple{Int, Int}, Vector{Tuple{Int, Int}}}()
    for (z, p) in balls
        cols = [(zp, p) for zp in window(p)]
        Mb[z, p] = hcat([Ablock(z, zp) for (zp, _) in cols]...)
        gb[z, p] = p == z ? Ablock(z, z) * f_ref : zeros(C * P)
        colsb[z, p] = cols
    end
    scale = Dict(z => norm(gb[z, z]) for z in 0:Z - 1)
    # global stacked LS point via blockwise normal equations
    nstk = length(stalks)
    Gg = zeros(nstk * nf, nstk * nf)
    hg = zeros(nstk * nf)
    for ball in balls
        M, g, cols = Mb[ball], gb[ball], colsb[ball]
        for (i, ci) in enumerate(cols)
            ii = sidx[ci]
            Mi = view(M, :, (i - 1) * nf + 1:i * nf)
            hg[(ii - 1) * nf + 1:ii * nf] .+= Mi' * g
            for (j, cj) in enumerate(cols)
                jj = sidx[cj]
                Gg[(ii - 1) * nf + 1:ii * nf, (jj - 1) * nf + 1:jj * nf] .+=
                    Mi' * view(M, :, (j - 1) * nf + 1:j * nf)
            end
        end
    end
    fhat = cholesky(Symmetric(Gg)) \ hg
    gather(f, ball) = vcat([f[(sidx[c] - 1) * nf + 1:sidx[c] * nf]
                            for c in colsb[ball]]...)
    floors = Dict(ball => norm(Mb[ball] * gather(fhat, ball) - gb[ball])
                  for ball in balls)
    # per-ball whitening
    LbT = Dict{Tuple{Int, Int}, Matrix{Float64}}()   # L'
    bb = Dict{Tuple{Int, Int}, Vector{Float64}}()
    cb = Dict{Tuple{Int, Int}, Float64}()
    for ball in balls
        M, g = Mb[ball], gb[ball]
        ch = cholesky(Symmetric(M' * M))
        b = ch.L \ (M' * g)
        LbT[ball] = Matrix(ch.U)
        bb[ball] = b
        cb[ball] = sqrt(max(dot(g, g) - dot(b, b), 0.0))
    end
    # acoustic-output budgets from f̂'s per-array radiated powers
    powers(f) = begin
        pw = zeros(Z)
        for (i, (zp, p)) in enumerate(stalks)
            x = f[(i - 1) * nf + 1:i * nf]
            pw[zp + 1] += dot(x, Qs[p + 1], x)
        end
        pw
    end
    betas = inst.theta .* sqrt.(powers(fhat))
    return (; Z, S, C, Lf, Lh, P, nf, gam, phis, Qs, RhatT, window,
            programs_at, stalks, sidx, balls, Mb, gb, colsb, scale, fhat,
            floors, LbT, bb, cb, betas, gather, powers, f_ref)
end

eps_ball(sys, inst, ball; eta = inst.eta) =
    hypot(sys.floors[ball], eta * sys.scale[ball[1]])

# -----------------------------------------------------------------------------
# Build
# -----------------------------------------------------------------------------
#
# Stalk order: filter stalks (array, program) | per reproduction ball:
# k = |W(p)| slice SOCs (s_j, w_j), dim nf+1, then hub (h, ŝ_1..k, σ), dim
# k+2 | per array: |P_a| budget slices (u_p, v_p), dim nf+1, then hub
# (h, û), dim |P_a|+1. Edges: slice ties (dense triangular L' blocks,
# nonzero rhs — E06's pattern), budget ties (dense R̂' blocks), s-links,
# pins (h = ε, σ = c; h = β). --joint replaces each ball's slices + hub by
# ONE SOC stalk (s, w, σ) of dim k·nf + 2 (the monolithic ablation).

function build_soundzone(sys, inst::SoundZoneInstance;
                         joint = JOINT, budget = !NOBUDGET, nocoop = false,
                         eta = inst.eta)
    nf, nstk = sys.nf, length(sys.stalks)
    used_balls = nocoop ? [b for b in sys.balls if b[1] == b[2]] : sys.balls
    use_budget = budget && !nocoop
    row_ids, col_ids, blocks = Int[], Int[], Matrix{Float64}[]
    rhs_val = Dict{Int, Vector{Float64}}()
    Inf_ = Matrix{Float64}(I, nf, nf)
    stalkdims = Tuple{Int, Symbol}[(nf, :free) for _ in 1:nstk]
    v = nstk
    rslice = Dict{Tuple{Tuple{Int, Int}, Int}, Int}()
    rhub = Dict{Tuple{Int, Int}, Int}()
    rjoint = Dict{Tuple{Int, Int}, Int}()
    for ball in used_balls
        k = length(sys.colsb[ball])
        if joint
            v += 1; rjoint[ball] = v
            push!(stalkdims, (k * nf + 2, :soc))
        else
            for j in 1:k
                v += 1; rslice[ball, j] = v
                push!(stalkdims, (nf + 1, :soc))
            end
            v += 1; rhub[ball] = v
            push!(stalkdims, (k + 2, :soc))
        end
    end
    bslice = Dict{Tuple{Int, Int}, Int}()
    bhub = Dict{Int, Int}()
    if use_budget
        for a in 0:sys.Z - 1
            Pa = sys.programs_at(a)
            for p in Pa
                v += 1; bslice[a, p] = v
                push!(stalkdims, (nf + 1, :soc))
            end
            v += 1; bhub[a] = v
            push!(stalkdims, (length(Pa) + 1, :soc))
        end
    end
    e = 0
    for ball in used_balls
        k = length(sys.colsb[ball])
        LT = sys.LbT[ball]
        if joint
            d = k * nf + 2
            e += 1                             # ties w = L'φ − b (all rows)
            Aw = zeros(k * nf, d); Aw[:, 2:k * nf + 1] .= I(k * nf)
            push!(row_ids, e); push!(col_ids, rjoint[ball]); push!(blocks, Aw)
            for m in 1:k
                cm = sys.sidx[sys.colsb[ball][m]]
                push!(row_ids, e); push!(col_ids, cm)
                push!(blocks, -LT[:, (m - 1) * nf + 1:m * nf])
            end
            rhs_val[e] = -sys.bb[ball]
            e += 1                             # pin s = ε
            Ap = zeros(1, d); Ap[1, 1] = 1.0
            push!(row_ids, e); push!(col_ids, rjoint[ball]); push!(blocks, Ap)
            rhs_val[e] = [eps_ball(sys, inst, ball; eta)]
            e += 1                             # pin σ = c
            As = zeros(1, d); As[1, d] = 1.0
            push!(row_ids, e); push!(col_ids, rjoint[ball]); push!(blocks, As)
            rhs_val[e] = [sys.cb[ball]]
        else
            for j in 1:k                       # slice ties (triangular, dense)
                e += 1
                Aw = zeros(nf, nf + 1); Aw[:, 2:nf + 1] .= Inf_
                push!(row_ids, e); push!(col_ids, rslice[ball, j])
                push!(blocks, Aw)
                for m in j:k
                    cm = sys.sidx[sys.colsb[ball][m]]
                    push!(row_ids, e); push!(col_ids, cm)
                    push!(blocks, -LT[(j - 1) * nf + 1:j * nf,
                                      (m - 1) * nf + 1:m * nf])
                end
                rhs_val[e] = -sys.bb[ball][(j - 1) * nf + 1:j * nf]
            end
            for j in 1:k                       # s-links ŝ_j = s_j
                e += 1
                Ah = zeros(1, k + 2); Ah[1, 1 + j] = 1.0
                As = zeros(1, nf + 1); As[1, 1] = -1.0
                push!(row_ids, e); push!(col_ids, rhub[ball]); push!(blocks, Ah)
                push!(row_ids, e); push!(col_ids, rslice[ball, j])
                push!(blocks, As)
            end
            e += 1                             # hub pin h = ε
            Ap = zeros(1, k + 2); Ap[1, 1] = 1.0
            push!(row_ids, e); push!(col_ids, rhub[ball]); push!(blocks, Ap)
            rhs_val[e] = [eps_ball(sys, inst, ball; eta)]
            e += 1                             # σ pin = c
            As = zeros(1, k + 2); As[1, k + 2] = 1.0
            push!(row_ids, e); push!(col_ids, rhub[ball]); push!(blocks, As)
            rhs_val[e] = [sys.cb[ball]]
        end
    end
    if use_budget
        for a in 0:sys.Z - 1
            Pa = sys.programs_at(a)
            for p in Pa                        # budget ties v = R̂'f
                e += 1
                Aw = zeros(nf, nf + 1); Aw[:, 2:nf + 1] .= Inf_
                push!(row_ids, e); push!(col_ids, bslice[a, p])
                push!(blocks, Aw)
                push!(row_ids, e); push!(col_ids, sys.sidx[a, p])
                push!(blocks, -sys.RhatT[p + 1])
            end
            for (jj, p) in enumerate(Pa)       # u-links
                e += 1
                Ah = zeros(1, length(Pa) + 1); Ah[1, 1 + jj] = 1.0
                As = zeros(1, nf + 1); As[1, 1] = -1.0
                push!(row_ids, e); push!(col_ids, bhub[a]); push!(blocks, Ah)
                push!(row_ids, e); push!(col_ids, bslice[a, p])
                push!(blocks, As)
            end
            e += 1                             # hub pin h = β
            Ap = zeros(1, length(Pa) + 1); Ap[1, 1] = 1.0
            push!(row_ids, e); push!(col_ids, bhub[a]); push!(blocks, Ap)
            rhs_val[e] = [sys.betas[a + 1]]
        end
    end
    B = blocksparse(row_ids, col_ids, blocks)
    g_rhs = zeros(size(B, 1))
    for (edge, val) in rhs_val
        g_rhs[rowrange(B, edge)] .= val
    end
    Q = IPM.allocblockdiag(B); fill!(Q, 0)
    for (i, (_, p)) in enumerate(sys.stalks)
        block(Q, i, i, i) .= sys.Qs[p + 1]
    end
    c = zeros(size(B, 2))
    K_cones = IPM.AbstractCone[kind == :free ? CofreeCone() :
                               SecondOrderCone()
                               for (_, kind) in stalkdims]
    prob = IPMProblem(Q, B, c, g_rhs, K_cones)
    ctx = (; sys, inst, stalks = stalkdims, joint, budget = use_budget,
           nocoop, eta)
    return prob, ctx
end

"Filter bank (nstk × nf rows) from the solution vector."
function extract_f(prob, ctx, p::AbstractVector)
    nstk, nf = length(ctx.sys.stalks), ctx.sys.nf
    F = zeros(nstk * nf)
    for i in 1:nstk
        F[(i - 1) * nf + 1:i * nf] .= p[colrange(prob.B, i)]
    end
    return F
end

soundzone_objective(ctx, F) =
    0.5 * sum(dot(F[(i - 1) * ctx.sys.nf + 1:i * ctx.sys.nf],
                  ctx.sys.Qs[p + 1],
                  F[(i - 1) * ctx.sys.nf + 1:i * ctx.sys.nf])
              for (i, (_, p)) in enumerate(ctx.sys.stalks))

# -----------------------------------------------------------------------------
# Native TRS comparator (Z = 1, ρ = 0, no budget) — Gay 1981, MS 1983
# -----------------------------------------------------------------------------

function trs_single(sys, inst)
    ball = sys.balls[1]
    LT, b = sys.LbT[ball], sys.bb[ball]
    Q = sys.Qs[1]
    ex = sqrt(max(eps_ball(sys, inst, ball)^2 - sys.cb[ball]^2, 0.0))
    H = Symmetric((LT' \ Q) / LT)
    F = eigen(H)
    β = F.vectors' * (H * b)
    @assert norm(b) > ex "constraint inactive"
    ynorm(ν) = sqrt(sum(abs2, β ./ (F.values .+ ν)))
    lo, hi = 0.0, 1.0
    while ynorm(hi) > ex
        hi *= 4
    end
    for _ in 1:200
        mid = 0.5 * (lo + hi)
        if ynorm(mid) > ex
            lo = mid
        else
            hi = mid
        end
    end
    ν = 0.5 * (lo + hi)
    y = F.vectors * (-β ./ (F.values .+ ν))
    f = LT \ (y + b)
    return f, 0.5 * dot(f, Q, f),
        (; ν, boundary = abs(norm(y) - ex),
         stat = norm(H * (y + b) + ν * y))
end

# -----------------------------------------------------------------------------
# Baseline (same conic program via JuMP)
# -----------------------------------------------------------------------------

function build_jump_soundzone(prob, ctx, optimizer)
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
    vals = [abs(U[i, j]) for j in 1:n for i in 1:j]
    count(vals .> tol * maximum(vals)) / length(vals)
end

function test_dense_by_nature(sys)
    dQ = minimum(density(Q) for Q in sys.Qs)
    dG = minimum(density(Matrix(sys.Mb[b]' * sys.Mb[b])) for b in sys.balls)
    dL = minimum(tridensity(LT) for LT in values(sys.LbT))
    dR = minimum(tridensity(R) for R in sys.RhatT)
    @assert min(dQ, dG, dL, dR) > 0.95 "dense-by-nature: $((dQ, dG, dL, dR))"
    println("  [PASS] dense by nature: density@1e-8 Q_p ≥ ",
        "$(round(dQ, digits = 3)), per-ball G ≥ $(round(dG, digits = 3)), ",
        "reproduction ties ≥ $(round(dL, digits = 3)), budget ties ≥ ",
        "$(round(dR, digits = 3))")
end

"Exclusion (i): one array per program ⟹ block-diagonal Gram, tri-density
1/(2ρ+1) — the structural boundary cooperative rendering crosses."
function test_exclusion(sys, inst)
    z = sys.Z ÷ 2
    W = sys.window(z)
    k = length(W)
    nf = sys.nf
    rows = Matrix{Float64}[]
    for (i, zp) in enumerate(W)
        j0 = findfirst(==((zp, z)), sys.colsb[z, z])
        Az = sys.Mb[z, z][:, (j0 - 1) * nf + 1:j0 * nf]
        push!(rows, hcat([q == i ? Az : zeros(size(Az))
                          for q in 1:k]...))
    end
    Gx = (V = vcat(rows...); V' * V)
    offmax = 0.0
    for i in 1:k, j in 1:k
        i == j && continue
        offmax = max(offmax, maximum(abs.(Gx[(i - 1) * nf + 1:i * nf,
                                             (j - 1) * nf + 1:j * nf])))
    end
    td = tridensity(Matrix(cholesky(Symmetric(Gx +
        1e-12 * maximum(Gx) * I)).U))
    @assert offmax == 0.0 && td < 1.2 / k "exclusion: $offmax, $td"
    println("  [PASS] exclusion (i) control: one array per program gives ",
        "cross-stalk Gram blocks EXACTLY 0, tri-density ",
        "$(round(td, digits = 3)) ~ 1/(2ρ+1) — cooperative rendering ",
        "restores 1.000")
end

function test_whitening(sys)
    rng = MersenneTwister(3)
    worst = 0.0
    for ball in sys.balls
        φ = randn(rng, size(sys.Mb[ball], 2))
        lhs = norm(sys.Mb[ball] * φ - sys.gb[ball])^2
        rhs = norm(sys.LbT[ball] * φ - sys.bb[ball])^2 + sys.cb[ball]^2
        worst = max(worst, abs(lhs - rhs) / lhs)
        @assert sys.cb[ball] <= sys.floors[ball] + 1e-10
    end
    a = sys.Z ÷ 2
    lhs = rhs = 0.0
    for p in sys.programs_at(a)
        x = randn(rng, sys.nf)
        lhs += dot(x, sys.Qs[p + 1], x)
        rhs += norm(sys.RhatT[p + 1] * x)^2
    end
    wb = abs(lhs - rhs) / lhs
    @assert worst < 1e-10 && wb < 1e-12 "whitening: $worst, $wb"
    println("  [PASS] whitening identities at random points: reproduction ",
        "$(round(worst, sigdigits = 2)), budget partition ",
        "$(round(wb, sigdigits = 2)); local floors c ≤ global floors")
end

"Exclusion (ii): connectivity of the stalk-constraint graph."
function test_connectivity(sys)
    comp(budget) = begin
        parent = collect(1:length(sys.stalks))
        function find(i)
            while parent[i] != i
                parent[i] = parent[parent[i]]
                i = parent[i]
            end
            return i
        end
        un!(i, j) = (parent[find(i)] = find(j))
        for ball in sys.balls
            ids = [sys.sidx[c] for c in sys.colsb[ball]]
            for i in ids[2:end]
                un!(ids[1], i)
            end
        end
        if budget
            for a in 0:sys.Z - 1
                ids = [sys.sidx[a, p] for p in sys.programs_at(a)]
                for i in ids[2:end]
                    un!(ids[1], i)
                end
            end
        end
        length(Set(find(i) for i in 1:length(sys.stalks)))
    end
    c0, c1 = comp(false), comp(true)
    @assert c0 == sys.Z && c1 == 1 "connectivity: $c0, $c1"
    println("  [PASS] connectivity certificate: $c0 components without ",
        "budgets — exact decomposition by program — vs $c1 with (exclusion ",
        "(ii): non-decomposability, not reach)")
end

function test_objective_identity(prob, ctx, res)
    F = extract_f(prob, ctx, res.p)
    o_direct = soundzone_objective(ctx, F)
    o_ipm = ipm_objective(prob, res)
    rel = abs(o_direct - o_ipm) / abs(o_direct)
    @assert rel < 1e-5 "assembly mismatch: $o_direct vs $o_ipm"
    println("  [PASS] objective identity: F(f) = ",
        "$(round(o_direct, digits = 5)) (rel $(round(rel, sigdigits = 2)))")
    return F
end

function test_endpoint(inst, Z, settings)
    errs = Float64[]
    ηs = (2e-2, 5e-3)
    sys = make_system(inst, Z)
    for η in ηs
        prob, ctx = build_soundzone(sys, inst; budget = false, eta = η)
        res = solve(prob, settings)
        @assert res.status in (OPTIMAL, NEAR_OPTIMAL) "η=$η: $(res.status)"
        push!(errs, maximum(abs.(extract_f(prob, ctx, res.p) .- sys.fhat)))
    end
    rate = log(errs[2] / errs[1]) / log(ηs[2] / ηs[1])
    @assert errs[2] < errs[1] && 0.6 < rate < 1.4 "endpoint: $errs, $rate (oracle: η^0.98)"
    println("  [PASS] endpoint (no-budget): ‖f − f̂‖∞ ",
        "$(round(errs[1], sigdigits = 2)) → $(round(errs[2], sigdigits = 2)) ",
        "over η = $(ηs[1]) → $(ηs[2]) (rate η^$(round(rate, digits = 2)); ",
        "oracle: 0.98)")
end

function test_trs(inst, settings)
    inst1 = soundzone_instance(; S = inst.S, C = inst.C, Lf = inst.Lf,
        lhfac = inst.lhfac, rho = 0, gam = inst.gam, eta = inst.eta,
        theta = inst.theta, seed = inst.seed)
    sys1 = make_system(inst1, 1)
    prob, ctx = build_soundzone(sys1, inst1; budget = false, joint = false)
    res = solve(prob, settings)
    @assert res.status in (OPTIMAL, NEAR_OPTIMAL) "TRS build: $(res.status)"
    v_ipm = ipm_objective(prob, res)
    _, v_trs, ms = trs_single(sys1, inst1)
    dv = abs(v_ipm - v_trs) / abs(v_trs)
    @assert ms.boundary < 1e-10 && ms.stat < 1e-6 "MS: $ms"
    @assert dv < 1e-5 "IPM vs TRS: $dv (oracle: 4.0e-10 conic-vs-secular)"
    println("  [PASS] TRS classical certificate (Z = 1, ρ = 0 reduction): ",
        "IPM == secular solve, rel Δv $(round(dv, sigdigits = 2)), ",
        "ν = $(round(ms.ν, digits = 3))")
end

function test_load_bearing(sys, inst, settings, F)
    prob, ctx = build_soundzone(sys, inst; nocoop = true)
    res = solve(prob, settings)
    @assert res.status in (OPTIMAL, NEAR_OPTIMAL) "nocoop: $(res.status)"
    Fnc = extract_f(prob, ctx, res.p)
    viol = maximum(norm(sys.Mb[b] * sys.gather(Fnc, b) - sys.gb[b]) /
                   eps_ball(sys, inst, b)
                   for b in sys.balls if b[1] != b[2])
    util = sqrt.(sys.powers(F)) ./ sys.betas
    @assert viol > 1.3 "load-bearing: ×$viol (oracle: ×1.8)"
    @assert maximum(util) < 1.0 + 1e-6 "budget violated: $util"
    println("  [PASS] coupling load-bearing: non-cooperative baseline ",
        "violates the dark balls by ×$(round(viol, digits = 1)) (oracle: ",
        "×1.8); budget utilization in [$(round(minimum(util), digits = 2)), ",
        "$(round(maximum(util), digits = 2))] (oracle: [0.83, 0.94])")
end

function test_ipm_vs_clarabel(prob, ctx, F_ipm)
    model, xv = build_jump_soundzone(prob, ctx, clarabel_opt(; tol = TOL))
    optimize!(model)
    F_cla = extract_f(prob, ctx, value.(xv))
    df = maximum(abs.(F_ipm .- F_cla))
    @assert df < 1e-3 "IPM vs Clarabel: $df"
    println("  [PASS] IPM vs Clarabel (same conic program): ‖Δf‖∞ = ",
        "$(round(df, sigdigits = 2))")
end

# -----------------------------------------------------------------------------
# Sweep
# -----------------------------------------------------------------------------

zsizes() = OPTS.tiny ? [4] : OPTS.quick ? [4, 8, 16] : [4, 8, 16, 32]
lsizes() = OPTS.tiny ? [16] : OPTS.quick ? [16, 32] : [16, 32, 48]

function run()
    println("\n", "=" ^ 78)
    println("  E12: distributed sound zones, cooperative rendering ",
        LDIAL ? "(dial: block Lf; Z = 8)" : "(dial: zones Z)",
        JOINT ? " [JOINT ablation]" : "", NOBUDGET ? " [NO-BUDGET]" : "")
    println("=" ^ 78)

    ipm_settings = IPMSettings{Float64}(
        feas_tol = TOL, gap_tol = TOL, itmax = 200)
    hsd_settings = HSDSettings{Float64}(
        feas_tol = TOL, gap_tol = TOL, itmax = 200)

    inst = soundzone_instance()
    sys = make_system(inst, 6)

    println("\n  Gate tests (Z = 6, S = 4, C = 3, Lf = 16, Lh = 96):")
    test_dense_by_nature(sys)
    test_exclusion(sys, inst)
    test_whitening(sys)
    test_connectivity(sys)
    prob, ctx = build_soundzone(sys, inst; joint = false, budget = true)
    res = solve(prob, ipm_settings)
    @assert res.status in (OPTIMAL, NEAR_OPTIMAL) "IPM status $(res.status)"
    F = test_objective_identity(prob, ctx, res)
    test_endpoint(inst, 6, ipm_settings)
    test_trs(inst, ipm_settings)
    test_load_bearing(sys, inst, ipm_settings, F)
    test_ipm_vs_clarabel(prob, ctx, F)
    println()

    cla_opt = clarabel_opt(; tol = TOL)
    msk_opt = OPTS.mosek ? mosek_opt(; tol = TOL) : nothing

    rows = []
    for sz in (LDIAL ? lsizes() : zsizes())
        Z = LDIAL ? 8 : sz
        instk = LDIAL ? soundzone_instance(; Lf = sz) : inst
        sysk = make_system(instk, Z)
        probk, ctxk = build_soundzone(sysk, instk)
        stats = problem_stats(probk)

        @printf("  %s=%-4d dof=%-6d n1=%-6d blk=%-4.0f  ",
            LDIAL ? "Lf" : "Z", sz, stats.N0, stats.N1, stats.med_block)

        m_ipm = measure_ipm(probk, ipm_settings; nruns = NRUNS)
        m_hsd = measure_ipm(probk, hsd_settings; nruns = NRUNS)
        m_cla = measure_jump(() -> first(build_jump_soundzone(probk, ctxk, cla_opt)); nruns = NRUNS)
        # Mosek benchmarked in DUAL form (Dualization.jl): ~1.1x faster than primal here.
        m_msk = msk_opt !== nothing ?
            measure_jump(() -> first(build_jump_soundzone(probk, ctxk, Dualization.dual_optimizer(msk_opt))); nruns = NRUNS) :
            (t = NaN, status = "", obj = NaN)

        ratio(b) = isfinite(b.t) && isfinite(m_ipm.t) ? b.t / m_ipm.t : NaN
        fmt_ratio(b) = isnan(ratio(b)) ? "—" : (@sprintf("%.2fx", ratio(b)))

        println("IPM ", fmt_time(m_ipm), "  HSD ", fmt_time(m_hsd), " (", fmt_ratio(m_hsd),
            ")  Cla ", fmt_time(m_cla), " (", fmt_ratio(m_cla),
            ")  Msk ", fmt_time(m_msk), " (", fmt_ratio(m_msk), ")")

        push!(rows, (size = sz, dof = stats.N0, ipm = m_ipm, hsd = m_hsd, cla = m_cla, msk = m_msk))
    end

    dofs = [r.dof for r in rows]
    println("\nFitted log-log slopes (time vs DOF):")
    for (name, get_t) in [("IPM", r -> r.ipm.t), ("HSD", r -> r.hsd.t), ("Clarabel", r -> r.cla.t), ("Mosek", r -> r.msk.t)]
        sl = loglog_slope(dofs, [get_t(r) for r in rows])
        print("  $name: ", isnan(sl) ? "n/a" : (@sprintf("DOF^%.2f", sl)))
    end
    println()
end

if abspath(PROGRAM_FILE) == @__FILE__
    run()
end

# =============================================================================
# Sample run: 2026-07-25 (--quick --mosek)  [Mosek benchmarked DUAL — ~1.1x faster than primal]
# -----------------------------------------------------------------------------
# Z=4    dof=3040   n1=2364   blk=65    IPM   65.6ms  HSD   72.6ms (1.11x)  Cla  200.3ms (3.05x)  Msk  188.9ms (2.88x)
# Z=8    dof=7004   n1=5512   blk=65    IPM  167.5ms  HSD  179.9ms (1.07x)  Cla  528.3ms (3.15x)  Msk  501.6ms (2.99x)
# Z=16   dof=14932  n1=11808  blk=65    IPM  381.8ms  HSD  444.0ms (1.16x)  Cla 1151.6ms (3.02x)  Msk 1282.6ms (3.36x)
# IPM: DOF^1.11  HSD: DOF^1.14  Clarabel: DOF^1.10  Mosek: DOF^1.20  (dual shaves ~12% off primal at Z=16; sheaf still leads)
# =============================================================================
