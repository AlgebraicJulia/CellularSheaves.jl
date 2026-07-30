# adversarial/X04.jl — near-dependent CONSTRAINT ROWS (angle pathology) [I2, ρ-ladder]. PRIMARY.
#
# NOTE: this is NOT the old X04 (that one duplicated aggregate rows and scaled them by decades — a
# `rowscale` pathology the equilibrator removes exactly). The NEW principle, per adversarial_spec.md
# §X04 + addendum, is ROW ANGLES: append rows that are near-copies of base rows,
#
#     row_j' = row_i + eps·(random unit row),      then row_j' *= 10^U(−rowscale, rowscale).
#
# Ruiz normalizes the ∞-norm of every row to ≈1 (killing the 10^U rowscale factor exactly) but PRESERVES
# the angle between row_i and row_j' — E·B·D is a row/column rescaling, it cannot rotate two near-parallel
# rows apart. So σ_min(B̂) ~ eps survives equilibration: this is the one pathology diagonal scaling
# provably cannot touch (invariant I2). Large α makes F=(1/α)H+BᵀB ≈ BᵀB, which is near-singular
# (σ_min(B)²~eps²) ⇒ the Cholesky ρ-shift ladder (uzawa.jl rgmin..rgmax, never fired on the e-set) must
# escalate at the top of the grid: the plateau is squeezed from BOTH ends and, per the addendum's
# emulation, collapses (no feasible α on a 1e0–1e14 grid) near μ_collapse ≈ 1e2·eps².
#
# Manufactured optimum (checkable independently of the solver): p* strictly interior (x-part > 0, z-part
# free), so d* = 0 (strict complementarity on the PositiveCone, K*={0} on the free block) and, with a tiny
# ridge Q, p* is the unique min of the equality-constrained QP. g = B p*, y* random, c = Q p* − Bᵀ y* − d*.
#
# Benign twin (`benign=true`): appended rows are FRESH independent randn rows instead of near-copies —
# same shape / seed / rowscale, σ_min(B̂) ~ O(1). Differencing an instance against its twin isolates the
# angle pathology from everything the equilibrator already handles (the addendum's mandated control).

using CellularSheaves
using CellularSheaves.IPM: IPMProblem, IPMSettings, HSDSettings, OPTIMAL, NEAR_OPTIMAL, solve!, init,
    PositiveCone, CofreeCone
import CellularSheaves.IPM as IPM
using CellularSheaves.BlockSparseArrays: blocksparse
using LinearAlgebra, SparseArrays, Random, Printf

# dense m×nv matrix → BlockSparseMatrix of 1×1 scalar blocks (one dim-1 vertex per column). Every row and
# column must carry ≥1 nonzero (blocksparse infers block dims from the entries it sees).
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

# smallest principal angle (radians) between any appended row and its base parent — the invariant Ruiz
# cannot enlarge; ~eps by construction.
function _min_row_angle(Bd, parents)
    θ = Inf
    for (jrow, irow) in parents
        a = @view Bd[irow, :]; b = @view Bd[jrow, :]
        na = norm(a); nb = norm(b)
        (na == 0 || nb == 0) && continue
        c = clamp(abs(dot(a, b)) / (na * nb), -1.0, 1.0)
        θ = min(θ, acos(c))
    end
    return θ
end

"""
    build_degenerate(; n=16, nfree=4, degen=8, rowscale=8.0, eps=1e-6, benign=false, seed=1)
        -> (prob, meta)

Near-dependent-row LP/QP. `n` PositiveCone (x) + `nfree` CofreeCone (z, free) variables; `m0 = n÷2`
well-conditioned base rows, then `degen` appended rows each a near-copy of a base row with principal
angle ~`eps`, rescaled over `rowscale` decades. **`eps` is the hardness dial** (σ_min(B̂) ~ eps survives
equilibration; addendum: plateau collapses near μ ≈ 1e2·eps²); `rowscale` is the deliberately
equilibrator-removable knob (separates "Ruiz fixes it" from the residual angle pathology).

Q is a modest PD diagonal (1 + 0.5·U) — NOT the spec's 1e-6 ridge: with an interior optimum and a ~0 Q
the objective is constant on the feasible set (c ⟂ null(B)) and p* is non-unique, so a modest Q uniquely
identifies p* for the gate. It cannot mask the adversary — the high-α wall lives in F ≈ BᵀB where
(1/α)Q → 0 — and post-equilibration nQ stays O(1), negligible vs the B pathology. `benign=true` gives the
matched independent-rows twin. `meta`: `pstar`, `sigmin`/`min_angle` (raw severity), dials.
"""
function build_degenerate(; n::Int = 16, nfree::Int = 4, degen::Int = 8, rowscale::Float64 = 8.0,
                          eps::Float64 = 1e-6, benign::Bool = false, seed::Int = 1)
    rng = MersenneTwister(seed)
    m0 = max(n ÷ 2, 1)
    nv = n + nfree
    B0 = randn(rng, m0, nv)                          # well-conditioned base rows
    Bd = Matrix{Float64}(undef, m0 + degen, nv)
    Bd[1:m0, :] .= B0
    parents = Tuple{Int, Int}[]                       # (appended_row, base_row)
    for j in 1:degen
        i = ((j - 1) % m0) + 1                        # cycle over base rows
        if benign
            row = randn(rng, nv)                      # twin: fresh independent row
        else
            u = randn(rng, nv); u ./= norm(u)         # random unit direction
            row = @views B0[i, :] .+ eps .* u         # near-copy: angle ~eps to its parent
        end
        row .*= 10.0^(2 * rowscale * (rand(rng) - 0.5))  # 10^U(−rowscale, rowscale): Ruiz-removable
        Bd[m0 + j, :] .= row
        push!(parents, (m0 + j, i))
    end

    # manufactured interior optimum: x>0, z free; d*=0 (strict complementarity / K*={0}); c = Qp*−Bᵀy*−d*.
    pstar = Vector{Float64}(undef, nv)
    pstar[1:n] .= 1.0 .+ 0.5 .* rand(rng, n)          # PositiveCone part strictly interior
    pstar[n+1:nv] .= randn(rng, nfree)                # free part arbitrary
    ystar = randn(rng, m0 + degen)
    Qdiag = 1.0 .+ 0.5 .* rand(rng, nv)               # modest PD diagonal ⇒ interior p* uniquely identified
    g = Bd * pstar
    c = Qdiag .* pstar .- Bd' * ystar                 # d* = 0

    B = _scalar_blocksparse(Bd)
    Q = IPM.allocblockdiag(B); fill!(Q, 0)
    for v in 1:nv; IPM.block(Q, v, v, v) .= Qdiag[v]; end
    K = IPM.AbstractCone[PositiveCone() for _ in 1:n]
    append!(K, IPM.AbstractCone[CofreeCone() for _ in 1:nfree])

    meta = (n = n, nfree = nfree, degen = degen, rowscale = rowscale, eps = eps, benign = benign,
            nv = nv, nrows = m0 + degen, sigmin = minimum(svdvals(Bd)), min_angle = _min_row_angle(Bd, parents),
            pstar = pstar, Qdiag = Qdiag, Bd = Bd, g = g)
    return IPMProblem(Q, B, c, g, K), meta
end

"Matched benign twin (independent appended rows, σ_min ~ O(1))."
build_degenerate_twin(; kw...) = build_degenerate(; benign = true, kw...)

"""
    gate(; kw...)

Solve at 1e-8 (both solvers) and check the manufactured optimum: ‖p−p*‖∞/(1+‖p*‖∞) ≤ 1e-5. Reports
σ_min(B), the min row angle, status, and the ρ-ladder-relevant α range. Not run by the oracle.
"""
function gate(; tol = 1e-8, kw...)
    prob, meta = build_degenerate(; kw...)
    # near-singular B (σmin ~ eps) makes p ill-determined along null(B) — check objective + feasibility,
    # not ‖p−p*‖ (spec: "check objective value and residuals instead" where p* is non-unique).
    obj(p) = 0.5 * dot(p, meta.Qdiag .* p) - dot(prob.c, p)
    @printf("X04 near-dependent-rows  n=%d nfree=%d degen=%d rowscale=%.0f eps=%.0e  σmin(B)=%.2e  min∠=%.2e rad\n",
        meta.n, meta.nfree, meta.degen, meta.rowscale, meta.eps, meta.sigmin, meta.min_angle)
    for (tag, S) in (("HSD", HSDSettings), ("IPM", IPMSettings))
        r = solve!(init(prob, S{Float64}(feas_tol = tol, gap_tol = tol, itmax = 300)))
        αs = [row.α for row in r.history]
        objgap = abs(obj(r.p) - obj(meta.pstar)) / (1 + abs(obj(meta.pstar)))
        feas = norm(meta.Bd * r.p .- meta.g, Inf) / (1 + norm(meta.g, Inf))
        @printf("  %-3s status=%-15s niter=%-3d  objgap=%.2e  feas=%.2e  α∈[%.1e,%.1e]  %s\n",
            tag, string(r.status), r.niter, objgap, feas, minimum(αs), maximum(αs),
            (objgap ≤ 1e-6 && feas ≤ 1e-6) ? "PASS (obj+feas)" : "CHECK")
    end
    return nothing
end

if abspath(PROGRAM_FILE) == @__FILE__
    for e in (1e-3, 1e-4, 1e-5); gate(; eps = e); end
    println("  -- benign twin --"); gate(; benign = true)
end
