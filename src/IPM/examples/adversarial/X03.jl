# adversarial/X03.jl — within-block COLUMN-norm spread under a single cone [I1]. Redesigned per addendum.
#
# NOTE: this is NOT the old X03 (free-cone nullity + decade-scaled coupling rows). The addendum flagged
# two prototype failures to avoid: (a) grading FREE columns makes F near-singular at every α (tiny B and
# tiny ridge-H on the same coordinate) — unsolvable, not controller-hostile; (b) pinning every row with
# the same full-scale column makes rows near-parallel (accidental X04 contamination). The mandated
# redesign, implemented here:
#
#   • Grade under the PositiveCone ONLY, where the log-barrier supplies H (so F stays PD at working α).
#   • A band of ~half the columns at full scale; each row is pinned by a DISTINCT band column (plus one
#     random second band entry), so the row's ∞-norm is set by an O(1) entry and E cannot inflate the
#     graded tail; the block-constant D (one scalar for the whole PositiveCone block) cannot reach inside.
#   • A graded tail (the other half of the columns) spanning `spread` decades: column norms 1 … 10^{−spread}.
#
# Invariant I1 (within-block column-norm ratio) then survives Ruiz exactly: post-equilibration ratio still
# ≈ 10^spread. The addendum's honest expectation (from the numpy emulation) is a SMALL plateau shift vs a
# spread=0 twin (1–2 decades) — unconstrained row scaling launders most within-block spread, because a row
# either carries a large entry (and does not need its small columns) or is boosted wholesale by E. A null
# result here (plateau barely moves) is confirmatory, not a failure of the run. Run with the extended grid
# (1e0–1e14).
#
# Manufactured optimum: p* strictly interior (all PositiveCone > 0) ⇒ d* = 0; ridge Q ⇒ p* the unique QP
# min. g = B p*, y* random, c = Q p* − Bᵀ y* − d*.

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

_colnorm_ratio(Bd) = (ns = [norm(@view Bd[:, j]) for j in axes(Bd, 2)]; maximum(ns) / minimum(ns))

"""
    build_narrow(; npos=24, spread=8.0, rows_per_col=3, ridge=1e-6, seed=1) -> (prob, meta)

Single PositiveCone block of dim `npos`, split into a band (`npos÷2` full-scale columns, one pinning each
row) and a graded tail (`npos÷2` columns spanning `spread` decades). Within-block column-norm ratio ≈
10^spread survives equilibration (invariant I1). `spread=0` is the matched benign twin. **`spread` is the
dial** (expected small effect — see header). `meta`: `pstar`, `colratio` (raw within-block spread), dials.
"""
function build_narrow(; npos::Int = 24, spread::Float64 = 8.0, rows_per_col::Int = 3,
                      ridge::Float64 = 1e-6, seed::Int = 1)
    rng = MersenneTwister(seed)
    nb = max(npos ÷ 2, 1)                 # band (full-scale) columns
    nt = npos - nb                        # graded-tail columns
    m = nb                                # one row pinned per band column
    Bd = zeros(m, npos)
    for i in 1:nb                         # row i pinned by distinct band column i, plus one random 2nd band entry
        Bd[i, i] = 1.0
        j = mod(i, nb) + 1                # a distinct second band column
        Bd[i, j] += 0.5 * randn(rng)
    end
    for t in 1:nt                         # graded tail: column norm 1 … 10^{−spread}
        k = nb + t
        sc = nt == 1 ? 1.0 : 10.0^(-spread * (t - 1) / (nt - 1))
        for _ in 1:min(rows_per_col, m)
            r = rand(rng, 1:m)
            Bd[r, k] += sc * (1 + 0.1 * randn(rng))
        end
        iszero(@view Bd[:, k]) && (Bd[rand(rng, 1:m), k] = sc)   # guarantee the column is nonempty
    end

    pstar = 1.0 .+ 0.5 .* rand(rng, npos)     # strictly interior
    ystar = randn(rng, m)
    g = Bd * pstar
    c = fill(ridge, npos) .* pstar .- Bd' * ystar     # d* = 0

    B = _scalar_blocksparse(Bd)
    Q = IPM.allocblockdiag(B); fill!(Q, 0)
    for v in 1:npos; IPM.block(Q, v, v, v) .= ridge; end
    K = IPM.AbstractCone[PositiveCone() for _ in 1:npos]

    meta = (npos = npos, nband = nb, ntail = nt, spread = spread, nv = npos, nrows = m,
            colratio = _colnorm_ratio(Bd), sigmin = minimum(svdvals(Bd)), pstar = pstar)
    return IPMProblem(Q, B, c, g, K), meta
end

"Matched benign twin (spread=0: tail columns at full scale)."
build_narrow_twin(; kw...) = build_narrow(; spread = 0.0, kw...)

function gate(; tol = 1e-8, kw...)
    prob, meta = build_narrow(; kw...)
    @printf("X03 within-block col-spread  npos=%d (band %d / tail %d)  spread=%.0f  colratio=%.2e  σmin(B)=%.2e\n",
        meta.npos, meta.nband, meta.ntail, meta.spread, meta.colratio, meta.sigmin)
    for (tag, S) in (("HSD", HSDSettings), ("IPM", IPMSettings))
        r = solve!(init(prob, S{Float64}(feas_tol = tol, gap_tol = tol, itmax = 300)))
        αs = [row.α for row in r.history]
        relerr = norm(r.p .- meta.pstar, Inf) / (1 + norm(meta.pstar, Inf))
        @printf("  %-3s status=%-15s niter=%-3d  relerr(p,p*)=%.2e  α∈[%.1e,%.1e]  %s\n",
            tag, string(r.status), r.niter, relerr, minimum(αs), maximum(αs), relerr ≤ 1e-5 ? "PASS" : "CHECK")
    end
    return nothing
end

if abspath(PROGRAM_FILE) == @__FILE__
    for s in (2.0, 4.0, 6.0, 8.0); gate(; spread = s); end
    println("  -- benign twin --"); gate(; spread = 0.0)
end
