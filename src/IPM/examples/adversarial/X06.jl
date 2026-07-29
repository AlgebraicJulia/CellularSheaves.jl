# adversarial/X06.jl — SOC corner optimum [cone-generality control]. Per adversarial_spec.md §X06.
#
# NOTE: not the old X06 (near-rank-deficient B). This is a CONTROL, not a scaling trick: the e-set oracle
# runs are PositiveCone-dominated, so the measured laws (q ≈ 1, floor constant range, N(μ) exponent) might
# be artifacts of the orthant's diagonal NT scaling. SecondOrderCone blocks whose optimum sits at the cone
# CORNER (apex coord x₁ = ‖x̄‖ with x̄ ≈ 0) maximize the NT Hessian's intrinsic within-block spread — a
# spread NO diagonal equilibrator touches, because it lives in the scaling point w, not in B or Q. If the
# floor/q laws survive here the controller design is cone-generic; if the N(μ) exponent shifts, the
# feedforward constant is cone-dependent and the r₀ controller's self-calibration becomes load-bearing.
#
# SOC convention (soc.jl): (x₁, x̄) ∈ K ⟺ x₁ ≥ ‖x̄‖; the cone is self-dual (K*=K); complementary pairs are
# J-reflections, x ∘ s = 0 ⟺ s = σ·(x₁, −x̄).
#
# Manufactured optimum (half corner, half interior):
#   • interior block: p*=(t, x̄), t = ‖x̄‖+1 (strictly inside) ⇒ d* = 0 (strict complementarity).
#   • corner block:   x̄ = r·û small, p*=(r, x̄) on the boundary near the apex; d* = σ·(r, −x̄) on the
#     reflected boundary ray (σ>0). Then p*∘d* = 0 exactly (verified in-header algebra).
# g = B p*, y* random, Q modest block-SPD (I + 0.1·MᵀM/dim), c = Q p* − Bᵀ y* − d*. B is a well-conditioned
# random Gaussian m×N with m = N÷2 (full row rank, σ_min O(1)) — no B pathology, the cone is the point.

using CellularSheaves
using CellularSheaves.IPM: IPMProblem, IPMSettings, HSDSettings, OPTIMAL, NEAR_OPTIMAL, solve!, init,
    SecondOrderCone
import CellularSheaves.IPM as IPM
using CellularSheaves.BlockSparseArrays: blocksparse
using LinearAlgebra, SparseArrays, Random, Printf

"""
    build_corner_soc(; nsoc=6, dim=10, corner_r=0.1, benign=false, seed=1) -> (prob, meta)

`nsoc` SecondOrderCone blocks of dimension `dim`. Half sit at the apex corner (radius `corner_r`), half
strictly interior; `benign=true` puts ALL blocks interior (the matched low-spread twin). B is well-
conditioned random, Q modest block-SPD. `meta`: `pstar`, `dstar`, `ncorner`, `Qd` (dense Q for the gate
objective check), dials.
"""
function build_corner_soc(; nsoc::Int = 6, dim::Int = 10, corner_r::Float64 = 0.1,
                          benign::Bool = false, seed::Int = 1)
    @assert dim ≥ 2
    rng = MersenneTwister(seed)
    N = nsoc * dim
    m = max(N ÷ 2, 1)
    Bd = randn(rng, m, N)                     # well-conditioned constraint block

    pstar = Vector{Float64}(undef, N)
    dstar = zeros(N)
    Qd = zeros(N, N)
    ncorner = 0
    for v in 1:nsoc
        idx = (v - 1) * dim + 1 : v * dim
        xbar = randn(rng, dim - 1)
        corner = !benign && isodd(v)          # half the blocks at the corner
        if corner
            ncorner += 1
            xbar .*= corner_r / max(norm(xbar), eps())      # ‖x̄‖ = corner_r (near apex)
            r = norm(xbar)
            pstar[idx] .= vcat(r, xbar)                     # on the boundary: apex coord = ‖x̄‖
            σ = 0.5 + rand(rng)
            dstar[idx] .= σ .* vcat(r, -xbar)               # reflected boundary ray ⇒ p*∘d* = 0
        else
            pstar[idx] .= vcat(norm(xbar) + 1.0, xbar)      # strictly interior ⇒ d* = 0
        end
        M = randn(rng, dim, dim)
        Qd[idx, idx] .= Matrix(I, dim, dim) .+ 0.1 .* (M' * M) ./ dim   # modest SPD
    end
    ystar = randn(rng, m)
    g = Bd * pstar
    c = Qd * pstar .- Bd' * ystar .- dstar

    rows = Int[]; cols = Int[]; blks = Matrix{Float64}[]
    for v in 1:nsoc
        idx = (v - 1) * dim + 1 : v * dim
        for i in 1:m
            push!(rows, i); push!(cols, v); push!(blks, reshape(Vector(@view Bd[i, idx]), 1, dim))
        end
    end
    B = blocksparse(rows, cols, blks)
    Q = IPM.allocblockdiag(B); fill!(Q, 0)
    for v in 1:nsoc
        idx = (v - 1) * dim + 1 : v * dim
        IPM.block(Q, v, v, v) .= @view Qd[idx, idx]
    end
    K = IPM.AbstractCone[SecondOrderCone() for _ in 1:nsoc]

    meta = (nsoc = nsoc, dim = dim, corner_r = corner_r, benign = benign, ncorner = ncorner,
            nv = N, nrows = m, pstar = pstar, dstar = dstar, Qd = Qd)
    return IPMProblem(Q, B, c, g, K), meta
end

"Matched benign twin (all blocks strictly interior — minimal NT within-block spread)."
build_corner_soc_twin(; kw...) = build_corner_soc(; benign = true, kw...)

function gate(; tol = 1e-8, kw...)
    prob, meta = build_corner_soc(; kw...)
    obj(p) = 0.5 * dot(p, meta.Qd * p) - dot(prob.c, p)
    @printf("X06 SOC corners  nsoc=%d dim=%d  corner=%d/%d  (N=%d, m=%d)\n",
        meta.nsoc, meta.dim, meta.ncorner, meta.nsoc, meta.nv, meta.nrows)
    for (tag, S) in (("HSD", HSDSettings), ("IPM", IPMSettings))
        r = solve!(init(prob, S{Float64}(feas_tol = tol, gap_tol = tol, itmax = 300)))
        αs = [row.α for row in r.history]
        relerr = norm(r.p .- meta.pstar, Inf) / (1 + norm(meta.pstar, Inf))
        objgap = abs(obj(r.p) - obj(meta.pstar)) / (1 + abs(obj(meta.pstar)))
        @printf("  %-3s status=%-15s niter=%-3d  relerr(p,p*)=%.2e  objgap=%.2e  α∈[%.1e,%.1e]  %s\n",
            tag, string(r.status), r.niter, relerr, objgap, minimum(αs), maximum(αs),
            objgap ≤ 1e-6 ? "PASS (obj)" : "CHECK")
    end
    return nothing
end

if abspath(PROGRAM_FILE) == @__FILE__
    gate()
    println("  -- benign twin (all interior) --"); gate(; benign = true)
end
