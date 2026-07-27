# =============================================================================
# e01/def.jl — Two-asset American basket put under Merton jump-diffusion
#         (dense-NOC habitat; implicit-Euler chain of QPs)
#
# Usage:
#   julia --project e01/run.jl              # Clarabel only
#   julia --project e01/run.jl --mosek      # + Mosek
#   julia --project e01/run.jl --quick      # small sweep
#
# Problem. Price an American basket put, payoff ψ = max(K − (S1+S2)/2, 0),
# under a two-asset Merton jump-diffusion (correlated Brownians, common
# Poisson jump times, iid Gaussian log-jump sizes). Backward in time this is
# a partial integro-differential complementarity problem; implicit Euler
# gives one LCP per timestep with a DENSE matrix (the jump integral couples
# every grid node to every other) — the classic dense-operator habitat
# (d'Halluin–Forsyth–Vetzal-style PIDE pricing).
#
# QP formulation. An LCP is the KKT system of a QP only for a SYMMETRIC
# matrix. The 2D integrating factor ρ_w = exp(2Σ⁻¹c·(x−lnK)) (Σ = diffusion
# covariance, c = risk-neutral log drifts) makes drift+diffusion self-adjoint
# in divergence form — the 2D analogue of ref/obstacle_option.md §3. The
# weighted JUMP kernel is symmetric only when 2μJ/σJ² = 2Σ⁻¹c, which generic
# parameters violate; the small antisymmetric residue Ga (‖Ga‖/‖Gs‖ ≈ 1e-4
# here), plus any jump coupling between nodes never co-resident in a patch
# (Gfar), is moved to the linear term and converged by Picard sweeps within
# each timestep (contraction ∝ λΔt ≈ 1e-3; ≤ ~4 sweeps to 1e-8·K). At the
# Picard fixed point the QP's KKT conditions are EXACTLY the full LCP: the
# splitting affects the iteration, not the answer (verified to 8e-13 against
# an independent full-matrix penalty solver — see b1_oracle.py).
#
# Sheaf decomposition. Q must stay SPD per patch. Entry-wise 1/#owners
# splitting can make overlap-corner blocks indefinite, so we element-split
# (the tensor_spline discipline): diffusion is assembled as PSD pieces —
# per-cell cross matrix grouped with a θ = |ρ|/2 share of the cell's four
# edge fluxes (θ* = |ρ|/2 exactly, independent of dx,dy,σ), plus PSD edge
# remainders — and mass/reaction/jump-diagonal are node-split. Each piece
# goes to the patches containing it, weight 1/#owners, so
# Σ_v embed(Q_v) + Gfar == Gs exactly and every block is a sum of PSD pieces
# plus a strictly positive mass diagonal (no ε regularization needed).
#
# Validation (b1_oracle.py, self-contained numpy):
#   * kernel mass (interior + payoff-extension boundary correction) = 1.000000
#   * European PIDE vs exact-simulation MC (4e6 paths): −0.124/−0.072/−0.039
#     at n=31/41/61 — clean O(h²) convergence to the MC value 5.6801 ± 0.0042
#   * weighted divergence form ≡ ρ_w·(standard FD PIDE): 1.1e-3 relative
#   * ‖u(sym-QP+Picard) − u(full-matrix penalty LCP)‖∞ = 8e-13 (same grid)
#   * American 5.8127 vs independent standard-FD chain 5.8109 at S0=(100,100)
# The gate tests below re-derive the fast subset of these in-process.
# =============================================================================

using AppleAccelerate
using LinearAlgebra, SparseArrays, Printf, Random

include("../utils.jl")

const OPTS = parse_args(ARGS)
const TOL = OPTS.tol
const NRUNS = 3

# -----------------------------------------------------------------------------
# Problem definition
# -----------------------------------------------------------------------------

struct MertonBasketInstance
    nx::Int            # grid nodes per dimension (incl. Dirichlet boundary)
    Px::Int; Py::Int   # patch grid
    olap::Int          # patch overlap (nodes) on the interior grid
    nsteps::Int        # implicit-Euler timesteps
    σ1::Float64; σ2::Float64; ρ::Float64; r::Float64; K::Float64
    λ::Float64; μJ::Float64; σJ::Float64; T::Float64
end

merton_instance(; nx = 21, Px = 2, Py = 2, olap = 5, nsteps = 10) =
    MertonBasketInstance(nx, Px, Py, olap, nsteps,
        0.3, 0.25, 0.5, 0.05, 100.0, 0.2, -0.1, 0.15, 0.5)

const LHALF = log(6.0)   # log-moneyness half-width: S ∈ [K/6, 6K]

# Flattening convention (column-major native): node (i,j) ↦ k = (j-1)*n + i,
# i = x-index (log S1), j = y-index (log S2). vec(M) matches this directly,
# and the separable jump kernel is kron(fy, fx).

gauss1(z, μJ, σJ) = exp(-0.5 * ((z - μJ) / σJ)^2) / (σJ * sqrt(2π))

"1D kernel slice: fx[i,p] = f(x[p] − x[i]) * dx (jump FROM node i TO node p)."
kernel1d(xfrom, xto, dx, μJ, σJ) =
    [gauss1(xt - xf, μJ, σJ) * dx for xf in xfrom, xt in xto]

# -----------------------------------------------------------------------------
# Build: grid, payoff, weighted symmetric operator, PSD element split
# -----------------------------------------------------------------------------

function build_merton_chain(inst::MertonBasketInstance)
    n = inst.nx
    N = n * n
    x = range(log(inst.K) - LHALF, log(inst.K) + LHALF; length = n)
    dx = step(x)
    dt = inst.T / inst.nsteps
    κ = exp(inst.μJ + 0.5 * inst.σJ^2) - 1.0
    c1 = inst.r - 0.5 * inst.σ1^2 - inst.λ * κ
    c2 = inst.r - 0.5 * inst.σ2^2 - inst.λ * κ
    Σ = [inst.σ1^2 inst.ρ*inst.σ1*inst.σ2; inst.ρ*inst.σ1*inst.σ2 inst.σ2^2]
    α = 2.0 * (Σ \ [c1, c2])

    ψmat = [max(inst.K - 0.5 * (exp(xi) + exp(yj)), 0.0) for xi in x, yj in x]
    wmat = [exp(α[1] * (xi - log(inst.K)) + α[2] * (yj - log(inst.K))) for xi in x, yj in x]
    ψ = vec(ψmat); w = vec(wmat)

    # ---- jump kernel (interior box) + payoff-extension boundary correction
    fx = kernel1d(x, x, dx, inst.μJ, inst.σJ)
    J = kron(fx, fx)                                # k=(j-1)n+i ⇒ kron(fy,fx); fy==fx here
    npad = ceil(Int, (5 * inst.σJ + abs(inst.μJ)) / dx)
    xe = range(first(x) - npad * dx, last(x) + npad * dx; length = n + 2npad)
    fxe = kernel1d(x, xe, dx, inst.μJ, inst.σJ)     # n × ne
    ψe = [max(inst.K - 0.5 * (exp(xi) + exp(yj)), 0.0) for xi in xe, yj in xe]
    # qb[i,j] = Σ_{ext nodes} f·ψ  =  (full extended sum) − (interior box sum)
    qbmat = fxe * ψe * fxe' .- fx * ψmat * fx'
    qb = vec(qbmat)
    kernel_mass = vec(sum(J, dims = 2)) .+ vec(fxe * ones(length(xe), length(xe)) * fxe') .-
                  vec(fx * ones(n, n) * fx')

    # ---- weighted symmetric step matrix Gw = diag(w)/dt + Aw − λ·diag(w)·J
    # (the −λu loss term lives inside Aw's (r+λ)ρ_w·u reaction — do not add it twice)
    Gw = zeros(N, N)
    idx(i, j) = (j - 1) * n + i
    # edge fluxes (xx along i, yy along j), half-point weights — PSD per edge
    px_e = zeros(n - 1, n); py_e = zeros(n, n - 1)
    for j in 1:n, i in 1:(n - 1)
        p = 0.5 * inst.σ1^2 * 0.5 * (wmat[i, j] + wmat[i + 1, j]) / dx^2
        px_e[i, j] = p
        a, b = idx(i, j), idx(i + 1, j)
        Gw[a, a] += p; Gw[b, b] += p; Gw[a, b] -= p; Gw[b, a] -= p
    end
    for j in 1:(n - 1), i in 1:n
        p = 0.5 * inst.σ2^2 * 0.5 * (wmat[i, j] + wmat[i, j + 1]) / dx^2
        py_e[i, j] = p
        a, b = idx(i, j), idx(i, j + 1)
        Gw[a, a] += p; Gw[b, b] += p; Gw[a, b] -= p; Gw[b, a] -= p
    end
    # cross term, cell-based one-point form, cell-center weights
    s12 = inst.ρ * inst.σ1 * inst.σ2
    gx = [-1.0, 1.0, -1.0, 1.0] ./ (2dx)
    gy = [-1.0, -1.0, 1.0, 1.0] ./ (2dx)
    C0 = 0.5 * s12 * (gx * gy' + gy * gx')          # 4×4, corner order (i,j),(i+1,j),(i,j+1),(i+1,j+1)
    bcell = zeros(n - 1, n - 1)
    for j in 1:(n - 1), i in 1:(n - 1)
        bc = 0.25 * (wmat[i, j] + wmat[i + 1, j] + wmat[i, j + 1] + wmat[i + 1, j + 1])
        bcell[i, j] = bc
        corners = (idx(i, j), idx(i + 1, j), idx(i, j + 1), idx(i + 1, j + 1))
        for (bq, kq) in enumerate(corners), (bp, kp) in enumerate(corners)
            Gw[kp, kq] += bc * C0[bp, bq]
        end
    end
    # reaction + mass (node diagonal), jump term
    for k in 1:N
        Gw[k, k] += (inst.r + inst.λ) * w[k] + w[k] / dt
    end
    Gw .-= inst.λ .* (w .* J)                       # row-scaled kernel: w[i]*J[i,:]

    # ---- Dirichlet elimination: interior unknowns, boundary pinned at ψ
    bnd = falses(n, n); bnd[1, :] .= true; bnd[n, :] .= true; bnd[:, 1] .= true; bnd[:, n] .= true
    bmask = vec(bnd); imask = .!bmask
    gidx = findall(imask)                            # global index of each interior unknown
    g2i = zeros(Int, N); for (t, k) in enumerate(gidx); g2i[k] = t; end
    Ni = length(gidx)
    Gii = Gw[gidx, findall(imask)]
    Gib = Gw[gidx, findall(bmask)]
    Gs = 0.5 .* (Gii .+ Gii')                        # symmetric part (QP quadratic form)
    Ga = Gii .- Gs                                   # antisymmetric jump residue → Picard
    Gw = nothing; J = nothing                        # drop full-grid dense N² matrices

    # ---- interior patches (rectangles on the (n−2)² interior grid)
    m = n - 2
    px = div(m - inst.olap, inst.Px) + inst.olap
    py = div(m - inst.olap, inst.Py) + inst.olap
    rects = NTuple{4, Int}[]
    for pj in 1:inst.Py, pi in 1:inst.Px
        # anchor the LAST patch to the end: guarantees full coverage of the
        # interior grid (plain strides leave a gap when m − olap is odd)
        i0 = pi == inst.Px ? m - px + 1 : (pi - 1) * (px - inst.olap) + 1
        j0 = pj == inst.Py ? m - py + 1 : (pj - 1) * (py - inst.olap) + 1
        push!(rects, (i0, i0 + px - 1, j0, j0 + py - 1))
    end
    iidx(ii, jj) = g2i[idx(ii + 1, jj + 1)]          # interior (ii,jj) → interior unknown id
    patches = [ [iidx(ii, jj) for jj in r[3]:r[4] for ii in r[1]:r[2]] for r in rects ]
    npatch = length(patches)
    loc = [zeros(Int, Ni) for _ in 1:npatch]         # interior id → local col in patch (0 = absent)
    for (v, p) in enumerate(patches), (t, k) in enumerate(p); loc[v][k] = t; end
    pmask = [UInt8(sum((1 << (v - 1) for v in 1:npatch if loc[v][k] > 0); init = 0)) for k in 1:Ni]
    node_owners = [count_ones(pmask[k]) for k in 1:Ni]
    @assert minimum(node_owners) ≥ 1 "interior node not covered by any patch"

    # ---- PSD element split of Gs into per-patch blocks + far remainder
    θ = abs(inst.ρ) / 2
    Qblocks = [zeros(length(p), length(p)) for p in patches]
    owners_of(ks::NTuple) = begin
        mask = ~UInt8(0)
        for k in ks; mask &= pmask[k]; end
        mask
    end
    addpiece!(ks, P) = begin
        mask = owners_of(ks)
        no = count_ones(mask)
        no ≥ 1 || return false
        for v in 1:npatch
            (mask >> (v - 1)) & 0x1 == 0x1 || continue
            Qv = Qblocks[v]; lv = loc[v]
            for (bq, kq) in enumerate(ks), (bp, kp) in enumerate(ks)
                Qv[lv[kp], lv[kq]] += P[bp, bq] / no
            end
        end
        true
    end
    # interior-cell pieces: cross + θ-share of the cell's 4 edge fluxes (PSD, weight b_c)
    F̂ = let e = [1.0 -1.0; -1.0 1.0] ./ dx^2
        Fs = [zeros(4, 4) for _ in 1:4]              # x-edges (1,2),(3,4); y-edges (1,3),(2,4)
        for (t, (a, b)) in enumerate(((1, 2), (3, 4), (1, 3), (2, 4)))
            Fs[t][a, a] = e[1, 1]; Fs[t][a, b] = e[1, 2]; Fs[t][b, a] = e[2, 1]; Fs[t][b, b] = e[2, 2]
        end
        Fs
    end
    edge_used_x = zeros(n - 1, n); edge_used_y = zeros(n, n - 1)
    for j in 1:(m - 1), i in 1:(m - 1)               # interior cells: corners all interior
        gi, gj = i + 1, j + 1                         # cell in full-grid coords
        bc = bcell[gi, gj]
        ks = (iidx(i, j), iidx(i + 1, j), iidx(i, j + 1), iidx(i + 1, j + 1))
        P = bc .* (C0 .+ θ .* (0.5 * inst.σ1^2 .* (F̂[1] .+ F̂[2]) .+ 0.5 * inst.σ2^2 .* (F̂[3] .+ F̂[4])))
        addpiece!(ks, P) || error("interior cell not covered by any patch (olap ≥ 1 required)")
        edge_used_x[gi, gj]     += θ * 0.5 * inst.σ1^2 * bc / dx^2
        edge_used_x[gi, gj + 1] += θ * 0.5 * inst.σ1^2 * bc / dx^2
        edge_used_y[gi, gj]     += θ * 0.5 * inst.σ2^2 * bc / dx^2
        edge_used_y[gi + 1, gj] += θ * 0.5 * inst.σ2^2 * bc / dx^2
    end
    # boundary-adjacent cells contribute only their interior-pair couplings to Gs;
    # those surviving pairs are edge fluxes and single cross couplings handled below
    # through the entry-wise remainder pass (they involve ≤ 2 interior nodes).
    # edge remainders (PSD: coefficient = assembled flux − θ-shares ≥ (1−|ρ|)·p > 0)
    for j in 2:(n - 1), i in 1:(n - 1)               # x-edges with both endpoints interior
        (i ≥ 2 && i + 1 ≤ n - 1) || continue
        pr = px_e[i, j] - edge_used_x[i, j]
        pr ≥ -1e-12 || error("negative edge remainder — θ grouping violated")
        ks = (iidx(i - 1, j - 1), iidx(i, j - 1))
        addpiece!(ks, pr .* [1.0 -1.0; -1.0 1.0]) || error("edge not covered")
    end
    for j in 1:(n - 1), i in 2:(n - 1)               # y-edges with both endpoints interior
        (j ≥ 2 && j + 1 ≤ n - 1) || continue
        pr = py_e[i, j] - edge_used_y[i, j]
        pr ≥ -1e-12 || error("negative edge remainder — θ grouping violated")
        ks = (iidx(i - 1, j - 1), iidx(i - 1, j))
        addpiece!(ks, pr .* [1.0 -1.0; -1.0 1.0]) || error("edge not covered")
    end
    # node pieces: everything on the diagonal of Gs not yet placed (mass, reaction,
    # jump diagonal, flux diagonals from boundary-adjacent edges/cells)
    placed_diag = zeros(Ni)
    for v in 1:npatch, (t, k) in enumerate(patches[v]); placed_diag[k] += Qblocks[v][t, t]; end
    for k in 1:Ni
        d = Gs[k, k] - placed_diag[k]
        for v in 1:npatch
            lv = loc[v][k]
            lv > 0 && (Qblocks[v][lv, lv] += d / node_owners[k])
        end
    end
    # remaining off-diagonal entries (jump kernel pairs + boundary-cell cross pairs):
    # pair-split where co-patched, else → far remainder (Picard RHS operator)
    placed = zeros(Ni, Ni)
    for v in 1:npatch
        p = patches[v]; Qv = Qblocks[v]
        for (tq, kq) in enumerate(p), (tp, kp) in enumerate(p)
            placed[kp, kq] += Qv[tp, tq]
        end
    end
    Gfar = zeros(Ni, Ni)
    for kq in 1:Ni, kp in 1:Ni
        kp == kq && continue
        d = Gs[kp, kq] - placed[kp, kq]
        abs(d) < 1e-300 && continue
        mask = pmask[kp] & pmask[kq]
        no = count_ones(mask)
        if no == 0
            Gfar[kp, kq] = d
        else
            for v in 1:npatch
                (mask >> (v - 1)) & 0x1 == 0x1 || continue
                Qblocks[v][loc[v][kp], loc[v][kq]] += d / no
            end
        end
    end
    Rop = Ga .+ Gfar                                  # Picard RHS operator

    # ---- agreement constraints on patch overlaps
    row_ids, col_ids, blocks = Int[], Int[], Matrix{Float64}[]
    edges = Tuple{Int, Int}[]
    for pj in 1:inst.Py, pi in 1:inst.Px
        p = (pj - 1) * inst.Px + pi
        pi < inst.Px && push!(edges, (p, p + 1))
        pj < inst.Py && push!(edges, (p, p + inst.Px))
    end
    for (e_idx, (p1, p2)) in enumerate(edges)
        overlap = intersect(patches[p1], patches[p2])
        isempty(overlap) && continue
        R1 = zeros(length(overlap), length(patches[p1]))
        R2 = zeros(length(overlap), length(patches[p2]))
        for (t, k) in enumerate(overlap)
            R1[t, loc[p1][k]] = 1; R2[t, loc[p2][k]] = 1
        end
        push!(row_ids, e_idx); push!(col_ids, p1); push!(blocks, R1)
        push!(row_ids, e_idx); push!(col_ids, p2); push!(blocks, -R2)
    end
    B = blocksparse(row_ids, col_ids, blocks)
    Q = IPM.allocblockdiag(B)
    for v in 1:npatch
        block(Q, v, v, v) .= Symmetric(0.5 .* (Qblocks[v] .+ Qblocks[v]'))
    end

    ψi = ψ[gidx]; ψb = ψ[findall(bmask)]; wi = w[gidx]; qbi = qb[gidx]
    h0 = Gii * ψi                                     # constant part of the linear term
    bpin = Gib * ψb
    c = zeros(sum(length.(patches)))
    g = zeros(size(B, 1))
    Kcones = [PositiveCone() for _ in 1:npatch]
    prob = IPMProblem(Q, B, c, g, Kcones)

    ctx = (; n, Ni, dt, gidx, patches, loc, node_owners, ψi, ψb, wi, qbi, h0, bpin,
             Gii, Gib, Gs, Ga, Gfar, Rop, kernel_mass, x, λ = inst.λ, K = inst.K,
             nsteps = inst.nsteps)
    prob, ctx
end

# -----------------------------------------------------------------------------
# Chain drivers (implicit Euler; Picard on the antisymmetric + far jump residue)
# -----------------------------------------------------------------------------

"Scatter a global interior linear term into per-patch c with 1/#owners weights."
function scatter_c!(c, ctx, cglob)
    off = 0
    for (v, p) in enumerate(ctx.patches)
        for (t, k) in enumerate(p)
            c[off + t] = cglob[k] / ctx.node_owners[k]
        end
        off += length(p)
    end
    c
end

"Gather a global interior vector from concatenated patch copies (owner-averaged)."
function gather_v(ctx, xsol)
    v = zeros(ctx.Ni)
    off = 0
    for (vp, p) in enumerate(ctx.patches)
        for (t, k) in enumerate(p)
            v[k] += xsol[off + t] / ctx.node_owners[k]
        end
        off += length(p)
    end
    v
end

function run_chain_ipm!(prob, ctx, settings; picard_tol = 1e-8, picard_max = 6)
    u = copy(ctx.ψi)
    sweeps = 0
    vbar = zeros(ctx.Ni)
    for _ in 1:ctx.nsteps
        b̃ = ctx.wi .* (u ./ ctx.dt .+ ctx.λ .* ctx.qbi) .- ctx.bpin
        vbar .= max.(u .- ctx.ψi, 0.0)
        local vres
        for _ in 1:picard_max
            cglob = b̃ .- ctx.h0 .- ctx.Rop * vbar  # negated for min ½p'Qp - c'p
            scatter_c!(prob.c, ctx, cglob)
            res = solve(prob, settings)
            res.status in (OPTIMAL, NEAR_OPTIMAL) || error("IPM status $(res.status) in chain")
            vres = gather_v(ctx, res.p)
            sweeps += 1
            δ = maximum(abs.(vres .- vbar))
            vbar .= vres
            δ < picard_tol * ctx.K && break
        end
        u .= vbar .+ ctx.ψi
    end
    u, sweeps
end

function build_jump_chain_model(prob, optimizer)
    nvar = size(prob.B, 2)
    model = Model(optimizer)
    set_silent(model)
    @variable(model, xv[1:nvar], lower_bound = 0)
    Qs = sparse(prob.Q); Bs = sparse(prob.B)
    @objective(model, Min, 0.5 * xv' * Qs * xv + prob.c' * xv)
    @constraint(model, Bs * xv .== prob.g)
    model, xv
end

"Chain with a prebuilt JuMP model (model construction excluded, as in measure_jump)."
function run_chain_jump!(model, xv, ctx; picard_tol = 1e-8, picard_max = 6)
    cvec = zeros(length(xv))
    u = copy(ctx.ψi)
    vbar = zeros(ctx.Ni)
    for _ in 1:ctx.nsteps
        b̃ = ctx.wi .* (u ./ ctx.dt .+ ctx.λ .* ctx.qbi) .- ctx.bpin
        vbar .= max.(u .- ctx.ψi, 0.0)
        for _ in 1:picard_max
            cglob = ctx.h0 .+ ctx.Rop * vbar .- b̃
            scatter_c!(cvec, ctx, cglob)
            if hasmethod(set_objective_coefficient,
                    Tuple{Model, Vector{VariableRef}, Vector{Float64}})
                set_objective_coefficient(model, xv, cvec)
            else
                set_objective_coefficient.(model, xv, cvec)
            end
            optimize!(model)
            vres = gather_v(ctx, value.(xv))
            δ = maximum(abs.(vres .- vbar))
            vbar .= vres
            δ < picard_tol * ctx.K && break
        end
        u .= vbar .+ ctx.ψi
    end
    u
end

function run_chain_jump(prob, ctx, optimizer; kwargs...)
    model, xv = build_jump_chain_model(prob, optimizer)
    run_chain_jump!(model, xv, ctx; kwargs...), model
end

# -----------------------------------------------------------------------------
# Monte Carlo oracle (exact simulation of the terminal law — no time stepping)
# -----------------------------------------------------------------------------

function mc_european(inst; npaths = 400_000, seed = 20260710)
    rng = MersenneTwister(seed)
    κ = exp(inst.μJ + 0.5 * inst.σJ^2) - 1.0
    c1 = inst.r - 0.5 * inst.σ1^2 - inst.λ * κ
    c2 = inst.r - 0.5 * inst.σ2^2 - inst.λ * κ
    ρc = sqrt(1 - inst.ρ^2)
    tot = 0.0; tot2 = 0.0
    for _ in 1:npaths
        z1 = randn(rng); z2 = inst.ρ * z1 + ρc * randn(rng)
        Np = rand_poisson(rng, inst.λ * inst.T)
        j1 = Np * inst.μJ + inst.σJ * sqrt(Np) * randn(rng)
        j2 = Np * inst.μJ + inst.σJ * sqrt(Np) * randn(rng)
        S1 = inst.K * exp(c1 * inst.T + inst.σ1 * sqrt(inst.T) * z1 + j1)
        S2 = inst.K * exp(c2 * inst.T + inst.σ2 * sqrt(inst.T) * z2 + j2)
        p = exp(-inst.r * inst.T) * max(inst.K - 0.5 * (S1 + S2), 0.0)
        tot += p; tot2 += p^2
    end
    μ = tot / npaths
    μ, sqrt(max(tot2 / npaths - μ^2, 0.0) / npaths)
end

function rand_poisson(rng, λ)
    L = exp(-λ); k = 0; p = 1.0
    while true
        p *= rand(rng)
        p ≤ L && return k
        k += 1
    end
end

# -----------------------------------------------------------------------------
# Gate tests
# -----------------------------------------------------------------------------

function test_kernel_and_split()
    inst = merton_instance(; nx = 21, olap = 5, nsteps = 10)
    prob, ctx = build_merton_chain(inst)
    merr = maximum(abs.(ctx.kernel_mass .- 1.0))
    @assert merr < 5e-3 "kernel mass off by $merr"
    # Σ_v embed(Q_v) + Gfar == Gs exactly
    S = copy(ctx.Gfar)
    off = 0
    for (v, p) in enumerate(ctx.patches)
        Qv = block(prob.Q, v, v, v)
        for (tq, kq) in enumerate(p), (tp, kp) in enumerate(p)
            S[kp, kq] += Qv[tp, tq]
        end
        off += length(p)
    end
    serr = maximum(abs.(S .- ctx.Gs)) / maximum(abs.(ctx.Gs))
    @assert serr < 1e-10 "split inconsistency: $serr"
    for v in 1:length(ctx.patches)
        cholesky(Symmetric(Matrix(block(prob.Q, v, v, v))))   # throws if not SPD
    end
    ratio = norm(ctx.Ga) / norm(ctx.Gs)
    println("  [PASS] kernel mass (err=$(round(merr, sigdigits=2))); ",
        "Σembed(Q)+Gfar==Gs (rel $(round(serr, sigdigits=2))); all blocks SPD; ",
        "‖Ga‖/‖Gs‖=$(round(ratio, sigdigits=2))")
end

function test_european_vs_mc()
    inst = merton_instance(; nx = 31, olap = 8, nsteps = 50)
    prob, ctx = build_merton_chain(inst)
    # direct dense chain on the full weighted system (obstacle off)
    n = inst.nx
    F = lu(ctx.Gii)
    u = copy(ctx.ψi)
    bglob = setdiff(1:(n * n), ctx.gidx)
    for mstep in 1:inst.nsteps
        τ = mstep * ctx.dt
        ub = [begin
            i = mod1(k, n); j = div(k - 1, n) + 1
            max(inst.K * exp(-inst.r * τ) - 0.5 * (exp(ctx.x[i]) + exp(ctx.x[j])), 0.0)
        end for k in bglob]
        rhs = ctx.wi .* (u ./ ctx.dt .+ ctx.λ .* ctx.qbi) .- ctx.Gib * ub
        u = F \ rhs
    end
    i0 = argmin(abs.(ctx.x .- log(inst.K)))
    kglob = (i0 - 1) * n + i0                       # k = (j-1)n+i with i=j=i0
    ti = findfirst(==(kglob), ctx.gidx)
    v_pide = u[ti]
    v_mc, se = mc_european(inst)
    rel = abs(v_pide - v_mc) / v_mc
    @assert rel < 0.04 "European PIDE $(v_pide) vs MC $(v_mc) (rel err $(rel))"
    println("  [PASS] European vs exact-sim MC: PIDE=$(round(v_pide, digits=4)) ",
        "MC=$(round(v_mc, digits=4))±$(round(se, digits=4)) (rel $(round(rel, sigdigits=2)))")
    v_mc
end

function test_american_chain(settings, v_euro_mc)
    inst = merton_instance(; nx = 21, olap = 5, nsteps = 10)
    prob, ctx = build_merton_chain(inst)
    u, sweeps = run_chain_ipm!(prob, ctx, settings)
    @assert minimum(u .- ctx.ψi) > -1e-6 * inst.K "obstacle violated"
    contact = count(abs.(u .- ctx.ψi) .< 1e-5 * inst.K) / ctx.Ni
    @assert contact > 0.2 "contact set suspiciously small: $contact"
    @assert sweeps ≤ inst.nsteps * 6 "Picard failed to converge"
    i0 = argmin(abs.(ctx.x .- log(inst.K)))
    ti = findfirst(==((i0 - 1) * inst.nx + i0), ctx.gidx)
    v_am = u[ti]
    @assert v_am > 0.8 * v_euro_mc "American value implausibly low: $v_am"
    # cross-check the whole chain against Clarabel solving the same QPs
    prob2, ctx2 = build_merton_chain(inst)
    u2, _ = run_chain_jump(prob2, ctx2, clarabel_opt(; tol = TOL))
    derr = maximum(abs.(u .- u2))
    @assert derr < 1e-4 * inst.K "IPM vs Clarabel chain mismatch: $derr"
    println("  [PASS] American chain: value=$(round(v_am, digits=4)) ",
        "(euro MC $(round(v_euro_mc, digits=4))), contact=$(round(contact, digits=2)), ",
        "picard sweeps/step=$(round(sweeps / inst.nsteps, digits=1)), ",
        "‖IPM−Clarabel‖∞=$(round(derr, sigdigits=2))")
end

# -----------------------------------------------------------------------------
# Sweep
# -----------------------------------------------------------------------------

sizes() = OPTS.tiny ? [(nx = 11, olap = 3)] :
    OPTS.quick ? [(nx = 11, olap = 3), (nx = 15, olap = 4), (nx = 21, olap = 5)] :
    [(nx = 11, olap = 3), (nx = 15, olap = 4), (nx = 21, olap = 5), (nx = 31, olap = 8),
     (nx = 41, olap = 10), (nx = 55, olap = 14), (nx = 71, olap = 17)]

const NSTEPS_BENCH = 10

function run()
    println("\n", "=" ^ 78)
    println("  E01: Two-asset American basket put, Merton jumps ",
        "(implicit-Euler chain × $NSTEPS_BENCH steps; dial: nx)")
    println("=" ^ 78)

    ipm_settings = IPMSettings{Float64}(
        feas_tol = TOL, gap_tol = TOL, itmax = 200)

    hsd_settings = HSDSettings{Float64}(
        feas_tol = TOL, gap_tol = TOL, itmax = 200)

    println("\n  Gate tests:")
    test_kernel_and_split()
    v_mc = test_european_vs_mc()
    test_american_chain(ipm_settings, v_mc)
    println()

    cla_opt = clarabel_opt(; tol = TOL)
    msk_opt = OPTS.mosek ? mosek_opt(; tol = TOL) : nothing

    rows = []
    for cfg in sizes()
        inst = merton_instance(; nx = cfg.nx, olap = cfg.olap, nsteps = NSTEPS_BENCH)
        prob, ctx = build_merton_chain(inst)
        stats = problem_stats(prob)

        @printf("  nx=%-4d dof=%-6d n1=%-5d blk=%-4.0f  ",
            cfg.nx, stats.N0, stats.N1, stats.med_block)

        t_ipm = try
            timed_min(() -> run_chain_ipm!(prob, ctx, ipm_settings); nruns = NRUNS)
        catch; NaN end
        m_ipm = (t = t_ipm, status = isfinite(t_ipm) ? "OPTIMAL" : "ERROR", obj = NaN)

        t_hsd = try
            timed_min(() -> run_chain_ipm!(prob, ctx, hsd_settings); nruns = NRUNS)
        catch; NaN end
        m_hsd = (t = t_hsd, status = isfinite(t_hsd) ? "OPTIMAL" : "ERROR", obj = NaN)

        t_cla = try
            model, xv = build_jump_chain_model(prob, cla_opt)
            timed_min(() -> run_chain_jump!(model, xv, ctx); nruns = NRUNS)
        catch; NaN end
        m_cla = (t = t_cla, status = isfinite(t_cla) ? "OPTIMAL" : "ERROR", obj = NaN)

        m_msk = if msk_opt !== nothing
            t = try
                model, xv = build_jump_chain_model(prob, msk_opt)
                timed_min(() -> run_chain_jump!(model, xv, ctx); nruns = NRUNS)
            catch; NaN end
            (t = t, status = isfinite(t) ? "OPTIMAL" : "ERROR", obj = NaN)
        else
            (t = NaN, status = "", obj = NaN)
        end
        # Mosek benchmarked in PRIMAL form (dual ~3.4x slower here).

        ratio(b, base) = isfinite(b.t) && isfinite(base.t) ? b.t / base.t : NaN
        fmt_ratio(b, base) = isnan(ratio(b, base)) ? "—" : @sprintf("%.2fx", ratio(b, base))

        println("IPM ", fmt_time(m_ipm), "  HSD ", fmt_time(m_hsd), " (", fmt_ratio(m_hsd, m_ipm),
            ")  Cla ", fmt_time(m_cla), " (", fmt_ratio(m_cla, m_ipm),
            ")  Msk ", fmt_time(m_msk), " (", fmt_ratio(m_msk, m_ipm), ")")

        push!(rows, (size = cfg.nx, dof = stats.N0, ipm = m_ipm, hsd = m_hsd, cla = m_cla, msk = m_msk))
    end

    dofs = [r.dof for r in rows]
    println("\nFitted log-log slopes (chain time vs DOF):")
    for (name, get_t) in [("IPM", r -> r.ipm.t), ("HSD", r -> r.hsd.t), ("Clarabel", r -> r.cla.t), ("Mosek", r -> r.msk.t)]
        sl = loglog_slope(dofs, [get_t(r) for r in rows])
        print("  $name: ", isnan(sl) ? "n/a" : @sprintf("DOF^%.2f", sl))
    end
    println()
end
