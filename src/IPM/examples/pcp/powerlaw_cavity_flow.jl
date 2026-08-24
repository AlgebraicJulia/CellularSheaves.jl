# =============================================================================
#
# # Power-Law Fluid Flow
#
# Most non-Newtonian liquids (e.g. polymer melts, slurries, blood, ice) have
# a that depends on how fast they are sheared. The simplest and most widely
# used model is the power-law fluid, whose stress grows like the (p - 1)st
# power of the strain rate. Different values of p model different types
# of fluids:
#
#  - p < 2: shear-thinning fluids (e.g. polymers)
#  - p = 2: Newtonian fluids
#  - p > 2: sheaf-thickening fluids (e.g. dense suspensions)
#
# Let Ω ⊆ ℝ² be a fluid domain, u: Ω → ℝ² a velocity field, and
#
#   d = 1/2 (∇u + ∇uᵀ)          ‖ d ‖ = √(2 d:d)
#
# the strain rate tensor and its norm. Given a consistency K > 0, an exponent
# p > 1, a body force f: Ω → ℝ², and boundary data g, the creeping flow
# minimizes the viscous dissipation
#
#   min  ∫ K/p ‖ d ‖ᵖ dx - ∫ f ⋅ u dx        (1)
#    u   Ω                  Ω
#
#   s.t. div u = 0  in  Ω
#            u = g  on ∂Ω
#
# Discretizing exactly as in viscoplastic_flow.jl (P2 elements, deviatoric
# components qⱼ = Gⱼ v and divergence rows hⱼᵀ v at quadrature points), and
# epigraphing the p-th power through a second order cone chained to a power
# cone, (1) becomes
#
#             m
#   min      ∑ ωⱼ K/p tⱼ - cᵀv               (2)
#  v,q,s,t  j=1
#
#   s.t.  qⱼ - Gⱼ v = 0,     ρⱼ - sⱼ = 0,       j = 1, …, m
#             hⱼᵀ v = 0,          σⱼ = 1,       j = 1, …, m
#
#         (sⱼ, qⱼ) ∈ Q³,   (tⱼ, σⱼ, ρⱼ) ∈ ℙ,    j = 1, …, m
#
# where
#
#   ℚ³ = { (s, q) : ‖ q ‖ ≤ s }
#
# is the second-order cone and
#
#    ℙ = { (a, b, c) ∈ ℝ³ : a ≥ 0, b ≥ 0, aᵗ b¹⁻ᵗ ≥ |c| }
#
# is the three-dimensional power cone with parameter t = 1/p.
#
# =============================================================================

using LinearAlgebra, Printf

using CellularSheaves.IPM
using CellularSheaves.IPM: CofreeCone, SecondOrderCone, PowerCone,
                           IPMSettings, HSDSettings
using CellularSheaves.BlockSparseArrays: colrange, block, blocksparse

include("../utils.jl")

# -----------------------------------------------------------------------------
# P2 reference element, mesh, element quantities
# (shared with socp/viscoplastic_flow.jl)
# -----------------------------------------------------------------------------

const GAUSS = [((1 / 6, 1 / 6), 1 / 6), ((2 / 3, 1 / 6), 1 / 6), ((1 / 6, 2 / 3), 1 / 6)]

# Gradients of the six P2 shape functions on the reference triangle, as a 2 × 6
# matrix. Nodes 1:3 are the vertices, 4:6 the midpoints of (1,2), (2,3), (3,1).
function p2grads(xi::Float64, eta::Float64)
    lam = (1 - xi - eta, xi, eta)
    dlam = [-1.0 1.0 0.0; -1.0 0.0 1.0]          # d(λ)/d(ξ, η), 2 × 3
    g = zeros(2, 6)

    for k in 1:3
        g[:, k] .= (4 * lam[k] - 1) .* view(dlam, :, k)
    end

    for (k, (a, b)) in enumerate(((1, 2), (2, 3), (3, 1)))
        g[:, 3 + k] .= 4 .* (view(dlam, :, a) .* lam[b] .+ lam[a] .* view(dlam, :, b))
    end

    return g
end

struct P2Mesh
    pts::Matrix{Float64}          # 2 × nnode
    elems::Matrix{Int}            # 6 × nelem
end

nnodes(m::P2Mesh) = size(m.pts, 2)
nelems(m::P2Mesh) = size(m.elems, 2)

# The unit square, N × N cells, each split into two triangles, with midside
# nodes added.
function unitsquare(N::Int)
    vid(i, j) = j * (N + 1) + i + 1                       # 0-based (i, j)
    nv = (N + 1)^2

    xy = zeros(2, nv)

    for j in 0:N, i in 0:N
        xy[1, vid(i, j)] = i / N
        xy[2, vid(i, j)] = j / N
    end

    tris = NTuple{3, Int}[]

    for j in 0:(N - 1), i in 0:(N - 1)
        a, b = vid(i, j), vid(i + 1, j)
        c, d = vid(i + 1, j + 1), vid(i, j + 1)
        push!(tris, (a, b, c))
        push!(tris, (a, c, d))
    end

    mid = Dict{Tuple{Int, Int}, Int}()
    extra = Vector{Float64}[]

    function midnode(u, v)
        key = (min(u, v), max(u, v))

        if !haskey(mid, key)
            push!(extra, 0.5 .* (xy[:, u] .+ xy[:, v]))
            mid[key] = nv + length(extra)
        end

        return mid[key]
    end

    elems = zeros(Int, 6, length(tris))

    for (e, (a, b, c)) in enumerate(tris)
        elems[:, e] .= (a, b, c, midnode(a, b), midnode(b, c), midnode(c, a))
    end

    pts = hcat(xy, reduce(hcat, extra))
    return P2Mesh(pts, elems)
end

# Which velocity components of each node are free, and the prescribed values of
# the rest.
struct Layout
    stalk::Vector{Int}
    comps::Vector{Vector{Int}}
    fixed::Matrix{Float64}        # 2 × nnode, prescribed values (0 where free)
    nvstalk::Int
end

# Lid-driven cavity: no slip everywhere, except the top edge which moves at
# (V, 0).
function cavity_layout(m::P2Mesh; V::Float64 = 1.0)
    n = nnodes(m)
    stalk = zeros(Int, n)
    comps = [Int[] for _ in 1:n]
    fixed = zeros(2, n)
    k = 0

    for i in 1:n
        x, y = m.pts[1, i], m.pts[2, i]
        onbnd = x ≈ 0 || x ≈ 1 || y ≈ 0 || y ≈ 1

        if onbnd
            fixed[1, i] = (y ≈ 1) ? V : 0.0
            fixed[2, i] = 0.0
        else
            k += 1
            stalk[i] = k
            comps[i] = [1, 2]
        end
    end

    return Layout(stalk, comps, fixed, k)
end

localcols(k::Int) = (2k - 1):(2k)

# Returns (weight, Ĝ, ĥ) per Gauss point, with Ĝ (2 × 12) the deviatoric strain
# operator and ĥ (1 × 12) the divergence; see viscoplastic_flow.jl.
function element_ops(m::P2Mesh, e::Int)
    el = view(m.elems, :, e)
    p1, p2, p3 = m.pts[:, el[1]], m.pts[:, el[2]], m.pts[:, el[3]]
    J = hcat(p2 .- p1, p3 .- p1)
    detJ = abs(det(J))
    Jit = transpose(inv(J))

    ops = Tuple{Float64, Matrix{Float64}, Vector{Float64}}[]

    for ((xi, eta), w) in GAUSS
        gr = Jit * p2grads(xi, eta)
        DN = zeros(3, 12)
        DN[1, 1:2:end] .= view(gr, 1, :)
        DN[2, 2:2:end] .= view(gr, 2, :)
        DN[3, 1:2:end] .= view(gr, 2, :)
        DN[3, 2:2:end] .= view(gr, 1, :)

        Ghat = vcat(transpose(DN[1, :] .- DN[2, :]), transpose(DN[3, :]))
        hhat = DN[1, :] .+ DN[2, :]
        push!(ops, (w * detJ, Ghat, hhat))
    end

    return ops
end

# -----------------------------------------------------------------------------
# Build
# -----------------------------------------------------------------------------

"""
    build_powerlaw(N; K, p, V)

Power-law fluid lid-driven cavity on an N × N structured mesh (lid speed V,
no body force).

Stalks: 1..nv velocity nodes (Cofree), then per Gauss point q the pair
nv + 2q - 1 holding (s, q₁, q₂) (SecondOrderCone) and nv + 2q holding
(t, σ, ρ) (PowerCone(1/p)). Row blocks per Gauss point: a 3-row physics block
(strain rows q - Ĝv = 0 and the incompressibility row ĥv = 0) and a 2-row
chain block (ρ - s = 0 and σ = 1).
"""
function build_powerlaw(N::Int; K::Float64 = 1.0, p::Float64 = 1.5, V::Float64 = 1.0)
    m = unitsquare(N)
    lay = cavity_layout(m; V)
    ng, NE = length(GAUSS), nelems(m)
    nv = lay.nvstalk
    nq = ng * NE

    socstalk(qi) = nv + 2 * qi - 1
    powstalk(qi) = nv + 2 * qi

    coldims = vcat(fill(0, nv), fill(3, 2 * nq))
    for i in 1:nnodes(m)
        lay.stalk[i] > 0 && (coldims[lay.stalk[i]] = length(lay.comps[i]))
    end

    rowdims = Vector{Int}(undef, 2 * nq)
    rowdims[1:2:end] .= 3
    rowdims[2:2:end] .= 2

    Brow, Bcol, Bblk = Int[], Int[], Matrix{Float64}[]
    g = zeros(sum(rowdims))
    c = zeros(sum(coldims))
    coff = cumsum(vcat(0, coldims))

    # the q-side of the strain rows; the divergence row has no cone entries
    Ssoc = zeros(3, 3)
    Ssoc[1:2, 2:3] .= Matrix{Float64}(I, 2, 2)

    # chain block selectors: ρ - s = 0 and σ = 1
    Csoc = [-1.0 0.0 0.0; 0.0 0.0 0.0]
    Cpow = [0.0 0.0 1.0; 0.0 1.0 0.0]

    for e in 1:NE
        el = view(m.elems, :, e)

        for (gi, (wt, Ghat, hhat)) in enumerate(element_ops(m, e))
            qi = (e - 1) * ng + gi
            r = 2 * qi - 1                     # physics block row
            goff = 5 * (qi - 1)                # scalar offset of this point's rows

            push!(Brow, r); push!(Bcol, socstalk(qi)); push!(Bblk, copy(Ssoc))
            push!(Brow, r + 1); push!(Bcol, socstalk(qi)); push!(Bblk, copy(Csoc))
            push!(Brow, r + 1); push!(Bcol, powstalk(qi)); push!(Bblk, copy(Cpow))

            c[coff[powstalk(qi)] + 1] = -wt * K / p
            g[goff + 5] = 1.0                  # σ = 1

            for (k, i) in enumerate(el)
                cols = localcols(k)

                if lay.stalk[i] > 0
                    blk = zeros(3, length(lay.comps[i]))
                    blk[1:2, :] .= .-view(Ghat, :, cols[lay.comps[i]])
                    blk[3, :] .= view(hhat, cols[lay.comps[i]])
                    push!(Brow, r); push!(Bcol, lay.stalk[i]); push!(Bblk, blk)
                end

                for j in 1:2
                    j in lay.comps[i] && continue
                    g[goff .+ (1:2)] .+= view(Ghat, :, cols[j]) .* lay.fixed[j, i]
                    g[goff + 3] -= hhat[cols[j]] * lay.fixed[j, i]
                end
            end
        end
    end

    B = blocksparse(Brow, Bcol, Bblk, rowdims, coldims)

    # no quadratic term: the whole rheology lives in the cones
    Qrow, Qcol, Qblk = Int[], Int[], Matrix{Float64}[]

    for s in 1:(nv + 2 * nq)
        push!(Qrow, s); push!(Qcol, s); push!(Qblk, zeros(coldims[s], coldims[s]))
    end

    Q = blocksparse(Qrow, Qcol, Qblk, coldims, coldims)

    cones = IPM.AbstractCone[CofreeCone() for _ in 1:nv]

    for _ in 1:nq
        push!(cones, SecondOrderCone())
        push!(cones, PowerCone(1.0 / p))
    end

    return IPMProblem(Q, B, c, g, cones)
end

"""
    e1plus_point(prob)

The E1+ warm start for a `build_powerlaw` cavity problem, returned as
`(p0, y0)` in original problem coordinates for the `p0`/`y0` solver API.

The t-coordinate of each power cone appears in no constraint row, so its
dual is pinned to the objective coefficient, dₜ = a = |fₜ|, at every
dual-feasible point. E1+ places each coupled Gauss point on that pin and
presets the chain-row multipliers so the initial dual residual is zero
and the primal residual is supported only on the physics rows. Per Gauss
point, with α the cone exponent and a the |t-coordinate objective coeff|:

    pow (t, σ, ρ):   t = (1 + 3α) / a,   σ = 1,   ρ = tᵃ / √2
    soc (s, q):      s = ρ,              q = 0
    y chain rows:    ρ − s row = 2/ρ,    σ row = −(4 − 3α)

The dual dₜ = a, d_ρ = −2/ρ, d_σ = 4 − 3α come from the barrier gradient;
the multipliers are their negatives, so d + f + Bᵀy = 0 exactly.
"""
function e1plus_point(prob)
    K = prob.K
    nv = count(c -> c isa CofreeCone, K)
    nq = (length(K) - nv) ÷ 2
    p0 = zeros(size(prob.B, 2))
    y0 = zeros(size(prob.B, 1))

    for qi in 1:nq
        soc = nv + 2 * qi - 1
        pow = nv + 2 * qi
        rpow = colrange(prob.B, pow)
        rsoc = colrange(prob.B, soc)
        α = K[pow].α
        a = abs(prob.f[first(rpow)])

        t = (1 + 3α) / a
        ρ = t^α / sqrt(2)

        p0[rpow[1]] = t                        # t
        p0[rpow[2]] = 1.0                      # σ
        p0[rpow[3]] = ρ                        # ρ
        p0[rsoc[1]] = ρ                        # s = ρ, q = 0

        goff = 5 * (qi - 1)
        y0[goff + 4] = 2 / ρ                   # ρ − s row
        y0[goff + 5] = -(4 - 3α)               # σ row
    end

    return p0, y0
end

# time an E1+-warm-started IPM solve, mirroring measure_ipm
function measure_ipm_e1plus(prob, settings)
    p0, y0 = e1plus_point(prob)
    res = solve(prob, settings; p0, y0)

    if res.status != OPTIMAL && res.status != NEAR_OPTIMAL
        return (t = NaN, status = string(res.status), obj = NaN)
    end

    t = @belapsed solve($prob, $settings; p0 = $p0, y0 = $y0)
    return (t = t, status = string(res.status), obj = res.pobj)
end

function build_powerlaw_jump(N::Int; K::Float64 = 1.0, p::Float64 = 1.5,
                             V::Float64 = 1.0, optimizer)
    m = unitsquare(N)
    lay = cavity_layout(m; V)
    ng, NE = length(GAUSS), nelems(m)
    nq = ng * NE

    model = Model(optimizer)
    @variable(model, v[1:2, 1:nnodes(m)])
    @variable(model, s[1:nq])
    @variable(model, t[1:nq])

    for i in 1:nnodes(m), j in 1:2
        j in lay.comps[i] || fix(v[j, i], lay.fixed[j, i]; force = true)
    end

    for e in 1:NE
        el = view(m.elems, :, e)

        for (gi, (wt, Ghat, hhat)) in enumerate(element_ops(m, e))
            qi = (e - 1) * ng + gi
            vel = [v[j, el[k]] for k in 1:6 for j in 1:2]

            q1 = sum(Ghat[1, a] * vel[a] for a in 1:12)
            q2 = sum(Ghat[2, a] * vel[a] for a in 1:12)

            @constraint(model, sum(hhat[a] * vel[a] for a in 1:12) == 0)
            @constraint(model, [s[qi], q1, q2] in JuMP.SecondOrderCone())
            @constraint(model, [t[qi], 1.0, s[qi]] in MOI.PowerCone(1.0 / p))

            set_objective_coefficient(model, t[qi], wt * K / p)
        end
    end

    @objective(model, Min, objective_function(model))
    return model
end

# -----------------------------------------------------------------------------
# Benchmark — lid-driven cavity
# -----------------------------------------------------------------------------

function benchmark(; ps = (1.5, 3.0), Ns = (8, 16, 32), K = 1.0, tol = 1.0e-5,
                     clarabel = false, mosek = false)
    @printf("%-6s %6s %6s %4s %9s %9s %9s %9s %8s %8s\n",
            "N", "nv", "cones", "p", "IPM+E1", "HSD", "Clarabel", "Mosek", "Cla/IPM", "Mos/IPM")

    for N in Ns, p in ps
        prob = build_powerlaw(N; K, p)

        mi = measure_ipm_e1plus(prob, IPMSettings{Float64}(feas_tol = tol, gap_tol = tol))
        mh = measure_ipm(prob, HSDSettings{Float64}(feas_tol = tol, gap_tol = tol))

        if clarabel
            mc = measure_jump(() -> build_powerlaw_jump(N; K, p, optimizer = clarabel_opt(; tol)))
        else
            mc = (t = NaN, status = "—", obj = NaN)
        end

        if mosek
            mk = measure_jump(() -> build_powerlaw_jump(N; K, p, optimizer = mosek_opt(; tol)))
        else
            mk = (t = NaN, status = "—", obj = NaN)
        end

        cratio = (isfinite(mc.t) && isfinite(mi.t)) ? @sprintf("%.2f", mc.t / mi.t) : "—"
        mratio = (isfinite(mk.t) && isfinite(mi.t)) ? @sprintf("%.2f", mk.t / mi.t) : "—"
        ncones = 2 * 3 * 2 * N * N

        @printf("%-6d %6d %6d %4.1f %s %s %s %s %8s %8s\n",
                N, 2 * (N - 1)^2, ncones, p,
                fmt_time(mi), fmt_time(mh), fmt_time(mc), fmt_time(mk), cratio, mratio)
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    opts = parse_args()
    tol = any(a -> startswith(a, "--tol="), ARGS) ? opts.tol : 1.0e-5
    benchmark(; clarabel = opts.clarabel, mosek = opts.mosek, tol)
end

# =============================================================================
# Sample run: 2026-08-24 (Apple M-series, -t 8, --clarabel --mosek), tol = 1e-5.
# IPM+E1 = E1+ warm start; HSD "—" = not converged; * = non-optimal solve.
#
#   N          nv  cones    p    IPM+E1       HSD  Clarabel     Mosek  Cla/IPM  Mos/IPM
#   8          98    768  1.5   26.9ms   34.0ms   26.2ms   33.6ms     0.97     1.25
#   8          98    768  3.0   31.9ms   58.5ms   34.7ms   38.3ms     1.08     1.20
#   16        450   3072  1.5  100.7ms —         197.1ms  261.5ms     1.96     2.60
#   16        450   3072  3.0  112.4ms  279.9ms 103.1ms*  281.3ms     0.92     2.50
#   32       1922  12288  1.5  448.7ms  620.2ms 3518.6ms* 3248.0ms     7.84     7.24
#   32       1922  12288  3.0  670.7ms 1712.3ms 676.5ms* 3481.6ms     1.01     5.19
# =============================================================================
