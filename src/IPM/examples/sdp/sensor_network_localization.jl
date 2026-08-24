# =============================================================================
#
# # Sensor Network Localization
#
# The sensor network localizaton problem seeks to discover the locations
#
#   X = [ x₁, …, xₙ ]                                    (1)
#
# of n sensors given the locations
#
#   A = [ a₁, …, aₘ ]                                    (2)
#
# of m anchors. Radio links certain sensors to sensors and sensors to
# anchors. If sensors i, j are in radio range, then we write (i, j) ∈ Nₓ;
# if a sensor i and an anchor k are in radio range, then we write (i, k) ∈ Nₐ.
# If two objects are in radio range, then we know their parwise distance.
#                
#   ‖ xᵢ - xⱼ ‖² = dᵢⱼ²   (i, j) ∈ Nx                    (3)
#   ‖ xᵢ - aₖ ‖² = eᵢₖ²   (i, k) ∈ Na
#                
# Recovering the exact locations (1) from (2) and (3) is NP-hard, but they
# can be approximated by solving an easier problem called the Biswas-Ye
# semidefinite relaxation:
#
#   find X, Y
#
#   s.t. Yᵢᵢ - 2Yᵢⱼ + Yⱼⱼ       = dᵢⱼ²,    (i, j) ∈ Nₓ   (5)
#        Yᵢᵢ - 2aₖᵀxᵢ + ‖ aₖ ‖² = eᵢₖ²,    (i, k) ∈ Nₐ
#
#        [ I  X ] ≥ 0
#        [ Xᵀ Y ]
#
# The ESDP relaxation is even easier:
#
#   find X, Y
#
#   s.t. Yᵢᵢ - 2Yᵢⱼ + Yⱼⱼ       = dᵢⱼ²,    (i, j) ∈ Nₓ   (6)
#        Yᵢᵢ - 2aₖᵀxᵢ + ‖ aₖ ‖² = eᵢₖ²,    (i, k) ∈ Nₐ
#
#        [ I    xᵢ   xⱼ  ]
#        [ xᵢᵀ  Yᵢᵢ  Yᵢⱼ ] ≥ 0,            (i, j) ∈ Nₓ
#        [ xⱼᵀ  Yᵢⱼ  Yⱼⱼ ]
#
# If the distances (3) are noisy, then the equalities in (6) may be relaxed
# by introducing slack variables uᵢⱼ and uⱼₖ for all (i, j) ∈ Nₓ and
# (i, k) ∈ Nₐ:
#
#   min  ∑ (uᵢⱼ + vᵢⱼ) + ∑ (uᵢₖ + vᵢₖ)                   (7)
#
#   s.t. Yᵢᵢ - 2Yᵢⱼ + Yⱼⱼ       - uᵢⱼ + vᵢⱼ = dᵢⱼ²,   (i, j) ∈ Nₓ
#        Yᵢᵢ - 2aₖᵀxᵢ + ‖ aₖ ‖² - uᵢₖ + vᵢₖ = eᵢₖ²,   (i, k) ∈ Nₐ
#
#        [ I    xᵢ   xⱼ  ]
#        [ xᵢᵀ  Yᵢᵢ  Yᵢⱼ ] ≥ 0,                       (i, j) ∈ Nₓ
#        [ xⱼᵀ  Yᵢⱼ  Yⱼⱼ ]
#
#        u, v ≥ 0
#
# # References
#
#   - Wang, Zheng, Ye & Boyd,
#     Further relaxations of the semidefinite programming approach to sensor
#     network localization,
#     SIAM J. Optim. 19(2):655-673, 2008.
#
# =============================================================================

using AppleAccelerate
using LinearAlgebra, Printf, Random
using Dualization

using CellularSheaves.IPM
using CellularSheaves.IPM: CofreeCone, PositiveCone, SemidefiniteCone,
                           IPMSettings, HSDSettings
using CellularSheaves.BlockSparseArrays: colrange, blocksparse

include("../utils.jl")

# Section 5.1 of the reference caps the number of selected edges at each sensor
# or anchor to keep the graph sparse. The paper does not say whether the cap is
# applied globally or per node; we select up to DEG_CAP candidates at each point
# and take the union, which leaves degrees near the cap without starving a point
# whose neighbours happen to saturate first.
const DEG_CAP = 7

const R2 = sqrt(2.0)

# -----------------------------------------------------------------------------
# Instances
# -----------------------------------------------------------------------------

struct SNLInstance
    xs::Matrix{Float64}                 # 2 × ns  true sensor positions
    as::Matrix{Float64}                 # 2 × na  anchor positions
    Nx::Vector{Tuple{Int, Int}}         # sensor-sensor edges
    Na::Vector{Tuple{Int, Int}}         # sensor-anchor edges
    dx::Vector{Float64}                 # measured distances, sensor-sensor
    da::Vector{Float64}                 # measured distances, sensor-anchor
end

nsensors(inst::SNLInstance) = size(inst.xs, 2)

# Points within radio range of each other, found by bucketing into rd-sized
# cells so that generation stays linear in the number of points.
function candidates(pts::Matrix{Float64}, rd::Float64)
    n = size(pts, 2)
    cell = max(rd, eps())
    buckets = Dict{Tuple{Int, Int}, Vector{Int}}()

    for i in 1:n
        key = (floor(Int, (pts[1, i] + 0.5) / cell), floor(Int, (pts[2, i] + 0.5) / cell))
        push!(get!(buckets, key, Int[]), i)
    end

    cand = [Int[] for _ in 1:n]

    for (key, members) in buckets
        neigh = Int[]

        for di in -1:1, dj in -1:1
            append!(neigh, get(buckets, (key[1] + di, key[2] + dj), Int[]))
        end

        for i in members, j in neigh
            i == j && continue

            if hypot(pts[1, i] - pts[1, j], pts[2, i] - pts[2, j]) <= rd
                push!(cand[i], j)
            end
        end
    end

    return cand
end

# n points uniform on [-1/2, 1/2]², m of them anchors, radio range rd, and a
# multiplicative measurement noise d = d̄ (1 + nf · randn), as in section 5.1.
# Points with no neighbour within radio range carry no measurement at all: their
# stalk would have an empty column in B and hence a singular block in Q + BᵀB,
# so they are dropped rather than left to the regularization.
function snl_instance(n::Int, m::Int, rd::Float64; nf::Float64 = 0.1, seed::Int = 0)
    rng = Xoshiro(seed)

    pts = rand(rng, 2, n) .- 0.5
    isanchor = falses(n)
    isanchor[randperm(rng, n)[1:m]] .= true

    cand = candidates(pts, rd)
    sel = Set{Tuple{Int, Int}}()

    for i in 1:n
        ci = cand[i]
        isempty(ci) && continue

        for t in randperm(rng, length(ci))[1:min(DEG_CAP, length(ci))]
            j = ci[t]
            isanchor[i] && isanchor[j] && continue      # anchor-anchor is useless
            push!(sel, (min(i, j), max(i, j)))
        end
    end

    # keep the non-anchor points that carry at least one measurement
    touched = falses(n)

    for (i, j) in sel
        touched[i] = true
        touched[j] = true
    end

    sid = zeros(Int, n)
    aid = zeros(Int, n)
    ns = na = 0

    for i in 1:n
        if isanchor[i]
            aid[i] = (na += 1)
        elseif touched[i]
            sid[i] = (ns += 1)
        end
    end

    xs = zeros(2, ns)
    as = zeros(2, na)

    for i in 1:n
        if isanchor[i]
            as[:, aid[i]] .= @view pts[:, i]
        elseif touched[i]
            xs[:, sid[i]] .= @view pts[:, i]
        end
    end

    Nx = Tuple{Int, Int}[]
    Na = Tuple{Int, Int}[]
    dx = Float64[]
    da = Float64[]

    noisy(d) = nf > 0 ? d * abs(1 + nf * randn(rng)) : d

    for (i, j) in sort!(collect(sel))
        d = hypot(pts[1, i] - pts[1, j], pts[2, i] - pts[2, j])

        if isanchor[i] || isanchor[j]
            s, a = isanchor[i] ? (j, i) : (i, j)
            push!(Na, (sid[s], aid[a]))
            push!(da, noisy(d))
        else
            push!(Nx, (sid[i], sid[j]))
            push!(dx, noisy(d))
        end
    end

    return SNLInstance(xs, as, Nx, Na, dx, da)
end

# -----------------------------------------------------------------------------
# Build
# -----------------------------------------------------------------------------

# Index of the (i, j) entry of a d × d symmetric matrix in the svec layout used
# by SemidefiniteCone: the lower triangle stored column by column, the diagonal
# unscaled and the off-diagonal scaled by √2.
#
#   (1,1), (2,1), …, (d,1), (2,2), (3,2), …, (d,d)
#
function svecidx(d::Int, i::Int, j::Int)
    i, j = max(i, j), min(i, j)
    return (j - 1) * d - ((j - 1) * (j - 2)) ÷ 2 + (i - j + 1)
end

# The stalks are numbered
#
#   1 … ns                    sensor    (xᵢ, Yᵢᵢ)          cofree,       width 3
#   ns + 1 … ns + nE          edge      svec Z₍₁,₂,ᵢ,ⱼ₎    semidefinite, width 10
#   ns + nE + 1 … ns + 2nE    slack     (uᵢⱼ, vᵢⱼ)        positive,     width 2
#   ns + 2nE + 1 … + nA       slack     (uᵢₖ, vᵢₖ)        positive,     width 2
#
# and within the edge stalk the 4 × 4 matrix is indexed (1, 2, i, j), so row 3
# belongs to endpoint i and row 4 to endpoint j.
function build_snl(inst::SNLInstance)
    ns = nsensors(inst)
    nE = length(inst.Nx)
    nA = length(inst.Na)

    S(i) = i
    P(e) = ns + e
    T(e) = ns + nE + e                  # e indexes Nx, then Na

    row_ids = Int[]
    col_ids = Int[]
    blks = Matrix{Float64}[]
    g = Float64[]

    r = 0

    for (e, (i, j)) in enumerate(inst.Nx)
        # -- the identity block: Z(1,1) = 1, Z(2,2) = 1, Z(2,1) = 0
        r += 1
        Bp = zeros(3, 10)
        Bp[1, svecidx(4, 1, 1)] = 1
        Bp[2, svecidx(4, 2, 2)] = 1
        Bp[3, svecidx(4, 2, 1)] = 1
        push!(row_ids, r); push!(col_ids, P(e)); push!(blks, Bp)
        append!(g, (1.0, 1.0, 0.0))

        # -- endpoint i: Z(3,1) = xᵢ₁, Z(3,2) = xᵢ₂, Z(3,3) = Yᵢᵢ
        #    the off-diagonal svec entries carry a √2, the diagonal does not
        r += 1
        Bp = zeros(3, 10)
        Bp[1, svecidx(4, 3, 1)] = 1
        Bp[2, svecidx(4, 3, 2)] = 1
        Bp[3, svecidx(4, 3, 3)] = 1
        push!(row_ids, r); push!(col_ids, P(e)); push!(blks, Bp)
        push!(row_ids, r); push!(col_ids, S(i)); push!(blks, [-R2 0 0; 0 -R2 0; 0 0 -1.0])
        append!(g, (0.0, 0.0, 0.0))

        # -- endpoint j: Z(4,1) = xⱼ₁, Z(4,2) = xⱼ₂, Z(4,4) = Yⱼⱼ
        r += 1
        Bp = zeros(3, 10)
        Bp[1, svecidx(4, 4, 1)] = 1
        Bp[2, svecidx(4, 4, 2)] = 1
        Bp[3, svecidx(4, 4, 4)] = 1
        push!(row_ids, r); push!(col_ids, P(e)); push!(blks, Bp)
        push!(row_ids, r); push!(col_ids, S(j)); push!(blks, [-R2 0 0; 0 -R2 0; 0 0 -1.0])
        append!(g, (0.0, 0.0, 0.0))

        # -- distance: Yᵢᵢ - 2Yᵢⱼ + Yⱼⱼ - u + v = d²,  and Yᵢⱼ = Z(4,3) / √2
        r += 1
        Bp = zeros(1, 10)
        Bp[1, svecidx(4, 3, 3)] = 1
        Bp[1, svecidx(4, 4, 3)] = -R2
        Bp[1, svecidx(4, 4, 4)] = 1
        push!(row_ids, r); push!(col_ids, P(e)); push!(blks, Bp)
        push!(row_ids, r); push!(col_ids, T(e)); push!(blks, [-1.0 1.0])
        push!(g, inst.dx[e]^2)
    end

    for (a, (i, k)) in enumerate(inst.Na)
        # -- anchor distance: Yᵢᵢ - 2aₖᵀxᵢ - u + v = d̄² - ‖ aₖ ‖²
        ak = @view inst.as[:, k]

        r += 1
        push!(row_ids, r); push!(col_ids, S(i)); push!(blks, [-2ak[1] -2ak[2] 1.0])
        push!(row_ids, r); push!(col_ids, T(nE + a)); push!(blks, [-1.0 1.0])
        push!(g, inst.da[a]^2 - (ak[1]^2 + ak[2]^2))
    end

    B = blocksparse(row_ids, col_ids, blks)

    # Q = 0: the objective ∑ (u + v) is linear
    Q = IPM.allocblockdiag(B)
    fill!(Q, 0)

    # min ½ pᵀQp - cᵀp with c = -1 on every slack coordinate
    c = zeros(size(B, 2))

    for e in 1:(nE + nA)
        c[colrange(B, T(e))] .= -1.0
    end

    cones = IPM.AbstractCone[CofreeCone() for _ in 1:ns]
    append!(cones, IPM.AbstractCone[SemidefiniteCone() for _ in 1:nE])
    append!(cones, IPM.AbstractCone[PositiveCone() for _ in 1:(nE + nA)])

    return IPMProblem(Q, B, c, g, cones), B
end

function build_snl_jump(inst::SNLInstance; optimizer)
    ns = nsensors(inst)
    model = Model(optimizer)

    @variable(model, x[1:2, 1:ns])
    @variable(model, Y[1:ns])

    obj = AffExpr(0.0)

    for (e, (i, j)) in enumerate(inst.Nx)
        Z = @variable(model, [1:4, 1:4], PSD)
        u = @variable(model, lower_bound = 0.0)
        v = @variable(model, lower_bound = 0.0)

        @constraint(model, Z[1, 1] == 1)
        @constraint(model, Z[2, 2] == 1)
        @constraint(model, Z[1, 2] == 0)
        @constraint(model, Z[3, 1] == x[1, i])
        @constraint(model, Z[3, 2] == x[2, i])
        @constraint(model, Z[4, 1] == x[1, j])
        @constraint(model, Z[4, 2] == x[2, j])
        @constraint(model, Z[3, 3] == Y[i])
        @constraint(model, Z[4, 4] == Y[j])
        @constraint(model, Z[3, 3] - 2Z[4, 3] + Z[4, 4] - u + v == inst.dx[e]^2)

        add_to_expression!(obj, u)
        add_to_expression!(obj, v)
    end

    for (a, (i, k)) in enumerate(inst.Na)
        ak = @view inst.as[:, k]
        u = @variable(model, lower_bound = 0.0)
        v = @variable(model, lower_bound = 0.0)

        @constraint(model, Y[i] - 2 * (ak[1] * x[1, i] + ak[2] * x[2, i]) - u + v ==
                           inst.da[a]^2 - (ak[1]^2 + ak[2]^2))

        add_to_expression!(obj, u)
        add_to_expression!(obj, v)
    end

    @objective(model, Min, obj)
    return model
end

# -----------------------------------------------------------------------------
# Localization quality
# -----------------------------------------------------------------------------

# The estimated positions and the individual traces (Y - XᵀX)ᵢᵢ, which are the
# reference's per-sensor confidence measure: a sensor with zero individual trace
# is exactly positioned by the relaxation.
function positions(inst::SNLInstance, p::AbstractVector, B)
    ns = nsensors(inst)
    X = zeros(2, ns)
    tr = zeros(ns)

    for i in 1:ns
        cols = colrange(B, i)
        X[1, i] = p[cols[1]]
        X[2, i] = p[cols[2]]
        tr[i] = p[cols[3]] - (X[1, i]^2 + X[2, i]^2)
    end

    return X, tr
end

# RMSD as defined in (4.3) of the reference.
function rmsd(X::AbstractMatrix, Xtrue::AbstractMatrix)
    return sqrt(sum(abs2, X .- Xtrue) / size(Xtrue, 2))
end

function accuracy(n, m, rd; nf = 0.1, tol = 1.0e-6, seed = 0)
    inst = snl_instance(n, m, rd; nf, seed)
    prob, B = build_snl(inst)
    res = solve(prob, IPMSettings{Float64}(feas_tol = tol, gap_tol = tol))
    X, tr = positions(inst, res.p, B)

    return (err = rmsd(X, inst.xs),
            conf = count(<(1.0e-4), tr) / length(tr),
            status = res.status)
end

# -----------------------------------------------------------------------------
# Benchmark
# -----------------------------------------------------------------------------

# The (n, m, rd) triples are the noisy test problems of Table 5.1, followed by
# the two problems of Table 5.2. The radio range shrinks as n grows, which is
# what keeps the average degree — and so the separators of the network — from
# growing with the problem.
const CASES = [(200, 5, 0.25), (400, 10, 0.20), (800, 20, 0.12)]

function benchmark(; nf = 0.1, tol = 1.0e-6, clarabel = false, mosek = false, cases = CASES)
    @printf("%-6s %5s %6s %6s %6s %7s %9s %9s %9s %9s %8s %8s\n",
            "n", "m", "rd", "sens", "|Nx|", "cols",
            "IPM", "HSD", "Clarabel", "Mosek", "Cla/IPM", "Mos/IPM")

    for (n, m, rd) in cases
        inst = snl_instance(n, m, rd; nf)
        prob, _ = build_snl(inst)

        mi = measure_ipm(prob, IPMSettings{Float64}(feas_tol = tol, gap_tol = tol))
        mh = measure_ipm(prob, HSDSettings{Float64}(feas_tol = tol, gap_tol = tol))

        if clarabel
            mc = measure_jump(() -> build_snl_jump(inst; optimizer = clarabel_opt(; tol)))
        else
            mc = (t = NaN, status = "—", obj = NaN)
        end

        if mosek
            mk = measure_jump(() -> build_snl_jump(inst; optimizer = dual_optimizer(mosek_opt(; tol))))
        else
            mk = (t = NaN, status = "—", obj = NaN)
        end

        cratio = (isfinite(mc.t) && isfinite(mi.t)) ? @sprintf("%.2f", mc.t / mi.t) : "—"
        mratio = (isfinite(mk.t) && isfinite(mi.t)) ? @sprintf("%.2f", mk.t / mi.t) : "—"

        @printf("%-6d %5d %6.3f %6d %6d %7d %s %s %s %s %8s %8s\n",
                n, m, rd, nsensors(inst), length(inst.Nx), size(prob.B, 2),
                fmt_time(mi), fmt_time(mh), fmt_time(mc), fmt_time(mk), cratio, mratio)
        flush(stdout)
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    opts = parse_args()
    tol = any(a -> startswith(a, "--tol="), ARGS) ? opts.tol : 1.0e-6
    benchmark(; clarabel = opts.clarabel, mosek = opts.mosek, tol)

    println()
    @printf("%-6s %6s %9s %9s\n", "n", "nf", "RMSD", "trace ok")

    for (n, m, rd) in ((1000, 100, 0.06), (4000, 400, 0.035))
        for nf in (0.0, 0.01, 0.1)
            a = accuracy(n, m, rd; nf, tol)
            @printf("%-6d %6.3f %9.2e %8.0f%%\n", n, nf, a.err, 100 * a.conf)
        end
    end
end

# =============================================================================
# Sample run: 2026-08-24 (Apple M-series, --clarabel --mosek, tol = 1e-6; Mosek via Dualization):
# * = time from a NEAR_OPTIMAL (not fully OPTIMAL) solve.
#
#   n          m     rd   sens   |Nx|    cols       IPM       HSD  Clarabel     Mosek  Cla/IPM  Mos/IPM
#   200        5  0.250    195   1153   14555 136.1ms* 177.0ms*  353.0ms  345.1ms     2.59     2.54
#   400       10  0.200    390   2432   30600  471.7ms 750.8ms* 1048.1ms  868.2ms     2.22     1.84
#   800       20  0.120    780   4722   59472 874.9ms* 1656.2ms* 2235.5ms 2152.3ms     2.56     2.46
#
# =============================================================================
