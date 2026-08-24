# =============================================================================
#
# # The Elastic Obstacle Problem
#
# A 3-D obstacle lies inside of a 2-D region. We model the region as
# a domain Ω ⊆ ℝ² and the height of the obstacle as a function ψ: Ω → ℝ.
# A membrane is stretched over the object and pinned above the boundary
# of the domain at height g: ∂Ω → ℝ. Given a body force f: Ω → ℝ, the
# displacement u: Ω → ℝ of the membrane at equilibrium minimizes the
# Dirichlet energy subject to a one-sided constraint:
#
#  min  1/2 ∫ | ∇u |² - ∫ f u    (1)
#   u
#
#  s.t. u ≥ ψ in  Ω
#       u = g on ∂Ω
#
# This is a variational problem. We turn it into a quadratic program
# using the Galerkin method, partitioning Ω into a mesh with nodes
# {x₁, …, xₙ} and searching over displacements that are flat on each
# element and continuous across shared edges. The set of all such
# displacements is a vector space Vₕ ⊆ H¹(Ω) with basis {φ₁, …, φₙ},
# where φᵢ: Ω → ℝ is the tent function
#
#   φᵢ(xⱼ) = { 1 i = j
#            { 0 i ≠ j
#
# of the ith node xᵢ. Every candidate displacement function uₕ ∈ Vₕ is
# of the form
#
#        n
#   uₕ = ∑ uᵢ φᵢ
#       i=1
#
# for some coefficients u₁, …, uₙ. If we relax the requirement that u ≥ ψ
# in Ω and impose it only on the nodes {x₁, …, xₙ}, then the problem (1)
# turns into an n-dimensional quadratic program
#
#   min  1/2 uᵀ K u - cᵀ u       (2)
#    u
#
#   s.t. uᵢ ≥ ψᵢ     for all nodes xᵢ
#        uᵢ = gᵢ     for all nodes xᵢ ∈ ∂Ω
#
# where
#
#   Kᵢⱼ = ∫ ∇φᵢᵀ ∇φⱼ dx        cᵢ = ∫ f φᵢ dx
#         Ω                          Ω
#
#   ψᵢ  = ψ(xᵢ)                gᵢ = g(xᵢ)
#
# # References
#
#   - Keith & Surowiec,
#     Proximal Galerkin: a structure-preserving finite element method for pointwise bound constraints,
#     Found. Comput. Math. 24(5):1511–1607, 2024.
#
# =============================================================================

using LinearAlgebra, Printf
using Dualization

using CellularSheaves.IPM
using CellularSheaves.IPM: PositiveCone, IPMSettings, HSDSettings
using CellularSheaves.BlockSparseArrays: colrange, blocksparse

include("../utils.jl")

# -----------------------------------------------------------------------------
# Exact solution
# -----------------------------------------------------------------------------

# free boundary radius: the root of  -a² ln a - (1/4 - a²)  on (0, 1/2)
function contact_radius()
    f(r) = -r^2 * log(r) - (0.25 - r^2)
    lo, hi = 1.0e-6, 0.5 - 1.0e-12

    for _ in 1:200
        mid = (lo + hi) / 2
        f(lo) * f(mid) <= 0 ? (hi = mid) : (lo = mid)
    end

    return (lo + hi) / 2
end

const RSTAR = contact_radius()                 # ≈ 0.3489826
const ASTAR = -RSTAR^2 / sqrt(0.25 - RSTAR^2)  # ≈ -0.3401297

# The cap of a sphere of radius r₀ = 1/2, continued along its tangent line
# beyond r = β r₀.  Cutting at β = 0.9 avoids the cap's vertical slope at
# r = r₀, so ψ is C¹ everywhere; matching value and derivative at the join
# gives ψ(r) = B + C r with B = √(r₀²-b²) + b²/√(r₀²-b²) and C = -b/√(r₀²-b²).
# This is `spherical_obstacle` from MFEM ex36 verbatim.  The cone descends
# through zero at r = 5/9 and reaches ≈ -1.77 at the corners of Ω, so it
# never binds outside the contact set r ≤ a ≈ 0.349 — where it agrees with
# the plain cap, which is why the exact solution below is unaffected.
function obstacle(x, y)
    r = sqrt(x^2 + y^2)
    r0, beta = 0.5, 0.9

    b = r0 * beta
    tmp = sqrt(r0^2 - b^2)
    B = tmp + b^2 / tmp
    C = -b / tmp

    return r > b ? B + r * C : sqrt(r0^2 - r^2)
end

function exact(x, y)
    r = sqrt(x^2 + y^2)
    return r > RSTAR ? ASTAR * log(r) : sqrt(max(0.25 - r^2, 0.0))
end

# -----------------------------------------------------------------------------
# Build
# -----------------------------------------------------------------------------

# L × L nodes on [-1, 1]², of which the (L-2)² interior ones are unknowns.
# Returns the interior coordinates, the obstacle, and the linear term c
# obtained by lifting the Dirichlet data into the right hand side.
function grid_data(L::Int)
    h = 2 / (L - 1)
    coord(i) = -1 + (i - 1) * h

    n = (L - 2)^2
    id = zeros(Int, L, L)                     # 0 on the boundary
    k = 0

    for j in 2:(L - 1), i in 2:(L - 1)
        id[i, j] = (k += 1)
    end

    xs = zeros(n)
    ys = zeros(n)
    psi = zeros(n)
    c = zeros(n)

    for j in 2:(L - 1), i in 2:(L - 1)
        v = id[i, j]
        x, y = coord(i), coord(j)
        xs[v], ys[v] = x, y
        psi[v] = obstacle(x, y)

        # cᵢ = ∫ f φᵢ - (K_∂ g)ᵢ.  Here f = 0, so only the boundary lift
        # survives; K_∂ has entry -1 per interior-boundary edge.  With f ≠ 0
        # this would need an extra h^2 * f(x, y), since ∫ φᵢ = h².
        for (di, dj) in ((-1, 0), (1, 0), (0, -1), (0, 1))
            if id[i + di, j + dj] == 0
                c[v] += exact(coord(i + di), coord(j + dj))
            end
        end
    end

    return (; n, id, xs, ys, psi, c)
end

# The eliminated form (3): the slacks sᵢ = uᵢ - ψᵢ ≥ 0 are the only variables,
# Q is the stiffness K, and there is no equality (B is empty). We assemble K
# once — full-symmetric, a diagonal block on every stalk as the union-pattern
# ordering requires — and use its action on ψ for the linear term (Kψ - c).
function build_obstacle(L::Int)
    d = grid_data(L)
    n = d.n

    neg_ = fill(-1.0, 1, 1)

    Qrow = Int[]
    Qcol = Int[]
    Qblk = Matrix{Float64}[]
    Kpsi = zeros(n)                            # (Kψ)ᵢ = 4ψᵢ - Σ_{interior nbrs} ψⱼ

    for v in 1:n
        push!(Qrow, v); push!(Qcol, v); push!(Qblk, fill(4.0, 1, 1))
        Kpsi[v] += 4.0 * d.psi[v]
    end

    for j in 2:(L - 1), i in 2:(L - 1)
        v = d.id[i, j]

        for (di, dj) in ((1, 0), (0, 1))       # each edge once, then mirrored
            w = d.id[i + di, j + dj]
            w == 0 && continue
            push!(Qrow, v); push!(Qcol, w); push!(Qblk, copy(neg_))
            push!(Qrow, w); push!(Qcol, v); push!(Qblk, copy(neg_))
            Kpsi[v] -= d.psi[w]; Kpsi[w] -= d.psi[v]
        end
    end

    Q = blocksparse(Qrow, Qcol, Qblk)

    # (3) is  min ½ sᵀKs + (Kψ - c)ᵀs, s ≥ 0.  In the solver's ½pᵀQp - cₚᵀp
    # form that means cₚ = c - Kψ.
    c = d.c .- Kpsi
    B = blocksparse(Int[], Int[], Matrix{Float64}[], Int[], ones(Int, n))   # no equality

    cones = IPM.AbstractCone[PositiveCone() for _ in 1:n]

    return IPMProblem(Q, B, c, Float64[], cones)
end

function build_obstacle_jump(L::Int; optimizer)
    d = grid_data(L)
    n = d.n

    model = Model(optimizer)
    @variable(model, u[1:n])
    @constraint(model, [v = 1:n], u[v] >= d.psi[v])

    # 1/2 uᵀKu - cᵀu, with K's diagonal 4 and each edge counted once
    obj = 2.0 * sum(u[v]^2 for v in 1:n) - sum(d.c[v] * u[v] for v in 1:n)

    for j in 2:(L - 1), i in 2:(L - 1)
        v = d.id[i, j]

        for (di, dj) in ((1, 0), (0, 1))
            w = d.id[i + di, j + dj]
            w == 0 && continue
            obj += -1.0 * u[v] * u[w]
        end
    end

    @objective(model, Min, obj)
    return model
end

# -----------------------------------------------------------------------------
# Accuracy — the point of using a benchmark with a closed form
# -----------------------------------------------------------------------------

function accuracy(L::Int; tol = 1.0e-8)
    d = grid_data(L)
    prob = build_obstacle(L)
    res = solve(prob, IPMSettings{Float64}(feas_tol = tol, gap_tol = tol))

    err = 0.0
    rad = 0.0

    for v in 1:d.n
        s = res.p[first(colrange(prob.B, v))]   # the eliminated variable is the slack
        u = d.psi[v] + s                          # recover the displacement u = ψ + s
        err = max(err, abs(u - exact(d.xs[v], d.ys[v])))
        # Contact set (s ≈ 0). The contact slacks bottom out at the solve's
        # complementarity floor, which grows with n (≈1e-7 at L=33 up to ≈1e-5 at
        # L=257), so a 1e-6 cutoff catches nothing on fine grids. 1e-3 sits cleanly
        # between the contact slacks and the O(1e-2) free-region slacks.
        s < 1.0e-3 && (rad = max(rad, hypot(d.xs[v], d.ys[v])))
    end

    return (; err, rad, status = res.status)
end

# -----------------------------------------------------------------------------
# Benchmark
# -----------------------------------------------------------------------------

function benchmark(; sides = (33, 65, 129), tol = 1.0e-8,
                     clarabel = false, mosek = false)
    @printf("%-6s %7s %9s %9s %9s %9s %8s %8s\n",
            "L", "n", "IPM", "HSD", "Clarabel", "Mosek", "Cla/IPM", "Mos/IPM")

    for L in sides
        n = (L - 2)^2
        prob = build_obstacle(L)

        mi = measure_ipm(prob, IPMSettings{Float64}(feas_tol = tol, gap_tol = tol))
        mh = measure_ipm(prob, HSDSettings{Float64}(feas_tol = tol, gap_tol = tol))

        if clarabel
            mc = measure_jump(() -> build_obstacle_jump(L; optimizer = clarabel_opt(; tol)))
        else
            mc = (t = NaN, status = "—", obj = NaN)
        end

        if mosek
            mk = measure_jump(() -> build_obstacle_jump(L; optimizer = dual_optimizer(mosek_opt(; tol))))
        else
            mk = (t = NaN, status = "—", obj = NaN)
        end

        cratio = (isfinite(mc.t) && isfinite(mi.t)) ? @sprintf("%.2f", mc.t / mi.t) : "—"
        mratio = (isfinite(mk.t) && isfinite(mi.t)) ? @sprintf("%.2f", mk.t / mi.t) : "—"

        @printf("%-6d %7d %s %s %s %s %8s %8s\n",
                L, n, fmt_time(mi), fmt_time(mh), fmt_time(mc), fmt_time(mk), cratio, mratio)
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    opts = parse_args()
    tol = any(a -> startswith(a, "--tol="), ARGS) ? opts.tol : 1.0e-8
    benchmark(; clarabel = opts.clarabel, mosek = opts.mosek, tol)

    println()
    @printf("%-6s %13s %13s\n", "L", "‖u-uₕ‖∞", "free bdry")

    for L in (33, 65, 129)
        a = accuracy(L; tol)
        @printf("%-6d %13.3e %13.5f\n", L, a.err, a.rad)
    end

    @printf("%-6s %13s %13.5f\n", "exact", "", RSTAR)
end

# =============================================================================
# Sample run: 2026-08-24 (--clarabel --mosek), tol = 1e-8.
#
#   L            n       IPM       HSD  Clarabel     Mosek  Cla/IPM  Mos/IPM
#   33         961    5.1ms    5.7ms    6.3ms   32.3ms     1.25     6.34
#   65        3969   26.5ms   26.1ms   36.7ms  191.5ms     1.38     7.23
#   129      16129  137.9ms  123.3ms  252.6ms 1121.5ms     1.83     8.13
# =============================================================================
