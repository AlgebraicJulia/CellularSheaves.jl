# =============================================================================
#
# # Constrained Kalman Smoothing
#
# Let
#
#   V = [ v₁, …, vₙ ]
#   W = [ w₁, …, wₙ ]
#   X = [ x₁, …, xₙ ]
#   Z = [ z₁, …, zₙ ]
#
# be random matrices distributed according to
#
#   xᵢ = Gᵢ xᵢ₋₁ + wᵢ,    wᵢ ~ N(0, Qᵢ)                                               (1)
#   zᵢ = Hᵢ xᵢ   + vᵢ,    vᵢ ~ N(0, Rᵢ)
#
# with w₁, …, wₙ and v₁, …, vₙ all mutually independent.
# This is called a Kalman-Bucy model; given a matrix Z
# of observations, and the maximum-likelihood estimate for
# the state matrix X can be obtained using the Kalman-Bucy
# smoother, which solves the unconstrained quadratic programming
# problem
#
#            n
#   min  1/2 ∑ (zᵢ - Hᵢ xᵢ)ᵀ Rᵢ⁻¹(zᵢ - Hᵢ xᵢ) + (xᵢ - Gᵢ xᵢ₋₁)ᵀ Qᵢ⁻¹(xᵢ - Gᵢ xᵢ₋₁)    (2)
#    X      i=1
#
# Appending inequality constraints to (2) yields the constrained
# Kalman-Bucy smoother
#
#            n
#   min  1/2 ∑ (zᵢ - Hᵢ xᵢ)ᵀ Rᵢ⁻¹(zᵢ - Hᵢ xᵢ) + (xᵢ - Gᵢ xᵢ₋₁)ᵀ Qᵢ⁻¹(xᵢ - Gᵢ xᵢ₋₁)    (3)
#    X      i=1
#
#   s.t. Fᵢ xᵢ + sᵢ = bᵢ,  i = 1, …, n
#                sᵢ ≥ 0,   i = 1, …, n
#
# which can be motivated by modeling the variables w₁, …, wₙ as *truncated*
# normal vectors.
#
# # References
#
#  - Bell, Burke & Pillonetto,
#    An inequality constrained nonlinear Kalman-Bucy smoother by interior
#    point likelihood maximization,
#    Automatica 45(1), 2009, 25-33.
#
# =============================================================================

using AppleAccelerate
using LinearAlgebra, Printf, Random

using CellularSheaves.IPM
using CellularSheaves.IPM: CofreeCone, PositiveCone, IPMSettings, HSDSettings
using CellularSheaves.BlockSparseArrays: colrange, blocksparse

include("../utils.jl")

# -----------------------------------------------------------------------------
# Model
# -----------------------------------------------------------------------------

# The measurement functional is scalar here (Hₖ is 1 × 2), and Gₖ, Wₖ, Rₖ are
# constant in k for k ≥ 2. The initial covariance W₁ is taken large so that the
# initial state estimate has no appreciable effect on the fit.
struct SplineModel
    N::Int
    dt::Float64
    sigma::Float64
    G::Matrix{Float64}          # Gₖ, k ≥ 2
    W::Matrix{Float64}          # Wₖ, k ≥ 2
    W1::Matrix{Float64}         # W₁
    H::Matrix{Float64}          # Hₖ
    x0::Vector{Float64}         # g₁(x₀), the initial state estimate
    z::Vector{Float64}          # measurements
    xtrue::Matrix{Float64}      # 2 × N, for reporting only
end

function spline_model(N::Int; sigma::Float64 = 0.5, seed::Int = 0)
    rng = Xoshiro(seed)
    dt = 2pi / N
    t = [k * dt for k in 1:N]

    G = [1.0 0.0; dt 1.0]
    W = [dt dt^2/2; dt^2/2 dt^3/3]
    W1 = 100.0 * Matrix{Float64}(I, 2, 2)
    H = [0.0 1.0]

    xtrue = vcat((-cos.(t))', (-sin.(t))')
    x0 = xtrue[:, 1]
    z = xtrue[2, :] .+ sigma .* randn(rng, N)

    return SplineModel(N, dt, sigma, G, W, W1, H, x0, z, xtrue)
end

# -----------------------------------------------------------------------------
# Build
# -----------------------------------------------------------------------------

function build_kalman_smoother(m::SplineModel; bound::Float64 = 1.0)
    N = m.N
    n = size(m.G, 1)
    l = 2n                                  # two bounds per state component

    Winv = inv(m.W)
    W1inv = inv(m.W1)
    Rinv = 1 / m.sigma^2
    HtRH = m.H' * m.H * Rinv

    # stalks: 1..N are the states xₖ, then N+k is the slack vector sₖ

    row_ids = Int[]
    col_ids = Int[]
    blks = Matrix{Float64}[]

    # Q: the block tridiagonal C of equation (2), stored in full (both Aₖ and
    # Aₖᵀ), zero on the slack stalks
    for k in 1:N
        Dk = (k == 1 ? W1inv : Winv) + HtRH

        if k < N
            Dk += m.G' * Winv * m.G         # the G_{k+1} term; G_{N+1} = 0
        end

        push!(row_ids, k); push!(col_ids, k); push!(blks, Dk)

        if k > 1
            Ak = -Winv * m.G
            push!(row_ids, k);     push!(col_ids, k - 1); push!(blks, Matrix(Ak))
            push!(row_ids, k - 1); push!(col_ids, k);     push!(blks, Matrix(Ak'))
        end
    end

    # every stalk needs a diagonal block, including the slacks
    for k in 1:N
        push!(row_ids, N + k); push!(col_ids, N + k); push!(blks, zeros(l, l))
    end

    Q = blocksparse(row_ids, col_ids, blks)

    # B: one row block per time point, Fₖ xₖ + sₖ = bₖ with the bound rows
    #
    #   -x₁ₖ ≤ 1,  x₁ₖ ≤ 1,  -x₂ₖ ≤ 1,  x₂ₖ ≤ 1
    F = zeros(l, n)

    for i in 1:n
        F[2i - 1, i] = -1.0
        F[2i, i] = 1.0
    end

    Il = Matrix{Float64}(I, l, l)

    brow_ids = Int[]
    bcol_ids = Int[]
    bblks = Matrix{Float64}[]

    for k in 1:N
        push!(brow_ids, k); push!(bcol_ids, k);     push!(bblks, copy(F))
        push!(brow_ids, k); push!(bcol_ids, N + k); push!(bblks, copy(Il))
    end

    B = blocksparse(brow_ids, bcol_ids, bblks)
    g = fill(bound, size(B, 1))

    # c: rₖ on the state stalks, zero on the slacks
    c = zeros(size(B, 2))

    for k in 1:N
        rk = (m.H' * [m.z[k]]) * Rinv

        if k == 1
            rk += W1inv * m.x0
        end

        c[colrange(B, k)] .= rk
    end

    cones = vcat(IPM.AbstractCone[CofreeCone() for _ in 1:N],
                 IPM.AbstractCone[PositiveCone() for _ in 1:N])

    return IPMProblem(Q, B, c, g, cones)
end

function build_ks_jump(m::SplineModel; bound::Float64, optimizer)
    N = m.N
    n = size(m.G, 1)

    Winv = inv(m.W)
    W1inv = inv(m.W1)
    Rinv = 1 / m.sigma^2

    model = Model(optimizer)
    @variable(model, -bound <= x[1:n, 1:N] <= bound)

    obj = QuadExpr()

    for k in 1:N
        resid = m.z[k] - sum(m.H[1, a] * x[a, k] for a in 1:n)
        add_to_expression!(obj, (0.5 * Rinv) * resid, resid)

        if k == 1
            r = [x[a, 1] - m.x0[a] for a in 1:n]
            M = W1inv
        else
            r = [x[a, k] - sum(m.G[a, b] * x[b, k - 1] for b in 1:n) for a in 1:n]
            M = Winv
        end

        for a in 1:n, b in 1:n
            add_to_expression!(obj, (0.5 * M[a, b]) * r[a], r[b])
        end
    end

    @objective(model, Min, obj)
    return model
end

# -----------------------------------------------------------------------------
# Benchmark
# -----------------------------------------------------------------------------

function benchmark(; sizes = (50, 200, 800, 3200), sigma = 0.5, bound = 1.0,
                     tol = 1e-6, clarabel = false, mosek = false)
    @printf("%-9s %6s %9s %9s %9s %9s %8s %8s\n",
            "model", "N", "IPM", "HSD", "Clarabel", "Mosek", "Cla/IPM", "Mos/IPM")

    for N in sizes
        m = spline_model(N; sigma)
        prob = build_kalman_smoother(m; bound)

        mi = measure_ipm(prob, IPMSettings{Float64}(feas_tol = tol, gap_tol = tol))
        mh = measure_ipm(prob, HSDSettings{Float64}(feas_tol = tol, gap_tol = tol))

        if clarabel
            mc = measure_jump(() -> build_ks_jump(m; bound, optimizer = clarabel_opt(; tol)))
        else
            mc = (t = NaN, status = "—", obj = NaN)
        end

        if mosek
            mk = measure_jump(() -> build_ks_jump(m; bound, optimizer = mosek_opt(; tol)))
        else
            mk = (t = NaN, status = "—", obj = NaN)
        end

        cratio = (isfinite(mc.t) && isfinite(mi.t)) ? @sprintf("%.2f", mc.t / mi.t) : "—"
        mratio = (isfinite(mk.t) && isfinite(mi.t)) ? @sprintf("%.2f", mk.t / mi.t) : "—"

        @printf("%-9s %6d %s %s %s %s %8s %8s\n",
                "spline", N, fmt_time(mi), fmt_time(mh), fmt_time(mc), fmt_time(mk), cratio, mratio)
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    opts = parse_args()
    tol = any(a -> startswith(a, "--tol="), ARGS) ? opts.tol : 1e-6
    benchmark(; clarabel = opts.clarabel, mosek = opts.mosek, tol)
end

# =============================================================================
# Sample run: 2026-08-24 (--clarabel --mosek), tol = 1e-6.
#
#   model          N       IPM       HSD  Clarabel     Mosek  Cla/IPM  Mos/IPM
#   spline        50     0.7ms    0.8ms    0.3ms    0.9ms     0.51     1.29
#   spline       200     3.5ms    3.7ms    1.8ms    4.4ms     0.52     1.26
#   spline       800    12.1ms   14.9ms    9.4ms   18.5ms     0.78     1.53
#   spline      3200    67.4ms   69.4ms  147.2ms   76.9ms     2.18     1.14
# =============================================================================
