# =============================================================================
#
# # Private Marginal Reconstruction
#
# A sensitive dataset
#
#   x: {1, …, m} → X₁ × ⋯ × Xₖ,   |Xᵢ| = nᵢ
#
# has records x₁, …, xₘ, Writing
#
#   n := n₁ ⋯ nₖ,
#
# the data vector p ∈ ℝⁿ counts the occurrences of each possible record.
# The vector p is enormous, so it is impractical for a database to
# deliver it to users. Instead, users submit queries, requesting smaller
# vectors called marginals. For a subset of attributes
#
#   r ⊆ {1, …, k},
#
# the r-marginal of p is a vector of length
#
#   nᵣ = ∏  nᵢ,
#       i∈r
#
# computed by summing out the indices i ∉ r. This is a linear
# transformation, and we write μᵣ = Mᵣ p.
#
# ## Residuals
#
# The vector space ℝⁿ splits as a sum of orthogonal subspaces
#
#   ℝⁿ = ⊕ Vᵣ,     r ⊆ {1, …, k},
#
# indexed by attribute subsets. Writing
#
#   Rᵣ: ℝⁿ → Vᵣ
#
# for the orthogonal projection, the vector αᵣ = Rᵣ p is called the
# r-residual of p. Marginals can be expressed as sums of residuals:
#
#   μᵣ = ∑ Sᵣₜᵀ αₜ.                               (1)
#       t⊆r
#
# where
#
#   Sᵣₜ: ℝⁿʳ → Vₜ
#
# is the projection onto Vₜ. 
#
# ## Privacy
#
# A database that answers queries exactly leaks its recors. Differentially
# private databases return noisy marginals, computed using (1) from noisy
# residuals
#
#   βᵣ = αᵣ + ϵ,   ϵ ~ N(0, I).
#
# Unlike exact marginals, noisy marginals can have negative entries. The
# weighted-least squares problem
#
#   min  ∑ (κₜ - βₜ)ᵀ Kₜ⁻¹ (κₜ - βₜ)               (2)
#    κ   t
#
#   s.t.  ∑ Sᵣₜᵀ κₜ ≥ 0        for all r ∈ W
#        t⊆r
#
#  
# solves for a similar set κᵣ of residuals which ensure that every marginal
# computed by (1) is entrywise non-negative.
#
# # Data
#
# The benchmarks use the datasets that ship with the reference implementation.
# To download them, run
#
#   git clone https://github.com/bcmullins/efficient-marginal-reconstruction \
#       src/IPM/data/efficient-marginal-reconstruction
#
# and point EMR_DATA at its datasets/ directory (the default below assumes that
# checkout).
#
# # References
#
#   - Mullins, Fuentes, Xiao, Kifer, Musco & Sheldon,
#     Efficient and Private Marginal Reconstruction with Local Non-Negativity,
#     NeurIPS 2024. arXiv:2410.01091.
#     Code and data: github.com/bcmullins/efficient-marginal-reconstruction
#
# =============================================================================

using AppleAccelerate
using LinearAlgebra, Random, Printf, Combinatorics, DelimitedFiles, JSON

using CellularSheaves.IPM
using CellularSheaves.IPM: CofreeCone, PositiveCone, IPMSettings, HSDSettings
using CellularSheaves.BlockSparseArrays: colrange, block, blocksparse
using Dualization: dual_optimizer

include("../utils.jl")

# EMR_DATA overrides the dataset directory; the default matches the clone in the
# header's Data section (gitignored under src/IPM/data/).
const DATADIR = get(ENV, "EMR_DATA",
                    joinpath(@__DIR__, "..", "..", "data",
                             "efficient-marginal-reconstruction", "datasets"))

# -----------------------------------------------------------------------------
# Residual algebra
# -----------------------------------------------------------------------------

# D₍ₖ₎ ∈ ℝ^{(n-1) × n}, successive differences.
function diffop(n::Int)
    D = zeros(n - 1, n)

    for i in 1:(n - 1)
        D[i, i] = 1.0
        D[i, i + 1] = -1.0
    end

    return D
end

# A_{γ,τ}, dense n_γ × m_τ (Proposition 1 of the reference).
function reconstruction_block(card::Vector{Int}, γ, τ)
    A = ones(1, 1)

    for k in γ
        nk = card[k]
        F = k in τ ? pinv(diffop(nk)) : fill(1 / nk, nk, 1)
        A = kron(A, F)
    end

    return A
end

# K_τ = 2^{|τ|} D_τ D_τᵀ, the loss weighting.
function residual_weight(card::Vector{Int}, τ)
    K = ones(1, 1)

    for k in τ
        D = diffop(card[k])
        K = kron(K, D * D')
    end

    return 2.0^length(τ) .* K
end

mtau(card, τ) = isempty(τ) ? 1 : prod(card[k] - 1 for k in τ)
ngam(card, γ) = isempty(γ) ? 1 : prod(card[k] for k in γ)

# -----------------------------------------------------------------------------
# Data + mechanism
# -----------------------------------------------------------------------------

struct MarginalInstance
    name::Symbol
    card::Vector{Int}
    W::Vector{Vector{Int}}
    Wd::Vector{Vector{Int}}
    z::Dict{Vector{Int},Vector{Float64}}
end

# ρ-zCDP budget matching (ε, δ), inverting δ = exp(-(ε-ρ)²/4ρ) by bisection.
function rho_from(ε::Float64, δ::Float64)
    lo, hi = 1.0e-12, 10.0

    for _ in 1:200
        mid = (lo + hi) / 2

        if exp(-((ε - mid)^2) / (4 * mid)) > δ
            hi = mid
        else
            lo = mid
        end
    end

    return lo
end

# ResidualPlanner's closed-form optimal noise allocation (getOptimalSigmasCF).
function optimal_sigmas(card, W, Wd, ρ)
    c = 2ρ
    attrmulti(t) = isempty(t) ? 1.0 : prod((card[k] - 1) / card[k] for k in t)

    function varsum(t)
        s = 0.0

        for γ in W
            issubset(t, γ) || continue
            q = prod(k in t ? (card[k] - 1) / card[k] : float(card[k])^-2 for k in γ)
            s += ngam(card, γ) * q
        end

        return s
    end

    p = Dict(t => attrmulti(t) for t in Wd)
    v = Dict(t => varsum(t) for t in Wd)
    T = sum(sqrt(v[t] * p[t]) for t in Wd)^2 / c

    return Dict(t => ((T * p[t]) / (c * v[t]))^0.25 for t in Wd)
end

# Count the τ-marginal directly from the records; p is never formed.
function marginal_counts(rows::Matrix{Int}, card, τ)
    μ = zeros(ngam(card, τ))

    for r in axes(rows, 1)
        idx = 0

        for k in τ
            idx = idx * card[k] + rows[r, k]
        end

        μ[idx + 1] += 1
    end

    return μ
end

# α_τ = D_τ μ_τ, applying D₍ₖ₎ along each mode of the τ-marginal.
function residual_from_marginal(μ::Vector{Float64}, card, τ)
    isempty(τ) && return copy(μ)

    dims = [card[k] for k in τ]
    A = reshape(μ, reverse(dims)...)

    for (j, k) in enumerate(reverse(τ))
        D = diffop(card[k])
        A = mapslices(v -> D * v, A; dims = j)
    end

    return vec(A)
end

function load_instance(name::Symbol; order::Int = 3, ε::Float64 = 1.0,
                       δ::Float64 = 1.0e-9, attrs = nothing, seed::Int = 0)
    rng = Xoshiro(seed)

    dom = JSON.parsefile(joinpath(DATADIR, "$(name)-domain.json"))

    raw, hdr = readdlm(joinpath(DATADIR, "$(name).csv"), ','; header = true)
    header = [strip(String(h)) for h in vec(hdr)]
    colof = Dict(h => i for (i, h) in enumerate(header))

    # Attributes are referenced by name (the CSV header / domain keys), never by
    # position — JSON.parsefile returns an unordered Dict. Default: CSV order.
    keep = attrs === nothing ? header : attrs
    card = Int[dom[a] for a in keep]
    rows = Matrix{Int}(undef, size(raw, 1), length(keep))

    for (j, a) in enumerate(keep)
        rows[:, j] = Int.(raw[:, colof[a]])
    end

    d = length(keep)
    W = [collect(g) for g in combinations(1:d, order)]
    Wd = Vector{Int}[]

    for k in 0:order, t in combinations(1:d, k)
        push!(Wd, collect(t))
    end

    ρ = rho_from(ε, δ)
    σ = optimal_sigmas(card, W, Wd, ρ)

    # z_τ = D_τ (μ_τ + σ_τ ξ),  ξ ~ N(0, I) on the marginal
    z = Dict{Vector{Int},Vector{Float64}}()

    for τ in Wd
        μ = marginal_counts(rows, card, τ)
        μnoisy = μ .+ σ[τ] .* randn(rng, length(μ))
        z[τ] = residual_from_marginal(μnoisy, card, τ)
    end

    return MarginalInstance(name, card, W, Wd, z)
end

# -----------------------------------------------------------------------------
# Conic model
# -----------------------------------------------------------------------------

function build_marginal_reconstruction(inst::MarginalInstance)
    card, W, Wd = inst.card, inst.W, inst.Wd
    nWd = length(Wd)

    aid = Dict(τ => i for (i, τ) in enumerate(Wd))
    sid = Dict(γ => nWd + i for (i, γ) in enumerate(W))

    coldims = vcat([mtau(card, τ) for τ in Wd], [ngam(card, γ) for γ in W])
    rowdims = [ngam(card, γ) for γ in W]

    row_ids = Int[]
    col_ids = Int[]
    blks = Matrix{Float64}[]
    gs = Float64[]

    for (r, γ) in enumerate(W)
        ng = ngam(card, γ)

        for k in 0:length(γ), τ in combinations(γ, k)
            t = collect(τ)
            push!(row_ids, r)
            push!(col_ids, aid[t])
            push!(blks, reconstruction_block(card, γ, t))
        end

        push!(row_ids, r)
        push!(col_ids, sid[γ])
        push!(blks, Matrix(-1.0I, ng, ng))

        append!(gs, zeros(ng))
    end

    B = blocksparse(row_ids, col_ids, blks, rowdims, coldims)

    Q = IPM.allocblockdiag(B)
    fill!(Q, 0)
    f = zeros(size(B, 2))

    for τ in Wd
        i = aid[τ]
        Kinv = inv(residual_weight(card, τ))
        block(Q, i, i, i) .= Kinv
        f[colrange(B, i)] .= Kinv * inst.z[τ]
    end

    cones = IPM.AbstractCone[CofreeCone() for _ in Wd]
    append!(cones, IPM.AbstractCone[PositiveCone() for _ in W])

    return IPMProblem(Q, B, f, gs, cones)
end

reconstructed(res, prob, inst, i) = res.p[colrange(prob.B, length(inst.Wd) + i)]

function build_mr_jump(inst::MarginalInstance; optimizer)
    card, W, Wd = inst.card, inst.W, inst.Wd

    model = Model(optimizer)
    α = Dict(τ => @variable(model, [1:mtau(card, τ)]) for τ in Wd)

    for γ in W
        expr = zeros(AffExpr, ngam(card, γ))

        for k in 0:length(γ), τ in combinations(γ, k)
            t = collect(τ)
            expr .+= reconstruction_block(card, γ, t) * α[t]
        end

        @constraint(model, expr .>= 0)
    end

    quad = zero(QuadExpr)

    for τ in Wd
        Kinv = inv(residual_weight(card, τ))
        r = α[τ] .- inst.z[τ]
        quad += 0.5 * sum(Kinv[i, j] * r[i] * r[j]
                          for i in eachindex(r), j in eachindex(r))
    end

    @objective(model, Min, quad)
    return model
end

# -----------------------------------------------------------------------------
# Accuracy
# -----------------------------------------------------------------------------

# Every reconstructed marginal must be non-negative, and must still sum to the
# same total as the unconstrained reconstruction (the total query is exact).
function accuracy(; tol = 1.0e-8, name = :titanic,
                  attrs = ["Survived", "Pclass", "Sex", "SibSp", "Parch", "Embarked"])
    inst = load_instance(name; attrs)
    prob = build_marginal_reconstruction(inst)
    res = solve(prob, IPMSettings{Float64}(feas_tol = tol, gap_tol = tol))

    worst = Inf
    spread = 0.0

    for i in eachindex(inst.W)
        μ̂ = reconstructed(res, prob, inst, i)
        worst = min(worst, minimum(μ̂))
        spread = max(spread, abs(sum(μ̂) - inst.z[Int[]][1]))
    end

    @printf("%-9s d=%d marginals=%d  obj %.8f  min μ̂ %.3e  |∑μ̂ - N̂| %.3e\n",
            name, length(inst.card), length(inst.W), res.pobj, worst, spread)
end

# -----------------------------------------------------------------------------
# Benchmark
# -----------------------------------------------------------------------------

# A Titanic ladder in the attribute count d, holding out the highest-cardinality
# attributes (Age = 91, Fare = 32) so the explicit conic form stays in memory:
# d = 5, then + Embarked, then + Cabin.
const PROBLEMS = [
    (:titanic, ["Survived", "Pclass", "Sex", "SibSp", "Parch"]),                       # d=5
    (:titanic, ["Survived", "Pclass", "Sex", "SibSp", "Parch", "Embarked"]),           # d=6, + Embarked
    (:titanic, ["Survived", "Pclass", "Sex", "SibSp", "Parch", "Embarked", "Cabin"]),  # d=7, + Cabin
]

function benchmark(; tol = 1.0e-6, clarabel = false, mosek = false,
                   problems = PROBLEMS, ε = 1.0)
    @printf("%-9s %4s %6s %9s %9s %9s %9s %9s %8s\n",
            "dataset", "d", "marg", "cols", "IPM", "HSD", "Clarabel", "Mosek", "Cla/IPM")

    for (name, attrs) in problems
        inst = load_instance(name; attrs, ε)
        prob = build_marginal_reconstruction(inst)

        mi = measure_ipm(prob, IPMSettings{Float64}(feas_tol = tol, gap_tol = tol))
        mh = measure_ipm(prob, HSDSettings{Float64}(feas_tol = tol, gap_tol = tol))

        # Fastest formulation per solver (measured): Clarabel dualizes faster on
        # this QP, Mosek is faster on the primal.
        mc = clarabel ?
            measure_jump(() -> build_mr_jump(inst; optimizer = dual_optimizer(clarabel_opt(; tol)))) :
            (t = NaN, status = "—", obj = NaN)

        mk = mosek ?
            measure_jump(() -> build_mr_jump(inst; optimizer = mosek_opt(; tol))) :
            (t = NaN, status = "—", obj = NaN)

        ratio = (isfinite(mc.t) && isfinite(mi.t)) ? @sprintf("%.2f", mc.t / mi.t) : "—"

        @printf("%-9s %4d %6d %9d %s %s %s %s %8s\n",
                name, length(inst.card), length(inst.W), size(prob.B, 2),
                fmt_time(mi), fmt_time(mh), fmt_time(mc), fmt_time(mk), ratio)
        flush(stdout)
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    opts = parse_args()
    tol = any(a -> startswith(a, "--tol="), ARGS) ? opts.tol : 1.0e-6
    benchmark(; clarabel = opts.clarabel, mosek = opts.mosek, tol)

    println()
    accuracy()
end

# =============================================================================
# Sample run: 2026-08-18 (Apple M1 Pro, -t 8, --clarabel --mosek), tol = 1e-6,
# real Titanic data (EMR_DATA), ε = 1, seed 0. Clarabel dualized, Mosek primal —
# each its faster formulation on this QP (measured).
#
#   dataset      d   marg      cols       IPM       HSD  Clarabel     Mosek  Cla/IPM
#   titanic      5     10      1360   31.9ms   37.6ms  117.2ms   57.2ms     3.68
#   titanic      6     20      2288   55.2ms   55.8ms  199.0ms   88.9ms     3.61
#   titanic      7     35      6332  385.9ms  468.5ms 2454.3ms  561.8ms     6.36
# =============================================================================
