# =============================================================================
#
# # Graph p-Laplacian Semi-Supervised Learning
#
# Let
#
#   X = [ x₁, …, xₙ ]
#
# be an m × n data matrix, and let G = (V, E, W) be an edge weighted
# graph with vertices V = {1, …, n}. The edge weights wᵢⱼ ∈ W quantify
# the similarity of the samples xᵢ and xⱼ.
#
# A small subset U ⊆ V of vertices carries labels g: U → ℝ. The task of
# semi-supervised learning is to propagate these labels to the rest of V.
#
#     g
#   U → ℝ
#   ↓ ↗ h
#   V
#
# which we can do by writing uᵢ = h(i) and minimizing the p-Dirichlet energy
#
#   min  1/p  ∑  wᵢⱼᵖ | uᵢ - uⱼ |ᵖ       (1)
#    u  {i,j} ∈ E
#
#   s.t. uᵢ = gᵢ,  i ∈ U
#
# When p = 2, this is Laplacian learning and the problem is linear.
# When the number |U| of labelled vertices is small, the Laplacian learning
# problem is degenerate and the labels fail to propagate. The solution is
# to choose p > d, where d is the intrinsic dimension of the data; then (1)
# is well-posed with arbitrarily few labels, and can be reformulated as a
# power-cone program
#
#   min    1/p ∑ τᵢⱼ                     (2)
#  u,τ,r {i,j} ∈ E
#
#   s.t.  rᵢⱼ - wᵢⱼ uᵢ + wᵢⱼ uⱼ = 0,   {i,j} ∈ E
#                           sᵢⱼ = 1,   {i,j} ∈ E
#               (τᵢⱼ, sᵢⱼ, rᵢⱼ) ∈ ℙ,   {i,j} ∈ E
#
# where the set
#
#   ℙ = { (a, b, c) ∈ ℝ³ : a ≥ 0, b ≥ 0, aᵗ b¹⁻ᵗ ≥ |c| }
#
# is the 3-dimensional power cone with parameter t = 1/p.
#
# # References
#
#   - Flores, Calder & Lerman,
#     Analysis and algorithms for ℓp-based semi-supervised learning on graphs,
#     Appl. Comput. Harmon. Anal. 60:77–122, 2022. arXiv:1901.05031.
#
# =============================================================================

using LinearAlgebra, Random, Statistics, Printf
using NearestNeighbors
using Dualization
import MLDatasets

using CellularSheaves.IPM
using CellularSheaves.IPM: CofreeCone, PowerCone, IPMSettings, HSDSettings
using CellularSheaves.BlockSparseArrays: colrange, block, blocksparse

include("../utils.jl")

# -----------------------------------------------------------------------------
# Data
# -----------------------------------------------------------------------------

function two_moons(n::Int; noise::Float64 = 0.1, seed::Int = 0)
    rng = Xoshiro(seed)
    m = n ÷ 2
    t1 = rand(rng, m) .* π
    t2 = rand(rng, n - m) .* π

    X = hcat(vcat(cos.(t1)', sin.(t1)'),
             vcat(1 .- cos.(t2)', 0.5 .- sin.(t2)'))
    X .+= noise .* randn(rng, 2, n)

    y = vcat(fill(-1.0, m), fill(1.0, n - m))
    return X, y
end

function pca(A::AbstractMatrix, d::Int)
    Ac = A .- mean(A, dims = 2)
    U, _, _ = svd(Ac)
    return U[:, 1:d]' * Ac
end

# Two digit classes of MNIST, projected to `dim` principal components.
function mnist_pair(a::Int, b::Int; nmax::Int = typemax(Int), dim::Int = 20, seed::Int = 0)
    ENV["DATADEPS_ALWAYS_ACCEPT"] = "true"
    data = MLDatasets.MNIST(split = :train)
    F = reshape(Float64.(data.features), 28 * 28, :)
    T = data.targets

    idx = findall(t -> t == a || t == b, T)

    if length(idx) > nmax
        idx = idx[randperm(Xoshiro(seed), length(idx))[1:nmax]]
    end

    X = pca(F[:, idx], dim)
    y = [t == a ? -1.0 : 1.0 for t in T[idx]]
    return X, y
end

# k-NN graph with self-tuning Gaussian weights: σᵢ is the distance from xᵢ to
# its kth neighbour, and wᵢⱼ = exp(-4 ‖xᵢ - xⱼ‖² / σᵢσⱼ). The graph is the
# union symmetrization of the directed k-NN graph.
function knn_weighted_edges(X::AbstractMatrix, k::Int)
    n = size(X, 2)
    idxs, dsts = knn(KDTree(X), X, k + 1, true)     # sorted; first neighbour is the point itself
    sigma = [max(dsts[i][end], eps()) for i in 1:n]

    seen = Set{Tuple{Int, Int}}()
    E = Tuple{Int, Int}[]
    w = Float64[]

    for i in 1:n
        for j in idxs[i]
            i == j && continue
            key = minmax(i, j)
            key in seen && continue
            push!(seen, key)
            d2 = sum(abs2, @view(X[:, i]) .- @view(X[:, j]))
            push!(E, key)
            push!(w, exp(-4 * d2 / (sigma[i] * sigma[j])))
        end
    end

    return E, w
end

function pick_labels(y::AbstractVector, nlab::Int; seed::Int = 0)
    rng = Xoshiro(seed)
    labels = Dict{Int, Float64}()

    for cls in (-1.0, 1.0)
        cand = findall(==(cls), y)
        for v in cand[randperm(rng, length(cand))[1:nlab]]
            labels[v] = cls
        end
    end

    return labels
end

# -----------------------------------------------------------------------------
# Build
# -----------------------------------------------------------------------------

# stalks: 1..nu are the unlabeled uᵥ, then nu + e is the power stalk
# (τₑ, sₑ, rₑ) of the eth active edge. Rows come in blocks of two per edge:
#
#   rₑ - wₑ uᵢ + wₑ uⱼ = wₑ (labeled terms)
#                   sₑ = 1
#
# Requires every unlabeled vertex to be connected to a labeled one through
# the graph; otherwise the minimizer is not unique.
function build_plap(nv::Int, E::Vector{Tuple{Int, Int}}, w::Vector{Float64},
                    labels::Dict{Int, Float64}, p::Float64)
    umap = Dict{Int, Int}()                 # vertex -> unlabeled stalk id

    for v in 1:nv
        haskey(labels, v) || (umap[v] = length(umap) + 1)
    end

    nu = length(umap)

    # [r-row; s-row] selectors into the (τ, s, r) stalk
    Sp = [0.0 0.0 1.0; 0.0 1.0 0.0]

    row_ids = Int[]
    col_ids = Int[]
    blks = Matrix{Float64}[]
    gs = Float64[]

    e = 0

    for (t, (i, j)) in enumerate(E)
        li = haskey(labels, i)
        lj = haskey(labels, j)
        li && lj && continue                # constant contribution; dropped

        e += 1
        rhs = 0.0

        push!(row_ids, e); push!(col_ids, nu + e); push!(blks, copy(Sp))

        if li
            rhs += w[t] * labels[i]
        else
            push!(row_ids, e); push!(col_ids, umap[i]); push!(blks, reshape([-w[t], 0.0], 2, 1))
        end

        if lj
            rhs -= w[t] * labels[j]
        else
            push!(row_ids, e); push!(col_ids, umap[j]); push!(blks, reshape([w[t], 0.0], 2, 1))
        end

        push!(gs, rhs); push!(gs, 1.0)
    end

    # every unlabeled stalk must be touched, or blocksparse infers zero width
    touched = falses(nu)

    for j in col_ids
        j <= nu && (touched[j] = true)
    end

    @assert all(touched) "every unlabeled vertex needs an incident active edge"

    B = blocksparse(row_ids, col_ids, blks)

    Q = IPM.allocblockdiag(B)
    fill!(Q, 0)

    # c: -1/p on each τₑ, so that -c'p = 1/p ∑ τₑ
    c = zeros(size(B, 2))

    for k in 1:e
        c[first(colrange(B, nu + k))] = -1.0 / p
    end

    cones = vcat(IPM.AbstractCone[CofreeCone() for _ in 1:nu],
                 IPM.AbstractCone[PowerCone(1.0 / p) for _ in 1:e])

    return IPMProblem(Q, B, c, gs, cones), umap
end

function build_plap_jump(nv::Int, E::Vector{Tuple{Int, Int}}, w::Vector{Float64},
                         labels::Dict{Int, Float64}, p::Float64; optimizer)
    umap = Dict{Int, Int}()

    for v in 1:nv
        haskey(labels, v) || (umap[v] = length(umap) + 1)
    end

    act = [(t, i, j) for (t, (i, j)) in enumerate(E)
           if !(haskey(labels, i) && haskey(labels, j))]

    model = Model(optimizer)
    @variable(model, u[1:length(umap)])
    @variable(model, tau[1:length(act)])
    @variable(model, r[1:length(act)])

    uexpr(v) = haskey(labels, v) ? labels[v] : u[umap[v]]

    for (e, (t, i, j)) in enumerate(act)
        @constraint(model, r[e] == w[t] * (uexpr(i) - uexpr(j)))
        @constraint(model, [tau[e], 1.0, r[e]] in MOI.PowerCone(1.0 / p))
    end

    @objective(model, Min, sum(tau) / p)
    return model
end

# -----------------------------------------------------------------------------
# Recovery and optimality residual
# -----------------------------------------------------------------------------

function recover_u(res, prob, nv::Int, labels::Dict{Int, Float64}, umap::Dict{Int, Int})
    u = zeros(nv)

    for v in 1:nv
        if haskey(labels, v)
            u[v] = labels[v]
        else
            u[v] = res.p[first(colrange(prob.B, umap[v]))]
        end
    end

    return u
end

# max over unlabeled vertices of the p-Laplace equation (3)
function plap_residual(u::Vector{Float64}, E::Vector{Tuple{Int, Int}},
                       w::Vector{Float64}, labels::Dict{Int, Float64}, p::Float64)
    r = zeros(length(u))

    for (t, (i, j)) in enumerate(E)
        d = u[j] - u[i]
        f = w[t]^p * abs(d)^(p - 2) * d
        r[i] += f
        r[j] -= f
    end

    return maximum(abs(r[v]) for v in 1:length(u) if !haskey(labels, v))
end

# -----------------------------------------------------------------------------
# Accuracy — the weighted path graph has a closed-form minimizer
# -----------------------------------------------------------------------------

# Path v₁ — v₂ — … — vₙ with weights wᵢ on edge (i, i+1) and labels u₁ = 0,
# uₙ = 1. Stationarity of ∑ wᵢᵖ Δᵢᵖ under ∑ Δᵢ = 1 gives increments
# Δᵢ ∝ wᵢ^(-q), q = p/(p-1), so uₖ = ∑_{i<k} wᵢ^(-q) / ∑ᵢ wᵢ^(-q).
function path_accuracy(n::Int, p::Float64; tol = 1.0e-8, seed::Int = 0)
    rng = Xoshiro(seed)
    w = exp.(rand(rng, n - 1) .* 2 .- 1)

    E = [(i, i + 1) for i in 1:(n - 1)]
    labels = Dict(1 => 0.0, n => 1.0)

    prob, umap = build_plap(n, E, w, labels, p)
    res = solve(prob, IPMSettings{Float64}(feas_tol = tol, gap_tol = tol))
    u = recover_u(res, prob, n, labels, umap)

    q = p / (p - 1)
    s = w .^ (-q)
    exact = vcat(0.0, cumsum(s) ./ sum(s))

    err = maximum(abs.(u .- exact))
    return (; err, status = res.status)
end

# -----------------------------------------------------------------------------
# Benchmark
# -----------------------------------------------------------------------------

function load_instances()
    insts = []

    X, y = two_moons(2000; noise = 0.1)
    E, w = knn_weighted_edges(X, 10)
    push!(insts, ("Moons-2k", size(X, 2), E, w, y))

    X, y = mnist_pair(3, 8; nmax = 2000)
    E, w = knn_weighted_edges(X, 10)
    push!(insts, ("MNIST38-2k", size(X, 2), E, w, y))

    return insts
end

function benchmark(; ps = (3.0, 10.0), nlab = 5, tol = 1.0e-6,
                     clarabel = false, mosek = false)
    @printf("%-12s %6s %7s %4s %9s %9s %9s %9s %8s %8s\n",
            "dataset", "n", "|E|", "p", "IPM", "HSD", "Clarabel", "Mosek", "Cla/IPM", "Mos/IPM")

    for (name, n, E, w, y) in load_instances()
        labels = pick_labels(y, nlab)

        for p in ps
            prob, _ = build_plap(n, E, w, labels, p)

            mi = measure_ipm(prob, IPMSettings{Float64}(feas_tol = tol, gap_tol = tol))
            mh = measure_ipm(prob, HSDSettings{Float64}(feas_tol = tol, gap_tol = tol))

            if clarabel
                mc = measure_jump(() -> build_plap_jump(n, E, w, labels, p; optimizer = clarabel_opt(; tol)))
            else
                mc = (t = NaN, status = "—", obj = NaN)
            end

            if mosek
                mk = measure_jump(() -> build_plap_jump(n, E, w, labels, p; optimizer = dual_optimizer(mosek_opt(; tol))))
            else
                mk = (t = NaN, status = "—", obj = NaN)
            end

            cratio = (isfinite(mc.t) && isfinite(mi.t)) ? @sprintf("%.2f", mc.t / mi.t) : "—"
            mratio = (isfinite(mk.t) && isfinite(mi.t)) ? @sprintf("%.2f", mk.t / mi.t) : "—"

            @printf("%-12s %6d %7d %4.0f %s %s %s %s %8s %8s\n",
                    name, n, length(E), p, fmt_time(mi), fmt_time(mh), fmt_time(mc), fmt_time(mk), cratio, mratio)
        end
    end
end

function classification(; ps = (3.0, 10.0), nlab = 5, tol = 1.0e-6)
    @printf("%-12s %4s %10s %13s\n", "dataset", "p", "class acc", "p-Lap resid")

    for (name, n, E, w, y) in load_instances()
        labels = pick_labels(y, nlab)

        for p in ps
            prob, umap = build_plap(n, E, w, labels, p)
            res = solve(prob, IPMSettings{Float64}(feas_tol = tol, gap_tol = tol))
            u = recover_u(res, prob, n, labels, umap)

            unl = [v for v in 1:n if !haskey(labels, v)]
            acc = mean(sign(u[v]) == y[v] for v in unl)
            rsd = plap_residual(u, E, w, labels, p)

            @printf("%-12s %4.0f %9.1f%% %13.3e\n", name, p, 100 * acc, rsd)
        end
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    opts = parse_args()
    tol = any(a -> startswith(a, "--tol="), ARGS) ? opts.tol : 1.0e-6
    benchmark(; clarabel = opts.clarabel, mosek = opts.mosek, tol)

    println()
    @printf("%-6s %4s %13s\n", "n", "p", "‖u-u*‖∞")

    for p in (3.0, 5.0, 10.0)
        a = path_accuracy(1000, p; tol = 1.0e-8)
        @printf("%-6d %4.0f %13.3e\n", 1000, p, a.err)
    end

    println()
    classification(; tol)
end

# =============================================================================
# Sample run: 2026-08-24 (-t 8, --clarabel --mosek), tol = 1e-6. Mosek via Dualization.
# * = non-optimal solve (Clarabel SLOW_PROGRESS on Moons-2k).
#
#   dataset           n     |E|    p       IPM       HSD  Clarabel     Mosek  Cla/IPM  Mos/IPM
#   Moons-2k       2000   12138    3  197.2ms  272.5ms 485.6ms*  467.9ms     2.46     2.37
#   Moons-2k       2000   12138   10  195.3ms  281.3ms 809.7ms*  434.9ms     4.15     2.23
#   MNIST38-2k     2000   13802    3  341.7ms  415.2ms 2272.6ms  570.5ms     6.65     1.67
#   MNIST38-2k     2000   13802   10  667.8ms  413.2ms 2288.7ms  526.1ms     3.43     0.79
# =============================================================================
