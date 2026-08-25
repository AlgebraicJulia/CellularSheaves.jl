# =============================================================================
#
# # Convex Clustering
#
# Let
#
#   A = [ a₁ a₂ … aₙ ]
#
# be a m × n data matrix with n observations and m features. Given an
# edge-weighted graph G = (V, E, W) with vertices V = {1, …, n},
#
#            w₁₂
#       1 ●────────● 2
#
# the weighted convex clustering model for A solves the convex optimization problem
#
#           n
#   min 1/2 ∑ ‖ xᵢ - aᵢ ‖² + γ ∑ wᵢⱼ ‖ xᵢ - xⱼ ‖   (1)
#    X     i=1               {i,j} ∈ E
#
# over all m × n centroid matrices
#
#   X = [ x₁ x₂ … xₙ ].
#
# Each distinct centroid xᵢ represents a cluster centered at xᵢ. Observations
# are assigned to the same cluster if their centroids are the same. The number
# γ > 0 is a tuning parameter. The problem (1) can be reformulated as a conic
# quadratic program
#
#               n
#     min   1/2 ∑ ‖ xᵢ ‖² - aᵢᵀ xᵢ + γ ∑ wᵢⱼ tᵢⱼ   (2)
#   X, U, T    i=1                   {i,j} ∈ E
#
#     s.t.  uᵢⱼ - xᵢ + xⱼ = 0,    {i,j} ∈ E
#                 ‖ uᵢⱼ ‖ ≤ tᵢⱼ,  {i,j} ∈ E
#
# The set { (t, u) : ‖ u ‖ ≤ t } ⊆ ℝᵐ⁺¹ is a convex cone called the second order
# cone.
#
# # References
#
#  - Yuan, Sun & Toh,
#    An Efficient Semismooth Newton Based Algorithm for Convex Clustering,
#    ICML 2018, PMLR 80, 5718–5726.
#
# =============================================================================

using LinearAlgebra, Statistics, Printf
using Graphs, NearestNeighbors
import MLDatasets, DataFrames

using CellularSheaves.IPM
using CellularSheaves.IPM: CofreeCone, SecondOrderCone, IPMSettings, HSDSettings
using CellularSheaves.BlockSparseArrays: colrange, block, blocksparse

include("../utils.jl")

# -----------------------------------------------------------------------------
# Data
# -----------------------------------------------------------------------------

zscore(A) = (A .- mean(A, dims = 2)) ./ (std(A, dims = 2) .+ eps())

function pca(A::AbstractMatrix, d::Int)
    Ac = A .- mean(A, dims = 2)
    U, _, _ = svd(Ac)
    return U[:, 1:d]' * Ac
end

function load_datasets()
    ENV["DATADEPS_ALWAYS_ACCEPT"] = "true"
    iris  = zscore(Matrix{Float64}(Matrix(MLDatasets.Iris().features)'))
    wine  = zscore(Matrix{Float64}(Matrix(MLDatasets.Wine().features)'))
    mnist = pca(reshape(Float64.(MLDatasets.MNIST(split = :train).features), 28 * 28, :)[:, 1:1000], 10)
    return [("Iris", iris), ("WINE", wine), ("MNIST-1k", mnist)]
end

function knn_edges(A::AbstractMatrix, k::Int)
    idxs, _ = knn(KDTree(A), A, k + 1)     # k+1: a point's nearest neighbour is itself
    G = Graph(size(A, 2))

    for i in eachindex(idxs)
        for j in idxs[i]
            i == j || add_edge!(G, i, j)
        end
    end

    return G
end

# -----------------------------------------------------------------------------
# Build
# -----------------------------------------------------------------------------

function build_convex_clustering(A::AbstractMatrix, G::AbstractGraph;
                                 gamma::Float64 = 1.0, phi::Float64 = 0.5)
    d, n = size(A)
    # stalks: 1..n are the x_i, then n+e is s_e for edge e

    Id = Matrix{Float64}(I, d, d)
    # [0 | I_d] : the u_e slot of s_e, skipping the leading t_e column
    Su = hcat(zeros(d, 1), Id)

    row_ids = Int[]
    col_ids = Int[]
    blks = Matrix{Float64}[]

    for (e, edge) in enumerate(edges(G))
        i = src(edge)
        j = dst(edge)
        push!(row_ids, e); push!(col_ids, n + e); push!(blks, copy(Su))
        push!(row_ids, e); push!(col_ids, i);     push!(blks, -copy(Id))
        push!(row_ids, e); push!(col_ids, j);     push!(blks, copy(Id))
    end

    B = blocksparse(row_ids, col_ids, blks)
    g = zeros(size(B, 1))                       # all rows are homogeneous

    # Q: identity on each x_i stalk (the fidelity term), zero on the SOC stalks
    Q = IPM.allocblockdiag(B)
    fill!(Q, 0)

    for i in 1:n
        copyto!(block(Q, i, i, i), I)
    end

    # c: a_i on x_i (so -c'p = -sum a_i' x_i), -gamma*w_ij on t_e
    c = zeros(size(B, 2))

    for i in 1:n
        c[colrange(B, i)] .= @view A[:, i]
    end

    for (e, edge) in enumerate(edges(G))
        i = src(edge)
        j = dst(edge)
        w = exp(-phi * sum(abs2, @view(A[:, i]) .- @view(A[:, j])))
        c[first(colrange(B, n + e))] = -gamma * w
    end

    cones = vcat(IPM.AbstractCone[CofreeCone() for _ in 1:n],
                 IPM.AbstractCone[SecondOrderCone() for _ in 1:ne(G)])

    return IPMProblem(Q, B, c, g, cones)
end

function build_cc_jump(A::AbstractMatrix, G::AbstractGraph;
                       gamma::Float64, phi::Float64, optimizer)
    d, n = size(A)
    model = Model(optimizer)
    @variable(model, x[1:d, 1:n])
    @variable(model, t[1:ne(G)])

    for (e, edge) in enumerate(edges(G))
        i = src(edge)
        j = dst(edge)
        @constraint(model, vcat(t[e], x[:, i] .- x[:, j]) in JuMP.SecondOrderCone())
    end

    obj = 0.5 * sum((x[a, i] - A[a, i])^2 for a in 1:d, i in 1:n)

    for (e, edge) in enumerate(edges(G))
        i = src(edge)
        j = dst(edge)
        w = exp(-phi * sum(abs2, @view(A[:, i]) .- @view(A[:, j])))
        obj += gamma * w * t[e]
    end

    @objective(model, Min, obj)
    return model
end

# -----------------------------------------------------------------------------
# Benchmark
# -----------------------------------------------------------------------------

function benchmark(; k = 10, phi = 0.5, gamma = 1.0, tol = 1e-6,
                     clarabel = false, mosek = false)
    @printf("%-9s %3s %6s %9s %9s %9s %9s %8s %8s\n",
            "dataset", "d", "n", "IPM", "HSD", "Clarabel", "Mosek", "Cla/IPM", "Mos/IPM")

    for (name, A) in load_datasets()
        d, n = size(A)
        G = knn_edges(A, k)
        prob = build_convex_clustering(A, G; gamma, phi)

        mi = measure_ipm(prob, IPMSettings{Float64}(feas_tol = tol, gap_tol = tol))
        mh = measure_ipm(prob, HSDSettings{Float64}(feas_tol = tol, gap_tol = tol))

        if clarabel
            mc = measure_jump(() -> build_cc_jump(A, G; gamma, phi, optimizer = clarabel_opt(; tol)))
        else
            mc = (t = NaN, status = "—", obj = NaN)
        end

        if mosek
            mk = measure_jump(() -> build_cc_jump(A, G; gamma, phi, optimizer = mosek_opt(; tol)))
        else
            mk = (t = NaN, status = "—", obj = NaN)
        end

        cratio = (isfinite(mc.t) && isfinite(mi.t)) ? @sprintf("%.2f", mc.t / mi.t) : "—"
        mratio = (isfinite(mk.t) && isfinite(mi.t)) ? @sprintf("%.2f", mk.t / mi.t) : "—"

        @printf("%-9s %3d %6d %s %s %s %s %8s %8s\n",
                name, d, n, fmt_time(mi), fmt_time(mh), fmt_time(mc), fmt_time(mk), cratio, mratio)
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    opts = parse_args()
    tol = any(a -> startswith(a, "--tol="), ARGS) ? opts.tol : 1e-6
    benchmark(; clarabel = opts.clarabel, mosek = opts.mosek, tol)
end

# =============================================================================
# Sample run: 2026-08-24 (--clarabel --mosek), tol = 1e-6, γ = 1, k = 10.
#
#   dataset     d      n       IPM       HSD  Clarabel     Mosek  Cla/IPM  Mos/IPM
#   Iris        4    150   19.3ms   22.8ms   27.5ms   21.8ms     1.42     1.13
#   WINE       13    178  140.9ms  177.5ms  704.2ms  144.4ms     5.00     1.02
#   MNIST-1k   10   1000 1540.4ms 1732.2ms 24435.3ms 2478.5ms    15.86     1.61
# =============================================================================
