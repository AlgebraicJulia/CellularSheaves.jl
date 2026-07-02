######################################################################
# wasserstein_graph.jl
#
# GRAPH-STRUCTURED OPTIMAL TRANSPORT as a cellular sheaf — the pairwise
# (local-consistency) relaxation of multimarginal OT with a cost that
# decouples over a graph G = (V, E):
#
#     min  Σ_{e=(i,j)∈E} ⟨C_e, π_e⟩  (+ ε/2 ‖·‖²)
#     s.t. rowmarg(π_e) = μ_i,  colmarg(π_e) = μ_j
#          μ_k = ν_k for pinned nodes (folded into g),
#          π_e, μ_free ≥ 0                       (PositiveCone cells)
#
# One coupling cell per COST-GRAPH edge, marginal cells only for free
# nodes; pinned marginals ride in g. Three readings of one assembly:
#
#   • all nodes pinned            → the pairwise RELAXATION of graphical
#     MMOT (Haasler–Singh–Zhang–Karlsson–Chen, IEEE-TIT 67:4647).
#     EXACT on trees (Haasler–Ringh–Chen–Karlsson, SIAM JCO 59:2428,
#     §4: tree MMOT ≡ sum of bimarginal costs), and NECESSARILY loose
#     on cyclic cost graphs (Altschuler–Boix-Adserà: high-treewidth
#     pairwise MMOT is NP-hard — same P≠NP logic that makes MAP-LP
#     loose; this file is map_lp.jl wearing transport clothes).
#   • star, center free           → the WASSERSTEIN BARYCENTER of the
#     pinned leaves (Agueh–Carlier); its famous tractability is now a
#     corollary: star ⊂ tree ⊂ tight.
#   • general graph, some pins    → Wasserstein propagation (Solomon
#     et al. 2014) — couplings are the honest variables, exact on any
#     topology.
#
# THE GAP IS HOLONOMY (measured, wasserstein_graph_oracle.py): on a
# cycle with quadratic costs where `nanti` edges want y = 1−x and the
# rest y = x, the pairwise optimum is 0 for every configuration (each
# coupling independently a Monge map, all marginals uniform — locally
# consistent!), while exact MMOT is 0 iff the COMPOSED Monge map around
# the cycle is the identity, i.e. iff nanti is even:
#     K=3: gap 0.3265 / 0.0000    K=4: 0.2551 / 0.0000
#     K=5: gap 0.2400 / 0.0000        (nanti = 1 / 2)
# — cost-graph frustration, the transport sibling of map_lp's odd-cycle
# story, whose binary agreement-cost triangle this file reproduces
# verbatim (pairwise 0, exact 1). Locally consistent couplings that are
# marginals of NO joint law: the positive-global-section obstruction.
#
# ENTROPIC DIAL. Adding (1/η)Σ π log π per coupling (ExponentialCone
# leaves, one per coupling coordinate) gives Sinkhorn-comparable
# problems — and a SECOND gap, on trees, orthogonal to topology: the
# pairwise and multimarginal regularizations differ even where the
# unregularized problems coincide (HRCK's smoothing observation;
# measured: the pairwise middle marginal is strictly more diffuse).
#
# H¹, measured then explained: ONE KANTOROVICH GAUGE PER COUPLING —
# both marginal blocks of a coupling imply the same total mass (the
# classical transportation-polytope redundancy; potentials (f,g) ↦
# (f+c, g−c)) — glued at free nodes, whose μ columns force the incident
# gauge constants to balance:
#     H¹ = |E| − rank(signed free-node incidence).
# Measured: path/triangle all-pinned 3, star center-free 2, 4-cycle
# one-free 3. `wg_h1_predict` implements the law. Rings of couplings
# are thus consistently rank-deficient by design — another Uzawa probe.
#
# PROPER BLOCK-DIAGONAL Q. Beyond the ε floor, the file carries genuine
# per-cell Grams — the vertex blocks of the KKT system stop being scalar:
# first-class η (quadratically regularized OT: sparse where entropic is
# dense), Courty-style Laplacian Grams on couplings ("transported
# neighbors stay neighbors", POT-cross-checked), evidence Grams ΦᵀΦ
# (inverse-OT: fit observed matching statistics), and H¹ smoothness
# priors on free marginals. See the "proper block-diagonal Q" section
# below and the md note.
#
# Written against the CellularSheaves.IPM PR-67 API; not executed here.
# wasserstein_graph_oracle.py is the numerical ground truth.
######################################################################

using CellularSheaves.IPM
using CellularSheaves.BlockSparseArrays: colrange, rowrange, blocksparse, block, nvtxs
using LinearAlgebra
using Printf
using Random

# ---- instance ----------------------------------------------------------

struct WGInstance
    supports::Vector{Vector{Float64}}
    edges::Vector{Tuple{Int, Int}}       # cost-graph edges (i, j)
    C::Vector{Matrix{Float64}}           # C[t] is n_i × n_j
    pins::Dict{Int, Vector{Float64}}
    free::Vector{Int}
    ε::Float64
    η::Float64                           # first-class ½η‖π‖² (quad-reg OT)
    Qc::Dict{Int, Matrix{Float64}}       # proper Gram per coupling t
    cx::Dict{Int, Vector{Float64}}       # extra linear term per coupling
    Qμ::Dict{Int, Matrix{Float64}}       # proper Gram per free node k
end

function wg_instance(supports, edges, C, pins; ε::Float64 = 1e-9, η::Float64 = 0.0,
                     Qc = Dict{Int, Matrix{Float64}}(), cx = Dict{Int, Vector{Float64}}(),
                     Qμ = Dict{Int, Matrix{Float64}}())
    free = [k for k in 1:length(supports) if !haskey(pins, k)]
    @assert !isempty(pins) "at least one pinned node (normalization enters through pins)"
    WGInstance(supports, edges, C, pins, free, ε, η, Qc, cx, Qμ)
end

quad_cost(xi, xj; anti::Bool = false) =
    [anti ? (x + y - 1)^2 : (x - y)^2 for x in xi, y in xj]

"""Cycle of K uniform marginals on an n-grid; `nanti` edges want
y = 1−x, the rest y = x. Gap iff nanti is odd (holonomy)."""
function wg_cycle_instance(; K::Int = 3, n::Int = 8, nanti::Int = 1, ε = 1e-9)
    xs = collect(range(0, 1, length = n))
    edges = [(k, mod1(k + 1, K)) for k in 1:K]
    C = [quad_cost(xs, xs; anti = k ≤ nanti) for k in 1:K]
    pins = Dict(k => fill(1.0 / n, n) for k in 1:K)
    wg_instance(fill(xs, K), edges, C, pins; ε)
end

"""map_lp's frustrated triangle, transport reading: binary uniform
marginals, cost 1 iff the pair agrees. Pairwise 0, exact MMOT 1."""
function wg_binary_triangle(; ε = 1e-9)
    xs = [0.0, 1.0]
    pins = Dict(k => [0.5, 0.5] for k in 1:3)
    wg_instance(fill(xs, 3), [(1, 2), (2, 3), (3, 1)],
                [Matrix{Float64}(I, 2, 2) for _ in 1:3], pins; ε)
end

function wg_path_instance(; K::Int = 4, n::Int = 10, seed::Int = 0, ε = 1e-9)
    rng = MersenneTwister(seed)
    xs = collect(range(0, 1, length = n))
    edges = [(k, k + 1) for k in 1:(K - 1)]
    C = [quad_cost(xs, xs) for _ in edges]
    pins = Dict{Int, Vector{Float64}}()
    for k in 1:K
        v = abs.(randn(rng, n)) .+ 0.05
        pins[k] = v ./ sum(v)
    end
    wg_instance(fill(xs, K), edges, C, pins; ε)
end

"""Star: node 1 free (the barycenter), leaves pinned to νs, edge costs
λ_k · quadratic. Tightness is the tree theorem's corollary."""
function wg_star_barycenter(νs::Vector{Vector{Float64}}, xs::Vector{Float64},
                            λs::Vector{Float64}; ε = 1e-9)
    L = length(νs)
    edges = [(1, k + 1) for k in 1:L]
    C = [λs[k] .* quad_cost(xs, xs) for k in 1:L]
    pins = Dict(k + 1 => νs[k] for k in 1:L)
    wg_instance(fill(xs, L + 1), edges, C, pins; ε)
end

# ---- proper block-diagonal Q: Laplacian / evidence / smoothness ---------
#
# Three literature-real dense Grams (see wasserstein_graph.md §"Proper Q"):
#   • Laplacian-regularized OT (Courty–Flamary–Tuia–Rakotomamonjy, TPAMI
#     2017; POT's emd_laplace): "transported neighbors stay neighbors".
#     Q_e = 2η[α (X_tX_tᵀ ⊗ L_s) + (1−α)(L_t ⊗ X_sX_sᵀ)], matching POT's
#     implemented f(G) exactly (oracle: Gram self-test 1.7e-13, our QP ≤
#     POT's CG value, Δ = 2.6e-5).
#   • Evidence Gram (inverse-OT / matching econometrics): ½λ‖Φπ − y‖²
#     from observed plan reads: Q_e += λΦᵀΦ, cx −= λΦᵀy.
#   • Smoothness prior on a free marginal: Qμ = λ·(path Laplacian).
# Plus the first-class η (quadratically regularized OT, Blondel et al.
# 2018 / Lorenz–Manns–Meyer): sparse plans where entropic is dense
# (oracle: support 83 vs 519 of 625 at comparable smoothing).

wg_laplacian(S::Matrix{Float64}) = Diagonal(vec(sum(S, dims = 2))) - S

function wg_gauss_sim(X::Matrix{Float64}, σ::Float64)
    n = size(X, 1)
    S = [exp(-sum(abs2, X[i, :] .- X[j, :]) / (2σ^2)) for i in 1:n, j in 1:n]
    for i in 1:n; S[i, i] = 0.0; end
    return S
end

"""Courty Gram: ½vec(G)ᵀQ vec(G) = η[α tr(XtᵀGᵀLsGXt) + (1−α)tr(XsᵀGLtGᵀXs)]."""
wg_laplace_gram(Xs, Xt, Ls, Lt; α::Float64 = 0.5, η::Float64 = 0.4) =
    2η .* (α .* kron(Xt * Xt', Matrix(Ls)) .+ (1 - α) .* kron(Matrix(Lt), Xs * Xs'))

function wg_path_laplacian(n::Int)
    L = 2.0 .* Matrix{Float64}(I, n, n)
    L[1, 1] = L[n, n] = 1.0
    for i in 1:(n - 1)
        L[i, i + 1] = L[i + 1, i] = -1.0
    end
    return L
end

"""Domain-adaptation-flavored instance: two 2D point clouds, squared cost,
Gaussian-similarity Laplacian Gram on the single coupling."""
function wg_laplace_instance(; ns::Int = 20, nt::Int = 20, seed::Int = 7,
                             α::Float64 = 0.5, η::Float64 = 0.4, σ::Float64 = 0.8)
    rng = MersenneTwister(seed)
    Xs = randn(rng, ns, 2)
    Xt = randn(rng, nt, 2) .+ [1.5 0.5]
    M = [sum(abs2, Xs[i, :] .- Xt[j, :]) for i in 1:ns, j in 1:nt]
    Q = wg_laplace_gram(Xs, Xt, wg_laplacian(wg_gauss_sim(Xs, σ)),
                        wg_laplacian(wg_gauss_sim(Xt, σ)); α, η)
    wg_instance([collect(1.0:ns), collect(1.0:nt)], [(1, 2)], [M],
                Dict(1 => fill(1 / ns, ns), 2 => fill(1 / nt, nt)); Qc = Dict(1 => Q))
end

"""Attach evidence ½λ‖Φ·vec(π_t) − y‖² to coupling t (inverse-OT reading:
observed matching statistics pull the plan)."""
function wg_with_evidence(inst::WGInstance, t::Int, Φ::Matrix{Float64},
                          y::Vector{Float64}, λ::Float64)
    Qc = copy(inst.Qc); cx = copy(inst.cx)
    Qc[t] = get(Qc, t, zeros(size(Φ, 2), size(Φ, 2))) .+ λ .* (Φ' * Φ)
    cx[t] = get(cx, t, zeros(size(Φ, 2))) .- λ .* (Φ' * y)
    WGInstance(inst.supports, inst.edges, inst.C, inst.pins, inst.free,
               inst.ε, inst.η, Qc, cx, inst.Qμ)
end

"""Barycenter star with an H¹ smoothness prior on the free center."""
function wg_smooth_barycenter(νs, xs, λs; λsmooth::Float64 = 0.02, ε = 1e-9)
    inst = wg_star_barycenter(νs, xs, λs; ε)
    WGInstance(inst.supports, inst.edges, inst.C, inst.pins, inst.free,
               inst.ε, inst.η, inst.Qc, inst.cx,
               Dict(1 => λsmooth .* wg_path_laplacian(length(xs))))
end

# ---- LIFTED Laplacian formulation ----------------------------------------
# Instead of one cell with dense ns*nt × ns*nt Gram, we lift:
#   • ns*nt vertex cells (one per coupling entry π_ij)
#   • edge cells for source/target Laplacian differences
# Q becomes diagonal on edge cells. Dense Kronecker → sparse lifted.

"""Extract edges from a similarity-based Laplacian (threshold small weights)."""
function laplacian_edges(S::Matrix{Float64}; thresh::Float64 = 1e-6)
    n = size(S, 1)
    edges = Tuple{Int, Int, Float64}[]
    for i in 1:n, j in (i + 1):n
        S[i, j] > thresh && push!(edges, (i, j, S[i, j]))
    end
    return edges
end

"""LIFTED Laplacian-regularized OT: each coupling entry is a vertex cell,
edge cells capture (π_ij - π_kj) for source edges (i,k) and (π_ij - π_il)
for target edges (j,l). Q is DIAGONAL on edge cells.

Variables:
  • v_{ij} = π_ij  for i=1..ns, j=1..nt  (vertex cells, PositiveCone)
  • y^s_{(ik)j} for each source edge (i,k) and target j  (edge cells, CofreeCone)
  • y^t_{i(jl)} for each target edge (j,l) and source i  (edge cells, CofreeCone)

Agreements:
  • y^s_{(ik)j} = π_ij - π_kj  (source edge × target)
  • y^t_{i(jl)} = π_ij - π_il  (source × target edge)
  • marginal constraints: Σ_j π_ij = μ_i, Σ_i π_ij = ν_j

Objective:
  • Σ_{ij} M_ij π_ij  (linear cost)
  • + α·η · Σ_{(ik),j} ||Xt_j||² (y^s_{(ik)j})²  (source smoothness)
  • + (1-α)·η · Σ_{i,(jl)} ||Xs_i||² (y^t_{i(jl)})²  (target smoothness)
"""
function build_lifted_laplace_problem(; ns::Int = 20, nt::Int = 20, seed::Int = 7,
                                      α::Float64 = 0.5, η::Float64 = 0.4, σ::Float64 = 0.8,
                                      edge_thresh::Float64 = 0.1, ε::Float64 = 1e-9)
    rng = MersenneTwister(seed)
    Xs = randn(rng, ns, 2)
    Xt = randn(rng, nt, 2) .+ [1.5 0.5]
    M = [sum(abs2, Xs[i, :] .- Xt[j, :]) for i in 1:ns, j in 1:nt]

    # Similarity matrices and edges
    Ss = wg_gauss_sim(Xs, σ)
    St = wg_gauss_sim(Xt, σ)
    Es = laplacian_edges(Ss; thresh = edge_thresh)
    Et = laplacian_edges(St; thresh = edge_thresh)

    # Feature norms squared
    Xs_norm2 = [sum(abs2, Xs[i, :]) for i in 1:ns]
    Xt_norm2 = [sum(abs2, Xt[j, :]) for j in 1:nt]

    # Indexing: v_{ij} at cell (i-1)*nt + j for i=1..ns, j=1..nt
    # Then source edge cells, then target edge cells
    nv_coupling = ns * nt
    nv_source_edges = length(Es) * nt
    nv_target_edges = length(Et) * ns
    nv = nv_coupling + nv_source_edges + nv_target_edges

    vdims = ones(Int, nv)  # all 1D cells

    # Build agreement matrix B and RHS g
    row_ids = Int[]; col_ids = Int[]; blks = Matrix{Float64}[]
    place!(e, v, m) = (push!(row_ids, e); push!(col_ids, v); push!(blks, reshape([m], 1, 1)))
    gval = Dict{Int, Vector{Float64}}()
    ec = Ref(0)

    # Source edge agreements: y^s_{(ik)j} = π_ij - π_kj
    for (t, (i, k, _)) in enumerate(Es)
        for j in 1:nt
            v_coup_ij = (i - 1) * nt + j
            v_coup_kj = (k - 1) * nt + j
            v_edge = nv_coupling + (t - 1) * nt + j
            id = (ec[] += 1)
            place!(id, v_coup_ij, 1.0)
            place!(id, v_coup_kj, -1.0)
            place!(id, v_edge, -1.0)
            gval[id] = [0.0]
        end
    end

    # Target edge agreements: y^t_{i(jl)} = π_ij - π_il
    for (t, (j, l, _)) in enumerate(Et)
        for i in 1:ns
            v_coup_ij = (i - 1) * nt + j
            v_coup_il = (i - 1) * nt + l
            v_edge = nv_coupling + nv_source_edges + (t - 1) * ns + i
            id = (ec[] += 1)
            place!(id, v_coup_ij, 1.0)
            place!(id, v_coup_il, -1.0)
            place!(id, v_edge, -1.0)
            gval[id] = [0.0]
        end
    end

    # Row marginal: Σ_j π_ij = 1/ns
    for i in 1:ns
        id = (ec[] += 1)
        for j in 1:nt
            place!(id, (i - 1) * nt + j, 1.0)
        end
        gval[id] = [1.0 / ns]
    end

    # Column marginal: Σ_i π_ij = 1/nt
    for j in 1:nt
        id = (ec[] += 1)
        for i in 1:ns
            place!(id, (i - 1) * nt + j, 1.0)
        end
        gval[id] = [1.0 / nt]
    end

    B = blocksparse(row_ids, col_ids, blks)
    g = zeros(size(B, 1))
    for (e, v) in gval
        g[rowrange(B, e)] .= v
    end

    # Build Q (diagonal on edge cells) and c (linear cost)
    Q = IPM.allocblockdiag(B); fill!(Q, 0)
    c = zeros(size(B, 2))

    # Linear cost on coupling cells
    for i in 1:ns, j in 1:nt
        v = (i - 1) * nt + j
        c[colrange(B, v)] .= M[i, j]
        block(Q, v, v, v) .= reshape([ε], 1, 1)  # tiny regularization
    end

    # Source edge Q: 2 * α * η * ||Xt_j||² * w_{ik}
    for (t, (i, k, w_ik)) in enumerate(Es)
        for j in 1:nt
            v_edge = nv_coupling + (t - 1) * nt + j
            qval = 2.0 * α * η * Xt_norm2[j] * w_ik
            block(Q, v_edge, v_edge, v_edge) .= reshape([qval + ε], 1, 1)
        end
    end

    # Target edge Q: 2 * (1-α) * η * ||Xs_i||² * w_{jl}
    for (t, (j, l, w_jl)) in enumerate(Et)
        for i in 1:ns
            v_edge = nv_coupling + nv_source_edges + (t - 1) * ns + i
            qval = 2.0 * (1.0 - α) * η * Xs_norm2[i] * w_jl
            block(Q, v_edge, v_edge, v_edge) .= reshape([qval + ε], 1, 1)
        end
    end

    # Cones: coupling cells are PositiveCone, edge cells are CofreeCone
    cones = AbstractCone[]
    for v in 1:nv_coupling
        push!(cones, PositiveCone())
    end
    for v in 1:(nv_source_edges + nv_target_edges)
        push!(cones, CofreeCone())
    end

    prob = IPMProblem(Q, B, c, g, cones)
    info = (ns = ns, nt = nt, nv = nv, nv_coupling = nv_coupling,
            n_source_edges = length(Es), n_target_edges = length(Et),
            cupcols = [colrange(B, v) for v in 1:nv_coupling])
    return prob, info
end

"""Lifted Laplacian demo: compare dense vs lifted formulations."""
function wg_lifted_laplace_demo(; ns::Int = 20, nt::Int = 20, kwargs...)
    # Dense version
    inst_dense = wg_laplace_instance(; ns, nt, kwargs...)
    prob_dense, info_dense = build_wg_problem(inst_dense)
    settings = wg_settings()
    res_dense = solve(prob_dense, settings)
    t_dense = @elapsed solve(prob_dense, settings)

    # Lifted version
    prob_lifted, info_lifted = build_lifted_laplace_problem(; ns, nt, kwargs...)
    res_lifted = solve(prob_lifted, settings)
    t_lifted = @elapsed solve(prob_lifted, settings)

    @printf("Laplacian-OT %dx%d:\n", ns, nt)
    @printf("  Dense:  |V|=%d  DOF=%d  obj=%.6f  %.1fms  it=%d %s\n",
            1, size(prob_dense.B, 2), prob_dense.c' * res_dense.p,
            t_dense * 1000, res_dense.niter, res_dense.status)
    @printf("  Lifted: |V|=%d  DOF=%d  obj=%.6f  %.1fms  it=%d %s\n",
            info_lifted.nv, size(prob_lifted.B, 2), prob_lifted.c' * res_lifted.p,
            t_lifted * 1000, res_lifted.niter, res_lifted.status)
    @printf("  Source edges: %d, Target edges: %d\n",
            info_lifted.n_source_edges, info_lifted.n_target_edges)
    return (dense = (prob_dense, res_dense), lifted = (prob_lifted, res_lifted))
end

# ---- LIFTED quadreg: split diagonal Q into scalar cells ----------------------

"""LIFTED quadratic-regularized OT: each π_ij is its own 1D cell with Q=η.
Transforms diagonal Q in one big cell → scalar Q in ns*nt cells."""
function build_lifted_quadreg_problem(; n::Int = 25, η::Float64 = 0.5, seed::Int = 7, ε::Float64 = 1e-9)
    rng = MersenneTwister(seed)
    xs = collect(range(0, 1, length = n))
    ν1 = abs.(randn(rng, n)) .+ 0.05; ν1 ./= sum(ν1)
    ν2 = abs.(randn(rng, n)) .+ 0.05; ν2 ./= sum(ν2)
    M = quad_cost(xs, xs)

    # Each coupling entry is a 1D cell
    nv = n * n
    vdims = ones(Int, nv)

    row_ids = Int[]; col_ids = Int[]; blks = Matrix{Float64}[]
    place!(e, v, m) = (push!(row_ids, e); push!(col_ids, v); push!(blks, reshape([m], 1, 1)))
    gval = Dict{Int, Vector{Float64}}()
    ec = Ref(0)

    # Row marginals: Σ_j π_ij = ν1_i
    for i in 1:n
        id = (ec[] += 1)
        for j in 1:n
            place!(id, (i - 1) * n + j, 1.0)
        end
        gval[id] = [ν1[i]]
    end

    # Column marginals: Σ_i π_ij = ν2_j
    for j in 1:n
        id = (ec[] += 1)
        for i in 1:n
            place!(id, (i - 1) * n + j, 1.0)
        end
        gval[id] = [ν2[j]]
    end

    B = blocksparse(row_ids, col_ids, blks)
    g = zeros(size(B, 1))
    for (e, v) in gval
        g[rowrange(B, e)] .= v
    end

    # Q = η on each cell (diagonal!)
    Q = IPM.allocblockdiag(B); fill!(Q, 0)
    c = zeros(size(B, 2))
    for i in 1:n, j in 1:n
        v = (i - 1) * n + j
        c[colrange(B, v)] .= M[i, j]
        block(Q, v, v, v) .= reshape([η + ε], 1, 1)
    end

    cones = AbstractCone[PositiveCone() for _ in 1:nv]
    prob = IPMProblem(Q, B, c, g, cones)
    info = (n = n, nv = nv, h1 = 1)
    return prob, info
end

# ---- H¹ law (measured; see header) --------------------------------------

function wg_h1_predict(inst::WGInstance)
    isempty(inst.free) && return length(inst.edges)
    Ainc = zeros(length(inst.free), length(inst.edges))
    fi = Dict(k => t for (t, k) in enumerate(inst.free))
    for (t, (i, j)) in enumerate(inst.edges)
        haskey(fi, i) && (Ainc[fi[i], t] += 1.0)
        haskey(fi, j) && (Ainc[fi[j], t] -= 1.0)
    end
    return length(inst.edges) - rank(Ainc)
end

# ---- pairwise (sheaf) builder --------------------------------------------

rowmarg(ni, nj) = kron(ones(1, nj), Matrix{Float64}(I, ni, ni))
colmarg(ni, nj) = kron(Matrix{Float64}(I, nj, nj), ones(1, ni))

"""Coupling cell per cost-graph edge (PositiveCone, vec column-major),
free-marginal cells (PositiveCone), two agreement base-edges per
coupling; pinned marginals in g."""
function build_wg_problem(inst::WGInstance)
    dims = length.(inst.supports)
    row_ids = Int[]; col_ids = Int[]; blks = Matrix{Float64}[]
    place!(e, v, mat) = (push!(row_ids, e); push!(col_ids, v); push!(blks, Matrix{Float64}(mat)))

    nE = length(inst.edges)
    vμ = Dict(k => nE + t for (t, k) in enumerate(inst.free))
    ec = Ref(0)
    pinval = Dict{Int, Vector{Float64}}()          # base-edge id -> g block
    for (t, (i, j)) in enumerate(inst.edges)
        ni, nj = dims[i], dims[j]
        e = (ec[] += 1)
        place!(e, t, rowmarg(ni, nj))
        haskey(inst.pins, i) ? (pinval[e] = inst.pins[i]) : place!(e, vμ[i], -Matrix{Float64}(I, ni, ni))
        e = (ec[] += 1)
        place!(e, t, colmarg(ni, nj))
        haskey(inst.pins, j) ? (pinval[e] = inst.pins[j]) : place!(e, vμ[j], -Matrix{Float64}(I, nj, nj))
    end

    B = blocksparse(row_ids, col_ids, blks)
    g = zeros(size(B, 1))
    for (e, v) in pinval
        g[rowrange(B, e)] .= v
    end
    nv = nE + length(inst.free)
    @assert nvtxs(B) == nv
    Q = IPM.allocblockdiag(B); fill!(Q, 0)
    c = zeros(size(B, 2))
    for v in 1:nv
        d = length(colrange(B, v))
        block(Q, v, v, v) .= (inst.ε + (v ≤ nE ? inst.η : 0.0)) .* Matrix{Float64}(I, d, d)
        if v ≤ nE
            haskey(inst.Qc, v) && (block(Q, v, v, v) .+= inst.Qc[v])
            c[colrange(B, v)] .= vec(inst.C[v])
            haskey(inst.cx, v) && (c[colrange(B, v)] .+= inst.cx[v])
        end
    end
    for (k, v) in vμ
        haskey(inst.Qμ, k) && (block(Q, v, v, v) .+= inst.Qμ[k])
    end
    prob = IPMProblem(Q, B, c, g, AbstractCone[PositiveCone() for _ in 1:nv])
    info = (cupcols = Dict(t => colrange(B, t) for t in 1:nE),
            mucols = Dict(k => colrange(B, v) for (k, v) in vμ),
            h1 = wg_h1_predict(inst), nv = nv)
    return prob, info
end

# ---- exact graphical MMOT as the ONE-CELL degenerate sheaf ----------------

function cost_tensor(inst::WGInstance)
    dims = Tuple(length.(inst.supports))
    T = zeros(dims)
    for (t, (i, j)) in enumerate(inst.edges), idx in CartesianIndices(dims)
        T[idx] += inst.C[t][idx[i], idx[j]]
    end
    return T
end

"""Joint tensor as a single PositiveCone cell, one marginalization pin
edge per node — the transport analogue of map_sdp's single-clique Shor
block (one big region, nothing to glue). Sizes: Πn (≤ ~10⁴ sensible)."""
function build_wg_mmot_problem(inst::WGInstance)
    @assert isempty(inst.free)
    dims = Tuple(length.(inst.supports))
    N = prod(dims)
    row_ids = Int[]; col_ids = Int[]; blks = Matrix{Float64}[]
    for k in 1:length(dims)
        M = zeros(dims[k], N)
        for (lin, idx) in enumerate(CartesianIndices(dims))
            M[idx[k], lin] = 1.0
        end
        push!(row_ids, k); push!(col_ids, 1); push!(blks, M)
    end
    B = blocksparse(row_ids, col_ids, blks)
    g = vcat((inst.pins[k] for k in 1:length(dims))...)
    Q = IPM.allocblockdiag(B); fill!(Q, 0)
    block(Q, 1, 1, 1) .= inst.ε .* Matrix{Float64}(I, N, N)
    prob = IPMProblem(Q, B, vec(cost_tensor(inst)), g, AbstractCone[PositiveCone()])
    return prob, (N = N,)
end

# ---- entropic dial ---------------------------------------------------------

"""Pairwise + (1/η)Σ π log π: one ExponentialCone leaf (x,y,z) per
coupling coordinate — library cone {(x,y,z): x,y>0, y·log(x/y) ≥ z},
bound LAST — wired by a 2-row edge [x = 1 (g); y − e_abᵀπ = 0], with
objective −1/η on z: at the optimum z = π·log(1/π) = −π log π (epigraph
tight), the mle_spline.jl pattern with the roles of x and y swapped."""
function build_wg_entropic_problem(inst::WGInstance; η::Float64 = 10.0)
    dims = length.(inst.supports)
    row_ids = Int[]; col_ids = Int[]; blks = Matrix{Float64}[]
    place!(e, v, mat) = (push!(row_ids, e); push!(col_ids, v); push!(blks, Matrix{Float64}(mat)))

    nE = length(inst.edges)
    vμ = Dict(k => nE + t for (t, k) in enumerate(inst.free))
    vc = Ref(nE + length(inst.free))
    ec = Ref(0)
    pinval = Dict{Int, Vector{Float64}}()
    leafz = Int[]                                   # vertex ids of the leaves
    for (t, (i, j)) in enumerate(inst.edges)
        ni, nj = dims[i], dims[j]
        e = (ec[] += 1); place!(e, t, rowmarg(ni, nj))
        haskey(inst.pins, i) ? (pinval[e] = inst.pins[i]) : place!(e, vμ[i], -Matrix{Float64}(I, ni, ni))
        e = (ec[] += 1); place!(e, t, colmarg(ni, nj))
        haskey(inst.pins, j) ? (pinval[e] = inst.pins[j]) : place!(e, vμ[j], -Matrix{Float64}(I, nj, nj))
        for ab in 1:(ni * nj)                       # entropy leaves
            v = (vc[] += 1); push!(leafz, v)
            e = (ec[] += 1)
            sel = zeros(1, ni * nj); sel[ab] = 1.0
            place!(e, v, [1.0 0.0 0.0; 0.0 1.0 0.0])    # rows: x = 1 ; y − π_ab = 0
            place!(e, t, [zeros(1, ni * nj); -sel])
            pinval[e] = [1.0, 0.0]
        end
    end

    B = blocksparse(row_ids, col_ids, blks)
    g = zeros(size(B, 1))
    for (e, v) in pinval
        g[rowrange(B, e)] .= v
    end
    nv = vc[]
    Q = IPM.allocblockdiag(B); fill!(Q, 0)
    c = zeros(size(B, 2))
    cones = AbstractCone[]
    for v in 1:nv
        d = length(colrange(B, v))
        block(Q, v, v, v) .= inst.ε .* Matrix{Float64}(I, d, d)
        if v ≤ nE
            c[colrange(B, v)] .= vec(inst.C[v]); push!(cones, PositiveCone())
        elseif v ≤ nE + length(inst.free)
            push!(cones, PositiveCone())
        else
            c[colrange(B, v)[3]] = -1.0 / η; push!(cones, ExponentialCone())
        end
    end
    prob = IPMProblem(Q, B, c, g, cones)
    info = (cupcols = Dict(t => colrange(B, t) for t in 1:nE),
            mucols = Dict(k => colrange(B, v) for (k, v) in vμ), nv = nv)
    return prob, info
end

# ---- demos -------------------------------------------------------------------

wg_settings(; raug = 1e5) = IPMSettings{Float64}(
    kkt = UzawaSettings{Float64}(raug = raug),
    feas_tol = 1e-8, gap_tol = 1e-8, itmax = 200)

wg_objective(inst, prob, res) = prob.c' * res.p          # linear part (ε reported separately)

"""McCann midpoint: barycenter of δ(.25) and δ(.75) with λ = (½,½) on a
grid containing 0.5 is δ(.5)."""
function wg_barycenter_demo(; n::Int = 21)
    xs = collect(range(0, 1, length = n))
    d1 = zeros(n); d1[6] = 1.0
    d2 = zeros(n); d2[16] = 1.0
    inst = wg_star_barycenter([d1, d2], xs, [0.5, 0.5])
    prob, info = build_wg_problem(inst)
    res = solve(prob, wg_settings())
    μ = res.p[info.mucols[1]]
    @printf("barycenter star (L=2 point masses):  H1=%d  it=%d %s\n",
            info.h1, res.niter, res.status)
    @printf("  center mass at midpoint 0.5: %.6f  (McCann: 1.0)\n", μ[11])
    return inst, res, info
end

"""Pairwise vs exact MMOT: the binary triangle (0 vs 1, map_lp verbatim)
and the quadratic holonomy family (gap iff nanti odd)."""
function wg_gap_demo(; n::Int = 8)
    for (name, inst) in (("binary triangle       ", wg_binary_triangle()),
                         ("K=3 quad, nanti=1     ", wg_cycle_instance(; K = 3, n, nanti = 1)),
                         ("K=3 quad, nanti=2     ", wg_cycle_instance(; K = 3, n, nanti = 2)))
        prob, info = build_wg_problem(inst)
        res = solve(prob, wg_settings())
        probx, _ = build_wg_mmot_problem(inst)
        resx = solve(probx, wg_settings())
        vp, ve = wg_objective(inst, prob, res), wg_objective(inst, probx, resx)
        @printf("%s pairwise %+.6f   exact %+.6f   gap %.6f   (H1=%d)\n",
                name, vp, ve, ve - vp, info.h1)
    end
end

"""Entropic pairwise on a pinned-ends path; middle marginal entropy
should exceed the multimarginal value (oracle: 2.51 vs 2.24 at η=10) —
the regularization gap that exists ON TREES."""
function wg_entropic_demo(; n::Int = 15, η::Float64 = 10.0)
    xs = collect(range(0, 1, length = n))
    ν1 = exp.(-((xs .- 0.25) ./ 0.08) .^ 2); ν1 ./= sum(ν1)
    ν2 = exp.(-((xs .- 0.75) ./ 0.08) .^ 2); ν2 ./= sum(ν2)
    inst = wg_instance(fill(xs, 3), [(1, 2), (2, 3)],
                       [quad_cost(xs, xs), quad_cost(xs, xs)],
                       Dict(1 => ν1, 3 => ν2))
    prob, info = build_wg_entropic_problem(inst; η)
    res = solve(prob, wg_settings())
    μ = res.p[info.mucols[2]]
    H = -sum(v > 1e-12 ? v * log(v) : 0.0 for v in μ)
    @printf("entropic path η=%.0f: middle-marginal entropy %.4f  it=%d %s\n",
            η, H, res.niter, res.status)
    @printf("  (oracle: pairwise 2.5083 vs multimarginal 2.2390 at η=10)\n")
    return inst, res, info
end

"""Quadratically regularized OT (first-class η): sparse plans where
entropic is dense; support grows with η toward the min-norm interior."""
function wg_quadreg_demo(; n::Int = 25, seed::Int = 7)
    rng = MersenneTwister(seed)
    xs = collect(range(0, 1, length = n))
    ν = [abs.(randn(rng, n)) .+ 0.05 for _ in 1:2]
    pins = Dict(1 => ν[1] ./ sum(ν[1]), 2 => ν[2] ./ sum(ν[2]))
    for η in (0.0, 0.05, 0.5)
        inst = wg_instance([xs, xs], [(1, 2)], [quad_cost(xs, xs)], pins; η)
        prob, info = build_wg_problem(inst)
        res = solve(prob, wg_settings())
        s = count(>(1e-6), res.p[info.cupcols[1]])
        @printf("  η=%.2f: support %d / %d   it=%d %s\n", η, s, n^2, res.niter, res.status)
    end
    println("  (oracle: LP 49, quad η=0.5: 83, entropic η=20: 519 of 625)")
end

"""Laplacian-regularized OT (Courty et al.): dense Kronecker Gram on the
coupling cell — the corpus's first structured Q on an orthant cell."""
function wg_laplace_demo(; kwargs...)
    inst = wg_laplace_instance(; kwargs...)
    prob, info = build_wg_problem(inst)
    res = solve(prob, wg_settings())
    π = res.p[info.cupcols[1]]
    Ω = 0.5 * π' * (inst.Qc[1] * π)
    @printf("laplace-OT: ⟨M,π⟩ = %.6f   ½πᵀQπ = %.6f   it=%d %s\n",
            prob.c' * res.p, Ω, res.niter, res.status)
    @printf("  (oracle: our QP 1.80278887 ≤ POT CG 1.80281443)\n")
    return inst, res, info
end

"""Barycenter with H¹ center prior: roughness falls, support spreads."""
function wg_smooth_barycenter_demo(; n::Int = 21)
    xs = collect(range(0, 1, length = n))
    d1 = zeros(n); d1[4] = 0.5; d1[9] = 0.5
    d2 = zeros(n); d2[13] = 0.5; d2[18] = 0.5
    for λ in (0.0, 0.02, 0.2)
        inst = λ == 0 ? wg_star_barycenter([d1, d2], xs, [0.5, 0.5]) :
                        wg_smooth_barycenter([d1, d2], xs, [0.5, 0.5]; λsmooth = λ)
        prob, info = build_wg_problem(inst)
        res = solve(prob, wg_settings())
        μ = res.p[info.mucols[1]]
        @printf("  λ=%.2f: roughness Σ(Δμ)² = %.5f   support %d/%d\n",
                λ, sum(abs2, diff(μ)), count(>(1e-6), μ), n)
    end
end

# =====================================================================
# JuMP reference
# =====================================================================

using JuMP

function build_jump_wg(inst::WGInstance, optimizer)
    dims = length.(inst.supports)
    model = Model(optimizer); set_silent(model)
    μ = Dict{Int, Any}(k => inst.pins[k] for k in keys(inst.pins))
    for k in inst.free
        μ[k] = @variable(model, [1:dims[k]], lower_bound = 0)
    end
    π = [@variable(model, [1:dims[i], 1:dims[j]], lower_bound = 0)
         for (i, j) in inst.edges]
    for (t, (i, j)) in enumerate(inst.edges)
        @constraint(model, vec(sum(π[t], dims = 2)) .== μ[i])
        @constraint(model, vec(sum(π[t], dims = 1)) .== μ[j])
    end
    vc(t) = vec(π[t])
    @objective(model, Min,
               sum(sum(inst.C[t] .* π[t]) +
                   (haskey(inst.cx, t) ? inst.cx[t]' * vc(t) : 0.0) +
                   0.5 * (inst.ε + inst.η) * sum(π[t][a, b]^2 for a in 1:size(π[t], 1),
                                                                  b in 1:size(π[t], 2)) +
                   (haskey(inst.Qc, t) ? 0.5 * vc(t)' * inst.Qc[t] * vc(t) : 0.0)
                   for t in 1:length(inst.edges)) +
               sum(0.5 * inst.ε * sum(μ[k][a]^2 for a in 1:dims[k]) +
                   (haskey(inst.Qμ, k) ? 0.5 * μ[k]' * inst.Qμ[k] * μ[k] : 0.0)
                   for k in inst.free; init = 0.0))
    return model, π, μ
end

"""JuMP model for lifted Laplacian OT (diagonal Q on auxiliary edge variables)."""
function build_jump_lifted_laplace(; ns::Int = 20, nt::Int = 20, seed::Int = 7,
                                   α::Float64 = 0.5, η::Float64 = 0.4, σ::Float64 = 0.8,
                                   edge_thresh::Float64 = 0.1, ε::Float64 = 1e-9, optimizer = nothing)
    rng = MersenneTwister(seed)
    Xs = randn(rng, ns, 2)
    Xt = randn(rng, nt, 2) .+ [1.5 0.5]
    M = [sum(abs2, Xs[i, :] .- Xt[j, :]) for i in 1:ns, j in 1:nt]
    Ss = wg_gauss_sim(Xs, σ)
    St = wg_gauss_sim(Xt, σ)
    Es = laplacian_edges(Ss; thresh = edge_thresh)
    Et = laplacian_edges(St; thresh = edge_thresh)
    Xs_norm2 = [sum(abs2, Xs[i, :]) for i in 1:ns]
    Xt_norm2 = [sum(abs2, Xt[j, :]) for j in 1:nt]

    model = Model(optimizer); set_silent(model)
    @variable(model, π[1:ns, 1:nt] >= 0)
    # Source edge variables: y^s_{(ik),j} = π_ij - π_kj
    @variable(model, ys[1:length(Es), 1:nt])
    # Target edge variables: y^t_{(jl),i} = π_ij - π_il
    @variable(model, yt[1:length(Et), 1:ns])

    # Edge definitions
    for (t, (i, k, _)) in enumerate(Es)
        for j in 1:nt
            @constraint(model, ys[t, j] == π[i, j] - π[k, j])
        end
    end
    for (t, (j, l, _)) in enumerate(Et)
        for i in 1:ns
            @constraint(model, yt[t, i] == π[i, j] - π[i, l])
        end
    end

    # Marginals: uniform
    @constraint(model, [i in 1:ns], sum(π[i, :]) == 1.0 / ns)
    @constraint(model, [j in 1:nt], sum(π[:, j]) == 1.0 / nt)

    # Objective: linear cost + diagonal quadratics on edge variables
    @objective(model, Min,
        sum(M[i, j] * π[i, j] + 0.5 * ε * π[i, j]^2 for i in 1:ns, j in 1:nt) +
        sum(0.5 * (2.0 * α * η * Xt_norm2[j] * Es[t][3] + ε) * ys[t, j]^2 for t in 1:length(Es), j in 1:nt) +
        sum(0.5 * (2.0 * (1.0 - α) * η * Xs_norm2[i] * Et[t][3] + ε) * yt[t, i]^2 for t in 1:length(Et), i in 1:ns))

    info = (ns = ns, nt = nt)
    return model, info
end

"""JuMP model for lifted quadreg OT (diagonal Q = η on each coupling entry)."""
function build_jump_lifted_quadreg(; n::Int = 25, η::Float64 = 0.5, seed::Int = 7,
                                   ε::Float64 = 1e-9, optimizer = nothing)
    rng = MersenneTwister(seed)
    xs = collect(range(0, 1, length = n))
    ν1 = abs.(randn(rng, n)) .+ 0.05; ν1 ./= sum(ν1)
    ν2 = abs.(randn(rng, n)) .+ 0.05; ν2 ./= sum(ν2)
    M = quad_cost(xs, xs)

    model = Model(optimizer); set_silent(model)
    @variable(model, π[1:n, 1:n] >= 0)

    @constraint(model, [i in 1:n], sum(π[i, :]) == ν1[i])
    @constraint(model, [j in 1:n], sum(π[:, j]) == ν2[j])

    @objective(model, Min, sum(M[i, j] * π[i, j] + 0.5 * (η + ε) * π[i, j]^2 for i in 1:n, j in 1:n))

    info = (n = n,)
    return model, info
end

# ---- benchmark -----------------------------------------------------------

function run_wg_benchmark(; optimizer = nothing, dual_optimizer = nothing,
                          nwarmup::Int = 2, nruns::Int = 5, solver_name::String = "Ref")
    optimizer === nothing && error("Pass optimizer, e.g. run_wg_benchmark(optimizer=Mosek.Optimizer)")
    # Standard cases (inst_fn, label, raug) — only truly pairwise problems
    cases = [
        (() -> wg_binary_triangle(), "tri LP", 1e5),
        (() -> wg_cycle_instance(; K = 4, n = 15), "4-cyc LP", 1e4),
        (() -> wg_smooth_barycenter([fill(1/20,20), fill(1/20,20), fill(1/20,20)],
                                    collect(range(0,1,20)), [0.33,0.33,0.34]; λsmooth = 0.1),
         "smooth bary", 1e5),
    ]
    # Lifted cases: diagonal Q decomposed into scalar cells
    lifted_cases = [
        (() -> build_lifted_quadreg_problem(; n = 25, η = 0.5), "quadreg 25 (lifted)", 1e3),
        (() -> build_lifted_laplace_problem(; ns = 15, nt = 15, edge_thresh = 0.3), "laplace 15 (lifted)", 1e2),
        (() -> build_lifted_laplace_problem(; ns = 20, nt = 20, edge_thresh = 0.4), "laplace 20 (lifted)", 1e2),
    ]
    println("="^95)
    @printf("Wasserstein Graph Benchmark: Sheaf IPM vs %s\n", solver_name)
    println("="^95)
    if dual_optimizer !== nothing
        @printf("%-20s %6s %5s %5s %5s %9s %9s %9s %7s %7s\n",
                "Config", "raug", "DOF", "|V|", "H1", "IPM(ms)", solver_name, solver_name * "-D", "P/IPM", "D/IPM")
    else
        @printf("%-20s %6s %5s %5s %5s %9s %9s %8s\n",
                "Config", "raug", "DOF", "|V|", "H1", "IPM(ms)", solver_name, "speedup")
    end
    println("-"^95)
    for (inst_fn, label, raug) in cases
        inst = inst_fn()
        prob, info = build_wg_problem(inst)
        settings = wg_settings(; raug)
        for _ in 1:nwarmup; solve(prob, settings); end
        t_ipm = minimum([@elapsed solve(prob, settings) for _ in 1:nruns])

        for _ in 1:nwarmup
            m, _, _ = build_jump_wg(inst, optimizer); optimize!(m)
        end
        t_mosek = minimum([@elapsed begin
            m, _, _ = build_jump_wg(inst, optimizer); optimize!(m)
        end for _ in 1:nruns])

        if dual_optimizer !== nothing
            for _ in 1:nwarmup
                m, _, _ = build_jump_wg(inst, dual_optimizer); optimize!(m)
            end
            t_dual = minimum([@elapsed begin
                m, _, _ = build_jump_wg(inst, dual_optimizer); optimize!(m)
            end for _ in 1:nruns])
            @printf("%-20s %6.0e %5d %5d %5d %9.1f %9.1f %9.1f %6.2fx %6.2fx\n",
                    label, raug, size(prob.B, 2), info.nv, info.h1,
                    t_ipm * 1000, t_mosek * 1000, t_dual * 1000,
                    t_mosek / t_ipm, t_dual / t_ipm)
        else
            @printf("%-20s %6.0e %5d %5d %5d %9.1f %9.1f %7.2fx\n",
                    label, raug, size(prob.B, 2), info.nv, info.h1,
                    t_ipm * 1000, t_mosek * 1000, t_mosek / t_ipm)
        end
    end

    # Lifted cases: IPM and Mosek both on lifted formulation
    for (build_fn, label, raug) in lifted_cases
        prob, info = build_fn()
        nv = info.nv
        h1 = haskey(info, :h1) ? info.h1 : 1
        settings = wg_settings(; raug)
        for _ in 1:nwarmup; solve(prob, settings); end
        t_ipm = minimum([@elapsed solve(prob, settings) for _ in 1:nruns])

        # Build lifted JuMP model (same structure)
        if haskey(info, :ns)  # laplace case
            ns, nt = info.ns, info.nt
            for _ in 1:nwarmup
                m, _ = build_jump_lifted_laplace(; ns, nt, edge_thresh = 0.3, optimizer); optimize!(m)
            end
            t_mosek = minimum([@elapsed begin
                m, _ = build_jump_lifted_laplace(; ns, nt, edge_thresh = 0.3, optimizer); optimize!(m)
            end for _ in 1:nruns])
            if dual_optimizer !== nothing
                for _ in 1:nwarmup
                    m, _ = build_jump_lifted_laplace(; ns, nt, edge_thresh = 0.3, optimizer = dual_optimizer); optimize!(m)
                end
                t_dual = minimum([@elapsed begin
                    m, _ = build_jump_lifted_laplace(; ns, nt, edge_thresh = 0.3, optimizer = dual_optimizer); optimize!(m)
                end for _ in 1:nruns])
            else
                t_dual = Inf
            end
        else  # quadreg case
            n = info.n
            for _ in 1:nwarmup
                m, _ = build_jump_lifted_quadreg(; n, optimizer); optimize!(m)
            end
            t_mosek = minimum([@elapsed begin
                m, _ = build_jump_lifted_quadreg(; n, optimizer); optimize!(m)
            end for _ in 1:nruns])
            if dual_optimizer !== nothing
                for _ in 1:nwarmup
                    m, _ = build_jump_lifted_quadreg(; n, optimizer = dual_optimizer); optimize!(m)
                end
                t_dual = minimum([@elapsed begin
                    m, _ = build_jump_lifted_quadreg(; n, optimizer = dual_optimizer); optimize!(m)
                end for _ in 1:nruns])
            else
                t_dual = Inf
            end
        end

        if dual_optimizer !== nothing && t_dual < Inf
            @printf("%-20s %6.0e %5d %5d %5d %9.1f %9.1f %9.1f %6.2fx %6.2fx\n",
                    label, raug, size(prob.B, 2), nv, h1,
                    t_ipm * 1000, t_mosek * 1000, t_dual * 1000,
                    t_mosek / t_ipm, t_dual / t_ipm)
        else
            @printf("%-20s %6.0e %5d %5d %5d %9.1f %9.1f %7.2fx\n",
                    label, raug, size(prob.B, 2), nv, h1,
                    t_ipm * 1000, t_mosek * 1000, t_mosek / t_ipm)
        end
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    include(joinpath(@__DIR__, "..", "benchmark_utils.jl"))
    opts = parse_benchmark_args(ARGS)

    if opts.mosek
        using MosekTools
        optimizer = Mosek.Optimizer
        dual_optimizer = Dualization.dual_optimizer(Mosek.Optimizer)
        solver_name = "Mosek"
        println("Solver: Mosek + Mosek-D")
    else
        using OSQP
        optimizer = OSQP.Optimizer
        dual_optimizer = Dualization.dual_optimizer(OSQP.Optimizer)
        solver_name = "OSQP"
        println("Solver: OSQP + OSQP-D (open-source)")
    end
    println("Runs: $(opts.nruns), Warmup: $(opts.nwarmup)\n")

    run_wg_benchmark(; optimizer, dual_optimizer, solver_name,
                     nwarmup = opts.nwarmup, nruns = opts.nruns)
end
