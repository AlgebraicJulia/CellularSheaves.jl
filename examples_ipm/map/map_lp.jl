######################################################################
# map_lp.jl
#
# A STOCHASTIC-KERNEL CELLULAR SHEAF solved as a conic LP/QP in the
# standard form of CellularSheaves.IPM:
#
#       min   ½ μᵀ Q μ + cᵀ μ            (Q = 0 ⇒ LP;  Q = ρI ⇒ SparseMAP)
#       s.t.  B μ = g                    (B = coboundary δ of a sheaf)
#             μ ∈ K = ∏_v K_v            (here every K_v = PositiveCone)
#
# The base object is a cellular sheaf whose RESTRICTION MAPS are general
# COLUMN-STOCHASTIC matrices (Markov kernels / channels), and whose
# stalks are ordered by the nonnegative cone.  A "positive global
# section" is a family of per-cell distributions that agree, after each
# is pushed through its channel, on every shared edge stalk.
#
# ---------------------------------------------------------------------
# TWO INDEPENDENT DIALS
# ---------------------------------------------------------------------
#   • OBJECTIVE axis:  Q = 0  (linear, the MAP-LP / local polytope)   or
#                      Q = ρI (native quadratic — SparseMAP: the
#                      Euclidean projection of −c/ρ onto the polytope,
#                      giving a UNIQUE, SPARSE, differentiable optimum).
#                      Quadratic is first-class: it rides in the
#                      objective, it is NOT lifted to a rotated-SOC cone,
#                      so δ's bipartite sparsity is preserved (Q only
#                      thickens each vertex's (1,1) KKT block).
#   • CONE axis:       K_v.  PositiveCone here (simplex with the pins);
#                      swap for SemidefiniteCone/ExponentialCone to get
#                      the moment/quantum-marginal/entropic siblings in
#                      the SAME assembly (see the note at end of file).
#
# ---------------------------------------------------------------------
# RESTRICTION MAPS = MARKOV KERNELS
# ---------------------------------------------------------------------
# A base-edge e joins cells u,w with a shared edge stalk ℝ^{mₑ} and two
# COLUMN-STOCHASTIC maps  Pₑ : ℝ^{nᵤ}→ℝ^{mₑ},  Rₑ : ℝ^{n_w}→ℝ^{mₑ}.
# The coboundary block is  (δμ)ₑ = Pₑ μᵤ − Rₑ μ_w,  so δμ = g demands the
# two cells agree after each is pushed through its channel.  Column-
# stochastic = nonnegative + columns sum to 1 = it maps distributions to
# distributions (a positive map on the ordered stalk), so the pushforward
# of a positive section stays positive.  Special cases:
#   • MARGINALISATION  (0/1 stochastic): the MRF local polytope 𝓛(𝔾).
#     Restriction {i}⊂{i,j} is "sum out xⱼ"; the node side is the identity.
#   • SOFT COARSE-GRAINING / NOISY CHANNEL (general stochastic): cells
#     lump/observe fine states into a shared coarse observable and must
#     agree there — the empirical-model reading of Abramsky–Brandenburger.
#   • (PERMUTATION / doubly-stochastic relabelings give the synchronisation
#     flavour, where a cycle's holonomy can obstruct LINEARLY, before
#     positivity — a source of "no global section" that marginalisation,
#     being surjective and holonomy-free, never exhibits.)
#
# ---------------------------------------------------------------------
# THE PUNCHLINE (unchanged, now with a quadratic and general kernels)
# ---------------------------------------------------------------------
# On a frustrated triangle the SparseMAP optimum is the UNIQUE fractional
# vertex of 𝓛∖𝓜: it satisfies δμ = g exactly (agreement holds), it is
# nonnegative, yet it is the restriction of NO global distribution — the
# positive-global-section obstruction (contextuality) a chain can't show.
#
# References:
#   Wainwright, Jordan, "Graphical models, exponential families, and
#     variational inference", Found. Trends ML 2008.        (local polytope)
#   Niculae, Martins, Blondel, Cardie, "SparseMAP: Differentiable sparse
#     structured inference", ICML 2018.                     (native-quadratic MAP)
#   Abramsky, Brandenburger, "The sheaf-theoretic structure of non-locality
#     and contextuality", New J. Phys. 2011.        (no global section)
#   Hanks, Riess, Cohen, Gross, Hale, Fairbanks, arXiv:2504.02049 (2025).
######################################################################

using CellularSheaves.IPM
using CellularSheaves.BlockSparseArrays: colrange, rowrange, blocksparse, block, nvtxs
using LinearAlgebra
using Random
using Printf

# =====================================================================
# SECTION 1.  Markov-kernel restriction maps (column-stochastic)
# =====================================================================

"""
    is_col_stochastic(S; tol) :: Bool

A restriction map is admissible iff it is COLUMN-STOCHASTIC: nonnegative,
with each column summing to 1.  Such maps send distributions to
distributions (positive maps), so they carry positive sections to
positive sections — the requirement for an *ordered* cellular sheaf.
"""
function is_col_stochastic(S::AbstractMatrix; tol::Real = 1e-9)
    all(S .>= -tol) && all(abs.(vec(sum(S, dims = 1)) .- 1) .<= tol)
end

"""
    lin(a, b, cj)  ->  Int      # 1-based joint index (a-1)·cj + b
"""
lin(a::Integer, b::Integer, cj::Integer) = (a - 1) * cj + b

"""
    marg_out_j(ci, cj) :: Matrix   (cᵢ × cᵢcⱼ)   — DETERMINISTIC kernel

Marginalisation "sum out xⱼ": (Mᵢ μᵢⱼ)_a = Σ_b μᵢⱼ(a,b).  A 0/1
column-stochastic map — the deterministic special case of a Markov kernel.
"""
function marg_out_j(ci::Integer, cj::Integer)
    M = zeros(ci, ci * cj)
    for a in 1:ci, b in 1:cj
        M[a, lin(a, b, cj)] = 1.0
    end
    return M
end

"""
    marg_out_i(ci, cj) :: Matrix   (cⱼ × cᵢcⱼ)   — DETERMINISTIC kernel

Marginalisation "sum out xᵢ": (Mⱼ μᵢⱼ)_b = Σ_a μᵢⱼ(a,b).
"""
function marg_out_i(ci::Integer, cj::Integer)
    M = zeros(cj, ci * cj)
    for a in 1:ci, b in 1:cj
        M[b, lin(a, b, cj)] = 1.0
    end
    return M
end

"""
    random_channel(m, n, rng) :: Matrix   (m × n)   — GENERAL kernel

A random column-stochastic m×n matrix: a soft coarse-graining / noisy
channel that (probabilistically) maps n fine states to m coarse ones.
"""
function random_channel(m::Integer, n::Integer, rng::AbstractRNG)
    M = rand(rng, m, n)
    M ./= sum(M, dims = 1)
    return M
end

# =====================================================================
# SECTION 2.  Sheaf specification (cells, stochastic edges, pins, cost)
# =====================================================================

"""
A stochastic base-edge  Pμᵤ − Rμ_w = rhs  into a shared edge stalk.
Both P and R must be column-stochastic (validated on construction).
"""
struct StochEdge
    u::Int; Pu::Matrix{Float64}
    w::Int; Rw::Matrix{Float64}
    rhs::Vector{Float64}
    function StochEdge(u, Pu, w, Rw, rhs)
        @assert is_col_stochastic(Pu) "Pu is not column-stochastic (not a Markov kernel)"
        @assert is_col_stochastic(Rw) "Rw is not column-stochastic (not a Markov kernel)"
        @assert size(Pu, 1) == size(Rw, 1) == length(rhs) "edge-stalk dimensions disagree"
        new(u, Matrix{Float64}(Pu), w, Matrix{Float64}(Rw), Vector{Float64}(rhs))
    end
end

"""A degree-1 pin  A μᵥ = rhs  (e.g. normalisation 1ᵀμᵥ = 1)."""
struct Pin
    v::Int; A::Matrix{Float64}; rhs::Vector{Float64}
end

"""A conic sheaf problem: per-cell dim/cone/linear-cost, stochastic edges, pins."""
struct SheafSpec
    dim::Vector{Int}
    cone::Vector{AbstractCone}
    c::Vector{Vector{Float64}}          # linear cost per cell (length dim[v])
    edges::Vector{StochEdge}
    pins::Vector{Pin}
end

# ---- MRF local polytope (marginalisation kernels) -------------------

"""
A pairwise-MRF MAP instance.  `edges` are node pairs (i<j); `card[i]` the
cardinality of xᵢ; potentials `node_pot[i]∈ℝ^{card[i]}`,
`edge_pot[e]∈ℝ^{card[i]·card[j]}` (row-major in the first variable).
"""
struct MAPInstance
    V::Int
    edges::Vector{Tuple{Int,Int}}
    card::Vector{Int}
    node_pot::Vector{Vector{Float64}}
    edge_pot::Vector{Vector{Float64}}
end

"""
    frustrated_triangle(; penalty) :: MAPInstance

Three binary variables on a triangle; agreement on every edge costs
`penalty`.  The odd cycle is not 2-colourable ⇒ integer MAP pays ≥1;
the relaxation splits every edge onto the anti-diagonal (cost 0) — the
canonical 𝓛 ⊋ 𝓜 witness.
"""
function frustrated_triangle(; penalty::Float64 = 1.0)
    edges = [(1, 2), (2, 3), (1, 3)]; card = [2, 2, 2]
    ep = Float64[penalty, 0.0, 0.0, penalty]        # (1,1)(1,2)(2,1)(2,2)
    return MAPInstance(3, edges, card, [zeros(2) for _ in 1:3], [copy(ep) for _ in edges])
end

"""
    mrf_spec(inst) :: (SheafSpec, info)

Lower a pairwise MRF to a stochastic-kernel sheaf: node-belief and edge-
belief cells (PositiveCone), MARGINALISATION restriction maps (node side
= identity), and normalisation pins.  Edge normalisation is implied.
"""
function mrf_spec(inst::MAPInstance)
    V, edges, card = inst.V, inst.edges, inst.card
    dim = Int[]; cone = AbstractCone[]; cvec = Vector{Float64}[]
    nodev = Int[]; edgev = Int[]
    for i in 1:V
        push!(dim, card[i]); push!(cone, PositiveCone()); push!(cvec, inst.node_pot[i])
        push!(nodev, length(dim))
    end
    for (ei, (i, j)) in enumerate(edges)
        push!(dim, card[i] * card[j]); push!(cone, PositiveCone()); push!(cvec, inst.edge_pot[ei])
        push!(edgev, length(dim))
    end
    sedges = StochEdge[]; pins = Pin[]
    for (ei, (i, j)) in enumerate(edges)
        ci, cj = card[i], card[j]
        push!(sedges, StochEdge(edgev[ei], marg_out_j(ci, cj), nodev[i], Matrix{Float64}(I, ci, ci), zeros(ci)))
        push!(sedges, StochEdge(edgev[ei], marg_out_i(ci, cj), nodev[j], Matrix{Float64}(I, cj, cj), zeros(cj)))
    end
    for i in 1:V
        push!(pins, Pin(nodev[i], ones(1, card[i]), [1.0]))
    end
    return SheafSpec(dim, cone, cvec, sedges, pins), (nodev = nodev, edgev = edgev, edges = edges, card = card)
end

# ---- soft coarse-graining / channel sheaf (general kernels) ---------

"""
    coarse_graining_spec(; ncell, nfine, ncoarse, ring, seed) :: (SheafSpec, info)

`ncell` cells, each a distribution over `nfine` states, connected in a
path (or ring) by edges that demand agreement of the two cells' images
under their own random column-stochastic channel into a shared `ncoarse`
observable.  A genuinely stochastic (non-0/1) restriction map; feasible
(a positive global section exists) whenever the channels are surjective
onto the coarse simplex, which holds generically.
"""
function coarse_graining_spec(; ncell::Int = 3, nfine::Int = 3, ncoarse::Int = 2,
                              ring::Bool = true, seed::Int = 0)
    rng = MersenneTwister(seed)
    chan = [random_channel(ncoarse, nfine, rng) for _ in 1:ncell]
    dim  = fill(nfine, ncell)
    cone = AbstractCone[PositiveCone() for _ in 1:ncell]
    cvec = [zeros(nfine) for _ in 1:ncell]
    pairs = ring ? [(v, mod1(v + 1, ncell)) for v in 1:ncell] :
                   [(v, v + 1) for v in 1:(ncell - 1)]
    sedges = [StochEdge(u, chan[u], w, chan[w], zeros(ncoarse)) for (u, w) in pairs]
    pins   = [Pin(v, ones(1, nfine), [1.0]) for v in 1:ncell]
    return SheafSpec(dim, cone, cvec, sedges, pins), (channels = chan, pairs = pairs)
end

# =====================================================================
# SECTION 3.  Assemble the IPMProblem   (B = coboundary, Q = 0 or ρI)
# =====================================================================

"""
    build_ipm_problem(spec; rho) :: (IPMProblem, info)

Assemble a SheafSpec as a CellularSheaves.IPM problem.  `rho = 0` gives a
pure conic LP; `rho > 0` adds the native SparseMAP quadratic Q = ρI on
every cell (block-diagonal — no cone lift, δ's sparsity untouched).
"""
function build_ipm_problem(spec::SheafSpec; rho::Float64 = 0.0)
    nvtx = length(spec.dim)
    row_ids = Int[]; col_ids = Int[]; blocks = Matrix{Float64}[]
    rhs = Dict{Int, Vector{Float64}}()
    ecnt = Ref(0); newedge!() = (ecnt[] += 1; ecnt[])
    place!(e, v, M) = (push!(row_ids, e); push!(col_ids, v); push!(blocks, Matrix{Float64}(M)))

    # stochastic edges:  Pμᵤ − Rμ_w = rhs   (degree-2)
    for se in spec.edges
        e = newedge!()
        place!(e, se.u, se.Pu)
        place!(e, se.w, -se.Rw)
        any(!=(0.0), se.rhs) && (rhs[e] = copy(se.rhs))
    end
    # pins:  Aμᵥ = rhs   (degree-1)
    for pn in spec.pins
        e = newedge!()
        place!(e, pn.v, pn.A)
        rhs[e] = copy(pn.rhs)
    end

    B = blocksparse(row_ids, col_ids, blocks)
    @assert nvtxs(B) == nvtx "every cell must appear in B (no empty column block)"

    g = zeros(size(B, 1))
    for (e, v) in rhs
        g[rowrange(B, e)] .= v
    end

    Q = IPM.allocblockdiag(B)
    fill!(Q, 0)
    if rho != 0.0                                   # native SparseMAP Hessian ρI
        for v in 1:nvtx
            Qv = block(Q, v, v, v)
            Qv .= rho .* Matrix{Float64}(I, spec.dim[v], spec.dim[v])
        end
    end

    c = zeros(size(B, 2))
    for v in 1:nvtx
        c[colrange(B, v)] .= spec.c[v]
    end

    cones = AbstractCone[spec.cone[v] for v in 1:nvtx]
    prob  = IPMProblem(Q, B, c, g, cones)
    info  = (B = B, dim = spec.dim, rho = rho)
    return prob, info
end

"""    cellvals(res, info) — per-cell primal blocks read from res.p."""
cellvals(res, info) = [res.p[colrange(info.B, v)] for v in 1:length(info.dim)]

# =====================================================================
# SECTION 4.  Reference model — the SAME problem as a JuMP model
# (orthant cones ⇒ nonneg variables; quadratic objective for rho>0).
# =====================================================================

using JuMP

function build_jump_model(spec::SheafSpec, optimizer; rho::Float64 = 0.0)
    n = length(spec.dim)
    model = Model(optimizer); set_silent(model)
    x = [@variable(model, [1:spec.dim[v]], lower_bound = 0.0) for v in 1:n]
    lin_term = sum(spec.c[v]' * x[v] for v in 1:n)
    @objective(model, Min, rho == 0.0 ? lin_term :
                            lin_term + (rho / 2) * sum(x[v]' * x[v] for v in 1:n))
    for se in spec.edges
        @constraint(model, se.Pu * x[se.u] .- se.Rw * x[se.w] .== se.rhs)
    end
    for pn in spec.pins
        @constraint(model, pn.A * x[pn.v] .== pn.rhs)
    end
    return model, x
end

# =====================================================================
# SECTION 5.  Solve / verify
# =====================================================================

default_settings() = IPMSettings{Float64}(feas_tol = 1e-9, gap_tol = 1e-9, itmax = 200)

"""
    solve_spec(spec; rho, settings) :: (res, info)

Solve any SheafSpec and print a generic report (feasibility of the
coboundary, objective, active-set sparsity).
"""
function solve_spec(spec::SheafSpec; rho::Float64 = 0.0, settings = default_settings())
    prob, info = build_ipm_problem(spec; rho = rho)
    res = solve(prob, settings)
    p = res.p; B = info.B
    obj = 0.5 * rho * dot(p, p) + dot(prob.c, p)
    @printf("status %s  iters %d   ‖δμ−g‖ %.2e   obj %.6f   nnz %d/%d\n",
            res.status, res.niter, norm(B * p - prob.g), obj,
            count(>(1e-6), p), length(p))
    return res, info
end

"""
    integer_map(inst) :: Float64   — exact discrete MAP by brute force.
"""
function integer_map(inst::MAPInstance)
    V, edges, card = inst.V, inst.edges, inst.card
    @assert prod(BigInt.(card)) <= 10_000_000 "state space too large to brute-force"
    best = Inf; x = ones(Int, V)
    energy() = sum(inst.node_pot[i][x[i]] for i in 1:V) +
               sum(inst.edge_pot[ei][lin(x[i], x[j], card[j])] for (ei, (i, j)) in enumerate(edges))
    while true
        best = min(best, energy())
        k = 1
        while k <= V
            x[k] += 1; x[k] <= card[k] && break; x[k] = 1; k += 1
        end
        k > V && break
    end
    return best
end

"""
    verify_mrf(inst; rho, settings)

Solve the MRF as a sheaf section (LP if rho=0, SparseMAP if rho>0) and
print a corridor-§8-style report exposing the local-vs-global gap.
"""
function verify_mrf(inst::MAPInstance; rho::Float64 = 0.0, settings = default_settings())
    spec, meta = mrf_spec(inst)
    intopt = integer_map(inst)

    # The LP (ρ=0) carries the CERTIFIED 𝓛⊋𝓜 witness: the LP value is a
    # valid lower bound on the integer MAP, so a positive gap proves the LP
    # vertex is no global section.  (SparseMAP's linear part is NOT a bound.)
    probLP, infoLP = build_ipm_problem(spec; rho = 0.0)
    resLP = solve(probLP, settings)
    linLP = dot(probLP.c, resLP.p)

    println(rho == 0.0 ? "MAP-LP  (linear objective)" :
                         "SparseMAP  (native quadratic Q = ρI, ρ = $(rho))")
    println("="^60)
    @printf("LP status %s (%d it)   ‖δμ−g‖ %.1e   (agreement holds)\n",
            resLP.status, resLP.niter, norm(infoLP.B * resLP.p - probLP.g))
    @printf("variables / equalities  δμ=g          : %d / %d\n", length(resLP.p), size(infoLP.B, 1))
    @printf("LP ⟨θ,μ⟩  vs  integer MAP             : %.6f  vs  %.6f\n", linLP, intopt)
    @printf("LP integrality gap                    : %.6f  %s\n", intopt - linLP,
            intopt - linLP > 1e-6 ? "⇒ 𝓛 ⊋ 𝓜 (LP vertex has no global section)" : "(tight: 𝓛 = 𝓜)")

    res, info = resLP, infoLP
    if rho != 0.0                                   # report the SparseMAP POINT
        prob, info = build_ipm_problem(spec; rho = rho)
        res = solve(prob, settings)
        p = res.p
        @printf("SparseMAP status %s (%d it)  obj %.6f  nnz %d/%d  point %s\n",
                res.status, res.niter, 0.5 * rho * dot(p, p) + dot(prob.c, p),
                count(>(1e-6), p), length(p),
                norm(p - resLP.p) < 1e-5 ? "= LP vertex (sparse)" : "interior of 𝓛 (dense)")
    end

    vals = cellvals(res, info)                       # diagnostic: agreement masses
    for (ei, (i, j)) in enumerate(meta.edges)
        ci, cj = meta.card[i], meta.card[j]
        if ci == cj
            eb = vals[meta.edgev[ei]]
            @printf("edge (%d,%d) agreement mass Σₐμᵢⱼ(a,a)  : %.4f\n", i, j,
                    sum(eb[lin(a, a, cj)] for a in 1:ci))
        end
    end
    return res, info
end

# ---------------------------------------------------------------------
# The graded family (same assembly, richer cone):
#
#   cell object            cone (IPM)             overlap coboundary (kernel)
#   ---------------------  ---------------------  --------------------------
#   local marginal μ_R     PositiveCone           marginalise / channel  (here)
#   pseudo-moment  M_R     SemidefiniteCone       shared-block consistency
#   reduced density ρ_R    SemidefiniteCone       partial-trace consistency
#   entropic MAP           ExponentialCone        marginalise + −H(μ) cost
#
# Only the per-cell cone (and the "shared readout" kernel) changes; the
# coboundary-agreement backbone is invariant, and in each the difficulty
# is the POSITIVE-global-section obstruction, not the linear one.
# ---------------------------------------------------------------------

# =====================================================================
# SECTION 6.  Benchmark: SparseMAP on various graph topologies
# =====================================================================

# ---- Graph generators ------------------------------------------------

function grid_mrf(n::Int; card::Int = 2, penalty::Float64 = 1.0)
    V = n * n
    edges = Tuple{Int,Int}[]
    for i in 1:n, j in 1:n
        v = (i - 1) * n + j
        j < n && push!(edges, (v, v + 1))
        i < n && push!(edges, (v, v + n))
    end
    return make_mrf(V, edges, card, penalty)
end

function cycle_mrf(n::Int; card::Int = 2, penalty::Float64 = 1.0)
    edges = [(i, mod1(i + 1, n)) for i in 1:n]
    return make_mrf(n, edges, card, penalty)
end

function complete_mrf(n::Int; card::Int = 2, penalty::Float64 = 1.0)
    edges = [(i, j) for i in 1:n for j in (i + 1):n]
    return make_mrf(n, edges, card, penalty)
end

function ladder_mrf(n::Int; card::Int = 2, penalty::Float64 = 1.0)
    V = 2 * n
    edges = Tuple{Int,Int}[]
    for i in 1:n
        push!(edges, (2i - 1, 2i))                         # rung
        i < n && push!(edges, (2i - 1, 2i + 1))            # top rail
        i < n && push!(edges, (2i, 2i + 2))                # bottom rail
    end
    return make_mrf(V, edges, card, penalty)
end

function triangular_mrf(n::Int; card::Int = 2, penalty::Float64 = 1.0)
    V = n * n
    edges = Tuple{Int,Int}[]
    for i in 1:n, j in 1:n
        v = (i - 1) * n + j
        j < n && push!(edges, (v, v + 1))                  # horizontal
        i < n && push!(edges, (v, v + n))                  # vertical
        i < n && j < n && push!(edges, (v, v + n + 1))     # diagonal
    end
    return make_mrf(V, edges, card, penalty)
end

function make_mrf(V::Int, edges::Vector{Tuple{Int,Int}}, card::Int, penalty::Float64)
    cards = fill(card, V)
    ep = zeros(card * card)
    for a in 1:card
        ep[lin(a, a, card)] = penalty
    end
    return MAPInstance(V, edges, cards, [zeros(card) for _ in 1:V], [copy(ep) for _ in edges])
end

# ---- Benchmark -------------------------------------------------------

function run_benchmark(; optimizer = nothing, dual_optimizer = nothing, solver_name::String = "Mosek", rho::Float64 = 1.0, nwarmup::Int = 2, nruns::Int = 5)
    optimizer === nothing && error("Pass optimizer, e.g. run_benchmark(optimizer=Mosek.Optimizer)")

    cases = [
        ("grid 50×50",       grid_mrf(50),       1e6),
        ("cycle",            cycle_mrf(4000),    1e6),
        ("complete",         complete_mrf(95),   1e6),
        ("ladder",           ladder_mrf(1200),   1e6),
        ("triangular 40×40", triangular_mrf(40), 1e6),
    ]

    sname = rpad(solver_name, 9)
    sname_d = rpad(solver_name * "-D", 9)
    println("="^80)
    println("SparseMAP Benchmark: Sheaf IPM vs $solver_name  (rho = $rho)")
    println("="^80)
    println()
    if dual_optimizer !== nothing
        @printf("%-18s %6s %7s %7s %9s %9s %9s %7s %7s\n",
                "Graph", "raug", "|V|", "|E|", "IPM(ms)", sname, sname_d, "P/IPM", "D/IPM")
    else
        @printf("%-18s %6s %7s %7s %9s %9s %8s\n",
                "Graph", "raug", "|V|", "|E|", "IPM(ms)", sname, "speedup")
    end
    println("-"^80)

    for (name, inst, raug) in cases
        spec, _ = mrf_spec(inst)
        prob, info = build_ipm_problem(spec; rho = rho)
        nV, nE = inst.V, length(inst.edges)

        settings = IPMSettings{Float64}(
            kkt = UzawaSettings{Float64}(raug = raug),
            feas_tol = 1e-8, gap_tol = 1e-8, itmax = 200)
        for _ in 1:nwarmup; solve(prob, settings); end
        t_ipm = minimum([@elapsed solve(prob, settings) for _ in 1:nruns])

        for _ in 1:nwarmup
            m, _ = build_jump_model(spec, optimizer; rho = rho)
            optimize!(m)
        end
        t_mosek = minimum([begin
            m, _ = build_jump_model(spec, optimizer; rho = rho)
            @elapsed optimize!(m)
        end for _ in 1:nruns])

        if dual_optimizer !== nothing
            for _ in 1:nwarmup
                m, _ = build_jump_model(spec, dual_optimizer; rho = rho)
                optimize!(m)
            end
            t_dual = minimum([begin
                m, _ = build_jump_model(spec, dual_optimizer; rho = rho)
                @elapsed optimize!(m)
            end for _ in 1:nruns])

            @printf("%-18s %6.0e %7d %7d %9.1f %9.1f %9.1f %6.2fx %6.2fx\n",
                    name, raug, nV, nE, t_ipm * 1000, t_mosek * 1000, t_dual * 1000,
                    t_mosek / t_ipm, t_dual / t_ipm)
        else
            @printf("%-18s %6.0e %7d %7d %9.1f %9.1f %7.2fx\n",
                    name, raug, nV, nE, t_ipm * 1000, t_mosek * 1000,
                    t_mosek / t_ipm)
        end
    end
end

# ---- Dense stochastic channel graphs ---------------------------------

"""
    random_doubly_stochastic(n, rng) :: Matrix  (n × n)

A random doubly-stochastic matrix via Sinkhorn iteration. The uniform
distribution is a fixed point: P · (1/n,...,1/n)' = (1/n,...,1/n)', so
agreement constraints are always feasible (uniform section glues).
"""
function random_doubly_stochastic(n::Integer, rng::AbstractRNG; iters::Int = 20)
    M = rand(rng, n, n) .+ 0.1
    for _ in 1:iters
        M ./= sum(M, dims = 1)  # normalize columns
        M ./= sum(M, dims = 2)  # normalize rows
    end
    M ./= sum(M, dims = 1)      # final column normalization
    return M
end

function dense_channel_spec(edges::Vector{Tuple{Int,Int}}, V::Int;
                            nstates::Int = 8, seed::Int = 0)
    rng = MersenneTwister(seed)
    chan = [random_doubly_stochastic(nstates, rng) for _ in 1:V]
    dim  = fill(nstates, V)
    cone = AbstractCone[PositiveCone() for _ in 1:V]
    cvec = [randn(rng, nstates) for _ in 1:V]  # random linear cost
    sedges = [StochEdge(u, chan[u], w, chan[w], zeros(nstates)) for (u, w) in edges]
    pins   = [Pin(v, ones(1, nstates), [1.0]) for v in 1:V]
    return SheafSpec(dim, cone, cvec, sedges, pins)
end

function grid_dense(n::Int; nstates::Int = 8, seed::Int = 0)
    V = n * n
    edges = Tuple{Int,Int}[]
    for i in 1:n, j in 1:n
        v = (i - 1) * n + j
        j < n && push!(edges, (v, v + 1))
        i < n && push!(edges, (v, v + n))
    end
    return dense_channel_spec(edges, V; nstates, seed)
end

function cycle_dense(n::Int; nstates::Int = 8, seed::Int = 0)
    edges = [(i, mod1(i + 1, n)) for i in 1:n]
    return dense_channel_spec(edges, n; nstates, seed)
end

function complete_dense(n::Int; nstates::Int = 8, seed::Int = 0)
    edges = [(i, j) for i in 1:n for j in (i + 1):n]
    return dense_channel_spec(edges, n; nstates, seed)
end

function ladder_dense(n::Int; nstates::Int = 8, seed::Int = 0)
    V = 2 * n
    edges = Tuple{Int,Int}[]
    for i in 1:n
        push!(edges, (2i - 1, 2i))
        i < n && push!(edges, (2i - 1, 2i + 1))
        i < n && push!(edges, (2i, 2i + 2))
    end
    return dense_channel_spec(edges, V; nstates, seed)
end

function triangular_dense(n::Int; nstates::Int = 8, seed::Int = 0)
    V = n * n
    edges = Tuple{Int,Int}[]
    for i in 1:n, j in 1:n
        v = (i - 1) * n + j
        j < n && push!(edges, (v, v + 1))
        i < n && push!(edges, (v, v + n))
        i < n && j < n && push!(edges, (v, v + n + 1))
    end
    return dense_channel_spec(edges, V; nstates, seed)
end

function run_benchmark_dense(; optimizer = nothing, dual_optimizer = nothing, solver_name::String = "Mosek", rho::Float64 = 1.0, nwarmup::Int = 2, nruns::Int = 5)
    optimizer === nothing && error("Pass optimizer, e.g. run_benchmark_dense(optimizer=Mosek.Optimizer)")

    cases = [
        ("grid 20×20",       grid_dense(20),       20*20, 2*20*(20-1), 1e6),
        ("cycle",            cycle_dense(400),     400,   400,         1e6),
        ("complete",         complete_dense(40),   40,    40*39÷2,     1e6),
        ("ladder",           ladder_dense(200),    400,   3*200-2,     1e6),
        ("triangular 15×15", triangular_dense(15), 15*15, 3*15*(15-1), 1e6),
    ]

    sname = rpad(solver_name, 9)
    sname_d = rpad(solver_name * "-D", 9)
    println("="^80)
    println("Dense Stochastic Channels: Sheaf IPM vs $solver_name  (rho = $rho)")
    println("="^80)
    println()
    if dual_optimizer !== nothing
        @printf("%-18s %6s %7s %7s %9s %9s %9s %7s %7s\n",
                "Graph", "raug", "|V|", "|E|", "IPM(ms)", sname, sname_d, "P/IPM", "D/IPM")
    else
        @printf("%-18s %6s %7s %7s %9s %9s %8s\n",
                "Graph", "raug", "|V|", "|E|", "IPM(ms)", sname, "speedup")
    end
    println("-"^80)

    for (name, spec, nV, nE, raug) in cases
        prob, info = build_ipm_problem(spec; rho = rho)

        settings = IPMSettings{Float64}(
            kkt = UzawaSettings{Float64}(raug = raug),
            feas_tol = 1e-8, gap_tol = 1e-8, itmax = 200)
        for _ in 1:nwarmup; solve(prob, settings); end
        t_ipm = minimum([@elapsed solve(prob, settings) for _ in 1:nruns])

        for _ in 1:nwarmup
            m, _ = build_jump_model(spec, optimizer; rho = rho)
            optimize!(m)
        end
        t_mosek = minimum([begin
            m, _ = build_jump_model(spec, optimizer; rho = rho)
            @elapsed optimize!(m)
        end for _ in 1:nruns])

        if dual_optimizer !== nothing
            for _ in 1:nwarmup
                m, _ = build_jump_model(spec, dual_optimizer; rho = rho)
                optimize!(m)
            end
            t_dual = minimum([begin
                m, _ = build_jump_model(spec, dual_optimizer; rho = rho)
                @elapsed optimize!(m)
            end for _ in 1:nruns])

            @printf("%-18s %6.0e %7d %7d %9.1f %9.1f %9.1f %6.2fx %6.2fx\n",
                    name, raug, nV, nE, t_ipm * 1000, t_mosek * 1000, t_dual * 1000,
                    t_mosek / t_ipm, t_dual / t_ipm)
        else
            @printf("%-18s %6.0e %7d %7d %9.1f %9.1f %7.2fx\n",
                    name, raug, nV, nE, t_ipm * 1000, t_mosek * 1000,
                    t_mosek / t_ipm)
        end
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    include(joinpath(@__DIR__, "..", "benchmark_utils.jl"))
    opts = parse_benchmark_args(ARGS)

    if opts.mosek
        using MosekTools
    else
        using Clarabel
    end
    optimizer, dual_optimizer = get_optimizers(opts)
    solver_name = opts.mosek ? "Mosek" : "Clarabel"
    print_benchmark_config(opts; lp_only = true)

    run_benchmark(; optimizer, dual_optimizer, solver_name,
                  nwarmup = opts.nwarmup, nruns = opts.nruns)
    println()
    run_benchmark_dense(; optimizer, dual_optimizer, solver_name,
                        nwarmup = opts.nwarmup, nruns = opts.nruns)
end
