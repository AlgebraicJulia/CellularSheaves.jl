######################################################################
# map_sdp.jl
#
# The PSD (MOMENT) TIGHTENING of the MAP sheaf, as a conic program in
# CellularSheaves.IPM — the SemidefiniteCone sibling of map_lp.jl.
#
#       min   ⟨C, M⟩ (+ ρ/2‖M‖_F²)        M = (M_R)_R  clique moment matrices
#       s.t.  δM = g                       shared-block agreement + pins
#             M_R ∈ SemidefiniteCone       each region's moment matrix ⪰ 0
#
# Same assembly as the LP note; only the cone and the meaning of the
# restriction change (see local_consistency_sheaf.md §4.1–4.2):
#
#   LP  (map_lp.jl)                 SDP moment  (this file)
#   -----------------------------   -------------------------------------
#   region belief  μ_R (simplex)    region moment matrix M_R (PSD)
#   orthant cone                    SemidefiniteCone
#   normalisation Σμ_R = 1          M_R[1,1] = 1   (constant monomial)
#   —                               localising  M_R[x_i,x_i] = 1  (x_i²=1)
#   marginalise (restriction)       principal-submatrix / compression
#   overlap = agree on marginals    overlap = agree on the SHARED BLOCK
#   global joint distribution       global PSD moment matrix (completion)
#   tight iff junction tree         completes iff CHORDAL (Grone–Johnson–…)
#
# ±1 variables (state 1 ↦ +1, state 2 ↦ −1).  A region over variables
# R = {i,j,…} carries a moment matrix over the monomials (1, x_i, x_j, …)
# of size |R|+1; coordinate 1 is the constant.  Matrices are stored as the
# solver's scaled svec (column-major lower triangle, off-diagonals ×√2,
# so ⟨svec A, svec B⟩ = ⟨A,B⟩ = tr(AB)).
#
# TWO ORTHOGONAL AXES (both verified in map_sdp_oracle.py):
#   • VALUE (region/clique size).  LP ≤ pairwise-edge sheaf ≤ single-clique
#     Shor ≤ integer.  Frustrated triangle: 0 ≤ 0 ≤ 3/4 ≤ 1.  The three
#     edge blocks overlap only on node means, so they DON'T couple the
#     correlations — pairwise sheaf = LP; the {i,j,k} clique block does.
#   • COMPLETABILITY (region-graph topology).  A positive global section
#     exists for every consistent local spec iff the pattern is chordal.
#     Triangle (complete ⇒ chordal): always completes.  Chordless 4-cycle:
#     correlations (.8,.8,.8,−.8) are pairwise-valid yet DON'T complete —
#     positive local sections agreeing on overlaps, no positive global
#     section, sitting exactly where the LP was tight.
#
# NOTE: written against the CellularSheaves.IPM PR-67 API (ships
# SemidefiniteCone, svec!/smat! with the √2 convention).  Not executed
# here (no Julia toolchain); map_sdp_oracle.py is the numerical ground
# truth for every value quoted above.
######################################################################

using CellularSheaves.IPM
using CellularSheaves.BlockSparseArrays: colrange, rowrange, blocksparse, block, nvtxs
using LinearAlgebra
using Printf

# reuse MAPInstance and frustrated_triangle from map_lp.jl
include("map_lp.jl")

# =====================================================================
# SECTION 1.  Scaled svec  (matches CellularSheaves.IPM/cone/sdp.jl)
# =====================================================================

const RT2 = sqrt(2.0)

"""Lower-triangular svec index: dict (i,j)↦position, i≥j, column-major."""
function svec_index(k::Integer)
    pos = Dict{Tuple{Int,Int},Int}(); c = 0
    for j in 1:k
        c += 1; pos[(j, j)] = c
        for i in j+1:k
            c += 1; pos[(i, j)] = c
        end
    end
    return pos                      # length k(k+1)÷2
end

svec_len(k::Integer) = k * (k + 1) ÷ 2
loij(a, b) = a ≥ b ? (a, b) : (b, a)

"""svec of a symmetric matrix M (√2 on off-diagonals): ⟨svec A,svec B⟩=tr(AB)."""
function svec(M::AbstractMatrix)
    k = size(M, 1); v = zeros(svec_len(k)); pos = svec_index(k)
    for a in 1:k, b in 1:a
        v[pos[(a, b)]] = a == b ? M[a, a] : RT2 * M[a, b]
    end
    return v
end

"""smat: invert svec back to a symmetric matrix (for reading solutions)."""
function smat(v::AbstractVector, k::Integer)
    M = zeros(k, k); pos = svec_index(k)
    for a in 1:k, b in 1:a
        x = v[pos[(a, b)]]
        M[a, b] = M[b, a] = a == b ? x : x / RT2
    end
    return M
end

# =====================================================================
# SECTION 2.  A moment region and its selections
# =====================================================================

"""
A region over MRF variables `vars` (sorted). Monomial coordinates:
1 ↦ constant, 1+t ↦ x_{vars[t]}.  Matrix size k = |vars|+1.
"""
struct MomentRegion
    vars::Vector{Int}
    k::Int
    MomentRegion(vars) = new(sort(vars), length(vars) + 1)
end

"""Monomial coordinate (1-based, 1 = constant) of MRF variable `i` in region."""
varcoord(R::MomentRegion, i::Integer) = 1 + findfirst(==(i), R.vars)

"""
    shared_selection(R, S) :: Matrix   (svec_len(|S|+1) × svec_len(R.k))

0/1 map reading the moment block on the shared variable set `S` (⊆ R.vars)
out of R's svec.  The constant coordinate is always shared; the √2 scaling
matches on both sides, so a selection suffices (no rescaling).
"""
function shared_selection(R::MomentRegion, S::Vector{Int})
    Ssub = MomentRegion(S)                         # coords 1..|S|+1
    posR = svec_index(R.k); posS = svec_index(Ssub.k)
    E = zeros(svec_len(Ssub.k), svec_len(R.k))
    # map each shared coordinate (incl. constant) to its coordinate in R
    coordR = Dict(1 => 1); coordS = Dict(1 => 1)
    for (t, v) in enumerate(Ssub.vars)
        coordS[1 + t] = 1 + t
        coordR[1 + t] = varcoord(R, v)
    end
    for a in 1:Ssub.k, b in 1:a
        E[posS[(a, b)], posR[loij(coordR[a], coordR[b])]] = 1.0
    end
    return E
end

"""Selection of the diagonal svec positions of a size-k region (for the pins)."""
function diag_selection(k::Integer)
    pos = svec_index(k); D = zeros(k, svec_len(k))
    for a in 1:k
        D[a, pos[(a, a)]] = 1.0
    end
    return D                                        # rows: M[a,a]
end

# =====================================================================
# SECTION 3.  ±1 potentials → moment cost matrix
# =====================================================================

spin(a::Integer) = a == 1 ? 1.0 : -1.0             # state 1 ↦ +1, state 2 ↦ −1

"""Symmetric cost matrix C on region R so that ⟨C,M_R⟩ = Σ (assigned potentials),
including the affine constant, which is folded onto C[1,1] (M[1,1]=1 is pinned)."""
function region_cost(R::MomentRegion, inst::MAPInstance,
                     my_nodes::Vector{Int}, my_edges::Vector{Int})
    C = zeros(R.k, R.k)
    for i in my_nodes                              # θ_i = α + β x_i
        θ = inst.node_pot[i]
        α = (θ[1] + θ[2]) / 2
        β = (θ[1] - θ[2]) / 2                       # state 1 ↦ +1, state 2 ↦ −1
        ci = varcoord(R, i)
        C[1, 1] += α
        C[1, ci] += β / 2; C[ci, 1] += β / 2
    end
    for ei in my_edges                             # θ_ij = a + b x_i + c x_j + d x_i x_j
        (i, j) = inst.edges[ei]; θ = inst.edge_pot[ei]   # θ[(a-1)*2 + b], states 1/2
        f(a, b) = θ[(a - 1) * 2 + b]
        a_ = (f(1,1) + f(1,2) + f(2,1) + f(2,2)) / 4
        b_ = (f(1,1) + f(1,2) - f(2,1) - f(2,2)) / 4
        c_ = (f(1,1) - f(1,2) + f(2,1) - f(2,2)) / 4
        d_ = (f(1,1) - f(1,2) - f(2,1) + f(2,2)) / 4
        ci, cj = varcoord(R, i), varcoord(R, j)
        C[1, 1] += a_
        C[1, ci] += b_ / 2; C[ci, 1] += b_ / 2
        C[1, cj] += c_ / 2; C[cj, 1] += c_ / 2
        C[ci, cj] += d_ / 2; C[cj, ci] += d_ / 2
    end
    return C
end

# =====================================================================
# SECTION 4.  Build the moment sheaf spec
# =====================================================================

struct SDPEdge          # shared-block agreement  E_u svec(M_u) − E_w svec(M_w) = 0
    u::Int; Eu::Matrix{Float64}
    w::Int; Ew::Matrix{Float64}
end
struct SDPPin           # A svec(M_v) = rhs   (here: diagonal = 1)
    v::Int; A::Matrix{Float64}; rhs::Vector{Float64}
end
struct MomentSpec
    regions::Vector{MomentRegion}
    c::Vector{Vector{Float64}}        # svec(C_R) per region
    edges::Vector{SDPEdge}
    pins::Vector{SDPPin}
end

"""
    moment_spec(inst; cliques) :: (MomentSpec, info)

Build the moment sheaf over the given `cliques` (each a vector of MRF
variables).  Overlapping cliques get a shared-block agreement edge; each
region gets its diagonal pinned to 1 (normalisation + x_i²=1).  Each node
and edge potential is assigned to one region containing it.

`cliques` defaults to the graph's EDGES (the pairwise sheaf).  Pass the
maximal cliques (e.g. `[[1,2,3]]` for a triangle) for the Shor tightening.
"""
function moment_spec(inst::MAPInstance; cliques::Vector{Vector{Int}} =
                     [collect(e) for e in inst.edges])
    regions = [MomentRegion(cl) for cl in cliques]
    nR = length(regions)

    # assign each node / edge potential to the first region containing it
    node_of = Dict{Int,Int}(); edge_of = Dict{Int,Int}()
    for i in 1:inst.V, r in 1:nR
        i in regions[r].vars && (get!(node_of, i, r); )
    end
    for (ei, (i, j)) in enumerate(inst.edges), r in 1:nR
        (i in regions[r].vars && j in regions[r].vars) && (get!(edge_of, ei, r); )
    end
    @assert length(edge_of) == length(inst.edges) "some edge has no clique containing both endpoints"

    cvec = Vector{Float64}[]
    for r in 1:nR
        mine_n = [i for (i, rr) in node_of if rr == r]
        mine_e = [ei for (ei, rr) in edge_of if rr == r]
        push!(cvec, svec(region_cost(regions[r], inst, mine_n, mine_e)))
    end

    # shared-block agreement on every overlapping pair of regions
    edges = SDPEdge[]
    for u in 1:nR, w in u+1:nR
        S = sort(collect(intersect(Set(regions[u].vars), Set(regions[w].vars))))
        isempty(S) && continue
        push!(edges, SDPEdge(u, shared_selection(regions[u], S),
                             w, shared_selection(regions[w], S)))
    end

    # pins: diagonal of every region equals 1
    pins = [SDPPin(r, diag_selection(regions[r].k), ones(regions[r].k)) for r in 1:nR]

    return MomentSpec(regions, cvec, edges, pins),
           (regions = regions, node_of = node_of, edge_of = edge_of)
end

# =====================================================================
# SECTION 5.  Assemble the IPMProblem  (SemidefiniteCone cells)
# =====================================================================

"""
    build_sdp_problem(spec; rho) :: (IPMProblem, info)

Every cell is a SemidefiniteCone whose column block is the region's svec.
`rho > 0` adds the native Frobenius quadratic ρ/2‖M‖_F² (= ρ/2‖svec‖²).
"""
function build_sdp_problem(spec::MomentSpec; rho::Float64 = 0.0)
    nR = length(spec.regions)
    dims = [svec_len(R.k) for R in spec.regions]
    row_ids = Int[]; col_ids = Int[]; blocks = Matrix{Float64}[]
    rhs = Dict{Int,Vector{Float64}}()
    ecnt = Ref(0); newedge!() = (ecnt[] += 1; ecnt[])
    place!(e, v, M) = (push!(row_ids, e); push!(col_ids, v); push!(blocks, Matrix{Float64}(M)))

    for se in spec.edges                       # E_u svec_u − E_w svec_w = 0
        e = newedge!(); place!(e, se.u, se.Eu); place!(e, se.w, -se.Ew)
    end
    for pn in spec.pins                        # diag = 1
        e = newedge!(); place!(e, pn.v, pn.A); rhs[e] = copy(pn.rhs)
    end

    B = blocksparse(row_ids, col_ids, blocks)
    @assert nvtxs(B) == nR
    g = zeros(size(B, 1)); for (e, v) in rhs; g[rowrange(B, e)] .= v; end

    Q = IPM.allocblockdiag(B); fill!(Q, 0)
    if rho != 0.0
        for v in 1:nR
            Qv = block(Q, v, v, v)
            Qv .= rho .* Matrix{Float64}(I, dims[v], dims[v])
        end
    end

    c = zeros(size(B, 2)); for v in 1:nR; c[colrange(B, v)] .= spec.c[v]; end
    cones = AbstractCone[SemidefiniteCone() for _ in 1:nR]
    return IPMProblem(Q, B, c, g, cones), (B = B, regions = spec.regions)
end

"""Read region r's moment matrix back from a solution."""
read_moment(res, info, r::Integer) =
    smat(res.p[colrange(info.B, r)], info.regions[r].k)

# =====================================================================
# SECTION 6.  Instances
# =====================================================================

"""Triangle as ONE clique {1,2,3}: a single 4×4 moment matrix = Shor SDP (→ 3/4)."""
triangle_clique(; penalty = 1.0) =
    moment_spec(frustrated_triangle(; penalty = penalty); cliques = [[1, 2, 3]])

"""Triangle as three EDGE cliques: pairwise sheaf, node-only overlaps (→ LP value 0)."""
triangle_pairwise(; penalty = 1.0) =
    moment_spec(frustrated_triangle(; penalty = penalty))

"""
Even 4-cycle as four edge cliques {(1,2),(2,3),(3,4),(4,1)}: genuine overlaps
(single nodes), value-tight but the site of the chordless completion obstruction.
"""
function cycle4(; penalty = 1.0)
    edges = [(1, 2), (2, 3), (3, 4), (1, 4)]; card = [2, 2, 2, 2]
    ep = Float64[penalty, 0.0, 0.0, penalty]
    inst = MAPInstance(4, edges, card, [zeros(2) for _ in 1:4], [copy(ep) for _ in edges])
    return moment_spec(inst)
end

# =====================================================================
# SECTION 7.  PSD-completion witness (feasibility)
# =====================================================================

"""
    completion_spec(edge_corr) :: MomentSpec

Single 4×4 region over (1,x_1,x_2,x_3,x_4)... actually over {1,2,3,4}: a
5×5 moment matrix with unit diagonal, the four C4-edge correlations PINNED
to `edge_corr` (Dict (i,j)=>value), chords free.  `solve` returns OPTIMAL
if the local spec completes to a global PSD moment matrix, INFEASIBLE if
not — the chordless-C4 obstruction.  With (.8,.8,.8,−.8): INFEASIBLE.
"""
function completion_spec(edge_corr::Dict{Tuple{Int,Int},Float64})
    R = MomentRegion([1, 2, 3, 4])                 # 5×5, coords 1=const, 1+i=x_i
    pos = svec_index(R.k)
    rows = Vector{Float64}[]; rhsv = Float64[]
    # unit diagonal
    for a in 1:R.k
        e = zeros(svec_len(R.k)); e[pos[(a, a)]] = 1.0; push!(rows, e); push!(rhsv, 1.0)
    end
    # pinned edge correlations M[x_i,x_j] = v  (off-diagonal ⇒ svec entry ×√2)
    for ((i, j), v) in edge_corr
        ci, cj = varcoord(R, i), varcoord(R, j)
        e = zeros(svec_len(R.k)); e[pos[loij(ci, cj)]] = 1.0
        push!(rows, e); push!(rhsv, RT2 * v)
    end
    A = permutedims(hcat(rows...))
    spec = MomentSpec([R], [zeros(svec_len(R.k))],
                      SDPEdge[], [SDPPin(1, A, rhsv)])
    return spec
end

# =====================================================================
# SECTION 8.  Solve / report
# =====================================================================

default_settings() = IPMSettings{Float64}(feas_tol = 1e-9, gap_tol = 1e-9, itmax = 200)

function solve_sdp(spec::MomentSpec; rho::Float64 = 0.0, settings = default_settings())
    prob, info = build_sdp_problem(spec; rho = rho)
    res = solve(prob, settings)
    p = res.p
    obj = 0.5 * rho * dot(p, p) + dot(prob.c, p)
    @printf("status %s  iters %d   ‖δM−g‖ %.2e   ⟨C,M⟩ %.6f\n",
            res.status, res.niter, norm(info.B * p - prob.g), dot(prob.c, p))
    # min-eigenvalue of each region block, as a positivity read-out
    for r in 1:length(info.regions)
        M = read_moment(res, info, r)
        @printf("  region %d  (%d×%d)  λ_min = %+.4f\n", r, info.regions[r].k,
                info.regions[r].k, minimum(eigvals(Symmetric(M))))
    end
    return res, info
end

# =====================================================================
# SECTION 9.  JuMP reference model (for Mosek comparison)
# =====================================================================

using JuMP

function build_jump_sdp(spec::MomentSpec, optimizer; rho::Float64 = 0.0)
    nR = length(spec.regions)
    model = Model(optimizer); set_silent(model)

    # PSD matrix variables for each region
    Ms = [@variable(model, [1:spec.regions[r].k, 1:spec.regions[r].k], PSD) for r in 1:nR]

    # Objective: Σ ⟨C_r, M_r⟩ + ρ/2 ‖M‖_F²
    obj = sum(dot(smat(spec.c[r], spec.regions[r].k), Ms[r]) for r in 1:nR)
    if rho != 0.0
        obj += (rho / 2) * sum(dot(Ms[r], Ms[r]) for r in 1:nR)
    end
    @objective(model, Min, obj)

    # Shared-block agreement: E_u svec(M_u) = E_w svec(M_w)
    for se in spec.edges
        Ru, Rw = spec.regions[se.u], spec.regions[se.w]
        # Extract shared variables
        Su = Set(Ru.vars); Sw = Set(Rw.vars)
        shared = sort(collect(intersect(Su, Sw)))
        Ssub = MomentRegion(shared)
        # Constrain shared block equality
        for a in 1:Ssub.k, b in 1:a
            ca_u = a == 1 ? 1 : varcoord(Ru, shared[a-1])
            cb_u = b == 1 ? 1 : varcoord(Ru, shared[b-1])
            ca_w = a == 1 ? 1 : varcoord(Rw, shared[a-1])
            cb_w = b == 1 ? 1 : varcoord(Rw, shared[b-1])
            @constraint(model, Ms[se.u][ca_u, cb_u] == Ms[se.w][ca_w, cb_w])
        end
    end

    # Pins: diagonal = 1
    for pn in spec.pins
        k = spec.regions[pn.v].k
        for a in 1:k
            @constraint(model, Ms[pn.v][a, a] == 1.0)
        end
    end

    return model, Ms
end

# =====================================================================
# SECTION 10.  Graph generators for benchmarks
# =====================================================================

"""Pairwise MRF instance on an n×n grid (anti-ferromagnetic)."""
function grid_mrf_sdp(n::Int; penalty::Float64 = 1.0)
    V = n * n
    edges = Tuple{Int,Int}[]
    for i in 1:n, j in 1:n
        v = (i - 1) * n + j
        j < n && push!(edges, (v, v + 1))
        i < n && push!(edges, (v, v + n))
    end
    card = fill(2, V)
    ep = Float64[penalty, 0.0, 0.0, penalty]
    return MAPInstance(V, edges, card, [zeros(2) for _ in 1:V], [copy(ep) for _ in edges])
end

"""Pairwise MRF instance on a cycle of length n."""
function cycle_mrf_sdp(n::Int; penalty::Float64 = 1.0)
    edges = [(i, mod1(i + 1, n)) for i in 1:n]
    card = fill(2, n)
    ep = Float64[penalty, 0.0, 0.0, penalty]
    return MAPInstance(n, edges, card, [zeros(2) for _ in 1:n], [copy(ep) for _ in edges])
end

"""Pairwise MRF instance on a ladder graph."""
function ladder_mrf_sdp(n::Int; penalty::Float64 = 1.0)
    V = 2 * n
    edges = Tuple{Int,Int}[]
    for i in 1:n
        push!(edges, (2i - 1, 2i))
        i < n && push!(edges, (2i - 1, 2i + 1))
        i < n && push!(edges, (2i, 2i + 2))
    end
    card = fill(2, V)
    ep = Float64[penalty, 0.0, 0.0, penalty]
    return MAPInstance(V, edges, card, [zeros(2) for _ in 1:V], [copy(ep) for _ in edges])
end

"""Pairwise MRF instance on a complete graph."""
function complete_mrf_sdp(n::Int; penalty::Float64 = 1.0)
    edges = [(i, j) for i in 1:n for j in (i+1):n]
    card = fill(2, n)
    ep = Float64[penalty, 0.0, 0.0, penalty]
    return MAPInstance(n, edges, card, [zeros(2) for _ in 1:n], [copy(ep) for _ in edges])
end

# =====================================================================
# SECTION 11.  Unital CPTP channels (bistochastic = doubly-stochastic lift)
# =====================================================================

"""
    random_unitary(k, rng) :: Matrix

Random k×k unitary via QR decomposition of a random matrix.
"""
function random_unitary(k::Integer, rng::AbstractRNG)
    A = randn(rng, k, k)
    Q, R = qr(A)
    # Ensure determinant +1 (proper rotation)
    Q = Matrix(Q)
    Q[:, 1] .*= sign(R[1, 1])
    return Q
end

"""
    unitary_channel_svec(U) :: Matrix

The svec-form matrix T such that svec(U M U') = T · svec(M).
Uses the √2 scaling convention for off-diagonals.
"""
function unitary_channel_svec(U::AbstractMatrix)
    k = size(U, 1)
    pos = svec_index(k)
    n = svec_len(k)
    T = zeros(n, n)
    # Build by applying to basis matrices
    for (idx_in, (a, b)) in enumerate(sort(collect(keys(pos)), by = p -> pos[p]))
        # Basis matrix E_ab (symmetric, with √2 scaling accounted for)
        E = zeros(k, k)
        if a == b
            E[a, a] = 1.0
        else
            E[a, b] = 1.0 / RT2
            E[b, a] = 1.0 / RT2
        end
        # Apply channel
        Mout = U * E * U'
        # Extract svec
        T[:, pos[(a, b)]] = svec(Mout)
    end
    return T
end

"""
    random_unitary_channel_svec(k, nkraus, rng) :: Matrix

Random unital channel Φ(M) = Σₖ pₖ Uₖ M Uₖ' as a svec-form matrix.
With nkraus=1: single unitary (sharp, basis-changing).
With nkraus large: approaches depolarizing (erases to I).
Unital: Φ(I) = I, so M = I is always feasible.
"""
function random_unitary_channel_svec(k::Integer, nkraus::Integer, rng::AbstractRNG)
    # Random probabilities
    p = rand(rng, nkraus)
    p ./= sum(p)
    # Sum of weighted unitary channels
    T = zeros(svec_len(k), svec_len(k))
    for i in 1:nkraus
        U = random_unitary(k, rng)
        T .+= p[i] .* unitary_channel_svec(U)
    end
    return T
end

"""
    is_unital(T, k) :: Bool

Check that T (svec-form channel) satisfies Φ(I) = I.
"""
function is_unital(T::AbstractMatrix, k::Integer)
    Isvec = svec(Matrix{Float64}(I, k, k))
    return norm(T * Isvec - Isvec) < 1e-9
end

# ---- Unital channel graph specs --------------------------------------

struct UnitalSDPEdge
    u::Int; Tu::Matrix{Float64}
    w::Int; Tw::Matrix{Float64}
end

struct UnitalMomentSpec
    k::Int                              # all cells same size k×k
    ncells::Int
    c::Vector{Vector{Float64}}          # svec cost per cell
    edges::Vector{UnitalSDPEdge}
    pins::Vector{SDPPin}
end

"""
    unital_channel_spec(edges, V; k, nkraus, seed) :: UnitalMomentSpec

Build a moment sheaf over V cells of size k×k, connected by random unital
channels. Each edge demands Φ_u(M_u) = Φ_w(M_w) in a shared k×k edge stalk.
Unital ⇒ M = I on all cells is always feasible (the quantum uniform).
"""
function unital_channel_spec(edges::Vector{Tuple{Int,Int}}, V::Int;
                             k::Int = 3, nkraus::Int = 3, seed::Int = 0)
    rng = MersenneTwister(seed)
    # Random unital channel per cell
    chans = [random_unitary_channel_svec(k, nkraus, rng) for _ in 1:V]
    # Verify unitality
    for (v, T) in enumerate(chans)
        @assert is_unital(T, k) "Channel $v is not unital"
    end
    # Random cost (pull away from I)
    cvec = [randn(rng, svec_len(k)) for _ in 1:V]
    # Agreement edges: T_u svec(M_u) - T_w svec(M_w) = 0
    sedges = [UnitalSDPEdge(u, chans[u], w, chans[w]) for (u, w) in edges]
    # Pins: diagonal = 1 (unit diagonal for ±1 moments)
    D = diag_selection(k)
    pins = [SDPPin(v, D, ones(k)) for v in 1:V]
    return UnitalMomentSpec(k, V, cvec, sedges, pins)
end

function build_unital_sdp_problem(spec::UnitalMomentSpec; rho::Float64 = 0.0)
    V = spec.ncells
    k = spec.k
    dim = svec_len(k)
    row_ids = Int[]; col_ids = Int[]; blocks = Matrix{Float64}[]
    rhs = Dict{Int,Vector{Float64}}()
    ecnt = Ref(0); newedge!() = (ecnt[] += 1; ecnt[])
    place!(e, v, M) = (push!(row_ids, e); push!(col_ids, v); push!(blocks, Matrix{Float64}(M)))

    # Agreement: T_u svec(M_u) - T_w svec(M_w) = 0
    for se in spec.edges
        e = newedge!(); place!(e, se.u, se.Tu); place!(e, se.w, -se.Tw)
    end
    # Pins: diagonal = 1
    for pn in spec.pins
        e = newedge!(); place!(e, pn.v, pn.A); rhs[e] = copy(pn.rhs)
    end

    B = blocksparse(row_ids, col_ids, blocks)
    @assert nvtxs(B) == V
    g = zeros(size(B, 1)); for (e, v) in rhs; g[rowrange(B, e)] .= v; end

    Q = IPM.allocblockdiag(B); fill!(Q, 0)
    if rho != 0.0
        for v in 1:V
            Qv = block(Q, v, v, v)
            Qv .= rho .* Matrix{Float64}(I, dim, dim)
        end
    end

    c = zeros(size(B, 2)); for v in 1:V; c[colrange(B, v)] .= spec.c[v]; end
    cones = AbstractCone[SemidefiniteCone() for _ in 1:V]
    return IPMProblem(Q, B, c, g, cones), (B = B, k = k)
end

# Graph generators for unital channels
function grid_unital(n::Int; k::Int = 3, nkraus::Int = 3, seed::Int = 0)
    V = n * n
    edges = Tuple{Int,Int}[]
    for i in 1:n, j in 1:n
        v = (i - 1) * n + j
        j < n && push!(edges, (v, v + 1))
        i < n && push!(edges, (v, v + n))
    end
    return unital_channel_spec(edges, V; k, nkraus, seed)
end

function cycle_unital(n::Int; k::Int = 3, nkraus::Int = 3, seed::Int = 0)
    edges = [(i, mod1(i + 1, n)) for i in 1:n]
    return unital_channel_spec(edges, n; k, nkraus, seed)
end

function ladder_unital(n::Int; k::Int = 3, nkraus::Int = 3, seed::Int = 0)
    V = 2 * n
    edges = Tuple{Int,Int}[]
    for i in 1:n
        push!(edges, (2i - 1, 2i))
        i < n && push!(edges, (2i - 1, 2i + 1))
        i < n && push!(edges, (2i, 2i + 2))
    end
    return unital_channel_spec(edges, V; k, nkraus, seed)
end

function complete_unital(n::Int; k::Int = 3, nkraus::Int = 3, seed::Int = 0)
    edges = [(i, j) for i in 1:n for j in (i+1):n]
    return unital_channel_spec(edges, n; k, nkraus, seed)
end

# JuMP builder for unital channel spec
function build_jump_unital_sdp(spec::UnitalMomentSpec, optimizer; rho::Float64 = 0.0)
    V = spec.ncells; k = spec.k
    model = Model(optimizer); set_silent(model)
    Ms = [@variable(model, [1:k, 1:k], PSD) for _ in 1:V]

    obj = sum(dot(smat(spec.c[v], k), Ms[v]) for v in 1:V)
    if rho != 0.0
        obj += (rho / 2) * sum(dot(Ms[v], Ms[v]) for v in 1:V)
    end
    @objective(model, Min, obj)

    # Agreement: T_u svec(M_u) = T_w svec(M_w)
    for se in spec.edges
        # Compute T_u * svec(M_u) and T_w * svec(M_w) symbolically
        # This is complex in JuMP, so we use explicit element constraints
        # For each output coordinate...
        for out_idx in 1:svec_len(k)
            lhs = sum(se.Tu[out_idx, in_idx] * (
                let (a, b) = collect(keys(svec_index(k)))[findfirst(==(in_idx), [svec_index(k)[p] for p in keys(svec_index(k))])]
                    a == b ? Ms[se.u][a, a] : RT2 * Ms[se.u][a, b]
                end
            ) for in_idx in 1:svec_len(k))
            rhs = sum(se.Tw[out_idx, in_idx] * (
                let (a, b) = collect(keys(svec_index(k)))[findfirst(==(in_idx), [svec_index(k)[p] for p in keys(svec_index(k))])]
                    a == b ? Ms[se.w][a, a] : RT2 * Ms[se.w][a, b]
                end
            ) for in_idx in 1:svec_len(k))
            @constraint(model, lhs == rhs)
        end
    end

    # Pins: diagonal = 1
    for pn in spec.pins
        for a in 1:k
            @constraint(model, Ms[pn.v][a, a] == 1.0)
        end
    end

    return model, Ms
end

# =====================================================================
# SECTION 13.  Benchmarks
# =====================================================================

function run_sdp_benchmark(; optimizer = nothing, dual_optimizer = nothing, solver_name::String = "Mosek", rho::Float64 = 1.0, nwarmup::Int = 2, nruns::Int = 5)
    optimizer === nothing && error("Pass optimizer, e.g. run_sdp_benchmark(optimizer=Mosek.Optimizer)")

    # Pairwise moment sheaf (edge cliques, 3×3 moment matrices)
    inst_grid = grid_mrf_sdp(20)
    inst_cycle = cycle_mrf_sdp(800)
    inst_complete = complete_mrf_sdp(12)
    inst_ladder = ladder_mrf_sdp(250)
    cases = [
        ("grid 20×20",  moment_spec(inst_grid),     inst_grid.V,     length(inst_grid.edges),     1e6),
        ("cycle",       moment_spec(inst_cycle),    inst_cycle.V,    length(inst_cycle.edges),    1e6),
        ("complete",    moment_spec(inst_complete), inst_complete.V, length(inst_complete.edges), 1e6),
        ("ladder",      moment_spec(inst_ladder),   inst_ladder.V,   length(inst_ladder.edges),   1e6),
    ]

    sname = rpad(solver_name, 9)
    sname_d = rpad(solver_name * "-D", 9)
    println("="^80)
    println("Moment SDP Benchmark: Sheaf IPM vs $solver_name  (rho = $rho)")
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

    for (name, (spec, _), nV, nE, raug) in cases
        prob, info = build_sdp_problem(spec; rho = rho)

        settings = IPMSettings{Float64}(
            kkt = UzawaSettings{Float64}(raug = raug),
            feas_tol = 1e-8, gap_tol = 1e-8, itmax = 200)
        for _ in 1:nwarmup; solve(prob, settings); end
        t_ipm = minimum([@elapsed solve(prob, settings) for _ in 1:nruns])

        for _ in 1:nwarmup
            m, _ = build_jump_sdp(spec, optimizer; rho = rho)
            optimize!(m)
        end
        t_mosek = minimum([@elapsed begin
            m, _ = build_jump_sdp(spec, optimizer; rho = rho)
            optimize!(m)
        end for _ in 1:nruns])

        if dual_optimizer !== nothing
            for _ in 1:nwarmup
                m, _ = build_jump_sdp(spec, dual_optimizer; rho = rho)
                optimize!(m)
            end
            t_dual = minimum([@elapsed begin
                m, _ = build_jump_sdp(spec, dual_optimizer; rho = rho)
                optimize!(m)
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

function run_unital_benchmark(; optimizer = nothing, dual_optimizer = nothing, solver_name::String = "Mosek", rho::Float64 = 1.0, nwarmup::Int = 2, nruns::Int = 5)
    optimizer === nothing && error("Pass optimizer")

    # Unital CPTP channels (3×3 moment matrices, 3 Kraus operators)
    cases = [
        ("grid 20×20",  grid_unital(20; k=3, nkraus=3),   400, 2*20*(20-1), 1e6),
        ("cycle",       cycle_unital(400; k=3, nkraus=3), 400, 400,         1e6),
        ("complete",    complete_unital(18; k=3, nkraus=3), 18, 18*17÷2,    1e6),
        ("ladder",      ladder_unital(200; k=3, nkraus=3), 400, 3*200-2,    1e6),
    ]

    sname = rpad(solver_name, 9)
    sname_d = rpad(solver_name * "-D", 9)
    println("="^80)
    println("Unital CPTP Channels: Sheaf IPM vs $solver_name  (rho = $rho)")
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
        prob, info = build_unital_sdp_problem(spec; rho = rho)

        settings = IPMSettings{Float64}(
            kkt = UzawaSettings{Float64}(raug = raug),
            feas_tol = 1e-8, gap_tol = 1e-8, itmax = 200)
        for _ in 1:nwarmup; solve(prob, settings); end
        t_ipm = minimum([@elapsed solve(prob, settings) for _ in 1:nruns])

        for _ in 1:nwarmup
            m, _ = build_jump_unital_sdp(spec, optimizer; rho = rho)
            optimize!(m)
        end
        t_mosek = minimum([@elapsed begin
            m, _ = build_jump_unital_sdp(spec, optimizer; rho = rho)
            optimize!(m)
        end for _ in 1:nruns])

        if dual_optimizer !== nothing
            for _ in 1:nwarmup
                m, _ = build_jump_unital_sdp(spec, dual_optimizer; rho = rho)
                optimize!(m)
            end
            t_dual = minimum([@elapsed begin
                m, _ = build_jump_unital_sdp(spec, dual_optimizer; rho = rho)
                optimize!(m)
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
        optimizer = Mosek.Optimizer
        dual_optimizer = Dualization.dual_optimizer(Mosek.Optimizer)
        solver_name = "Mosek"
    else
        using Clarabel
        optimizer = Clarabel.Optimizer
        dual_optimizer = Dualization.dual_optimizer(Clarabel.Optimizer)
        solver_name = "Clarabel"
    end
    println("Solver: $solver_name")
    println("Runs: $(opts.nruns), Warmup: $(opts.nwarmup)\n")

    run_sdp_benchmark(; optimizer, dual_optimizer, solver_name,
                      nwarmup = opts.nwarmup, nruns = opts.nruns)
    println()
    run_unital_benchmark(; optimizer, dual_optimizer, solver_name,
                         nwarmup = opts.nwarmup, nruns = opts.nruns)
end
