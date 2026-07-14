#
# Conic (SOC) receding-horizon MPC via SheafSDP's primal-dual interior point method.
#
# Non-trivial conic program (minimum-fuel rendezvous):
#   N planar double-integrator agents must meet (terminal position consensus on a
#   graph) while minimizing total fuel measured as Σ_t ‖u_t‖₂ — a genuinely conic,
#   *non-quadratic* objective. The ℓ₂ (un-squared) norm switches whole control
#   vectors on/off, producing sparse, bang-off-bang thrusts rather than the smooth
#   spread-out controls a quadratic ½‖u‖² penalty gives. This is the group-sparse
#   SOC formulation from SheafSDP.jl/test/small/soc.jl, here driven in closed loop.
#
# SOC epigraph bundle (per agent i, step t):
#   ζ_i^t = (s_i^t; u_i^t) ∈ Q^{1+nu},  objective Σ s,  hard box |u| ≤ ū.
#
# The sheaf coboundary B encodes per-agent dynamics + box slacks (POS cones) and
# the cross-agent terminal-consensus rows (the sheaf restriction maps). The IPM
# solves the resulting block-structured KKT systems with a chordal factorization.
#
# This module is solver-API faithful to SheafSDP; the cone scaling (1/√2 on the
# control-extraction maps and the cost head) matches SheafSDP's isometric SOC
# coordinates, where the cone identity is √2·e₁.
#

module SOCMPC

using LinearAlgebra
using SparseArrays
using BlockSparseArrays: colrange, rowrange, blocksparse, block

export SOCParams, build_soc_problem, SOCLayout, extract_u, extract_x, extract_s,
       feasibility_residual, col_x, col_ζ, col_sp, col_sm

# ---------------------------------------------------------------------------
# Problem parameters
# ---------------------------------------------------------------------------

"""
    SOCParams

Planar double-integrator fleet for the minimum-fuel rendezvous conic program.

Fields
- `N`        : number of agents
- `T`        : planning-horizon length (states x^1..x^T, controls u^1..u^{T-1})
- `h`        : timestep
- `ubar`     : per-component control box bound |u| ≤ ubar
- `edges`    : consensus edges (default complete graph Kₙ)
- `A`,`B`    : discrete double-integrator matrices (nx=4, nu=2)
- `Pproj`    : terminal-position projection (picks x,y out of [x,y,vx,vy])
"""
struct SOCParams
    N::Int
    T::Int
    h::Float64
    ubar::Float64
    edges::Vector{Tuple{Int,Int}}
    A::Matrix{Float64}
    B::Matrix{Float64}
    Pproj::Matrix{Float64}
end

function SOCParams(; N=5, T=12, h=0.1, ubar=4.0, edges=nothing)
    A = [1.0 0.0 h   0.0;
         0.0 1.0 0.0 h;
         0.0 0.0 1.0 0.0;
         0.0 0.0 0.0 1.0]
    B = [0.0 0.0;
         0.0 0.0;
         h   0.0;
         0.0 h]
    Pproj = [1.0 0.0 0.0 0.0;
             0.0 1.0 0.0 0.0]
    e = edges === nothing ? [(i, j) for i in 1:N for j in i+1:N] : edges
    return SOCParams(N, T, h, ubar, e, A, B, Pproj)
end

nx(::SOCParams) = 4
nu(::SOCParams) = 2
npos(p::SOCParams) = size(p.Pproj, 1)

# ---------------------------------------------------------------------------
# Block layout — column/row index maps (mirrors SheafSDP test/small/soc.jl)
# ---------------------------------------------------------------------------

"""
    SOCLayout

Column/row block-index bookkeeping for one `build_soc_problem` instance, plus the
constructed coboundary `B` so solution blocks can be sliced back out.
"""
struct SOCLayout
    p::SOCParams
    blocks_per_agent::Int
    rows_per_agent::Int
    ne::Int
    B  # BlockSparseMatrix
end

# Column blocks per agent: T states, then (T-1) × (ζ, box+slack, box-slack)
col_x(L::SOCLayout, i, t) = (i - 1) * L.blocks_per_agent + t
col_ζ(L::SOCLayout, i, t) = (i - 1) * L.blocks_per_agent + L.p.T + 3 * (t - 1) + 1
col_sp(L::SOCLayout, i, t) = (i - 1) * L.blocks_per_agent + L.p.T + 3 * (t - 1) + 2
col_sm(L::SOCLayout, i, t) = (i - 1) * L.blocks_per_agent + L.p.T + 3 * (t - 1) + 3

# Row blocks per agent: init, then (T-1) × (dynamics, box+, box-); then edge rows
row_init(L::SOCLayout, i) = (i - 1) * L.rows_per_agent + 1
row_dyn(L::SOCLayout, i, t) = (i - 1) * L.rows_per_agent + 1 + t
row_boxp(L::SOCLayout, i, t) = (i - 1) * L.rows_per_agent + L.p.T + (t - 1) + 1
row_boxm(L::SOCLayout, i, t) = (i - 1) * L.rows_per_agent + L.p.T + (L.p.T - 1) + (t - 1) + 1
row_coord(L::SOCLayout, e) = L.p.N * L.rows_per_agent + e

# ---------------------------------------------------------------------------
# Problem construction
# ---------------------------------------------------------------------------

"""
    build_soc_problem(p::SOCParams, x0; consensus=:terminal)

Build the SOC minimum-fuel rendezvous program from current states `x0`
(`Vector` of length-`nx` agent states). Returns `(c, g, B, Q, cones, layout)`
in the exact form SheafSDP's `IPMProblem(c, g, B, Q, cones)` expects.

Cone scaling: ζ lives in SheafSDP's isometric SOC coordinates, so the stored
block is √2·(s; u). Every map that reads `u` (or `s`) out of ζ — the dynamics'
B-application, the box-slack extraction, and the cost head — is scaled by 1/√2.
"""
function build_soc_problem(p::SOCParams, x0::Vector{<:AbstractVector})
    @assert length(x0) == p.N "need one initial state per agent"
    nxv, nuv = nx(p), nu(p)
    T = p.T
    edges = p.edges
    ne = length(edges)

    blocks_per_agent = T + 3 * (T - 1)
    rows_per_agent = 1 + 3 * (T - 1)
    L = SOCLayout(p, blocks_per_agent, rows_per_agent, ne, nothing)

    invrt2 = 1 / sqrt(2.0)
    # [0 B] reads u from ζ=(s;u) and applies the input matrix; scaled to isometric coords
    B_on_ζ = hcat(zeros(nxv, 1), p.B) .* invrt2          # nx × (1+nu)
    extract_u_box = hcat(zeros(nuv, 1), Matrix(1.0I, nuv, nuv)) .* invrt2  # nu × (1+nu)

    row_ids, col_ids, blocks = Int[], Int[], Matrix{Float64}[]

    for i in 1:p.N
        # init: x_i^1 = x0_i
        push!(row_ids, row_init(L, i)); push!(col_ids, col_x(L, i, 1))
        push!(blocks, Matrix(1.0I, nxv, nxv))

        for t in 1:T-1
            # dynamics: x^{t+1} - A x^t - B u^t = 0  (u from ζ tail)
            push!(row_ids, row_dyn(L, i, t)); push!(col_ids, col_x(L, i, t)); push!(blocks, -p.A)
            push!(row_ids, row_dyn(L, i, t)); push!(col_ids, col_x(L, i, t + 1)); push!(blocks, Matrix(1.0I, nxv, nxv))
            push!(row_ids, row_dyn(L, i, t)); push!(col_ids, col_ζ(L, i, t)); push!(blocks, -B_on_ζ)

            # box+: u + s+ = ū
            push!(row_ids, row_boxp(L, i, t)); push!(col_ids, col_ζ(L, i, t)); push!(blocks, extract_u_box)
            push!(row_ids, row_boxp(L, i, t)); push!(col_ids, col_sp(L, i, t)); push!(blocks, Matrix(1.0I, nuv, nuv))

            # box-: -u + s- = ū
            push!(row_ids, row_boxm(L, i, t)); push!(col_ids, col_ζ(L, i, t)); push!(blocks, -extract_u_box)
            push!(row_ids, row_boxm(L, i, t)); push!(col_ids, col_sm(L, i, t)); push!(blocks, Matrix(1.0I, nuv, nuv))
        end
    end

    # Consensus (sheaf restriction maps): P x_i^T - P x_j^T = 0
    for (e, (i, j)) in enumerate(edges)
        push!(row_ids, row_coord(L, e)); push!(col_ids, col_x(L, i, T)); push!(blocks, -p.Pproj)
        push!(row_ids, row_coord(L, e)); push!(col_ids, col_x(L, j, T)); push!(blocks, p.Pproj)
    end

    B = blocksparse(row_ids, col_ids, blocks)
    L = SOCLayout(p, blocks_per_agent, rows_per_agent, ne, B)

    # Cost: c = 1/√2 on the head (s) of each ζ
    c = zeros(size(B, 2))
    for i in 1:p.N, t in 1:T-1
        c[colrange(B, col_ζ(L, i, t))[1]] = invrt2
    end

    # RHS
    g = zeros(size(B, 1))
    for i in 1:p.N
        g[rowrange(B, row_init(L, i))] .= x0[i]
        for t in 1:T-1
            g[rowrange(B, row_boxp(L, i, t))] .= p.ubar
            g[rowrange(B, row_boxm(L, i, t))] .= p.ubar
        end
    end

    # Q = 0 (objective is purely the conic Σ‖u‖₂, no quadratic term)
    nv = p.N * blocks_per_agent
    cones = Vector{Symbol}(undef, nv)
    for i in 1:p.N
        for t in 1:T
            cones[col_x(L, i, t)] = :NOC
        end
        for t in 1:T-1
            cones[col_ζ(L, i, t)] = :SOC
            cones[col_sp(L, i, t)] = :POS
            cones[col_sm(L, i, t)] = :POS
        end
    end

    return c, g, B, cones, L
end

# ---------------------------------------------------------------------------
# Solution extraction (undo the 1/√2 isometric scaling)
# ---------------------------------------------------------------------------

const INVRT2 = 1 / sqrt(2.0)

"True control u_i^t recovered from the ζ block of a primal solution `pvec`."
function extract_u(L::SOCLayout, pvec::AbstractVector, i, t)
    r = colrange(L.B, col_ζ(L, i, t))
    return pvec[r[2:end]] .* INVRT2
end

"State x_i^t from a primal solution."
function extract_x(L::SOCLayout, pvec::AbstractVector, i, t)
    return pvec[colrange(L.B, col_x(L, i, t))]
end

"Epigraph value s_i^t (≈ ‖u_i^t‖) from a primal solution."
function extract_s(L::SOCLayout, pvec::AbstractVector, i, t)
    r = colrange(L.B, col_ζ(L, i, t))
    return pvec[r[1]] * INVRT2
end

# ---------------------------------------------------------------------------
# Structural feasibility check (runs WITHOUT the IPM solver)
# ---------------------------------------------------------------------------

"""
    feasibility_residual(p, L, B, g; controls) -> Float64

Construct an explicit feasible primal point from a chosen control schedule
(roll the dynamics forward, set box slacks and the SOC epigraph s=‖u‖), and
return ‖B p − g‖∞. A value ≈ 0 certifies the affine layout (dynamics + box +
consensus rows, cone scaling, RHS) is correctly wired — independent of the
solver. Consensus rows are only satisfiable if the chosen controls actually
rendezvous, so pass `check_consensus=false` to skip those rows.
"""
function feasibility_residual(p::SOCParams, L::SOCLayout, B, g;
                              controls=nothing, check_consensus=false, x0)
    nxv, nuv = nx(p), nu(p)
    T = p.T
    pvec = zeros(size(B, 2))

    if controls === nothing
        controls = [[zeros(nuv) for _ in 1:T-1] for _ in 1:p.N]
    end

    for i in 1:p.N
        x = copy(x0[i])
        pvec[colrange(B, col_x(L, i, 1))] .= x
        for t in 1:T-1
            u = controls[i][t]
            # ζ stored as √2·(s; u)
            r = colrange(B, col_ζ(L, i, t))
            pvec[r[1]] = sqrt(2.0) * norm(u)
            pvec[r[2:end]] .= sqrt(2.0) .* u
            # box slacks: s± = ū ∓ u
            pvec[colrange(B, col_sp(L, i, t))] .= p.ubar .- u
            pvec[colrange(B, col_sm(L, i, t))] .= p.ubar .+ u
            # propagate
            x = p.A * x + p.B * u
            pvec[colrange(B, col_x(L, i, t + 1))] .= x
        end
    end

    r = B * pvec - g
    if !check_consensus
        # zero out consensus rows (only feasible if controls rendezvous)
        for e in 1:L.ne
            r[rowrange(B, row_coord(L, e))] .= 0.0
        end
    end
    return norm(r, Inf), pvec
end

end # module SOCMPC
