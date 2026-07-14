# ipm_harmonic.jl — the CONIC harmonic extension, solved by Richard's sheaf interior-point method.
#
# The ECC-paper coordination step is the harmonic extension
#       q* = argmin_z  ½‖δz‖²   s.t.  z_b = p_b  (pinned boundary),
# i.e. minimize sheaf disagreement with the targets fixed. `CellularSheaves.harmonic_extension`
# does this as a linear solve. Here we solve the SAME problem as a sheaf conic QP through the IPM,
# in the primal form (P): min ½pᵀQp + cᵀp  s.t. Bp = g, p ∈ K.
#
# Variables (IPM "cells"):  a vertex cell z_v (free) per sheaf vertex,
#                           an edge cell e_k  per sheaf edge, with identity cost  ½‖e_k‖².
# Constraints (IPM rows):   e_k − ρ_u z_u + ρ_v z_v = 0   (coboundary consistency, one per edge),
#                           z_b = p_b                     (boundary pin, one per target).
# With every vertex cone FREE this reproduces `harmonic_extension` exactly (validation step 1);
# swap a vertex's free cone for an SOC / positive cone to get a FEASIBILITY-CONSTRAINED
# coordination target (step 2) — that is the "conic harmonic extension".
#
#   julia --project=. -e 'include("examples/RL/lib/ipm_harmonic.jl")'   # loads (root has IPM)
using CellularSheaves
using CellularSheaves: vertex_stalks, edge_stalks, underlying_graph, get_restriction_map, get_edge_stalk
import CellularSheaves.NetworkSheaves.EuclideanSheaves: harmonic_extension
const IPM = CellularSheaves.IPM
using CellularSheaves.IPM: IPMProblem, IPMSettings, UzawaSettings, CofreeCone, SecondOrderCone,
                           PositiveCone, AbstractCone, solve
using CellularSheaves.BlockSparseArrays: blocksparse, block, colrange, rowrange
using Graphs: edges, src, dst
using LinearAlgebra

# raug=1e7 / itmax=200: robust across the multi-cone families (a subset SOC cap is ill-conditioned at
# raug=1e5 → NUMERICAL_FAILURE, though the answer is still correct; the stronger augmentation converges
# it cleanly). Easy problems (free, full ball) still converge early, so the itmax bump is nearly free.
default_ipm_settings() = IPMSettings{Float64}(itmax = 200, kkt = UzawaSettings{Float64}(raug = 1e7))

"""
    ipm_harmonic_problem(s, boundary; vertex_cones=nothing, ball=nothing) -> (prob, info)

Assemble the sheaf conic QP once. In the control loop the sheaf STRUCTURE (Q, B, cones) is fixed and
only the boundary changes → build once, then `update_boundary!` + `solve` per step (see below).
`info` carries the per-vertex column ranges for reading `q*` out of `res.p`, and the boundary-row
ranges for cheap `g` updates.

`ball::Dict{Int,Float64}` (vertex → radius r) adds a REACHABILITY second-order cone `‖q*_v‖ ≤ r` on
each listed vertex. The IPM's `SecondOrderCone` on an n-vector is the ice-cream cone `x₁ ≥ ‖x₂..ₙ‖`,
not a ball, so each capped vertex gets an auxiliary cell `t_v = (r; q*_v) ∈ ℝ^{1+d}` pinned to radius
`r` in its first coordinate and coupled to the vertex in the rest; the SOC on `t_v` then reads
`r ≥ ‖q*_v‖`. With `ball=nothing` this is exactly the linear harmonic extension.
"""
function ipm_harmonic_problem(s, boundary::Dict{Int,<:AbstractVector}; vertex_cones=nothing, ball=nothing)
    vs   = vertex_stalks(s)
    NV   = length(vs)
    el   = collect(edges(underlying_graph(s)))
    NE   = length(el)
    edgecell(k) = NV + k                                    # edge cells follow the vertex cells

    row_ids = Int[]; col_ids = Int[]; blks = Matrix{Float64}[]
    place!(r, c, M) = (push!(row_ids, r); push!(col_ids, c); push!(blks, Matrix{Float64}(M)))
    gval = Dict{Int,Vector{Float64}}(); rc = Ref(0); edim = Int[]

    for (k, e) in enumerate(el)                             # coboundary rows: e_k − ρ_u z_u + ρ_v z_v = 0
        u, v = src(e), dst(e)
        d  = get_edge_stalk(s, u, v)
        ρu = get_restriction_map(s, u, v)                  # d × stalk(u)
        ρv = get_restriction_map(s, v, u)                  # d × stalk(v)
        r  = (rc[] += 1)
        place!(r, edgecell(k), Matrix{Float64}(I, d, d))
        place!(r, u, -ρu); place!(r, v, ρv)
        gval[r] = zeros(d); push!(edim, d)
    end
    brow = Dict{Int,Int}()                                  # boundary vertex → its constraint-row id
    for (b, val) in boundary                                # boundary rows: z_b = p_b
        r = (rc[] += 1); brow[b] = r
        place!(r, b, Matrix{Float64}(I, vs[b], vs[b]))
        gval[r] = Vector{Float64}(val)
    end
    NCAP = 0                                                # reachability caps: t_v − S z_v = (r;0), SOC on t_v
    if ball !== nothing
        for (v, radius) in ball
            d = vs[v]; capcell = NV + NE + (NCAP += 1)
            S = zeros(1 + d, d); for i in 1:d; S[i + 1, i] = 1.0; end
            r = (rc[] += 1)
            place!(r, capcell, Matrix{Float64}(I, 1 + d, 1 + d)); place!(r, v, -S)
            gval[r] = vcat(Float64(radius), zeros(d))
        end
    end

    B = blocksparse(row_ids, col_ids, blks)
    g = zeros(size(B, 1)); for (r, val) in gval; g[rowrange(B, r)] .= val; end
    Q = IPM.allocblockdiag(B); fill!(Q, 0)                  # ½‖e‖²: identity on edge cells, 0 elsewhere
    for k in 1:NE; block(Q, edgecell(k), edgecell(k), edgecell(k)) .= Matrix{Float64}(I, edim[k], edim[k]); end
    c = zeros(size(B, 2))

    cones = AbstractCone[]
    for cell in 1:(NV + NE + NCAP)
        push!(cones, cell ≤ NV ? (vertex_cones !== nothing ? vertex_cones[cell] : CofreeCone()) :
                     cell ≤ NV + NE ? CofreeCone() : SecondOrderCone())
    end

    prob = IPMProblem(Q, B, c, g, cones)
    info = (B = B, NV = NV, vs = vs, brow = brow, off = [0; cumsum(vs)])
    prob, info
end

"Update the pinned targets in-place (cheap: only `g` changes, structure is reused)."
function update_boundary!(prob, info, boundary::Dict{Int,<:AbstractVector})
    for (b, val) in boundary; prob.g[rowrange(info.B, info.brow[b])] .= val; end
    prob
end

"Read the vertex cochain q* out of an IPM result."
function extract_qstar(res, info)
    z = zeros(info.off[end])
    for v in 1:info.NV; z[info.off[v]+1:info.off[v+1]] = res.p[colrange(info.B, v)]; end
    z
end

"""One-shot conic harmonic extension (build + solve + extract). See `ipm_harmonic_problem` for the
loop-efficient path that reuses the assembled problem across steps."""
function build_ipm_harmonic(s, boundary::Dict{Int,<:AbstractVector};
                            vertex_cones=nothing, ball=nothing, settings=default_ipm_settings())
    prob, info = ipm_harmonic_problem(s, boundary; vertex_cones, ball)
    res = solve(prob, settings)
    extract_qstar(res, info), res
end
