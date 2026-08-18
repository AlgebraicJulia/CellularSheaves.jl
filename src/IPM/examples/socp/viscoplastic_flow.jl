# =============================================================================
#
# # Viscoplastic Flow
#
# A yield stress fluid flows only when subjected to a stress larger than a
# critical value; below it, the material deforms a finite amount and then
# behaves like a solid. Let Ω ⊆ ℝ² be a two-dimensional fluid domain and
#
#   u: Ω → ℝ²
#
# a velocity field. Write
#
#   d = 1/2 (∇u + ∇uᵀ)          ‖ d ‖ = √(2 d:d)
#
# for the strain rate tensor and its norm. Given a viscosity μ > 0, a yield
# shear stress τ₀ ≥ 0, a volume body force density f: Ω → ℝ², and a prescribed
# velocity g: ∂Ω → ℝ², the creeping flow of a Bingham fluid minimizes a sum of
# viscous and plastic dissipation over all incompressible velocity fields:
#
#   min  ∫ μ/2 ‖ d ‖² + τ₀ ‖ d ‖ dx - ∫ f ⋅ u dx        (1)
#    u   Ω                            Ω
#
#   s.t. div u = 0  in  Ω
#            u = g  on ∂Ω
#
# The second term is not differentiable where d = 0, and the solution divides Ω
# into a rigid region where d = 0 and a yielded region where d ≠ 0; the interface
# between them is not known in advance.
#
# We turn (1) into a conic program using the Galerkin method, partitioning Ω
# into a mesh of triangles and searching over velocities that are quadratic on
# each element and continuous across shared edges. Quadratic interpolation
# avoids the volumetric locking that incompressibility induces for linear
# elements. Each candidate uₕ is determined by its values v ∈ ℝⁿ at the nodes
# {x₁, …, xₙ}, one per vertex and one per edge midpoint. Its strain rate is
# linear on each element, so both integrals in (1) are approximated by a
# quadrature rule with points {y₁, …, y_m} and weights {ω₁, …, ω_m}.
#
# Incompressibility makes d traceless, so at each quadrature point the strain
# rate is carried by its two deviatoric components
#
#   q = (d₁₁ - d₂₂, 2d₁₂)       ‖ q ‖ = ‖ d ‖
#
# and q = Gⱼ v for a matrix Gⱼ determined by the element containing yⱼ. Writing
# hⱼ for the corresponding divergence row and epigraphing the plastic term with
# a variable tⱼ ≥ ‖ qⱼ ‖, the problem (1) becomes
#
#              m
#   min    1/2 ∑ ωⱼ μ ‖ qⱼ ‖² + ωⱼ τ₀ tⱼ - cᵀv          (2)
#   v,q,t   j=1
#
#   s.t.  qⱼ - Gⱼ v = 0,   j = 1, …, m
#             hⱼᵀ v = 0,   j = 1, …, m
#            ‖ qⱼ ‖ ≤ tⱼ,  j = 1, …, m
#
# where cᵢ = ∫ f ⋅ φᵢ dx for the vector basis function φᵢ of the ith node, and
# the nodes on ∂Ω are fixed to g. The set { (t, q) : ‖ q ‖ ≤ t } ⊆ ℝ³ is a
# convex cone called the second order cone.
#
# # References
#
#  - Bleyer, Maillard, de Buhan & Coussot,
#    Efficient numerical computations of yield stress fluid flows using
#    second-order cone programming,
#    Comput. Methods Appl. Mech. Engrg. 283:599-614, 2015.
#
# =============================================================================

using AppleAccelerate
using LinearAlgebra, Printf

using CellularSheaves.IPM
using CellularSheaves.IPM: CofreeCone, SecondOrderCone, IPMSettings, HSDSettings
using CellularSheaves.BlockSparseArrays: colrange, blocksparse

include("../utils.jl")

# 3-point Gauss rule on the reference triangle (weights sum to 1/2, its area)
const GAUSS = [((1 / 6, 1 / 6), 1 / 6), ((2 / 3, 1 / 6), 1 / 6), ((1 / 6, 2 / 3), 1 / 6)]

# -----------------------------------------------------------------------------
# P2 reference element
# -----------------------------------------------------------------------------

# Gradients of the six P2 shape functions on the reference triangle, as a 2 × 6
# matrix. Nodes 1:3 are the vertices, 4:6 the midpoints of (1,2), (2,3), (3,1).
function p2grads(xi::Float64, eta::Float64)
    lam = (1 - xi - eta, xi, eta)
    dlam = [-1.0 1.0 0.0; -1.0 0.0 1.0]          # d(λ)/d(ξ, η), 2 × 3
    g = zeros(2, 6)

    for k in 1:3
        g[:, k] .= (4 * lam[k] - 1) .* view(dlam, :, k)
    end

    for (k, (a, b)) in enumerate(((1, 2), (2, 3), (3, 1)))
        g[:, 3 + k] .= 4 .* (view(dlam, :, a) .* lam[b] .+ lam[a] .* view(dlam, :, b))
    end

    return g
end

# -----------------------------------------------------------------------------
# Mesh
# -----------------------------------------------------------------------------

struct P2Mesh
    pts::Matrix{Float64}          # 2 × nnode
    elems::Matrix{Int}            # 6 × nelem
end

nnodes(m::P2Mesh) = size(m.pts, 2)
nelems(m::P2Mesh) = size(m.elems, 2)

# The unit square, N × N cells, each split into two triangles, with midside
# nodes added.
function unitsquare(N::Int)
    vid(i, j) = j * (N + 1) + i + 1                       # 0-based (i, j)
    nv = (N + 1)^2

    xy = zeros(2, nv)

    for j in 0:N, i in 0:N
        xy[1, vid(i, j)] = i / N
        xy[2, vid(i, j)] = j / N
    end

    tris = NTuple{3, Int}[]

    for j in 0:(N - 1), i in 0:(N - 1)
        a, b = vid(i, j), vid(i + 1, j)
        c, d = vid(i + 1, j + 1), vid(i, j + 1)
        push!(tris, (a, b, c))
        push!(tris, (a, c, d))
    end

    mid = Dict{Tuple{Int, Int}, Int}()
    extra = Vector{Float64}[]

    function midnode(u, v)
        key = (min(u, v), max(u, v))

        if !haskey(mid, key)
            push!(extra, 0.5 .* (xy[:, u] .+ xy[:, v]))
            mid[key] = nv + length(extra)
        end

        return mid[key]
    end

    elems = zeros(Int, 6, length(tris))

    for (e, (a, b, c)) in enumerate(tris)
        elems[:, e] .= (a, b, c, midnode(a, b), midnode(b, c), midnode(c, a))
    end

    pts = hcat(xy, reduce(hcat, extra))
    return P2Mesh(pts, elems)
end

# -----------------------------------------------------------------------------
# Boundary conditions and stalk layout
# -----------------------------------------------------------------------------

# Which velocity components of each node are free, and the prescribed values of
# the rest. `stalk[i]` is the velocity stalk index of node i, or 0 if the node
# is fully constrained; `comps[i]` lists its free components.
struct Layout
    stalk::Vector{Int}
    comps::Vector{Vector{Int}}
    fixed::Matrix{Float64}        # 2 × nnode, prescribed values (0 where free)
    nvstalk::Int
end

# Lid-driven cavity: no slip everywhere, except the top edge which moves at
# (V, 0). Every boundary node is fully constrained.
function cavity_layout(m::P2Mesh; V::Float64 = 1.0)
    n = nnodes(m)
    stalk = zeros(Int, n)
    comps = [Int[] for _ in 1:n]
    fixed = zeros(2, n)
    k = 0

    for i in 1:n
        x, y = m.pts[1, i], m.pts[2, i]
        onbnd = x ≈ 0 || x ≈ 1 || y ≈ 0 || y ≈ 1

        if onbnd
            fixed[1, i] = (y ≈ 1) ? V : 0.0
            fixed[2, i] = 0.0
        else
            k += 1
            stalk[i] = k
            comps[i] = [1, 2]
        end
    end

    return Layout(stalk, comps, fixed, k)
end

# -----------------------------------------------------------------------------
# Element quantities
# -----------------------------------------------------------------------------


localcols(k::Int) = (2k - 1):(2k)

# Returns (weight, Ĝ, ĥ) per Gauss point, with Ĝ (2 × 12) the deviatoric strain
# operator and ĥ (1 × 12) the divergence. Incompressibility makes d traceless, so
# on the feasible set ‖d‖² = 2dxx² + 2dyy² + 4dxy² = (dxx - dyy)² + (2dxy)²: the
# strain rate is carried by the two deviatoric components q = (dxx - dyy, 2dxy),
# the cone is 3-dimensional (t, q) ∈ Q³, and the divergence becomes its own row
# on v — 3 columns and 3 rows per Gauss point.
function element_ops(m::P2Mesh, e::Int)
    el = view(m.elems, :, e)
    p1, p2, p3 = m.pts[:, el[1]], m.pts[:, el[2]], m.pts[:, el[3]]
    J = hcat(p2 .- p1, p3 .- p1)
    detJ = abs(det(J))
    Jit = transpose(inv(J))

    ops = Tuple{Float64, Matrix{Float64}, Vector{Float64}}[]

    for ((xi, eta), w) in GAUSS
        gr = Jit * p2grads(xi, eta)
        DN = zeros(3, 12)
        DN[1, 1:2:end] .= view(gr, 1, :)
        DN[2, 2:2:end] .= view(gr, 2, :)
        DN[3, 1:2:end] .= view(gr, 2, :)
        DN[3, 2:2:end] .= view(gr, 1, :)

        Ghat = vcat(transpose(DN[1, :] .- DN[2, :]), transpose(DN[3, :]))
        hhat = DN[1, :] .+ DN[2, :]
        push!(ops, (w * detJ, Ghat, hhat))
    end

    return ops
end

# -----------------------------------------------------------------------------
# Build
# -----------------------------------------------------------------------------

"""
    build_bingham(N; Bi, mu, V)

Lid-driven cavity for a Bingham fluid on an N × N structured mesh. The viscous
term is carried by the solver's quadratic objective and the plastic term by the
cone.

Stalks: 1..nv velocity nodes (Cofree), then nv + (e-1)*ng + g holding
(t, q₁, q₂) (SecondOrderCone). Row block per Gauss point: rows 1:2 are the
strain rows q - Ĝv = 0, row 3 is the incompressibility row ĥv = 0.
"""
function build_bingham(N::Int; Bi::Float64 = 2.0, mu::Float64 = 1.0, V::Float64 = 1.0)
    m = unitsquare(N)
    lay = cavity_layout(m; V)
    tau0 = Bi * mu
    ng, NE = length(GAUSS), nelems(m)
    nv = lay.nvstalk
    nsoc = ng * NE

    socstalk(e, g) = nv + (e - 1) * ng + g
    rowblk(e, g) = (e - 1) * ng + g

    coldims = vcat(fill(0, nv), fill(3, nsoc))
    for i in 1:nnodes(m)
        lay.stalk[i] > 0 && (coldims[lay.stalk[i]] = length(lay.comps[i]))
    end
    rowdims = fill(3, nsoc)

    Brow, Bcol, Bblk = Int[], Int[], Matrix{Float64}[]
    g = zeros(3 * nsoc)
    c = zeros(sum(coldims))
    coff = cumsum(vcat(0, coldims))
    Qacc = Dict{Tuple{Int, Int}, Matrix{Float64}}()

    # the q-side of the linking rows; the divergence row has no cone entries
    Ssoc = zeros(3, 3)
    Ssoc[1:2, 2:3] .= Matrix{Float64}(I, 2, 2)

    for e in 1:NE
        el = view(m.elems, :, e)

        for (gi, (wt, Ghat, hhat)) in enumerate(element_ops(m, e))
            r = rowblk(e, gi)
            s = socstalk(e, gi)

            push!(Brow, r); push!(Bcol, s); push!(Bblk, copy(Ssoc))
            c[coff[s] + 1] = -wt * tau0

            for (k, i) in enumerate(el)
                cols = localcols(k)

                if lay.stalk[i] > 0
                    blk = zeros(3, length(lay.comps[i]))
                    blk[1:2, :] .= .-view(Ghat, :, cols[lay.comps[i]])
                    blk[3, :] .= view(hhat, cols[lay.comps[i]])
                    push!(Brow, r); push!(Bcol, lay.stalk[i]); push!(Bblk, blk)
                end

                for j in 1:2
                    j in lay.comps[i] && continue
                    g[3 * (r - 1) .+ (1:2)] .+= view(Ghat, :, cols[j]) .* lay.fixed[j, i]
                    g[3 * (r - 1) + 3] -= hhat[cols[j]] * lay.fixed[j, i]
                end
            end

            KE = (wt * mu) .* (transpose(Ghat) * Ghat)

            for (k, i) in enumerate(el), (l, j) in enumerate(el)
                si = lay.stalk[i]
                si > 0 || continue

                blk = view(KE, localcols(k), localcols(l))
                sj = lay.stalk[j]

                if sj > 0
                    A = get!(() -> zeros(length(lay.comps[i]), length(lay.comps[j])),
                             Qacc, (si, sj))
                    A .+= view(blk, lay.comps[i], lay.comps[j])
                end

                pres = setdiff(1:2, lay.comps[j])

                if !isempty(pres)
                    c[coff[si] .+ (1:length(lay.comps[i]))] .-=
                        view(blk, lay.comps[i], pres) * view(lay.fixed, pres, j)
                end
            end
        end
    end

    Qrow, Qcol, Qblk = Int[], Int[], Matrix{Float64}[]

    for ((i, j), A) in Qacc
        push!(Qrow, i); push!(Qcol, j); push!(Qblk, A)
    end

    for s in (nv + 1):(nv + nsoc)
        push!(Qrow, s); push!(Qcol, s); push!(Qblk, zeros(3, 3))
    end

    Q = blocksparse(Qrow, Qcol, Qblk, coldims, coldims)
    B = blocksparse(Brow, Bcol, Bblk, rowdims, coldims)

    cones = vcat(IPM.AbstractCone[CofreeCone() for _ in 1:nv],
                 IPM.AbstractCone[SecondOrderCone() for _ in 1:nsoc])

    return IPMProblem(Q, B, c, g, cones)
end

# JuMP twin of `build_bingham`.
function build_bingham_jump(N::Int; Bi::Float64 = 2.0, mu::Float64 = 1.0,
                                V::Float64 = 1.0, optimizer)
    m = unitsquare(N)
    lay = cavity_layout(m; V)
    tau0 = Bi * mu
    ng, NE = length(GAUSS), nelems(m)

    model = Model(optimizer)
    @variable(model, v[1:2, 1:nnodes(m)])

    for i in 1:nnodes(m), j in 1:2
        j in lay.comps[i] || fix(v[j, i], lay.fixed[j, i]; force = true)
    end

    @variable(model, t[1:NE, 1:ng])
    obj = QuadExpr()

    for e in 1:NE
        el = view(m.elems, :, e)

        for (gi, (wt, Ghat, hhat)) in enumerate(element_ops(m, e))
            q = [sum(Ghat[a, localcols(k)[j]] * v[j, el[k]] for k in 1:6, j in 1:2)
                 for a in 1:2]

            @constraint(model, vcat(t[e, gi], q) in JuMP.SecondOrderCone())
            @constraint(model, sum(hhat[localcols(k)[j]] * v[j, el[k]]
                                   for k in 1:6, j in 1:2) == 0)

            add_to_expression!(obj, wt * tau0, t[e, gi])
            add_to_expression!(obj, wt * mu / 2, sum(qa^2 for qa in q))
        end
    end

    @objective(model, Min, obj)
    return model
end

# -----------------------------------------------------------------------------
# Reporting
# -----------------------------------------------------------------------------

# Horizontal velocity along the vertical mid-plane x = 1/2, for comparison with
# the reference profiles.
function midplane_profile(N::Int, prob, res)
    m = unitsquare(N)
    lay = cavity_layout(m)
    ys, us = Float64[], Float64[]

    for i in 1:nnodes(m)
        m.pts[1, i] ≈ 0.5 || continue
        push!(ys, m.pts[2, i])

        if lay.stalk[i] == 0
            push!(us, lay.fixed[1, i])
        else
            push!(us, res.p[first(colrange(prob.B, lay.stalk[i]))])
        end
    end

    ord = sortperm(ys)
    return ys[ord], us[ord]
end

# -----------------------------------------------------------------------------
# Benchmark
# -----------------------------------------------------------------------------

const BINGHAM_NUMBERS = [0.0, 2.0, 5.0, 20.0, 50.0, 200.0, 500.0]

const PROBLEMS = [
    ("Cavity-11", 11),      # NE =  242, comparable to the reference's mesh 1
    ("Cavity-32", 32),      # NE = 2048, ~ mesh 2
    ("Cavity-54", 54),      # NE = 5832, ~ mesh 3
]

function benchmark(; Bi = 500.0, tol = 1e-6, clarabel = false, mosek = false)
    @printf("%-11s %6s %7s %7s %9s %9s %9s %9s %8s\n",
            "problem", "NE", "n", "m", "IPM", "HSD", "Clarabel", "Mosek", "Cla/IPM")

    for (name, N) in PROBLEMS
        prob = build_bingham(N; Bi)

        mi = measure_ipm(prob, IPMSettings{Float64}(feas_tol = tol, gap_tol = tol))
        mh = measure_ipm(prob, HSDSettings{Float64}(feas_tol = tol, gap_tol = tol))

        if clarabel
            mc = measure_jump(() -> build_bingham_jump(N; Bi, optimizer = clarabel_opt(; tol)))
        else
            mc = (t = NaN, status = "—", obj = NaN)
        end

        if mosek
            mk = measure_jump(() -> build_bingham_jump(N; Bi, optimizer = mosek_opt(; tol)))
        else
            mk = (t = NaN, status = "—", obj = NaN)
        end

        ratio = (isfinite(mc.t) && isfinite(mi.t)) ? @sprintf("%.2f", mc.t / mi.t) : "—"

        @printf("%-11s %6d %7d %7d %s %s %s %s %8s\n",
                name, 2N^2, size(prob.B, 2), size(prob.B, 1),
                fmt_time(mi), fmt_time(mh), fmt_time(mc), fmt_time(mk), ratio)
    end
end

# Sweep the Bingham number on a fixed mesh: at Bi = 0 the problem is a Stokes
# flow and the cones are inactive; as Bi grows the yield term dominates and the
# problem approaches a limit analysis. The quadratic term goes from carrying the
# whole objective to carrying none of it.
function sweep(; N = 32, tol = 1e-6)
    @printf("%8s %9s %9s %12s\n", "Bi", "IPM", "HSD", "objective")

    for Bi in BINGHAM_NUMBERS
        prob = build_bingham(N; Bi)
        mi = measure_ipm(prob, IPMSettings{Float64}(feas_tol = tol, gap_tol = tol))
        mh = measure_ipm(prob, HSDSettings{Float64}(feas_tol = tol, gap_tol = tol))
        @printf("%8.1f %s %s %12.6f\n", Bi, fmt_time(mi), fmt_time(mh), mi.obj)
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    opts = parse_args()
    tol = any(a -> startswith(a, "--tol="), ARGS) ? opts.tol : 1e-6
    benchmark(; clarabel = opts.clarabel, mosek = opts.mosek, tol)

    if "--sweep" in ARGS
        println()
        sweep(; tol)
    end
end

# =============================================================================
# Sample run: 2026-08-14 (Apple M-series, --clarabel --mosek), tol = 1e-6, Bi = 500.
#
#   problem      NE       n       m       IPM       HSD  Clarabel     Mosek  Cla/IPM
#   Cavity-11   242    3060    2178   22.5ms   25.8ms   22.5ms   39.9ms     1.00
#   Cavity-32  2048   26370   18432  276.2ms  315.4ms  335.3ms  445.6ms     1.21
#   Cavity-54  5832   75386   52488 1095.2ms 1299.0ms 2187.9ms 2084.7ms     2.00
# =============================================================================
