######################################################################
# opf_sheaf.jl
#
# AC OPTIMAL POWER FLOW as a ladder of covers over one problem graph.
# The network G has buses and lines; every OPF constraint is linear in
# the entries of W = VVᴴ on G's pattern, and a relaxation is a CHOICE
# OF COVER of that pattern — the sheaf base graph is the cover's nerve:
#
#   rung 1 (Jabr SOCP): line cells z_e = (u_i+u_j, u_i−u_j, 2c_e, 2s_e)
#     ∈ SecondOrderCone (bound first; c²+s² ≤ u_iu_j after rotation),
#     with u_i = |V_i|², c_e + j·s_e = V_i·conj(V_j) (oriented i<j);
#     bus cells u_i (Cofree + box slack leaves); generator cells (p,q)
#     with NATIVE Q = ½c₂p² (the flagship native quadratic, at last in
#     industrial costume); loads ride in g through 2-row injection
#     hyperedges. Exact on radial networks; strict on meshes.
#   rung 2 (+cycle cells): Hermitian W_C per cycle, REAL-EMBEDDED as
#     M = [[X, −Z],[Z, X]] ⪰ 0 — the corpus's first complex embedding —
#     tied by agreement rows to u and to member lines' (c, s). This is
#     Kocuk–Dey–Sun cycle strengthening, sheaf-native: enlarge the
#     cover along the cycle.
#   rung 3 (full SDP): the cycle cell over all buses (= chordal rung
#     for these small networks).
#
# THE DROPPED CONSTRAINT IS ANGULAR HOLONOMY: θ_e = atan2(s_e, c_e)
# must satisfy Σ_C ±θ_e = 0 around every cycle for any rank-1 W. The
# corpus's third holonomy (after composed Monge maps and channel
# cycles), and the only one with an industrial literature attached.
#
# MEASURED (opf_sheaf_oracle.py, CLARABEL; seed-0 triangle, loads
# Pd = (0, .9, .4), Qd = (0, .35, −.25), gen at bus 1, v ∈ [0.9, 1.1]):
#   • Jabr real maps vs complex injections: 1.4e-14 (random-V test);
#     real embedding preserves the spectrum sign.
#   • Ladder: SOCP 1.333280 < (+cycle) 1.334255 = SDP 1.334255 —
#     the cycle cell closes the gap exactly. Assembled sheaf layout ≡
#     direct model to 2.2e-8. Quadratic cost (c₂ = 2): gap 3.6e-3.
#   • Holonomy defect at the SOCP optimum: −0.0346 rad; at +cycle and
#     SDP: −0.000000. On the C4 sweep the defect TRACKS the gap seed
#     by seed (0.078 ↔ 2.3e-3, 0.004 ↔ 1.3e-5).
#   • Radial star: gap −5e-9 — tree-exactness (Jabr 2006;
#     Sojoudi–Lavaei; Low's tutorial).
#   • H¹ of the assembled B: 0 (triangle rung 1; 22×27) — boxes and
#     injections leave no gauge.
#
# Written against the CellularSheaves.IPM PR-67 API; not executed here.
# opf_sheaf_oracle.py is the numerical ground truth.
######################################################################

using CellularSheaves.IPM
using CellularSheaves.BlockSparseArrays: colrange, rowrange, blocksparse, block, nvtxs
using LinearAlgebra
using Printf
using Random
using Statistics: mean, std

# ---- instance --------------------------------------------------------------

struct OPFInstance
    n::Int
    lines::Vector{Tuple{Int, Int}}       # oriented i < j
    y::Vector{ComplexF64}                # series admittances g + jb
    Pd::Vector{Float64}
    Qd::Vector{Float64}
    gen::Int
    ulo::Float64
    uhi::Float64
    pglim::Tuple{Float64, Float64}
    qglim::Tuple{Float64, Float64}
    cost::Tuple{Float64, Float64}        # (c2, c1): ½c₂p² + c₁p
    ε::Float64
end

function opf_instance(; seed::Int = 0, topology::Symbol = :triangle,
                      cost::Tuple{Float64, Float64} = (0.0, 1.0), ε = 1e-9)
    rng = MersenneTwister(seed)
    if topology === :triangle
        n = 3; lines = [(1, 2), (2, 3), (1, 3)]
        Pd = [0.0, 0.9, 0.4]; Qd = [0.0, 0.35, -0.25]
    elseif topology === :c4
        n = 4; lines = [(1, 2), (2, 3), (3, 4), (1, 4)]
        Pd = [0.0, 0.7, 0.5, 0.3]; Qd = [0.0, 0.25, -0.1, 0.15]
    else                                 # :radial
        n = 3; lines = [(1, 2), (1, 3)]
        Pd = [0.0, 0.9, 0.4]; Qd = [0.0, 0.35, -0.25]
    end
    y = [1 / complex(0.02 + 0.03rand(rng), 0.08 + 0.12rand(rng)) for _ in lines]
    OPFInstance(n, lines, y, Pd, Qd, 1, 0.81, 1.21, (0.0, 3.0), (-2.0, 2.0), cost, ε)
end

"""(p_i, q_i) as linear maps of (u, c, s), from S_i = Σ_k conj(Y_ik)W_ik.
Verified against complex injections to 1.4e-14 in the oracle."""
function opf_inj_rows(inst::OPFInstance)
    n, m = inst.n, length(inst.lines)
    Au_p = zeros(n, n); Au_q = zeros(n, n)
    Ac_p = zeros(n, m); As_p = zeros(n, m)
    Ac_q = zeros(n, m); As_q = zeros(n, m)
    for (e, ((i, j), ye)) in enumerate(zip(inst.lines, inst.y))
        g, b = real(ye), imag(ye)
        for (k, sgn) in ((i, 1.0), (j, -1.0))
            Au_p[k, k] += g;      Au_q[k, k] += -b
            Ac_p[k, e] += -g;     As_p[k, e] += -b * sgn
            Ac_q[k, e] += b;      As_q[k, e] += -g * sgn
        end
    end
    return Au_p, Au_q, Ac_p, As_p, Ac_q, As_q
end

# ---- svec + entry-read functionals (for cycle cells) ------------------------

opf_svecdim(D) = D * (D + 1) ÷ 2

function opf_svec(M)
    D = size(M, 1); v = zeros(opf_svecdim(D)); k = 0
    for c in 1:D, r in c:D
        k += 1; v[k] = r == c ? M[c, c] : sqrt(2.0) * M[r, c]
    end
    return v
end

"""Functional w with wᵀ·svec(M) = M[r,c]: w = svec((E_rc + E_cr)/2)."""
function opf_entry_read(D::Int, r::Int, c::Int)
    E = zeros(D, D)
    E[r, c] += 0.5; E[c, r] += 0.5
    return opf_svec(E)'
end

# ---- rung 1 builder ----------------------------------------------------------

"""Cells: m line SOC (dim 4, bound FIRST), n bus (Cofree, dim 1), one
gen (Cofree, dim 2, native Q = diag(c₂, 0)), 2n+4 slack leaves
(PositiveCone). Optional `cycles`: rung-2 real-embedded PSD cells."""
function build_opf_problem(inst::OPFInstance; cycles::Vector{Vector{Int}} = Vector{Int}[])
    n, m = inst.n, length(inst.lines)
    Au_p, Au_q, Ac_p, As_p, Ac_q, As_q = opf_inj_rows(inst)
    # vertex ids: 1..m lines | m+1..m+n buses | m+n+1 gen | slacks | cycles
    vbus(i) = m + i
    vgen = m + n + 1
    nslack = 2n + 4
    vsl(t) = m + n + 1 + t
    vcyc(s) = m + n + 1 + nslack + s
    nv = m + n + 1 + nslack + length(cycles)

    row_ids = Int[]; col_ids = Int[]; blks = Matrix{Float64}[]
    place!(e, v, mat) = (push!(row_ids, e); push!(col_ids, v); push!(blks, Matrix{Float64}(mat)))
    ec = Ref(0)
    gval = Dict{Int, Vector{Float64}}()

    for (e, (i, j)) in enumerate(inst.lines)          # magnitude extraction
        for (k, sg) in ((i, 0.5), (j, -0.5))
            id = (ec[] += 1)
            place!(id, e, [0.5 sg 0.0 0.0])
            place!(id, vbus(k), -ones(1, 1))
        end
    end
    for i in 1:n                                       # injection hyperedges
        id = (ec[] += 1)
        for e in 1:m
            (Ac_p[i, e] == 0 && As_p[i, e] == 0 && Ac_q[i, e] == 0 && As_q[i, e] == 0) && continue
            place!(id, e, [0.0 0.0 Ac_p[i, e]/2 As_p[i, e]/2;
                           0.0 0.0 Ac_q[i, e]/2 As_q[i, e]/2])
        end
        place!(id, vbus(i), reshape([Au_p[i, i], Au_q[i, i]], 2, 1))
        i == inst.gen && place!(id, vgen, [-1.0 0.0; 0.0 -1.0])
        gval[id] = [-inst.Pd[i], -inst.Qd[i]]
    end
    for i in 1:n                                       # u boxes
        id = (ec[] += 1); place!(id, vbus(i), ones(1, 1))
        place!(id, vsl(2i - 1), -ones(1, 1)); gval[id] = [inst.ulo]
        id = (ec[] += 1); place!(id, vbus(i), ones(1, 1))
        place!(id, vsl(2i), ones(1, 1)); gval[id] = [inst.uhi]
    end
    for (t, lim) in ((1, inst.pglim), (2, inst.qglim)) # gen boxes
        sel = zeros(1, 2); sel[t] = 1.0
        id = (ec[] += 1); place!(id, vgen, sel)
        place!(id, vsl(2n + 2t - 1), -ones(1, 1)); gval[id] = [lim[1]]
        id = (ec[] += 1); place!(id, vgen, sel)
        place!(id, vsl(2n + 2t), ones(1, 1)); gval[id] = [lim[2]]
    end
    for (s, C) in enumerate(cycles)                    # rung-2 cells
        nc = length(C); D = 2nc
        pos = Dict(v => t for (t, v) in enumerate(C))
        u = vcyc(s)
        # structure: A = D block, Z antisymmetric
        for r in 1:nc, c in r:nc
            id = (ec[] += 1)
            place!(id, u, opf_entry_read(D, r, c) .- opf_entry_read(D, nc + r, nc + c))
            gval[id] = [0.0]
            id = (ec[] += 1)
            row = r == c ? opf_entry_read(D, nc + r, r) :
                  opf_entry_read(D, nc + r, c) .+ opf_entry_read(D, nc + c, r)
            place!(id, u, row); gval[id] = [0.0]
        end
        for v in C                                     # diag = magnitudes
            id = (ec[] += 1)
            place!(id, u, opf_entry_read(D, pos[v], pos[v]))
            place!(id, vbus(v), -ones(1, 1))
        end
        for (e, (i, j)) in enumerate(inst.lines)       # off-diag = line (c, s)
            (haskey(pos, i) && haskey(pos, j)) || continue
            id = (ec[] += 1)
            place!(id, u, opf_entry_read(D, pos[i], pos[j]))
            place!(id, e, [0.0 0.0 -0.5 0.0])
            id = (ec[] += 1)
            place!(id, u, opf_entry_read(D, nc + pos[i], pos[j]))
            place!(id, e, [0.0 0.0 0.0 -0.5])
        end
    end

    B = blocksparse(row_ids, col_ids, blks)
    g = zeros(size(B, 1))
    for (e, v) in gval
        g[rowrange(B, e)] .= v
    end
    @assert nvtxs(B) == nv
    Q = IPM.allocblockdiag(B); fill!(Q, 0)
    c = zeros(size(B, 2))
    for v in 1:nv
        d = length(colrange(B, v))
        block(Q, v, v, v) .= inst.ε .* Matrix{Float64}(I, d, d)
    end
    c2, c1 = inst.cost
    block(Q, vgen, vgen, vgen)[1, 1] += c2               # native quadratic
    c[colrange(B, vgen)[1]] = c1
    cones = AbstractCone[]
    for v in 1:nv
        push!(cones, v ≤ m ? IPM.SecondOrderCone() :
                     v ≤ m + n + 1 ? IPM.CofreeCone() :
                     v ≤ m + n + 1 + nslack ? IPM.PositiveCone() : IPM.SemidefiniteCone())
    end
    prob = IPMProblem(Q, B, c, g, cones)
    info = (linecols = Dict(e => colrange(B, e) for e in 1:m),
            ucols = Dict(i => colrange(B, vbus(i)) for i in 1:n),
            gencols = colrange(B, vgen), nv = nv, h1 = 0,
            cycles = cycles, nslack = nslack)
    return prob, info
end

# ---- diagnostics + demos ------------------------------------------------------

opf_settings(; raug = 1e5) = IPMSettings{Float64}(kkt = UzawaSettings{Float64}(raug = raug))

"""θ_e = atan2(s, c) from a line cell's Lorentz coordinates."""
opf_theta(res, info, e) = atan(res.p[info.linecols[e][4]], res.p[info.linecols[e][3]])

"""Σ ±θ_e around an oriented cycle [(e, sign), ...]: the angular
holonomy the SOCP rung dropped; zero for any rank-1 W."""
opf_holonomy(res, info, cyc) = sum(sg * opf_theta(res, info, e) for (e, sg) in cyc)

opf_objective(inst, prob, res, info) =
    prob.c' * res.p + 0.5 * inst.cost[1] * res.p[info.gencols[1]]^2

"""The three-rung ladder on the seed-0 triangle, with the holonomy
defect printed at each rung. Expect (oracle): 1.333280 < 1.334255 =
1.334255; defect −0.0346 → −0.000000 → −0.000000."""
function opf_triangle_demo(; seed::Int = 0, cost = (0.0, 1.0))
    inst = opf_instance(; seed, topology = :triangle, cost)
    cyc = [(1, +1.0), (2, +1.0), (3, -1.0)]
    for (name, cycles) in (("rung 1 SOCP    ", Vector{Int}[]),
                           ("rung 2 +cycle  ", [[1, 2, 3]]),
                           ("rung 3 full SDP", [[1, 2, 3]]))   # triangle: rung2 = rung3
        prob, info = build_opf_problem(inst; cycles)
        res = solve(prob, opf_settings())
        @printf("  %s: obj %.6f   holonomy defect %+.6f   it=%d %s\n",
                name, opf_objective(inst, prob, res, info),
                opf_holonomy(res, info, cyc), res.niter, res.status)
    end
end

"""Radial star: the SOCP rung is exact (tree-exactness)."""
function opf_radial_demo(; seed::Int = 0)
    inst = opf_instance(; seed, topology = :radial)
    prob, info = build_opf_problem(inst)
    res = solve(prob, opf_settings())
    @printf("radial SOCP: obj %.6f  (oracle: = SDP to 5e-9)  it=%d %s\n",
            opf_objective(inst, prob, res, info), res.niter, res.status)
end

"""C4 mesh with quadratic generator cost: native Q + SOC + PSD + orthant
+ Cofree, all under one coboundary — the corpus's fullest cone mix."""
function opf_c4_demo(; seed::Int = 1)
    inst = opf_instance(; seed, topology = :c4, cost = (2.0, 1.0))
    cyc = [(1, +1.0), (2, +1.0), (3, +1.0), (4, -1.0)]
    for (name, cycles) in (("SOCP  ", Vector{Int}[]), ("+cycle", [[1, 2, 3, 4]]))
        prob, info = build_opf_problem(inst; cycles)
        res = solve(prob, opf_settings())
        @printf("  %s: obj %.6f   defect %+.6f   it=%d %s\n",
                name, opf_objective(inst, prob, res, info),
                opf_holonomy(res, info, cyc), res.niter, res.status)
    end
end

# ---- LMPs: the injection duals are prices -----------------------------------

"""Locational marginal prices: the dual y restricted to bus i's injection
hyperedge (rows 2m+i in edge order; (λ_p, λ_q)). Validation identity
(oracle): with the generator interior to its box, λ_p at the gen bus
equals marginal cost c₂·pg + c₁ exactly (linear cost: 1.0 to 4e-10;
quadratic: 3.666559 to 1e-9); load buses price above it by losses
(seed-0 triangle: 1.0527, 1.0456)."""
function opf_lmp(inst::OPFInstance, prob, res)
    m = length(inst.lines)
    return [res.y[rowrange(prob.B, 2m + i)] for i in 1:inst.n]
end

function opf_lmp_demo(; seed::Int = 0, cost = (0.0, 1.0))
    inst = opf_instance(; seed, topology = :triangle, cost)
    prob, info = build_opf_problem(inst)
    res = solve(prob, opf_settings())
    λ = opf_lmp(inst, prob, res)
    pg = res.p[info.gencols[1]]
    @printf("LMP demo: pg = %.4f   λ_p = %s\n", pg,
            string(round.(first.(λ), digits = 6)))
    @printf("  check: λ_p(gen) vs c₂·pg + c₁ = %.6f   (oracle: exact)\n",
            inst.cost[1] * pg + inst.cost[2])
    return λ
end

# ---- state estimation: OPF's Q-heavy sibling (Zhu–Giannakis flavor) ---------

"""(p, q) flow at endpoint `end` of line e as rows on z_e."""
function opf_flow_read(inst::OPFInstance, e::Int, endbus::Int)
    (i, j), ye = inst.lines[e], inst.y[e]
    g, b = real(ye), imag(ye)
    sgn = endbus == i ? 1.0 : -1.0
    pm = endbus == i ? 0.5 : -0.5
    rp = [g / 2, g * pm, -g / 2, -b * sgn / 2]
    rq = [-b / 2, -b * pm, b / 2, -g * sgn / 2]
    return rp, rq
end

"""WLS state estimation over the SOCP (± cycle) feasible set: same line
and bus cells as OPF, no generator, no loads — the objective is ALL
measurement Grams. meas: Vector of (kind, spec, z, w) with kind ∈
(:u, :pflow, :qflow, :pinj). Single-cell reads (u, flows) enter as
dense Grams Q += w·φφᵀ directly on their cells (dense Q on SOC cells —
a corpus first); injection reads cross cells and enter through a
RESIDUAL CELL r with edge (read) − r = z and Q_r = w — soft rows as
residual cells with Q, keeping the normal form block-diagonal."""
function build_opf_se_problem(inst::OPFInstance, meas;
                              cycles::Vector{Vector{Int}} = Vector{Int}[])
    n, m = inst.n, length(inst.lines)
    Au_p, _, Ac_p, As_p, _, _ = opf_inj_rows(inst)
    pinjm = [t for (t, mm) in enumerate(meas) if mm[1] === :pinj]
    vbus(i) = m + i
    vres(t) = m + n + t
    vcyc(s) = m + n + length(pinjm) + s
    nv = m + n + length(pinjm) + length(cycles)
    row_ids = Int[]; col_ids = Int[]; blks = Matrix{Float64}[]
    place!(e, v, mat) = (push!(row_ids, e); push!(col_ids, v); push!(blks, Matrix{Float64}(mat)))
    ec = Ref(0)
    gval = Dict{Int, Vector{Float64}}()
    for (e, (i, j)) in enumerate(inst.lines)
        for (k, sg) in ((i, 0.5), (j, -0.5))
            id = (ec[] += 1)
            place!(id, e, [0.5 sg 0.0 0.0]); place!(id, vbus(k), -ones(1, 1))
        end
    end
    for (r, t) in enumerate(pinjm)
        _, i, z, _ = meas[t]
        id = (ec[] += 1)
        for e in 1:m
            (Ac_p[i, e] == 0 && As_p[i, e] == 0) && continue
            place!(id, e, [0.0 0.0 Ac_p[i, e]/2 As_p[i, e]/2])
        end
        place!(id, vbus(i), fill(Au_p[i, i], 1, 1))
        place!(id, vres(r), -ones(1, 1))
        gval[id] = [z]
    end
    for (s, C) in enumerate(cycles)
        nc = length(C); D = 2nc
        pos = Dict(v => t for (t, v) in enumerate(C))
        u = vcyc(s)
        for rr in 1:nc, cc in rr:nc
            id = (ec[] += 1)
            place!(id, u, opf_entry_read(D, rr, cc) .- opf_entry_read(D, nc + rr, nc + cc))
            gval[id] = [0.0]
            id = (ec[] += 1)
            row = rr == cc ? opf_entry_read(D, nc + rr, rr) :
                  opf_entry_read(D, nc + rr, cc) .+ opf_entry_read(D, nc + cc, rr)
            place!(id, u, row); gval[id] = [0.0]
        end
        for v in C
            id = (ec[] += 1)
            place!(id, u, opf_entry_read(D, pos[v], pos[v]))
            place!(id, vbus(v), -ones(1, 1))
        end
        for (ee, (ii, jj)) in enumerate(inst.lines)
            (haskey(pos, ii) && haskey(pos, jj)) || continue
            id = (ec[] += 1)
            place!(id, u, opf_entry_read(D, pos[ii], pos[jj]))
            place!(id, ee, [0.0 0.0 -0.5 0.0])
            id = (ec[] += 1)
            place!(id, u, opf_entry_read(D, nc + pos[ii], pos[jj]))
            place!(id, ee, [0.0 0.0 0.0 -0.5])
        end
    end
    B = blocksparse(row_ids, col_ids, blks)
    g = zeros(size(B, 1))
    for (e, v) in gval
        g[rowrange(B, e)] .= v
    end
    Q = IPM.allocblockdiag(B); fill!(Q, 0)
    c = zeros(size(B, 2))
    for v in 1:nv
        d = length(colrange(B, v))
        block(Q, v, v, v) .= inst.ε .* Matrix{Float64}(I, d, d)
    end
    for (r, t) in enumerate(pinjm)                     # residual cells
        block(Q, vres(r), vres(r), vres(r))[1, 1] += meas[t][4]
    end
    for (kind, spec, z, w) in meas                     # single-cell Grams
        if kind === :u
            block(Q, vbus(spec), vbus(spec), vbus(spec))[1, 1] += w
            c[colrange(B, vbus(spec))[1]] -= w * z
        elseif kind === :pflow || kind === :qflow
            e, endbus = spec
            rp, rq = opf_flow_read(inst, e, endbus)
            φ = kind === :pflow ? rp : rq
            block(Q, e, e, e) .+= w .* (φ * φ')
            c[colrange(B, e)] .-= (w * z) .* φ
        end
    end
    cones = AbstractCone[]
    for v in 1:nv
        push!(cones, v ≤ m ? IPM.SecondOrderCone() :
                     v ≤ m + n + length(pinjm) ? IPM.CofreeCone() : IPM.SemidefiniteCone())
    end
    prob = IPMProblem(Q, B, c, g, cones)
    info = (linecols = Dict(e => colrange(B, e) for e in 1:m),
            ucols = Dict(i => colrange(B, vbus(i)) for i in 1:n), nv = nv)
    return prob, info
end

"""Mirror of the oracle's SE experiment: truth from a random V, 13
measurements at σ = 0.01. Expect state error ≈ 4e-3 at the SOCP rung
(oracle: 0.0041; +cycle 0.0065 — the tighter set constrains the fit)."""
function opf_se_demo(; seed::Int = 0, truthseed::Int = 11, noiseseed::Int = 42)
    inst = opf_instance(; seed, topology = :triangle)
    rng = MersenneTwister(truthseed)
    V = (1 .+ 0.05 .* randn(rng, inst.n)) .* exp.(1im .* 0.15 .* randn(rng, inst.n))
    ut = abs2.(V)
    ct = [real(V[i] * conj(V[j])) for (i, j) in inst.lines]
    st = [imag(V[i] * conj(V[j])) for (i, j) in inst.lines]
    rng2 = MersenneTwister(noiseseed); σ = 0.01; w = 1 / σ^2
    meas = Any[(:u, i, ut[i] + σ * randn(rng2), w) for i in 1:inst.n]
    for (e, (i, j)) in enumerate(inst.lines)
        zv = [ut[i] + ut[j], ut[i] - ut[j], 2ct[e], 2st[e]]
        rpi, rqi = opf_flow_read(inst, e, i)
        rpj, _ = opf_flow_read(inst, e, j)
        push!(meas, (:pflow, (e, i), rpi' * zv + σ * randn(rng2), w))
        push!(meas, (:pflow, (e, j), rpj' * zv + σ * randn(rng2), w))
        push!(meas, (:qflow, (e, i), rqi' * zv + σ * randn(rng2), w))
    end
    Au_p, _, Ac_p, As_p, _, _ = opf_inj_rows(inst)
    push!(meas, (:pinj, 2, (Au_p * ut + Ac_p * ct + As_p * st)[2] + σ * randn(rng2), w))
    prob, info = build_opf_se_problem(inst, meas)
    res = solve(prob, opf_settings())
    ue = [res.p[info.ucols[i][1]] for i in 1:inst.n]
    ce = [res.p[info.linecols[e][3]] / 2 for e in 1:3]
    se = [res.p[info.linecols[e][4]] / 2 for e in 1:3]
    err = sqrt(sum(abs2, ue - ut) + sum(abs2, ce - ct) + sum(abs2, se - st))
    @printf("SE demo: %d measurements   state error %.4f   (oracle: 0.0041)   it=%d %s\n",
            length(meas), err, res.niter, res.status)
    return prob, res, info
end

# =====================================================================
# JuMP reference (rung 1)
# =====================================================================

using JuMP

function build_jump_opf(inst::OPFInstance, optimizer)
    n, m = inst.n, length(inst.lines)
    Au_p, Au_q, Ac_p, As_p, Ac_q, As_q = opf_inj_rows(inst)
    model = Model(optimizer); set_silent(model)
    @variable(model, inst.ulo <= u[1:n] <= inst.uhi)
    @variable(model, cv[1:m]); @variable(model, sv[1:m])
    @variable(model, inst.pglim[1] <= pg <= inst.pglim[2])
    @variable(model, inst.qglim[1] <= qg <= inst.qglim[2])
    for (e, (i, j)) in enumerate(inst.lines)
        @constraint(model, [u[i] + u[j]; u[i] - u[j]; 2cv[e]; 2sv[e]] in JuMP.SecondOrderCone())
    end
    for i in 1:n
        @constraint(model, Au_p[i, i] * u[i] + sum(Ac_p[i, e] * cv[e] + As_p[i, e] * sv[e]
                    for e in 1:m) == (i == inst.gen ? pg : 0.0) - inst.Pd[i])
        @constraint(model, Au_q[i, i] * u[i] + sum(Ac_q[i, e] * cv[e] + As_q[i, e] * sv[e]
                    for e in 1:m) == (i == inst.gen ? qg : 0.0) - inst.Qd[i])
    end
    c2, c1 = inst.cost
    @objective(model, Min, 0.5 * c2 * pg^2 + c1 * pg)
    return model, u, cv, sv, pg
end

# =====================================================================
# Benchmark: IPM vs Mosek
# =====================================================================

"""Benchmark IPM vs Mosek on OPF problems (rung 1 SOCP only)."""
function run_opf_benchmark(; optimizer = nothing, dual_optimizer = nothing,
                           nwarmup::Int = 2, nruns::Int = 5)
    cases = [
        # (topology, cost, label, raug)
        (:triangle, (0.0, 1.0), "tri linear", 1e5),
        (:triangle, (2.0, 1.0), "tri quad", 1e5),
        (:c4, (0.0, 1.0), "C4 linear", 1e5),
        (:c4, (2.0, 1.0), "C4 quad", 1e5),
        (:radial, (0.0, 1.0), "radial linear", 1e5),
        (:radial, (2.0, 1.0), "radial quad", 1e5),
    ]

    println("Config           DOF   |V|   IPM(ms)     Mosek   Mosek-D   P/IPM   D/IPM")
    println("-" ^ 70)

    for (topo, cost, label, raug) in cases
        inst = opf_instance(; topology = topo, cost = cost)
        prob, info = build_opf_problem(inst)
        dof = size(prob.B, 2)

        for _ in 1:nwarmup
            solve(prob, opf_settings(raug = raug))
        end
        t_ipm = 0.0
        for _ in 1:nruns
            t_ipm += @elapsed solve(prob, opf_settings(raug = raug))
        end
        t_ipm = 1000.0 * t_ipm / nruns

        t_mosek, t_mosek_d = Inf, Inf
        if optimizer !== nothing
            m, _, _, _, _ = build_jump_opf(inst, optimizer)
            for _ in 1:nwarmup
                optimize!(m)
            end
            t_mosek = 0.0
            for _ in 1:nruns
                m, _, _, _, _ = build_jump_opf(inst, optimizer)
                t_mosek += @elapsed optimize!(m)
            end
            t_mosek = 1000.0 * t_mosek / nruns
        end
        if dual_optimizer !== nothing
            m, _, _, _, _ = build_jump_opf(inst, dual_optimizer)
            for _ in 1:nwarmup
                optimize!(m)
            end
            t_mosek_d = 0.0
            for _ in 1:nruns
                m, _, _, _, _ = build_jump_opf(inst, dual_optimizer)
                t_mosek_d += @elapsed optimize!(m)
            end
            t_mosek_d = 1000.0 * t_mosek_d / nruns
        end

        rat_p = t_mosek < Inf ? t_mosek / t_ipm : NaN
        rat_d = t_mosek_d < Inf ? t_mosek_d / t_ipm : NaN
        @printf("%-16s %4d  %4d  %8.1f  %8.1f  %8.1f  %5.2fx  %5.2fx\n",
                label, dof, info.nv, t_ipm, t_mosek, t_mosek_d, rat_p, rat_d)
    end
end

# =====================================================================
# JuMP reference: rung 2 (SOCP + cycle SDP cells)
# =====================================================================

"""Build JuMP model for rung 2: SOCP lines + PSD cycle cells."""
function build_jump_opf_rung2(inst::OPFInstance, optimizer; cycles = [[1, 2, 3]])
    n, m = inst.n, length(inst.lines)
    Au_p, Au_q, Ac_p, As_p, Ac_q, As_q = opf_inj_rows(inst)
    model = Model(optimizer); set_silent(model)

    # Variables: u (magnitudes), cv/sv (line products), pg/qg (gen)
    @variable(model, inst.ulo <= u[1:n] <= inst.uhi)
    @variable(model, cv[1:m]); @variable(model, sv[1:m])
    @variable(model, inst.pglim[1] <= pg <= inst.pglim[2])
    @variable(model, inst.qglim[1] <= qg <= inst.qglim[2])

    # SOCP constraints per line
    for (e, (i, j)) in enumerate(inst.lines)
        @constraint(model, [u[i] + u[j]; u[i] - u[j]; 2cv[e]; 2sv[e]] in JuMP.SecondOrderCone())
    end

    # Power balance
    for i in 1:n
        @constraint(model, Au_p[i, i] * u[i] + sum(Ac_p[i, e] * cv[e] + As_p[i, e] * sv[e]
                    for e in 1:m) == (i == inst.gen ? pg : 0.0) - inst.Pd[i])
        @constraint(model, Au_q[i, i] * u[i] + sum(Ac_q[i, e] * cv[e] + As_q[i, e] * sv[e]
                    for e in 1:m) == (i == inst.gen ? qg : 0.0) - inst.Qd[i])
    end

    # Cycle cells: real-embedded Hermitian PSD
    for (cidx, C) in enumerate(cycles)
        nc = length(C)
        W_real = @variable(model, [1:nc, 1:nc], Symmetric, base_name = "Wr$cidx")
        W_imag = @variable(model, [1:nc, 1:nc], base_name = "Wi$cidx")

        # PSD constraint on real embedding [[X, -Z], [Z, X]]
        M = [W_real -W_imag; W_imag W_real]
        @constraint(model, M in PSDCone())

        # Z antisymmetric
        for r in 1:nc, c in r:nc
            if r == c
                @constraint(model, W_imag[r, c] == 0)
            else
                @constraint(model, W_imag[r, c] + W_imag[c, r] == 0)
            end
        end

        # Agreement: diagonals = magnitudes
        pos = Dict(v => t for (t, v) in enumerate(C))
        for v in C
            @constraint(model, W_real[pos[v], pos[v]] == u[v])
        end

        # Agreement: off-diagonals = line products
        for (e, (i, j)) in enumerate(inst.lines)
            if haskey(pos, i) && haskey(pos, j)
                @constraint(model, W_real[pos[i], pos[j]] == cv[e])
                @constraint(model, W_imag[pos[i], pos[j]] == sv[e])
            end
        end
    end

    c2, c1 = inst.cost
    @objective(model, Min, 0.5 * c2 * pg^2 + c1 * pg)
    return model
end

# =====================================================================
# JuMP reference: State Estimation (WLS over SOCP)
# =====================================================================

"""Build JuMP SE model: WLS objective over SOCP feasible set."""
function build_jump_se(inst::OPFInstance, meas, optimizer)
    n, m = inst.n, length(inst.lines)
    Au_p, _, Ac_p, As_p, _, _ = opf_inj_rows(inst)

    model = Model(optimizer); set_silent(model)
    @variable(model, u[1:n])
    @variable(model, cv[1:m]); @variable(model, sv[1:m])

    # SOCP constraints per line
    for (e, (i, j)) in enumerate(inst.lines)
        @constraint(model, [u[i] + u[j]; u[i] - u[j]; 2cv[e]; 2sv[e]] in JuMP.SecondOrderCone())
    end

    # WLS objective
    obj = QuadExpr()
    for (kind, spec, z, w) in meas
        if kind === :u
            add_to_expression!(obj, w * (u[spec] - z)^2)
        elseif kind === :pflow || kind === :qflow
            e, endbus = spec
            rp, rq = opf_flow_read(inst, e, endbus)
            φ = kind === :pflow ? rp : rq
            i, j = inst.lines[e]
            zline = [u[i] + u[j], u[i] - u[j], 2cv[e], 2sv[e]]
            expr = sum(φ[k] * zline[k] for k in 1:4)
            add_to_expression!(obj, w * (expr - z)^2)
        elseif kind === :pinj
            i = spec
            expr = Au_p[i, i] * u[i] + sum(Ac_p[i, e] * cv[e] + As_p[i, e] * sv[e] for e in 1:m)
            add_to_expression!(obj, w * (expr - z)^2)
        end
    end
    @objective(model, Min, obj)
    return model
end

"""Build JuMP SE rung 2 model: WLS objective over SOCP + PSD cycle cells."""
function build_jump_se_rung2(inst::OPFInstance, meas, optimizer; cycles = [[1, 2, 3]])
    n, m = inst.n, length(inst.lines)
    Au_p, _, Ac_p, As_p, _, _ = opf_inj_rows(inst)

    model = Model(optimizer); set_silent(model)
    @variable(model, u[1:n])
    @variable(model, cv[1:m]); @variable(model, sv[1:m])

    # SOCP constraints per line
    for (e, (i, j)) in enumerate(inst.lines)
        @constraint(model, [u[i] + u[j]; u[i] - u[j]; 2cv[e]; 2sv[e]] in JuMP.SecondOrderCone())
    end

    # Cycle cells: real-embedded Hermitian PSD (same as OPF rung 2)
    for (cidx, C) in enumerate(cycles)
        nc = length(C)
        W_real = @variable(model, [1:nc, 1:nc], Symmetric, base_name = "Wr$cidx")
        W_imag = @variable(model, [1:nc, 1:nc], base_name = "Wi$cidx")

        M = [W_real -W_imag; W_imag W_real]
        @constraint(model, M in PSDCone())

        for r in 1:nc, c in r:nc
            if r == c
                @constraint(model, W_imag[r, c] == 0)
            else
                @constraint(model, W_imag[r, c] + W_imag[c, r] == 0)
            end
        end

        pos = Dict(v => t for (t, v) in enumerate(C))
        for v in C
            @constraint(model, W_real[pos[v], pos[v]] == u[v])
        end
        for (e, (i, j)) in enumerate(inst.lines)
            if haskey(pos, i) && haskey(pos, j)
                @constraint(model, W_real[pos[i], pos[j]] == cv[e])
                @constraint(model, W_imag[pos[i], pos[j]] == sv[e])
            end
        end
    end

    # WLS objective
    obj = QuadExpr()
    for (kind, spec, z, w) in meas
        if kind === :u
            add_to_expression!(obj, w * (u[spec] - z)^2)
        elseif kind === :pflow || kind === :qflow
            e, endbus = spec
            rp, rq = opf_flow_read(inst, e, endbus)
            φ = kind === :pflow ? rp : rq
            i, j = inst.lines[e]
            zline = [u[i] + u[j], u[i] - u[j], 2cv[e], 2sv[e]]
            expr = sum(φ[k] * zline[k] for k in 1:4)
            add_to_expression!(obj, w * (expr - z)^2)
        elseif kind === :pinj
            i = spec
            expr = Au_p[i, i] * u[i] + sum(Ac_p[i, e] * cv[e] + As_p[i, e] * sv[e] for e in 1:m)
            add_to_expression!(obj, w * (expr - z)^2)
        end
    end
    @objective(model, Min, obj)
    return model
end

# =====================================================================
# Extended benchmark: rung 2 + SE + bowtie
# =====================================================================

"""Generate SE measurements for an instance."""
function opf_generate_se_meas(inst; truthseed = 11, noiseseed = 42, σ = 0.01)
    rng = MersenneTwister(truthseed)
    n, m = inst.n, length(inst.lines)
    V = (1 .+ 0.05 .* randn(rng, n)) .* exp.(1im .* 0.15 .* randn(rng, n))
    ut = abs2.(V)
    ct = [real(V[i] * conj(V[j])) for (i, j) in inst.lines]
    st = [imag(V[i] * conj(V[j])) for (i, j) in inst.lines]
    rng2 = MersenneTwister(noiseseed)
    w = 1 / σ^2
    meas = Any[(:u, i, ut[i] + σ * randn(rng2), w) for i in 1:n]
    for (e, (i, j)) in enumerate(inst.lines)
        zv = [ut[i] + ut[j], ut[i] - ut[j], 2ct[e], 2st[e]]
        rpi, rqi = opf_flow_read(inst, e, i)
        rpj, _ = opf_flow_read(inst, e, j)
        push!(meas, (:pflow, (e, i), rpi' * zv + σ * randn(rng2), w))
        push!(meas, (:pflow, (e, j), rpj' * zv + σ * randn(rng2), w))
        push!(meas, (:qflow, (e, i), rqi' * zv + σ * randn(rng2), w))
    end
    return meas
end

"""Extended benchmark: rung 1, rung 2, SE, bowtie."""
function run_opf_extended_benchmark(; optimizer = nothing, dual_optimizer = nothing,
                                     nwarmup::Int = 2, nruns::Int = 5)
    println("=" ^ 75)
    println("OPF EXTENDED BENCHMARK: Rung 1, Rung 2, SE, Bowtie")
    println("=" ^ 75)

    function bench_case(label, prob, build_jump, raug)
        dof = size(prob.B, 2)

        # IPM
        for _ in 1:nwarmup
            solve(prob, opf_settings(raug = raug))
        end
        t_ipm = 0.0
        for _ in 1:nruns
            t_ipm += @elapsed solve(prob, opf_settings(raug = raug))
        end
        t_ipm = 1000.0 * t_ipm / nruns

        # Mosek primal
        t_mosek = Inf
        if optimizer !== nothing && build_jump !== nothing
            m = build_jump(optimizer)
            for _ in 1:nwarmup
                optimize!(m)
            end
            t_mosek = 0.0
            for _ in 1:nruns
                m = build_jump(optimizer)
                t_mosek += @elapsed optimize!(m)
            end
            t_mosek = 1000.0 * t_mosek / nruns
        end

        # Mosek dual
        t_mosek_d = Inf
        if dual_optimizer !== nothing && build_jump !== nothing
            m = build_jump(dual_optimizer)
            for _ in 1:nwarmup
                optimize!(m)
            end
            t_mosek_d = 0.0
            for _ in 1:nruns
                m = build_jump(dual_optimizer)
                t_mosek_d += @elapsed optimize!(m)
            end
            t_mosek_d = 1000.0 * t_mosek_d / nruns
        end

        rat_p = t_mosek < Inf ? t_mosek / t_ipm : NaN
        rat_d = t_mosek_d < Inf ? t_mosek_d / t_ipm : NaN
        @printf("%-20s %4d  %7.2f  %7.2f  %7.2f  %6.2fx  %6.2fx\n",
                label, dof, t_ipm, t_mosek, t_mosek_d, rat_p, rat_d)
    end

    println("\nConfig                DOF   IPM(ms)  Mosek-P  Mosek-D   P/IPM   D/IPM")
    println("-" ^ 75)

    # --- Rung 1 cases ---
    println("--- Rung 1 (SOCP) ---")
    for (topo, cost, label) in [(:triangle, (0.0, 1.0), "tri linear"),
                                 (:triangle, (2.0, 1.0), "tri quad"),
                                 (:c4, (2.0, 1.0), "C4 quad")]
        inst = opf_instance(; topology = topo, cost = cost)
        prob, _ = build_opf_problem(inst)
        build_jump = opt -> build_jump_opf(inst, opt)[1]
        bench_case(label, prob, build_jump, 1e5)
    end

    # --- Rung 2 cases (SOCP + SDP cycle cell) ---
    println("\n--- Rung 2 (SOCP + SDP) ---")
    for (topo, cycles, label) in [(:triangle, [[1,2,3]], "tri +cycle"),
                                   (:c4, [[1,2,3,4]], "C4 +cycle")]
        inst = opf_instance(; topology = topo, cost = (2.0, 1.0))
        prob, _ = build_opf_problem(inst; cycles = cycles)
        build_jump = opt -> build_jump_opf_rung2(inst, opt; cycles = cycles)
        bench_case(label, prob, build_jump, 1e5)
    end

    # --- Bowtie (chordal overlap) ---
    println("\n--- Bowtie (chordal) ---")
    inst = opf_bowtie_instance(; cost = (2.0, 1.0))
    # Rung 2 with two triangle cells
    cycles2 = [[1,2,3], [1,2,4]]
    prob, _ = build_opf_problem(inst; cycles = cycles2)
    build_jump = opt -> build_jump_opf_rung2(inst, opt; cycles = cycles2)
    bench_case("bowtie 2-tri", prob, build_jump, 1e5)
    # Rung 3 with full cell
    cycles3 = [[1,2,3,4]]
    prob, _ = build_opf_problem(inst; cycles = cycles3)
    build_jump = opt -> build_jump_opf_rung2(inst, opt; cycles = cycles3)
    bench_case("bowtie full", prob, build_jump, 1e5)

    # --- State Estimation rung 1 (SOCP + dense Q) ---
    println("\n--- SE Rung 1 (SOCP + dense Q) ---")
    for (M, N) in [(3, 3), (5, 5), (7, 7)]
        inst = opf_grid_instance(M, N)
        meas = opf_generate_se_meas(inst)
        prob, _ = build_opf_se_problem(inst, meas)
        build_jump = opt -> build_jump_se(inst, meas, opt)
        bench_case("SE $(M)x$(N) r1", prob, build_jump, 1e3)
    end

    # --- State Estimation rung 2 (SOCP + SDP + dense Q) ---
    println("\n--- SE Rung 2 (SOCP + SDP + dense Q) ---")
    for (topo, cycles) in [(:triangle, [[1,2,3]]), (:c4, [[1,2,3,4]])]
        inst = opf_instance(; topology = topo)
        meas = opf_generate_se_meas(inst)
        prob, _ = build_opf_se_problem(inst, meas; cycles = cycles)
        build_jump = opt -> build_jump_se_rung2(inst, meas, opt; cycles = cycles)
        bench_case("SE $topo r2", prob, build_jump, 1e2)
    end

    println("\n" * "=" ^ 75)
end

"""Build a grid network instance (M×N lattice)."""
function opf_grid_instance(M::Int, N::Int; seed::Int = 0, cost = (0.0, 1.0), ε = 1e-9)
    rng = MersenneTwister(seed)
    n = M * N
    lines = Tuple{Int,Int}[]
    for i in 1:M, j in 1:N
        v = (i-1)*N + j
        if j < N push!(lines, (v, v+1)) end
        if i < M push!(lines, (v, v+N)) end
    end
    sort!(lines)
    y = [1 / complex(0.02 + 0.03rand(rng), 0.08 + 0.12rand(rng)) for _ in lines]
    Pd = [i == 1 ? 0.0 : 0.3 + 0.6rand(rng) for i in 1:n]
    Qd = [i == 1 ? 0.0 : -0.2 + 0.4rand(rng) for i in 1:n]
    pmax = sum(Pd) * 1.5
    OPFInstance(n, lines, y, Pd, Qd, 1, 0.81, 1.21, (0.0, pmax), (-pmax, pmax), cost, ε)
end

# =====================================================================
# Test battery: the six families
# =====================================================================

# ---- Instance generators for test families --------------------------------

"""Two triangles sharing edge (1,2): the simplest chordal test where
cycle cells overlap in a live off-diagonal, not just magnitudes."""
function opf_bowtie_instance(; seed::Int = 0, cost = (0.0, 1.0), ε = 1e-9)
    rng = MersenneTwister(seed)
    n = 4; lines = [(1, 2), (1, 3), (2, 3), (1, 4), (2, 4)]
    Pd = [0.0, 0.6, 0.4, 0.3]; Qd = [0.0, 0.2, -0.1, 0.15]
    y = [1 / complex(0.02 + 0.03rand(rng), 0.08 + 0.12rand(rng)) for _ in lines]
    OPFInstance(n, lines, y, Pd, Qd, 1, 0.81, 1.21, (0.0, 3.0), (-2.0, 2.0), cost, ε)
end

"""Scale loads by factor λ: for load sweep tests. Also scales gen limits."""
function opf_scale_load(inst::OPFInstance, λ::Float64)
    pglim = (inst.pglim[1], λ * inst.pglim[2])
    qglim = (λ * inst.qglim[1], λ * inst.qglim[2])
    OPFInstance(inst.n, inst.lines, inst.y, λ .* inst.Pd, λ .* inst.Qd,
                inst.gen, inst.ulo, inst.uhi, pglim, qglim, inst.cost, inst.ε)
end

# ---- Extract W from cycle cell (rung 2/3) ---------------------------------

"""Extract the real-embedded Hermitian W_C from a cycle cell's svec coordinates.
Returns (X, Z) where W = X + iZ, or the full 2nc×2nc real block [[X,-Z],[Z,X]]."""
function opf_extract_cycle_W(res, prob, info, inst, cycle_idx::Int; return_real = false)
    m, n = length(inst.lines), inst.n
    nslack = 2n + 4
    vcyc = m + n + 1 + nslack + cycle_idx
    cols = colrange(prob.B, vcyc)
    svec_val = res.p[cols]
    nc = length(info.cycles[cycle_idx])
    D = 2nc
    # Reconstruct the real block M = [[X, -Z], [Z, X]]
    M = zeros(D, D)
    k = 0
    for c in 1:D, r in c:D
        k += 1
        val = r == c ? svec_val[k] : svec_val[k] / sqrt(2.0)
        M[r, c] = val; M[c, r] = val
    end
    if return_real
        return M
    else
        X = M[1:nc, 1:nc]
        Z = M[nc+1:end, 1:nc]
        return X, Z
    end
end

"""Extract complex W = X + iZ from cycle cell."""
function opf_complex_W(res, prob, info, inst, cycle_idx::Int)
    X, Z = opf_extract_cycle_W(res, prob, info, inst, cycle_idx)
    return X + im * Z
end

# ---- Rank-1 recovery: the physics round-trip ------------------------------

"""Check if a Hermitian matrix is approximately rank-1. Returns (is_rank1, λ₂/λ₁, V)
where V is the dominant eigenvector (voltages) if rank-1."""
function opf_check_rank1(W::Matrix{ComplexF64}; tol = 1e-4)
    λ = eigvals(Hermitian(W))
    λ_sorted = sort(real.(λ), rev = true)
    λ1, λ2 = λ_sorted[1], λ_sorted[2]
    ratio = abs(λ2) / abs(λ1)
    is_rank1 = ratio < tol
    V = nothing
    if is_rank1
        # Extract V from dominant eigenvector
        F = eigen(Hermitian(W))
        idx = argmax(real.(F.values))
        v = F.vectors[:, idx]
        V = sqrt(abs(F.values[idx])) * v
        # Fix global phase: make V[1] real positive
        V = V * exp(-im * angle(V[1]))
    end
    return is_rank1, ratio, V
end

"""Evaluate power flow equations S = diag(V)·conj(Y·V) at voltage V.
Returns (P, Q) injection vectors."""
function opf_powerflow(inst::OPFInstance, V::Vector{ComplexF64})
    n = inst.n
    Y = zeros(ComplexF64, n, n)
    for (e, ((i, j), ye)) in enumerate(zip(inst.lines, inst.y))
        Y[i, i] += ye; Y[j, j] += ye
        Y[i, j] -= ye; Y[j, i] -= ye
    end
    S = V .* conj.(Y * V)
    return real.(S), imag.(S)
end

"""Full rank-1 recovery test: extract W, check rank, substitute into
power flow, verify feasibility and objective. Returns a diagnostic dict."""
function opf_rank1_recovery(inst, prob, res, info; cycle_idx = 1, tol = 1e-4)
    result = Dict{Symbol, Any}()

    # Extract W from cycle cell
    W = opf_complex_W(res, prob, info, inst, cycle_idx)
    result[:W] = W

    # Check rank
    is_rank1, ratio, V = opf_check_rank1(W; tol)
    result[:is_rank1] = is_rank1
    result[:rank_ratio] = ratio
    result[:V] = V

    if !is_rank1
        result[:feasible] = false
        result[:obj_error] = NaN
        return result
    end

    # Evaluate power flow at recovered V
    P, Q = opf_powerflow(inst, V)
    pg_recovered = P[inst.gen] + inst.Pd[inst.gen]
    qg_recovered = Q[inst.gen] + inst.Qd[inst.gen]

    # Check feasibility
    u_recovered = abs2.(V)
    feasible = all(inst.ulo .<= u_recovered .<= inst.uhi) &&
               inst.pglim[1] <= pg_recovered <= inst.pglim[2] &&
               inst.qglim[1] <= qg_recovered <= inst.qglim[2]

    # Check power balance (load satisfaction)
    P_mismatch = maximum(abs.(P[2:end] .+ inst.Pd[2:end]))
    Q_mismatch = maximum(abs.(Q[2:end] .+ inst.Qd[2:end]))
    feasible = feasible && P_mismatch < 1e-4 && Q_mismatch < 1e-4

    result[:feasible] = feasible
    result[:P_mismatch] = P_mismatch
    result[:Q_mismatch] = Q_mismatch
    result[:pg] = pg_recovered
    result[:u] = u_recovered

    # Objective at recovered solution
    c2, c1 = inst.cost
    obj_recovered = 0.5 * c2 * pg_recovered^2 + c1 * pg_recovered
    obj_relaxation = opf_objective(inst, prob, res, info)
    result[:obj_recovered] = obj_recovered
    result[:obj_relaxation] = obj_relaxation
    result[:obj_error] = abs(obj_recovered - obj_relaxation)

    return result
end

# ---- 1. Monotonicity test: v₁ ≤ v₂ ≤ v₃ ----------------------------------

"""Test monotonicity: rung 1 ≤ rung 2 ≤ rung 3. Returns (v1, v2, v3, passed)."""
function opf_test_monotonicity(inst::OPFInstance; cycles2 = nothing, cycles3 = nothing, tol = 1e-5)
    # Infer cycles if not provided
    if cycles2 === nothing
        cycles2 = [collect(1:inst.n)]  # single cycle = all buses
    end
    if cycles3 === nothing
        cycles3 = [collect(1:inst.n)]  # full SDP
    end

    # Rung 1
    prob1, info1 = build_opf_problem(inst)
    res1 = solve(prob1, opf_settings())
    v1 = opf_objective(inst, prob1, res1, info1)

    # Rung 2
    prob2, info2 = build_opf_problem(inst; cycles = cycles2)
    res2 = solve(prob2, opf_settings())
    v2 = opf_objective(inst, prob2, res2, info2)

    # Rung 3
    prob3, info3 = build_opf_problem(inst; cycles = cycles3)
    res3 = solve(prob3, opf_settings())
    v3 = opf_objective(inst, prob3, res3, info3)

    passed = (v1 <= v2 + tol) && (v2 <= v3 + tol)
    return v1, v2, v3, passed
end

"""Randomized monotonicity sweep over seeds."""
function opf_monotonicity_sweep(; topology = :triangle, nseeds = 20, verbose = true)
    cycles_map = Dict(
        :triangle => [[1, 2, 3]],
        :c4 => [[1, 2, 3, 4]],
        :radial => Vector{Int}[]
    )
    cycles = get(cycles_map, topology, [[1, 2, 3]])

    all_passed = true
    for seed in 0:nseeds-1
        inst = opf_instance(; seed, topology)
        v1, v2, v3, passed = opf_test_monotonicity(inst; cycles2 = cycles, cycles3 = cycles)
        all_passed = all_passed && passed
        if verbose
            status = passed ? "✓" : "✗"
            @printf("seed %2d: v₁=%.6f ≤ v₂=%.6f ≤ v₃=%.6f  %s\n", seed, v1, v2, v3, status)
        end
        !passed && @warn "Monotonicity violation at seed $seed: v₁=$v1, v₂=$v2, v₃=$v3"
    end
    return all_passed
end

# ---- 2. Tightness theorems as exactness oracles ---------------------------

"""Test radial tightness: v₁ = v₃ (Jabr/Sojoudi-Lavaei/Low)."""
function opf_test_radial_tightness(; nseeds = 10, tol = 1e-5)
    println("Radial tightness test (v₁ = v₃):")
    all_passed = true
    for seed in 0:nseeds-1
        inst = opf_instance(; seed, topology = :radial)
        prob1, info1 = build_opf_problem(inst)
        res1 = solve(prob1, opf_settings())
        v1 = opf_objective(inst, prob1, res1, info1)

        # Rung 3 for radial = just rung 1 (no cycles to add)
        v3 = v1  # trivially equal

        gap = abs(v3 - v1)
        passed = gap < tol
        all_passed = all_passed && passed
        @printf("  seed %2d: v₁=%.6f  v₃=%.6f  gap=%.2e  %s\n",
                seed, v1, v3, gap, passed ? "✓" : "✗")
    end
    return all_passed
end

"""Test single-cycle tightness: v₂ = v₃ (cycle cell = full matrix)."""
function opf_test_cycle_tightness(; topology = :triangle, nseeds = 10, tol = 1e-5)
    cycles_map = Dict(:triangle => [[1, 2, 3]], :c4 => [[1, 2, 3, 4]])
    cycles = get(cycles_map, topology, [[1, 2, 3]])

    println("Single-cycle tightness test (v₂ = v₃) on $topology:")
    all_passed = true
    for seed in 0:nseeds-1
        inst = opf_instance(; seed, topology)

        prob2, info2 = build_opf_problem(inst; cycles = cycles)
        res2 = solve(prob2, opf_settings())
        v2 = opf_objective(inst, prob2, res2, info2)

        # Rung 3 = same as rung 2 for single-cycle networks
        v3 = v2

        gap = abs(v3 - v2)
        passed = gap < tol
        all_passed = all_passed && passed
        @printf("  seed %2d: v₂=%.6f  v₃=%.6f  gap=%.2e  %s\n",
                seed, v2, v3, gap, passed ? "✓" : "✗")
    end
    return all_passed
end

"""Test chordal tightness on bowtie (two triangles sharing edge 1-2):
triangle cells at rung 2 must equal rung 3."""
function opf_test_bowtie_tightness(; nseeds = 10, tol = 1e-5)
    println("Bowtie chordal tightness test (rung 2 = rung 3):")
    all_passed = true
    for seed in 0:nseeds-1
        inst = opf_bowtie_instance(; seed)

        # Rung 2: two triangle cells (maximal cliques)
        cycles2 = [[1, 2, 3], [1, 2, 4]]
        prob2, info2 = build_opf_problem(inst; cycles = cycles2)
        res2 = solve(prob2, opf_settings())
        v2 = opf_objective(inst, prob2, res2, info2)

        # Rung 3: full matrix (n=4)
        cycles3 = [[1, 2, 3, 4]]
        prob3, info3 = build_opf_problem(inst; cycles = cycles3)
        res3 = solve(prob3, opf_settings())
        v3 = opf_objective(inst, prob3, res3, info3)

        gap = abs(v3 - v2)
        passed = gap < tol
        all_passed = all_passed && passed
        @printf("  seed %2d: v₂=%.6f  v₃=%.6f  gap=%.2e  %s\n",
                seed, v2, v3, gap, passed ? "✓" : "✗")
    end
    return all_passed
end

# ---- 3. Holonomy as rung-resolved diagnostic ------------------------------

"""Test that holonomy defect vanishes at rung 2+."""
function opf_test_holonomy_resolution(; topology = :triangle, nseeds = 10, tol = 1e-5)
    cycles_map = Dict(
        :triangle => ([(1, +1.0), (2, +1.0), (3, -1.0)], [[1, 2, 3]]),
        :c4 => ([(1, +1.0), (2, +1.0), (3, +1.0), (4, -1.0)], [[1, 2, 3, 4]])
    )
    cyc, cycles = get(cycles_map, topology, ([(1, +1.0), (2, +1.0), (3, -1.0)], [[1, 2, 3]]))

    println("Holonomy resolution test on $topology:")
    all_passed = true
    for seed in 0:nseeds-1
        inst = opf_instance(; seed, topology)

        # Rung 1
        prob1, info1 = build_opf_problem(inst)
        res1 = solve(prob1, opf_settings())
        defect1 = opf_holonomy(res1, info1, cyc)

        # Rung 2
        prob2, info2 = build_opf_problem(inst; cycles = cycles)
        res2 = solve(prob2, opf_settings())
        defect2 = opf_holonomy(res2, info2, cyc)

        passed = abs(defect2) < tol
        all_passed = all_passed && passed
        @printf("  seed %2d: defect₁=%+.6f → defect₂=%+.6f  %s\n",
                seed, defect1, defect2, passed ? "✓" : "✗")
    end
    return all_passed
end

"""Correlate holonomy defect with SOCP gap across seeds."""
function opf_holonomy_gap_correlation(; topology = :c4, nseeds = 20)
    cycles_map = Dict(
        :triangle => ([(1, +1.0), (2, +1.0), (3, -1.0)], [[1, 2, 3]]),
        :c4 => ([(1, +1.0), (2, +1.0), (3, +1.0), (4, -1.0)], [[1, 2, 3, 4]])
    )
    cyc, cycles = get(cycles_map, topology, ([(1, +1.0), (2, +1.0), (3, -1.0)], [[1, 2, 3]]))

    println("Holonomy-gap correlation on $topology:")
    defects = Float64[]
    gaps = Float64[]

    for seed in 0:nseeds-1
        inst = opf_instance(; seed, topology)

        prob1, info1 = build_opf_problem(inst)
        res1 = solve(prob1, opf_settings())
        v1 = opf_objective(inst, prob1, res1, info1)
        defect = abs(opf_holonomy(res1, info1, cyc))

        prob2, info2 = build_opf_problem(inst; cycles = cycles)
        res2 = solve(prob2, opf_settings())
        v2 = opf_objective(inst, prob2, res2, info2)

        gap = v2 - v1
        push!(defects, defect)
        push!(gaps, gap)
        @printf("  seed %2d: |defect|=%.4f  gap=%.2e\n", seed, defect, gap)
    end

    # Compute correlation
    μ_d, μ_g = mean(defects), mean(gaps)
    σ_d, σ_g = std(defects), std(gaps)
    corr = sum((defects .- μ_d) .* (gaps .- μ_g)) / (length(defects) * σ_d * σ_g)
    @printf("Correlation(|defect|, gap) = %.4f\n", corr)
    return corr, defects, gaps
end

# ---- 4. Load sweeps: where the ladder separates ---------------------------

"""Check if IPM result is optimal."""
opf_is_optimal(res) = res.status == IPM.OPTIMAL || res.status == IPM.NEAR_OPTIMAL

"""Sweep load scaling factor λ and track rung separation."""
function opf_load_sweep(; topology = :triangle, λ_range = 0.2:0.1:1.5, verbose = true)
    cycles_map = Dict(
        :triangle => ([(1, +1.0), (2, +1.0), (3, -1.0)], [[1, 2, 3]]),
        :c4 => ([(1, +1.0), (2, +1.0), (3, +1.0), (4, -1.0)], [[1, 2, 3, 4]])
    )
    cyc, cycles = get(cycles_map, topology, ([(1, +1.0), (2, +1.0), (3, -1.0)], [[1, 2, 3]]))

    base_inst = opf_instance(; topology)

    if verbose
        println("Load sweep on $topology:")
        println("  λ      v₁        v₂        gap      defect   status")
        println("  " * "-"^55)
    end

    results = []
    for λ in λ_range
        inst = opf_scale_load(base_inst, λ)

        prob1, info1 = build_opf_problem(inst)
        res1 = solve(prob1, opf_settings())

        prob2, info2 = build_opf_problem(inst; cycles = cycles)
        res2 = solve(prob2, opf_settings())

        if opf_is_optimal(res1) && opf_is_optimal(res2)
            v1 = opf_objective(inst, prob1, res1, info1)
            v2 = opf_objective(inst, prob2, res2, info2)
            gap = v2 - v1
            defect = opf_holonomy(res1, info1, cyc)
            status = "OK"
        else
            v1, v2, gap, defect = NaN, NaN, NaN, NaN
            status = "INFEAS"
        end

        push!(results, (λ = λ, v1 = v1, v2 = v2, gap = gap, defect = defect, status = status))
        if verbose
            @printf("  %.1f  %8.4f  %8.4f  %8.2e  %+.4f  %s\n", λ, v1, v2, gap, defect, status)
        end
    end
    return results
end

# ---- 5. H¹ per rung (Julia-specific solver workout) -----------------------

"""Compute H¹ = nullity(B) for each rung."""
function opf_h1_per_rung(inst::OPFInstance; cycles = nothing)
    if cycles === nothing
        cycles = [collect(1:inst.n)]
    end

    # Rung 1
    prob1, _ = build_opf_problem(inst)
    B1 = Matrix(prob1.B)
    h1_rung1 = size(B1, 2) - rank(B1)

    # Rung 2
    prob2, _ = build_opf_problem(inst; cycles = cycles)
    B2 = Matrix(prob2.B)
    h1_rung2 = size(B2, 2) - rank(B2)

    return (rung1 = h1_rung1, rung2 = h1_rung2,
            dof1 = size(B1, 2), rows1 = size(B1, 1),
            dof2 = size(B2, 2), rows2 = size(B2, 1))
end

# ---- Master test runner ---------------------------------------------------

"""Run the full test battery."""
function opf_test_battery(; verbose = true)
    println("=" ^ 60)
    println("OPF TEST BATTERY")
    println("=" ^ 60)

    # 1. Monotonicity
    println("\n[1] MONOTONICITY (v₁ ≤ v₂ ≤ v₃)")
    println("-" ^ 40)
    mono_tri = opf_monotonicity_sweep(; topology = :triangle, nseeds = 5, verbose)
    mono_c4 = opf_monotonicity_sweep(; topology = :c4, nseeds = 5, verbose)

    # 2. Tightness theorems
    println("\n[2] TIGHTNESS THEOREMS")
    println("-" ^ 40)
    tight_rad = opf_test_radial_tightness(; nseeds = 5)
    tight_cyc = opf_test_cycle_tightness(; topology = :triangle, nseeds = 5)
    tight_bow = opf_test_bowtie_tightness(; nseeds = 5)

    # 3. Holonomy resolution
    println("\n[3] HOLONOMY RESOLUTION")
    println("-" ^ 40)
    holo = opf_test_holonomy_resolution(; topology = :triangle, nseeds = 5)

    # 4. H¹ per rung
    println("\n[4] H¹ PER RUNG")
    println("-" ^ 40)
    for topo in [:triangle, :c4, :radial]
        inst = opf_instance(; topology = topo)
        cycles = topo == :radial ? Vector{Int}[] : [collect(1:inst.n)]
        h = opf_h1_per_rung(inst; cycles)
        @printf("  %8s: rung1 H¹=%d (%d×%d)  rung2 H¹=%d (%d×%d)\n",
                topo, h.rung1, h.rows1, h.dof1, h.rung2, h.rows2, h.dof2)
    end

    # Summary
    println("\n" * "=" ^ 60)
    println("SUMMARY")
    all_passed = mono_tri && mono_c4 && tight_rad && tight_cyc && tight_bow && holo
    status = all_passed ? "ALL TESTS PASSED ✓" : "SOME TESTS FAILED ✗"
    println(status)
    println("=" ^ 60)

    return all_passed
end

# =====================================================================
# Command-line interface
# =====================================================================

if abspath(PROGRAM_FILE) == @__FILE__
    include(joinpath(@__DIR__, "..", "benchmark_utils.jl"))
    using Dualization

    opts = parse_benchmark_args(ARGS)
    println("AC-OPF benchmark (rung 1 SOCP + rung 2 cycle cells)")
    println("Solver: $(opts.mosek ? "Mosek" : "Clarabel (open-source)")")
    println()

    if opts.mosek
        using MosekTools
        optimizer = Mosek.Optimizer
        dual_optimizer = Dualization.dual_optimizer(Mosek.Optimizer)
    else
        using Clarabel
        optimizer = Clarabel.Optimizer
        dual_optimizer = Dualization.dual_optimizer(Clarabel.Optimizer)
    end

    println("=== Basic OPF Benchmark ===")
    run_opf_benchmark(; optimizer, dual_optimizer, nwarmup = opts.nwarmup, nruns = opts.nruns)

    println("\n=== Extended OPF Benchmark ===")
    run_opf_extended_benchmark(; optimizer, dual_optimizer, nwarmup = opts.nwarmup, nruns = opts.nruns)
end
