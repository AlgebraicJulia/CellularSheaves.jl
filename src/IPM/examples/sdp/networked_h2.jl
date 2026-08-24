# =============================================================================
#
# # H₂ Analysis of Networked Systems via Chordal Decomposition
#
# A network of n linear subsystems interacts over a graph G = (V, E).
# The dynamics of subsystem i are
#
#   ẋᵢ = Aᵢᵢ xᵢ + ∑ Aᵢⱼ xⱼ + Bᵢ wᵢ
#                j ∈ Nᵢ                                          (1)
#   yᵢ = Cᵢ xᵢ
#
# with local state xᵢ ∈ Xᵢ, disturbance wᵢ ∈ Wᵢ, and output yᵢ ∈ Yᵢ.
# The squared H₂ norm of the system is
#
#   ‖ C (sI - A)⁻¹ B ‖²  = min  tr(Bᵀ P B)                       (2)
#                     H₂    P
#                          s.t. Aᵀ P + P A + Cᵀ C ≤ 0
#                                               P ≥ 0
#
# It measures how much the system amplifies persistent random noise.
# The semidefinite constraints in (2) couple every subsystem, so for
# large networks the problem is intractable. A tractable problem may
# be produced by restricting P to be block diagonal
#
#   P = [ P₁      ]
#       [    ⋱    ]
#       [      Pₙ ]
#
# so that the matrix S = Aᵀ P + P A + Cᵀ C inherits the sparsity pattern
# of E ∪ Eᵀ. The generalized Agler theorem allows us to replace the 
# constraint S ≤ 0 with a collection of smaller ones.
#
#         o
#   S = - ∑ Eₖᵀ Zₖ Eₖ,  Zₖ ≥ 0, k = 1, …, o
#        k=1
#
# This yields the restricted problem
#
#        n
#   min  ∑ tr(Bᵢᵀ Pᵢ Bᵢ)                                         (3)
#       i=1
#
#                            o
#   s.t. Aᵀ P + P A + Cᵀ C + ∑ Eₖᵀ Zₖ Eₖ = -ε I
#                           k=1
#                                     Pᵢ ≥ δ I,  i = 1, …, n
#                                     Zₖ ≥ 0,    k = 1, …, o
#
# for some small, positive robustness margins ε and δ. The minimum of
# (3) is an upper bound to the squared H₂ norm (2).
#
# # References
#
#   - Zheng, Kamgarpour, Sootla & Papachristodoulou,
#     Scalable analysis of linear networked systems via chordal decomposition,
#     ECC 2018. arXiv:1803.05996.
#
# =============================================================================

using LinearAlgebra, Printf, Random
using SparseArrays
using Dualization

using CellularSheaves.IPM
using CellularSheaves.IPM: SemidefiniteCone, IPMSettings, HSDSettings
using CellularSheaves.BlockSparseArrays: colrange, blocksparse

include("../utils.jl")

const R2 = sqrt(2.0)

# -----------------------------------------------------------------------------
# Instances
# -----------------------------------------------------------------------------

# A chain of n subsystems: subsystem i interacts with subsystems i ± 1, so the
# graph is a path and its maximal cliques are Cₑ = {e, e+1}, e = 1, …, n-1
# (a path is already chordal). This is the benchmark family of section V of
# the reference.
struct ChainInstance
    alpha::Vector{Int}                  # local state dimensions
    Adiag::Vector{Matrix{Float64}}      # Aᵢᵢ
    Asub::Vector{Matrix{Float64}}       # A₍ₑ₊₁₎ₑ — influence of e on e+1
    Asup::Vector{Matrix{Float64}}       # Aₑ₍ₑ₊₁₎ — influence of e+1 on e
    Bs::Vector{Matrix{Float64}}         # Bᵢ
    Cs::Vector{Matrix{Float64}}         # Cᵢ
end

nsubs(inst::ChainInstance) = length(inst.alpha)
nstate(inst::ChainInstance) = sum(inst.alpha)

# Random instance following section V of the reference: state dimensions
# αᵢ ~ U{5, …, 10}, disturbance and output dimensions mᵢ, dᵢ ~ U{1, …, 5},
# and random Gaussian blocks. The reference stabilizes A by subtracting
# (λₘₐₓ + 5)I where λₘₐₓ is the largest real part of the eigenvalues of A,
# giving a stability margin of exactly 5.
function chain_instance(n::Int; seed::Int = 0)
    @assert n >= 2
    rng = Xoshiro(seed)

    alpha = rand(rng, 5:10, n)
    m = rand(rng, 1:5, n)
    d = rand(rng, 1:5, n)

    Adiag = [randn(rng, alpha[i], alpha[i]) for i in 1:n]
    Asub = [randn(rng, alpha[e + 1], alpha[e]) for e in 1:(n - 1)]
    Asup = [randn(rng, alpha[e], alpha[e + 1]) for e in 1:(n - 1)]

    Bs = [randn(rng, alpha[i], m[i]) for i in 1:n]
    Cs = [randn(rng, d[i], alpha[i]) for i in 1:n]

    inst = ChainInstance(alpha, Adiag, Asub, Asup, Bs, Cs)

    A, _ = assemble(inst)
    mu = maximum(real, eigvals(A))

    for i in 1:n
        Adiag[i] -= (mu + 5.0) * I
    end

    return inst
end

# The assembled N × N state matrix, for the ground truth and the JuMP model.
function assemble(inst::ChainInstance)
    n = nsubs(inst)
    off = cumsum([1; inst.alpha])
    A = zeros(nstate(inst), nstate(inst))

    for i in 1:n
        A[off[i]:(off[i + 1] - 1), off[i]:(off[i + 1] - 1)] = inst.Adiag[i]
    end

    for e in 1:(n - 1)
        A[off[e + 1]:(off[e + 2] - 1), off[e]:(off[e + 1] - 1)] = inst.Asub[e]
        A[off[e]:(off[e + 1] - 1), off[e + 1]:(off[e + 2] - 1)] = inst.Asup[e]
    end

    return A, off
end

# -----------------------------------------------------------------------------
# Build
# -----------------------------------------------------------------------------

# Index of the (i, j) entry of a d × d symmetric matrix in the svec layout used
# by SemidefiniteCone: the lower triangle stored column by column, the diagonal
# unscaled and the off-diagonal scaled by √2.
function svecidx(d::Int, i::Int, j::Int)
    i, j = max(i, j), min(i, j)
    return (j - 1) * d - ((j - 1) * (j - 2)) ÷ 2 + (i - j + 1)
end

svecw(d::Int) = d * (d + 1) ÷ 2

# Factor mapping the svec coordinate (a, b) to the matrix entry: the svec
# coordinate stores √2 Mₐᵦ off the diagonal.
svecsc(a::Int, b::Int) = a == b ? 1.0 : 1.0 / R2

# The svec vector of a symmetric matrix M, satisfying ⟨M, X⟩ = svec(M)ᵀ svec(X).
function svecm(M::AbstractMatrix)
    d = size(M, 1)
    v = zeros(svecw(d))

    for b in 1:d, a in b:d
        v[svecidx(d, a, b)] = M[a, b] * (a == b ? 1.0 : R2)
    end

    return v
end

# The stalks are numbered
#
#   1 … n               subsystem   svec P̃ᵢ    semidefinite, width αᵢ(αᵢ+1)/2
#   n + 1 … n + (n-1)   clique      svec Zₑ    semidefinite, width sₑ(sₑ+1)/2
#
# with sₑ = αₑ + αₑ₊₁, and within clique stalk e the aggregated state is
# ordered (xₑ, xₑ₊₁), so subsystem e occupies rows 1:αₑ and subsystem e+1
# occupies rows αₑ+1:sₑ.
#
# The equality constraints run over the blocks of the chordal pattern: one
# block row per diagonal block (i, i) — the lower triangle, αᵢ(αᵢ+1)/2
# scalar rows — and one per subdiagonal block (e+1, e) — αₑ₊₁ αₑ scalar
# rows. Scalar rows are ordered column-major within each block.
function build_h2(inst::ChainInstance; eps::Float64 = 1.0e-6, delta::Float64 = 1.0e-6)
    n = nsubs(inst)
    alpha = inst.alpha

    P(i) = i
    Z(e) = n + e

    row_ids = Int[]
    col_ids = Int[]
    blks = Matrix{Float64}[]
    g = Float64[]

    r = 0

    # -- diagonal blocks (i, i):
    #
    #      A₍ᵢᵢ₎ᵀ P̃ᵢ + P̃ᵢ A₍ᵢᵢ₎ + ∑ (Zₑ)ᵢᵢ = -(Cᵢᵀ Cᵢ + εI + δ(Aᵢᵢᵀ + Aᵢᵢ))
    #                             e ∋ i
    for i in 1:n
        ai = alpha[i]
        wi = svecw(ai)
        A = inst.Adiag[i]

        r += 1
        M = zeros(wi, wi)

        for b in 1:ai, a in b:ai
            row = svecidx(ai, a, b)

            for c in 1:ai
                M[row, svecidx(ai, a, c)] += A[c, b] * svecsc(a, c)
                M[row, svecidx(ai, c, b)] += A[c, a] * svecsc(c, b)
            end
        end

        push!(row_ids, r); push!(col_ids, P(i)); push!(blks, M)

        for e in (i - 1, i)
            1 <= e <= n - 1 || continue
            s = alpha[e] + alpha[e + 1]
            q = i == e ? 0 : alpha[e]           # position of subsystem i in clique e

            Zb = zeros(wi, svecw(s))

            for b in 1:ai, a in b:ai
                Zb[svecidx(ai, a, b), svecidx(s, q + a, q + b)] = svecsc(q + a, q + b)
            end

            push!(row_ids, r); push!(col_ids, Z(e)); push!(blks, Zb)
        end

        Gd = -(inst.Cs[i]' * inst.Cs[i] + eps * I + delta * (A' + A))

        for b in 1:ai, a in b:ai
            push!(g, Gd[a, b])
        end
    end

    # -- subdiagonal blocks (i, j) = (e+1, e):
    #
    #      A₍ⱼᵢ₎ᵀ P̃ⱼ + P̃ᵢ A₍ᵢⱼ₎ + (Zₑ)ᵢⱼ = -δ(A₍ⱼᵢ₎ᵀ + A₍ᵢⱼ₎)
    for e in 1:(n - 1)
        i, j = e + 1, e
        ai, aj = alpha[i], alpha[j]
        Aij = inst.Asub[e]                      # A₍ₑ₊₁₎ₑ
        Aji = inst.Asup[e]                      # Aₑ₍ₑ₊₁₎

        r += 1
        Mi = zeros(ai * aj, svecw(ai))
        Mj = zeros(ai * aj, svecw(aj))
        Zb = zeros(ai * aj, svecw(aj + ai))

        for b in 1:aj, a in 1:ai
            row = (b - 1) * ai + a

            for c in 1:ai
                Mi[row, svecidx(ai, a, c)] += Aij[c, b] * svecsc(a, c)
            end

            for c in 1:aj
                Mj[row, svecidx(aj, c, b)] += Aji[c, a] * svecsc(c, b)
            end

            Zb[row, svecidx(aj + ai, aj + a, b)] = svecsc(aj + a, b)
        end

        push!(row_ids, r); push!(col_ids, P(i)); push!(blks, Mi)
        push!(row_ids, r); push!(col_ids, P(j)); push!(blks, Mj)
        push!(row_ids, r); push!(col_ids, Z(e)); push!(blks, Zb)

        Gd = -delta * (Aji' + Aij)
        append!(g, vec(Gd))
    end

    B = blocksparse(row_ids, col_ids, blks)

    # Q = 0: the objective ∑ Tr(Bᵢᵀ P̃ᵢ Bᵢ) is linear
    Q = IPM.allocblockdiag(B)
    fill!(Q, 0)

    # min ½ pᵀQp - cᵀp with c = -svec(Bᵢ Bᵢᵀ) on subsystem stalk i, so that
    # -cᵀp = ∑ ⟨Bᵢ Bᵢᵀ, P̃ᵢ⟩ = ∑ Tr(Bᵢᵀ P̃ᵢ Bᵢ)
    c = zeros(size(B, 2))

    for i in 1:n
        c[colrange(B, P(i))] .= -svecm(inst.Bs[i] * inst.Bs[i]')
    end

    cones = IPM.AbstractCone[SemidefiniteCone() for _ in 1:(2n - 1)]

    return IPMProblem(Q, B, c, g, cones), B
end

# The H₂ norm recovered from the solver objective: the objective of (3) is
# Tr(Bᵀ P̃ B) = Tr(Bᵀ P B) - δ Tr(Bᵀ B).
function h2norm(pobj::Float64, inst::ChainInstance; delta::Float64 = 1.0e-6)
    return sqrt(pobj + delta * sum(sum(abs2, Bi) for Bi in inst.Bs))
end

# The dense method the reference compares against: the same block-diagonal
# Lyapunov function, but a single PSD constraint of side N = ∑ αᵢ instead of
# the clique-wise decomposition. S is assembled block by block — the AffExpr
# work is O(∑ αᵢ³), only the cone is big.
function build_h2_jump(inst::ChainInstance; eps::Float64 = 1.0e-6, delta::Float64 = 1.0e-6, optimizer)
    n = nsubs(inst)
    _, off = assemble(inst)
    N = nstate(inst)

    model = Model(optimizer)
    Pt = [@variable(model, [1:a, 1:a] in PSDCone()) for a in inst.alpha]

    S = zeros(AffExpr, N, N)
    rows(i) = off[i]:(off[i + 1] - 1)

    for i in 1:n
        Ad = inst.Adiag[i]
        S[rows(i), rows(i)] = Ad' * Pt[i] + Pt[i] * Ad +
                              inst.Cs[i]' * inst.Cs[i] + eps * I + delta * (Ad' + Ad)
    end

    for e in 1:(n - 1)
        Aij = inst.Asub[e]                      # A₍ₑ₊₁₎ₑ
        Aji = inst.Asup[e]                      # Aₑ₍ₑ₊₁₎
        Se = Aji' * Pt[e] + Pt[e + 1] * Aij + delta * (Aji' + Aij)
        S[rows(e + 1), rows(e)] = Se
        S[rows(e), rows(e + 1)] = permutedims(Se)
    end

    @constraint(model, LinearAlgebra.Symmetric(-S) in PSDCone())
    @objective(model, Min, sum(tr(inst.Bs[i]' * Pt[i] * inst.Bs[i]) for i in 1:n))

    return model
end

# -----------------------------------------------------------------------------
# Ground truth
# -----------------------------------------------------------------------------

# The exact H₂ norm via the Lyapunov equation Aᵀ P + P A + CᵀC = 0. Because
# B is block diagonal, Tr(Bᵀ P B) only reads the diagonal blocks of P.
function h2true(inst::ChainInstance)
    A, off = assemble(inst)
    n = nsubs(inst)

    CtC = zeros(nstate(inst), nstate(inst))

    for i in 1:n
        rows = off[i]:(off[i + 1] - 1)
        CtC[rows, rows] = inst.Cs[i]' * inst.Cs[i]
    end

    P = lyap(Matrix(A'), CtC)

    return sqrt(sum(tr(inst.Bs[i]' * P[off[i]:(off[i + 1] - 1), off[i]:(off[i + 1] - 1)] * inst.Bs[i]) for i in 1:n))
end

# -----------------------------------------------------------------------------
# Benchmark
# -----------------------------------------------------------------------------

# The chain lengths of section V of the reference.
const CASES = [20, 40, 80]

function benchmark(; tol = 1.0e-6, clarabel = false, mosek = false, cases = CASES)
    @printf("%-6s %6s %7s %7s %9s %9s %9s %9s %8s %8s\n",
            "n", "N", "cols", "rows",
            "IPM", "HSD", "Clarabel", "Mosek", "Cla/IPM", "Mos/IPM")

    for n in cases
        inst = chain_instance(n)
        prob, _ = build_h2(inst)

        mi = measure_ipm(prob, IPMSettings{Float64}(feas_tol = tol, gap_tol = tol))
        mh = measure_ipm(prob, HSDSettings{Float64}(feas_tol = tol, gap_tol = tol))

        if clarabel
            mc = measure_jump(() -> build_h2_jump(inst; optimizer = clarabel_opt(; tol)))
        else
            mc = (t = NaN, status = "—", obj = NaN)
        end

        if mosek
            mk = measure_jump(() -> build_h2_jump(inst; optimizer = dual_optimizer(mosek_opt(; tol))))
        else
            mk = (t = NaN, status = "—", obj = NaN)
        end

        cratio = (isfinite(mc.t) && isfinite(mi.t)) ? @sprintf("%.2f", mc.t / mi.t) : "—"
        mratio = (isfinite(mk.t) && isfinite(mi.t)) ? @sprintf("%.2f", mk.t / mi.t) : "—"

        @printf("%-6d %6d %7d %7d %s %s %s %s %8s %8s\n",
                n, nstate(inst), size(prob.B, 2), size(prob.B, 1),
                fmt_time(mi), fmt_time(mh), fmt_time(mc), fmt_time(mk), cratio, mratio)
        flush(stdout)
    end
end

# The computed H₂ upper bound against the Lyapunov ground truth: the gap is
# the conservatism of the block-diagonal Lyapunov function (cf. Table I of
# the reference, where the decomposed and dense solutions agree with each
# other but exceed the accurate value).
function accuracy(; tol = 1.0e-8, cases = [20, 40, 80])
    @printf("%-6s %9s %9s %9s\n", "n", "true", "IPM", "gap")

    for n in cases
        inst = chain_instance(n)
        prob, _ = build_h2(inst)
        res = solve(prob, IPMSettings{Float64}(feas_tol = tol, gap_tol = tol))

        hi = h2norm(res.pobj, inst)
        ht = h2true(inst)

        @printf("%-6d %9.4f %9.4f %+8.1e\n", n, ht, hi, hi - ht)
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
# Sample run: 2026-08-24 (--clarabel --mosek), tol = 1e-6; accuracy at 1e-8.
#
#   n           N    cols    rows       IPM       HSD  Clarabel     Mosek  Cla/IPM  Mos/IPM
#   20        142    2706    1548  119.0ms  166.5ms  424.2ms  388.8ms     3.56     3.27
#   40        288    5635    3234  311.8ms  355.8ms  822.1ms 1425.7ms     2.64     4.57
#   80        585   11803    6757  738.3ms 1066.3ms 2432.9ms 7961.1ms     3.30    10.78
#
#   n       true       IPM       gap
#   20    7.8201    9.0634  +1.2e+00
#   40   11.3446   13.2833  +1.9e+00
#   80   16.3517   19.4563  +3.1e+00
# =============================================================================
