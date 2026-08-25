# =============================================================================
#
# # Single-Facility Location with Mixed Norms
#
# Let
#
#   B = [ b₁, … bₘ ]
#
# be a n × m matrix, each column of which is the location of a facility.
# Let w ∈ ℝᵐ and p ∈ ℝᵐ be vectors of weights and exponents:
#
#   w = [ w₁, …, wₘ ]
#   p = [ p₁, …, pₘ ]
#
# We wish to build a new facility that is close to all the other ones.
# The single-facility location problem solves this problem by minimizing
# the sum
#
#
#   min  ∑ wᵢ ‖ x - bᵢ ‖                                             (1)
#    x  i=1             pᵢ
#
#
# Problem (1) is convex, and it can be reformulated as a power-cone
# program
#
#        m   n
#   min  ∑   ∑ wᵢ yᵢⱼ                                                (2)
#   x,y i=1 j=1
#
#   s.t. ( yᵢⱼ , ∑ yᵢₖ , xⱼ - (bᵢ)ⱼ ) ∈ ℙᵢ,   i = 1 … m,  j = 1 … n
#               k=1
#
# where
#
#    ℙᵢ = { (a, b, c) ∈ ℝ³ : a ≥ 0, b ≥ 0, aᵗ b¹⁻ᵗ ≥ |c| }
#
# is the three-dimensional power cone with parameter t = 1/pᵢ.
#
# # References
#
#   - Chares & Glineur,
#     An interior-point method for the single-facility location problem with
#     mixed norms using a conic formulation,
#     Math. Meth. Oper. Res. 68(3):383–405, 2008. CORE Discussion Paper 2007/71.
#
# =============================================================================

using LinearAlgebra, Random, Printf

using CellularSheaves.IPM
using CellularSheaves.IPM: CofreeCone, PowerCone, IPMSettings, HSDSettings
using CellularSheaves.BlockSparseArrays: colrange, block, blocksparse
using Dualization: dual_optimizer

include("../utils.jl")

# -----------------------------------------------------------------------------
# Data — the reference's generator: pᵢ ~ U[1,3], Bᵢ ~ U[0,1]ⁿ, cᵢ = 1
# -----------------------------------------------------------------------------

struct MixedNormInstance
    n::Int                      # dimension of the space
    m::Int                      # number of fixed points
    B::Matrix{Float64}          # m × n fixed points
    p::Vector{Float64}          # exponents
    c::Vector{Float64}          # weights
end

function mixed_norm_instance(n::Int, m::Int; seed::Int = 0)
    rng = Xoshiro(seed)
    B = rand(rng, m, n)
    p = 1.0 .+ 2.0 .* rand(rng, m)
    c = ones(m)
    return MixedNormInstance(n, m, B, p, c)
end

# The objective (1), evaluated directly — used to check the conic optimum.
function direct_objective(inst::MixedNormInstance, x::AbstractVector)
    total = 0.0

    for i in 1:inst.m
        s = 0.0

        for j in 1:inst.n
            s += abs(x[j] - inst.B[i, j])^inst.p[i]
        end

        total += inst.c[i] * s^(1 / inst.p[i])
    end

    return total
end

# -----------------------------------------------------------------------------
# Conic model (2)
# -----------------------------------------------------------------------------

# Stalk layout (§3). x is split into n width-1 stalks — one width-n stalk would
# force every row block touching it to store a dense h×n block for a single
# nonzero. Each yᵢⱼ lives in a width-3 PowerCone stalk (z₁,z₂,z₃) =
# (yᵢⱼ, tᵢ, xⱼ-Bᵢⱼ), and the accumulator tᵢ = Σₖ yᵢₖ is built from width-1
# CofreeCone partials. `branch` selects the accumulator:
#
#   branch = b (integer)   b-ary tree summing the yᵢₖ; b ≥ n is one fan-in row
#                          (a 3n+1 clique in BᵀB), b = 2 a binary tree
#   branch = :chain        sequential partials sᵢₖ = sᵢ,ₖ₋₁ + yᵢₖ
#
# Fan-in has fewest variables but a clique that is catastrophic at large n;
# chain has least fill but most variables; the tree interpolates. Rows: one
# height-2 block per (i,j) linking z₃-xⱼ and z₂-tᵢ, plus the accumulator rows.
function build_mixed_norm(inst::MixedNormInstance; branch = min(inst.n, 10))
    n, m, B, p, c = inst.n, inst.m, inst.B, inst.p, inst.c

    coldims = Int[]
    cones = IPM.AbstractCone[]
    alloc!(w, cone) = (push!(coldims, w); push!(cones, cone); length(coldims))

    xid = [alloc!(1, CofreeCone()) for _ in 1:n]
    zid = [alloc!(3, PowerCone(1 / p[i])) for i in 1:m, j in 1:n]

    row_ids = Int[]
    col_ids = Int[]
    blks = Matrix{Float64}[]
    gs = Float64[]
    rowdims = Int[]
    r = 0

    Z = [0.0 0.0 1.0; 0.0 1.0 0.0]     # selects (z₃, z₂)
    X = reshape([-1.0, 0.0], 2, 1)     # -xⱼ in row 1
    A = reshape([0.0, -1.0], 2, 1)     # -tᵢ in row 2
    e1 = [-1.0 0.0 0.0]                # -z₁ selector

    for i in 1:m
        # accumulator over (z₁ᵢ₁ … z₁ᵢₙ); root stalk is tᵢ
        if branch === :chain
            prev = 0
            for k in 1:n
                pid = alloc!(1, CofreeCone())
                r += 1; push!(rowdims, 1)
                push!(row_ids, r); push!(col_ids, pid);       push!(blks, ones(1, 1))
                push!(row_ids, r); push!(col_ids, zid[i, k]); push!(blks, copy(e1))
                if prev != 0
                    push!(row_ids, r); push!(col_ids, prev); push!(blks, -ones(1, 1))
                end
                push!(gs, 0.0)
                prev = pid
            end
            root = prev
        else
            terms = Tuple{Int, Matrix{Float64}}[(zid[i, k], e1) for k in 1:n]
            while length(terms) > 1
                nxt = Tuple{Int, Matrix{Float64}}[]
                for grp in Iterators.partition(terms, branch)
                    pid = alloc!(1, CofreeCone())
                    r += 1; push!(rowdims, 1)
                    push!(row_ids, r); push!(col_ids, pid); push!(blks, ones(1, 1))
                    for (tid, sel) in grp
                        push!(row_ids, r); push!(col_ids, tid); push!(blks, copy(sel))
                    end
                    push!(gs, 0.0)
                    push!(nxt, (pid, -ones(1, 1)))
                end
                terms = nxt
            end
            root = terms[1][1]
        end

        # height-2 linking block per (i, j): z₃ᵢⱼ - xⱼ = -Bᵢⱼ ; z₂ᵢⱼ - tᵢ = 0
        for j in 1:n
            r += 1; push!(rowdims, 2)
            push!(row_ids, r); push!(col_ids, zid[i, j]); push!(blks, copy(Z))
            push!(row_ids, r); push!(col_ids, xid[j]);    push!(blks, copy(X))
            push!(row_ids, r); push!(col_ids, root);      push!(blks, copy(A))
            push!(gs, -B[i, j]); push!(gs, 0.0)
        end
    end

    B_bs = blocksparse(row_ids, col_ids, blks, rowdims, coldims)

    Q = IPM.allocblockdiag(B_bs)
    fill!(Q, 0)

    # f: -cᵢ on each z₁ᵢⱼ, so that -fᵀp = ∑ᵢ cᵢ ∑ⱼ yᵢⱼ
    f = zeros(size(B_bs, 2))
    for i in 1:m, j in 1:n
        f[first(colrange(B_bs, zid[i, j]))] = -c[i]
    end

    return IPMProblem(Q, B_bs, f, gs, cones)
end

# Recover the facility position from the split width-1 x stalks.
function recover_x(res, prob, inst::MixedNormInstance)
    return [res.p[first(colrange(prob.B, j))] for j in 1:inst.n]
end

# JuMP model of (2). The accumulator tᵢ = Σₖ yᵢₖ is either one fan-in equality
# (chain = false) or a chain of partial sums (chain = true) — the direct analogue
# of `branch` in build_mixed_norm, exposed so each external solver can be given
# its faster formulation in the benchmark.
function build_mn_jump(inst::MixedNormInstance; optimizer, chain::Bool = false)
    n, m, B, p, c = inst.n, inst.m, inst.B, inst.p, inst.c

    model = Model(optimizer)
    @variable(model, x[1:n])
    @variable(model, y[1:m, 1:n] >= 0)

    if chain
        @variable(model, s[1:m, 1:n])
        for i in 1:m
            @constraint(model, s[i, 1] == y[i, 1])
            for k in 2:n
                @constraint(model, s[i, k] == s[i, k - 1] + y[i, k])
            end
            for j in 1:n
                @constraint(model, [y[i, j], s[i, n], x[j] - B[i, j]] in MOI.PowerCone(1 / p[i]))
            end
        end
        @objective(model, Min, sum(c[i] * s[i, n] for i in 1:m))
    else
        @variable(model, t[1:m])
        for i in 1:m
            @constraint(model, t[i] == sum(y[i, k] for k in 1:n))
            for j in 1:n
                @constraint(model, [y[i, j], t[i], x[j] - B[i, j]] in MOI.PowerCone(1 / p[i]))
            end
        end
        @objective(model, Min, sum(c[i] * t[i] for i in 1:m))
    end

    return model
end

# -----------------------------------------------------------------------------
# Accuracy — the conic optimum must equal the value of (1) at the recovered x
# -----------------------------------------------------------------------------

function accuracy(; tol = 1.0e-8, cases = [(2, 100), (3, 200), (10, 100)])
    @printf("%-4s %6s %14s %14s %10s\n", "n", "m", "conic", "direct", "rel.err")

    for (n, m) in cases
        inst = mixed_norm_instance(n, m)
        prob = build_mixed_norm(inst)
        res = solve(prob, IPMSettings{Float64}(feas_tol = tol, gap_tol = tol))

        x = recover_x(res, prob, inst)
        d = direct_objective(inst, x)
        err = abs(res.pobj - d) / max(1.0, abs(d))

        @printf("%-4d %6d %14.9f %14.9f %10.2e\n", n, m, res.pobj, d, err)
    end
end

# -----------------------------------------------------------------------------
# Benchmark — the (n, m) grid of Tables 1–3 of the reference
# -----------------------------------------------------------------------------

const PROBLEMS = [
    (2, 10), (2, 100), (2, 1000), (2, 10000),
    (10, 10), (10, 50), (10, 100), (10, 500), (10, 1000),
    (50, 10), (50, 50), (50, 70), (50, 100), (50, 200),
]

# Three representative sizes recorded in the sample run: the reference's hardest
# instance (small n, huge m), the mid grid, and the large-n regime (where the
# chain accumulator matters).
const SAMPLE = [(2, 10000), (10, 1000), (50, 200)]

# Each solver is given its fastest formulation, measured over the chain × fan-in
# × dualize grid (see sample run): our IPM uses fan-in for n ≤ 10 and chain for
# larger n; Clarabel the fan-in primal (chain and dualization both hurt or fail
# it — Dualization rejects the power cone); Mosek the dualized model, fan-in for
# n ≤ 10 and chain above. Cla/IPM is Clarabel's best over IPM's best.
function benchmark(; tol = 1.0e-6, clarabel = false, mosek = false, problems = PROBLEMS)
    @printf("%-4s %6s %7s %8s %9s %9s %9s %9s %8s %8s\n",
            "n", "m", "cones", "cols", "IPM", "HSD", "Clarabel", "Mosek", "Cla/IPM", "Mos/IPM")

    for (n, m) in problems
        inst = mixed_norm_instance(n, m)
        br = n <= 10 ? n : :chain
        prob = build_mixed_norm(inst; branch = br)

        mi = measure_ipm(prob, IPMSettings{Float64}(feas_tol = tol, gap_tol = tol))
        mh = measure_ipm(prob, HSDSettings{Float64}(feas_tol = tol, gap_tol = tol))

        if clarabel
            mc = measure_jump(() -> build_mn_jump(inst; optimizer = clarabel_opt(; tol), chain = false))
        else
            mc = (t = NaN, status = "—", obj = NaN)
        end

        if mosek
            mk = measure_jump(() -> build_mn_jump(inst;
                    optimizer = dual_optimizer(mosek_opt(; tol)), chain = n > 10))
        else
            mk = (t = NaN, status = "—", obj = NaN)
        end

        cratio = (isfinite(mc.t) && isfinite(mi.t)) ? @sprintf("%.2f", mc.t / mi.t) : "—"
        mratio = (isfinite(mk.t) && isfinite(mi.t)) ? @sprintf("%.2f", mk.t / mi.t) : "—"

        @printf("%-4d %6d %7d %8d %s %s %s %s %8s %8s\n",
                n, m, m * n, size(prob.B, 2),
                fmt_time(mi), fmt_time(mh), fmt_time(mc), fmt_time(mk), cratio, mratio)
        flush(stdout)
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    opts = parse_args()
    tol = any(a -> startswith(a, "--tol="), ARGS) ? opts.tol : 1.0e-6
    benchmark(; problems = SAMPLE, clarabel = opts.clarabel, mosek = opts.mosek, tol)

    println()
    accuracy()
end

# =============================================================================
# Sample run: 2026-08-24 (Apple M1 Pro, -t 8, --clarabel --mosek), tol = 1e-6,
# seed 0. Fastest formulation per solver: IPM fan-in for n ≤ 10 and chain above;
# Clarabel fan-in primal (dualization rejects the power cone); Mosek dualized,
# fan-in for n ≤ 10 and chain above.
#
#   n         m   cones     cols       IPM       HSD  Clarabel     Mosek  Cla/IPM  Mos/IPM
#   2     10000   20000    70002  440.3ms  516.9ms  439.1ms  714.1ms     1.00     1.62
#   10     1000   10000    31010  294.1ms  334.4ms  224.7ms  403.0ms     0.76     1.37
#   50      200   10000    40050  308.5ms  303.4ms  257.8ms  529.9ms     0.84     1.72
# =============================================================================
