######################################################################
# nonnegative_spline_exact.jl
#
# Shape-constrained Bézier-spline estimation — EXACT variant.
# Enforces the true pointwise shape condition on each piece via the
# Karlin–Studden representation of nonnegative univariate polynomials
# (Papp & Alizadeh 2014). The certificate is an auxiliary cone vertex
# per piece, tied to the coefficients by one equality edge — so this
# reifies (cone-on-auxiliary-vertex), unlike the orthant LP file.
#
# ONE construction, auto-selecting the cone by degree of the certified
# polynomial nd:
#     block size 1  -> PositiveCone         (trivial)
#     block size 2  -> SecondOrderCone(3)   ⟹ SOCP   (nd ≤ 3)
#     block size ≥3 -> SemidefiniteCone     ⟹ SDP    (nd ≥ 4)
# A 2×2 PSD block [[a,b],[b,c]] ⪰ 0  ⟺  (a+c, a−c, 2b) ∈ Q₃, which is why
# low-degree pieces stay second-order-cone. For min-crackle (n=9) the
# blocks are 5×5 and this is a genuine SDP.
#
#   shape = :nonneg    -> certify f      (degree nd = n)
#   shape = :monotone  -> certify f'     (degree nd = n−1, coeffs n·ΔP)
#   shape = :convex    -> certify f''    (degree nd = n−2, coeffs ∝ Δ²P)
#
# Reuses the instance generator, objective, and numerics from the LP file.
# See nonnegative_spline.md for the full derivation.
######################################################################

include("nonnegative_spline_lp.jl")

# ---- Karlin–Studden certificate machinery ---------------------------

"""Bernstein→monomial matrix on [0,1], (n+1)×(n+1): c = M·P."""
function bern2mono(n::Int)
    Mm = zeros(n + 1, n + 1)
    for ℓ in 0:n, i in 0:ℓ
        Mm[ℓ + 1, i + 1] = binomial(n, i) * binomial(n - i, ℓ - i) * (-1)^(ℓ - i)
    end
    return Mm
end

svecdim(m::Int) = m * (m + 1) ÷ 2

"""Antidiagonal-sum selector Σ_{i+j=ℓ} A_ij in the library's svec coords
(col-major LOWER triangle, diagonal first per column, off-diag ×√2)."""
function antidiag(m::Int, ℓ::Int)
    s = zeros(svecdim(m)); p = 0
    for j in 0:(m - 1)
        p += 1
        (2j == ℓ) && (s[p] = 1.0)
        for i in (j + 1):(m - 1)
            p += 1
            (i + j == ℓ) && (s[p] = sqrt(2))
        end
    end
    return s
end

"""
Cone + coeffmap for a certificate block of matrix-size m.
Returns (C, cone, conedim) where C maps the cone stalk to svec coords.
"""
function block_cone(m::Int)
    if m == 1
        return Matrix{Float64}(I, 1, 1), IPM.PositiveCone(), 1
    elseif m == 2
        # aux stalk z=(a+c, a−c, 2b) ∈ Q₃ (bound first);  svec₂=(a, √2 b, c)
        C = [0.5 0.5 0.0; 0.0 0.0 sqrt(2)/2; 0.5 -0.5 0.0]
        return C, IPM.SecondOrderCone(), 3
    else
        d = svecdim(m)
        # NOTE: confirm your SemidefiniteCone's stalk layout matches this svec
        # (col-major upper triangle, off-diagonals scaled by √2). Adjust here if not.
        return Matrix{Float64}(I, d, d), IPM.SemidefiniteCone(), d
    end
end

"""Contribution map E (rows ℓ=0..nd) of a block to the monomial coefficients, from KS roles."""
function block_E(m::Int, nd::Int, roles)
    C, cone, cdim = block_cone(m)
    E = zeros(nd + 1, cdim)
    for ℓ in 0:nd, (shift, sgn) in roles
        E[ℓ + 1, :] .+= sgn .* (C' * antidiag(m, ℓ - shift))
    end
    return E, cone, cdim
end

"""Block sizes and KS roles (shift, sign) for a degree-nd nonnegativity certificate on [0,1]."""
function ks_spec(nd::Int)
    if nd % 2 == 1                       # odd nd = 2κ+1: X,Y both (κ+1)×(κ+1)
        κ = (nd - 1) ÷ 2
        return κ + 1, κ + 1, [(1, 1.0)], [(0, 1.0), (1, -1.0)]
    else                                 # even nd = 2κ: X (κ+1)×(κ+1), Y κ×κ
        κ = nd ÷ 2
        return κ + 1, κ, [(0, 1.0)], [(1, 1.0), (2, -1.0)]
    end
end

# ---- native model ---------------------------------------------------

"""
    build_shape_exact_problem(inst) :: (IPMProblem, info)

Segment stalks are free and carry the objective; each piece gets a
Karlin–Studden certificate (SOC blocks for nd ≤ 3, SDP blocks for nd ≥ 4)
tied to its coefficients by one equality edge. C^k continuity on the spine.
"""
function build_shape_exact_problem(inst::ShapeSplineInstance)
    n, k, M = inst.n, inst.k, inst.M
    d = 1; ns = n + 1
    L = jet_operator(n, d, k, 1.0, :left)
    R = jet_operator(n, d, k, 1.0, :right)
    m = inst.shape === :monotone ? 1 : inst.shape === :convex ? 2 : 0
    nd = n - m
    @assert nd ≥ 0 "degree too low for this shape"

    # segment coefficients → certified-polynomial monomial coefficients
    Aseg = bern2mono(nd) * (falling(n, m) .* diffop(n, m))          # (nd+1)×(n+1)
    sx, sy, rolesX, rolesY = ks_spec(nd)
    EX, coneX, cdimX = block_E(sx, nd, rolesX)
    hasY = sy ≥ 1
    EY = coneY = cdimY = nothing
    hasY && ((EY, coneY, cdimY) = block_E(sy, nd, rolesY))
    max_block = max(sx, sy)

    row_ids = Int[]; col_ids = Int[]; blocks = Matrix{Float64}[]
    col_cone = AbstractCone[]; col_dim = Int[]
    rhs = Dict{Int, Vector{Float64}}()
    new_vertex!(dim, cone) = (push!(col_dim, dim); push!(col_cone, cone); length(col_dim))
    ec = Ref(0); new_edge!() = (ec[] += 1; ec[])
    place!(e, v, mat) = (push!(row_ids, e); push!(col_ids, v); push!(blocks, Matrix{Float64}(mat)))

    seg = [new_vertex!(ns, CofreeCone()) for _ in 1:M]

    # certificate per piece:  Aseg·P − EX·wX − EY·wY = 0
    for s in 1:M
        xv = new_vertex!(cdimX, coneX); e = new_edge!()
        place!(e, seg[s], Aseg); place!(e, xv, -EX)
        if hasY
            yv = new_vertex!(cdimY, coneY); place!(e, yv, -EY)
        end
    end

    # C^k continuity
    for s in 1:(M - 1)
        e = new_edge!(); place!(e, seg[s], R); place!(e, seg[s + 1], -L)
    end

    # optional ∫ f = 1
    if inst.density
        e = new_edge!()
        for s in 1:M
            w = (inst.knots[s + 1] - inst.knots[s]) / (n + 1)
            place!(e, seg[s], reshape(fill(w, ns), 1, ns))
        end
        rhs[e] = [1.0]
    end

    B = blocksparse(row_ids, col_ids, blocks)
    nvtx = length(col_dim)
    @assert nvtxs(B) == nvtx "column blocks not contiguous 1..nvtx"
    g = zeros(size(B, 1)); for (e, v) in rhs; g[rowrange(B, e)] .= v; end
    Q = IPM.allocblockdiag(B); fill!(Q, 0)
    for s in 1:M; block(Q, seg[s], seg[s], seg[s]) .= inst.Qpiece[s]; end
    c = zeros(size(B, 2))
    for s in 1:M; c[colrange(B, seg[s])] .= inst.cpiece[s]; end
    cones = AbstractCone[col_cone[v] for v in 1:nvtx]

    prob = IPMProblem(Q, B, c, g, cones)
    info = (segcols = [colrange(B, seg[s]) for s in 1:M], nvtx = nvtx,
            cert_degree = nd, max_block = max_block,
            cone_kind = max_block ≤ 2 ? :SOCP : :SDP)
    return prob, info
end

"""Quick build+solve demo; reports whether the certificate is SOCP or SDP."""
function shape_exact_demo(; shape = :nonneg, kwargs...)
    inst = generate_shape_instance(; shape = shape, kwargs...)
    prob, info = build_shape_exact_problem(inst)
    solve(prob, lp_settings())                         # JIT warmup
    t = @elapsed (res = solve(prob, lp_settings()))
    @printf("EXACT %-9s  n=%d nd=%d block=%d ⇒ %-4s  M=%d vtx=%d  %.1f ms  iters=%d  status=%s\n",
            string(shape), inst.n, info.cert_degree, info.max_block, string(info.cone_kind),
            inst.M, info.nvtx, 1e3t, res.niter, res.status)
    return inst, res, info
end

# =====================================================================
# JuMP reference model (for Mosek comparison)
# =====================================================================

"""Antidiagonal sum Σ_{i+j=ℓ} M[i+1,j+1] as a JuMP expression."""
function antidiag_expr(M, m::Int, ℓ::Int)
    return sum(M[i + 1, j + 1] for i in 0:(m - 1) for j in 0:(m - 1) if i + j == ℓ; init = 0.0)
end

"""Build the exact shape-spline as a JuMP model with PSD certificates."""
function build_jump_exact(inst::ShapeSplineInstance, optimizer)
    n, k, M = inst.n, inst.k, inst.M
    ns = n + 1
    L = jet_operator(n, 1, k, 1.0, :left)
    R = jet_operator(n, 1, k, 1.0, :right)
    m = inst.shape === :monotone ? 1 : inst.shape === :convex ? 2 : 0
    nd = n - m
    sx, sy, rolesX, rolesY = ks_spec(nd)
    hasY = sy ≥ 1

    model = Model(optimizer); set_silent(model)

    # segment variables (free)
    P = [@variable(model, [1:ns]) for _ in 1:M]

    # PSD certificate matrices per piece
    Xs = [@variable(model, [1:sx, 1:sx], PSD) for _ in 1:M]
    Ys = hasY ? [@variable(model, [1:sy, 1:sy], PSD) for _ in 1:M] : nothing

    # objective: Σ_s ½ P_s' Q_s P_s + c_s' P_s
    obj = sum(0.5 * P[s]' * inst.Qpiece[s] * P[s] + inst.cpiece[s]' * P[s] for s in 1:M)
    @objective(model, Min, obj)

    # Karlin–Studden certificate: monomial coefficients = antidiag sums
    Aseg = bern2mono(nd) * (falling(n, m) .* diffop(n, m))
    for s in 1:M
        mono = Aseg * P[s]   # monomial coefficients of certified polynomial
        for ℓ in 0:nd
            lhs = mono[ℓ + 1]
            rhs = 0.0
            for (shift, sgn) in rolesX
                rhs = rhs + sgn * antidiag_expr(Xs[s], sx, ℓ - shift)
            end
            if hasY
                for (shift, sgn) in rolesY
                    rhs = rhs + sgn * antidiag_expr(Ys[s], sy, ℓ - shift)
                end
            end
            @constraint(model, lhs == rhs)
        end
    end

    # C^k continuity
    for s in 1:(M - 1)
        @constraint(model, R * P[s] .== L * P[s + 1])
    end

    # optional ∫ f = 1
    if inst.density
        int_coeffs = [(inst.knots[s + 1] - inst.knots[s]) / (n + 1) for s in 1:M]
        @constraint(model, sum(int_coeffs[s] * sum(P[s]) for s in 1:M) == 1.0)
    end

    return model, P
end

# =====================================================================
# Benchmark: Sheaf IPM vs Mosek
# =====================================================================

function run_exact_benchmark(; optimizer = nothing, dual_optimizer = nothing, nwarmup::Int = 2, nruns::Int = 5)
    optimizer === nothing && error("Pass optimizer, e.g. run_exact_benchmark(optimizer=Mosek.Optimizer)")

    cases = [
        (:nonneg,   50, 6),   # M=50 pieces, degree 6 -> nd=6, block=4 (SDP)
        (:nonneg,   100, 6),
        (:monotone, 50, 6),   # nd=5, block=3 (SDP)
        (:monotone, 100, 6),
        (:convex,   50, 6),   # nd=4, block=3 (SDP)
        (:convex,   100, 6),
    ]

    println("="^80)
    println("Shape-Spline Exact (SDP) Benchmark: Sheaf IPM vs Mosek vs Mosek-D")
    println("="^80)
    println()
    if dual_optimizer !== nothing
        @printf("%-12s %7s %6s %5s %5s %9s %9s %9s %7s %7s\n",
                "Shape", "splines", "degree", "|V|", "|E|", "IPM(ms)", "Mosek", "Mosek-D", "P/IPM", "D/IPM")
    else
        @printf("%-12s %7s %6s %5s %5s %9s %9s %8s\n",
                "Shape", "splines", "degree", "|V|", "|E|", "IPM(ms)", "Mosek", "speedup")
    end
    println("-"^80)

    for (shape, M, n) in cases
        inst = generate_shape_instance(; shape, M, n, N = 20 * M)
        prob, info = build_shape_exact_problem(inst)
        nV = info.nvtx
        # |E|: certificate edges (M) + continuity edges (M-1)
        nE = M + (M - 1)

        # warmup
        for _ in 1:nwarmup; solve(prob, lp_settings()); end
        t_ipm = minimum([@elapsed solve(prob, lp_settings()) for _ in 1:nruns])

        for _ in 1:nwarmup
            m, _ = build_jump_exact(inst, optimizer)
            optimize!(m)
        end
        t_mosek = minimum([@elapsed begin
            m, _ = build_jump_exact(inst, optimizer)
            optimize!(m)
        end for _ in 1:nruns])

        if dual_optimizer !== nothing
            for _ in 1:nwarmup
                m, _ = build_jump_exact(inst, dual_optimizer)
                optimize!(m)
            end
            t_dual = minimum([@elapsed begin
                m, _ = build_jump_exact(inst, dual_optimizer)
                optimize!(m)
            end for _ in 1:nruns])

            @printf("%-12s %7d %6d %5d %5d %9.1f %9.1f %9.1f %6.2fx %6.2fx\n",
                    string(shape), M, n, nV, nE, t_ipm * 1000, t_mosek * 1000, t_dual * 1000,
                    t_mosek / t_ipm, t_dual / t_ipm)
        else
            @printf("%-12s %7d %6d %5d %5d %9.1f %9.1f %7.2fx\n",
                    string(shape), M, n, nV, nE, t_ipm * 1000, t_mosek * 1000,
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
        println("Solver: Mosek + Mosek-D")
    else
        using Clarabel
        optimizer = Clarabel.Optimizer
        dual_optimizer = Dualization.dual_optimizer(Clarabel.Optimizer)
        println("Solver: Clarabel + Clarabel-D (open-source)")
    end
    println("Runs: $(opts.nruns), Warmup: $(opts.nwarmup)\n")

    run_exact_benchmark(; optimizer, dual_optimizer,
                        nwarmup = opts.nwarmup, nruns = opts.nruns)
end
