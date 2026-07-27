# =============================================================================
# e03/def.jl — Tensor-spline shape-constrained regression (dense-QP habitat)
#
# Usage:
#   julia --project e03/run.jl              # Clarabel only
#   julia --project e03/run.jl --mosek      # + Mosek
#   julia --project e03/run.jl --quick      # small sweep
# =============================================================================

using AppleAccelerate
using LinearAlgebra, SparseArrays, Printf, Random

include("../utils.jl")

const OPTS = parse_args(ARGS)
const TOL = OPTS.tol
const NRUNS = 5

# -----------------------------------------------------------------------------
# Bernstein polynomial numerics
# -----------------------------------------------------------------------------

falling(n::Integer, r::Integer) = prod(n - j for j in 0:(r-1); init = 1.0)

function bernstein_gram(m::Integer)
    G = Matrix{Float64}(undef, m + 1, m + 1)
    for i in 0:m, k in 0:m
        G[i+1, k+1] = binomial(m, i) * binomial(m, k) / ((2m + 1) * binomial(2m, i + k))
    end
    G
end

function endpoint_jet_coeffs(n::Integer, r::Integer, T::Real, side::Symbol)
    w = zeros(n + 1); scale = falling(n, r) / T^r
    if side === :left
        for l in 0:r; w[l+1] = scale * (-1)^(r-l) * binomial(r, l); end
    else
        for l in 0:r; w[(n-r+l)+1] = scale * (-1)^(r-l) * binomial(r, l); end
    end
    w
end

function jet_operator(n::Integer, d::Integer, k::Integer, T::Real, side::Symbol)
    R = zeros((k+1) * d, (n+1) * d)
    for r in 0:k
        w = endpoint_jet_coeffs(n, r, T, side)
        for i in 0:n, a in 1:d
            R[r*d + a, i*d + a] = w[i+1]
        end
    end
    R
end

function diffop(n::Int, m::Int)
    D = zeros(n + 1 - m, n + 1)
    for i in 0:(n-m), j in 0:m
        D[i+1, i+j+1] = (-1)^(m-j) * binomial(m, j)
    end
    D
end

function bernstein_eval(n::Int, t::Float64)
    φ = zeros(n + 1)
    for i in 0:n; φ[i+1] = binomial(n, i) * t^i * (1 - t)^(n - i); end
    φ
end

deriv_gram(n::Int, a::Int) = a == 0 ? bernstein_gram(n) :
    (falling(n, a)^2) .* (diffop(n, a)' * bernstein_gram(n - a) * diffop(n, a))

function thinplate_gram(n::Int, hx::Float64, hy::Float64)
    G0, G1, G2 = deriv_gram(n, 0), deriv_gram(n, 1), deriv_gram(n, 2)
    G = hx * hy .* (kron(G0, G2) ./ hx^4 .+ 2 .* kron(G1, G1) ./ (hx^2 * hy^2) .+
                    kron(G2, G0) ./ hy^4)
    G ./= tr(G)
    G
end

truth2(x, y) = exp(-((x - 0.3)^2 + (y - 0.35)^2) / (2 * 0.12^2)) +
               0.6 * exp(-((x - 0.7)^2 + (y - 0.75)^2) / (2 * 0.10^2))

locate01(t, M) = (i = clamp(floor(Int, t * M) + 1, 1, M); (i, clamp(t * M - (i - 1), 0.0, 1.0)))

# -----------------------------------------------------------------------------
# Problem definition (Tensor-product Bézier surface regression)
# -----------------------------------------------------------------------------

struct TensorSplineInstance
    Mx::Int; My::Int; n::Int; k::Int
    active::BitMatrix
    kx::Matrix{Int}; ky::Matrix{Int}
    Qp::Dict{Tuple{Int,Int}, Matrix{Float64}}
    cp::Dict{Tuple{Int,Int}, Vector{Float64}}
    λr::Float64; ε::Float64
end

function sample_active(rng, active, Mx, My)
    while true
        x, y = rand(rng), rand(rng)
        a, _ = locate01(x, Mx); b, _ = locate01(y, My)
        active[a, b] && return (x, y)
    end
end

function generate_tensor_instance(; Mx::Int = 2, My::Int = 2, n::Int = 4, k::Int = 2,
                                   N::Int = 100 * Mx * My, σ::Float64 = 0.05,
                                   λr::Float64 = 1e-4, ε::Float64 = 1e-8, seed::Int = 1)
    rng = MersenneTwister(seed)
    ns = n + 1; d = ns * ns
    active = trues(Mx, My)
    kx = fill(k, Mx - 1, My); ky = fill(k, Mx, My - 1)

    Qp = Dict{Tuple{Int,Int}, Matrix{Float64}}()
    cp = Dict{Tuple{Int,Int}, Vector{Float64}}()
    for b in 1:My, a in 1:Mx
        active[a, b] || continue
        Qp[(a, b)] = zeros(d, d); cp[(a, b)] = zeros(d)
    end

    for _ in 1:N
        (x, y) = sample_active(rng, active, Mx, My)
        z = truth2(x, y) + σ * randn(rng)
        a, u = locate01(x, Mx); b, v = locate01(y, My)
        φ = kron(bernstein_eval(n, v), bernstein_eval(n, u))
        Qp[(a, b)] .+= 2.0 .* (φ * φ')
        cp[(a, b)] .-= 2.0 .* z .* φ
    end

    TensorSplineInstance(Mx, My, n, k, active, kx, ky, Qp, cp, λr, ε)
end

function patch_QC(inst::TensorSplineInstance)
    ns = inst.n + 1; d = ns * ns
    hx, hy = 1.0 / inst.Mx, 1.0 / inst.My
    Qbase = inst.λr .* thinplate_gram(inst.n, hx, hy) .+ inst.ε .* Matrix{Float64}(I, d, d)
    Q = Dict{Tuple{Int,Int}, Matrix{Float64}}()
    C = Dict{Tuple{Int,Int}, Vector{Float64}}()
    for b in 1:inst.My, a in 1:inst.Mx
        inst.active[a, b] || continue
        Qi = copy(Qbase) .+ inst.Qp[(a, b)]
        ci = copy(inst.cp[(a, b)])
        Q[(a, b)] = Qi; C[(a, b)] = ci
    end
    Q, C
end

function build_tensor_problem(inst::TensorSplineInstance)
    Mx, My, n, k = inst.Mx, inst.My, inst.n, inst.k
    ns = n + 1; d = ns * ns
    Ins = Matrix{Float64}(I, ns, ns)
    Rj(ke) = jet_operator(n, 1, ke, 1.0, :right)
    Lj(ke) = jet_operator(n, 1, ke, 1.0, :left)

    row_ids, col_ids, blocks = Int[], Int[], Matrix{Float64}[]
    col_cone = AbstractCone[]; col_dim = Int[]
    new_vertex!(dim, cone) = (push!(col_dim, dim); push!(col_cone, cone); length(col_dim))
    ec = Ref(0); new_edge!() = (ec[] += 1; ec[])
    place!(e, v, mat) = (push!(row_ids, e); push!(col_ids, v); push!(blocks, Matrix{Float64}(mat)))

    patch = Dict{Tuple{Int,Int}, Int}()
    for b in 1:My, a in 1:Mx
        inst.active[a, b] || continue
        patch[(a, b)] = new_vertex!(d, PositiveCone())
    end

    for b in 1:My, a in 1:Mx
        haskey(patch, (a, b)) || continue
        if a < Mx && haskey(patch, (a+1, b)) && inst.kx[a, b] >= 0
            ke = inst.kx[a, b]; e = new_edge!()
            place!(e, patch[(a, b)], kron(Ins, Rj(ke)))
            place!(e, patch[(a+1, b)], -kron(Ins, Lj(ke)))
        end
        if b < My && haskey(patch, (a, b+1)) && inst.ky[a, b] >= 0
            ke = inst.ky[a, b]; e = new_edge!()
            place!(e, patch[(a, b)], kron(Rj(ke), Ins))
            place!(e, patch[(a, b+1)], -kron(Lj(ke), Ins))
        end
    end

    B = blocksparse(row_ids, col_ids, blocks)
    nvtx = length(col_dim)
    g = zeros(size(B, 1))

    Qd, Cd = patch_QC(inst)
    Q = IPM.allocblockdiag(B); fill!(Q, 0)
    c = zeros(size(B, 2))
    for ((a, b), v) in patch
        block(Q, v, v, v) .= Qd[(a, b)]
        c[colrange(B, v)] .= .-Cd[(a, b)]
    end

    cones = AbstractCone[col_cone[v] for v in 1:nvtx]
    IPMProblem(Q, B, c, g, cones), (patchcols = Dict(cell => colrange(B, v) for (cell, v) in patch), nvtx = nvtx)
end

function build_jump_tensor(inst::TensorSplineInstance, optimizer)
    prob, _ = build_tensor_problem(inst)
    n = size(prob.B, 2)

    model = Model(optimizer)
    set_silent(model)

    @variable(model, x[1:n], lower_bound=0)

    Q_sparse = sparse(prob.Q)
    B_sparse = sparse(prob.B)

    @objective(model, Min, 0.5 * x' * Q_sparse * x - prob.c' * x)
    @constraint(model, B_sparse * x .== prob.g)

    model, x
end

# -----------------------------------------------------------------------------
# Sweep
# -----------------------------------------------------------------------------

sizes() = OPTS.tiny ? [2] : OPTS.quick ? [2, 4, 8] : [2, 4, 8, 12, 16, 20]

function run()
    println("\n", "=" ^ 78)
    println("  E03: Tensor-spline regression (dial: M)")
    println("=" ^ 78)

    ipm_settings = IPMSettings{Float64}(
        feas_tol = TOL, gap_tol = TOL, itmax = 200)
    hsd_settings = HSDSettings{Float64}(
        feas_tol = TOL, gap_tol = TOL, itmax = 200)

    cla_opt = clarabel_opt(; tol = TOL)
    msk_opt = OPTS.mosek ? mosek_opt(; tol = TOL) : nothing

    rows = []
    for M in sizes()
        inst = generate_tensor_instance(; Mx = M, My = M, n = 4, k = 2, N = 50 * M * M)
        prob, _ = build_tensor_problem(inst)
        stats = problem_stats(prob)

        @printf("  M=%-4d dof=%-6d n1=%-5d blk=%-4.0f  ",
            M, stats.N0, stats.N1, stats.med_block)

        m_ipm = measure_ipm(prob, ipm_settings; nruns = NRUNS)
        m_hsd = measure_ipm(prob, hsd_settings; nruns = NRUNS)
        m_cla = measure_jump(() -> first(build_jump_tensor(inst, cla_opt)); nruns = NRUNS)
        m_msk = msk_opt !== nothing ?
            measure_jump(() -> first(build_jump_tensor(inst, msk_opt)); nruns = NRUNS) :
            (t = NaN, status = "", obj = NaN)
        # Mosek benchmarked in PRIMAL form (faster than dualized here).

        ratio(b, base) = isfinite(b.t) && isfinite(base.t) ? b.t / base.t : NaN
        fmt_ratio(b, base) = isnan(ratio(b, base)) ? "—" : @sprintf("%.2fx", ratio(b, base))

        println("IPM ", fmt_time(m_ipm), "  HSD ", fmt_time(m_hsd), " (", fmt_ratio(m_hsd, m_ipm),
            ")  Cla ", fmt_time(m_cla), " (", fmt_ratio(m_cla, m_ipm),
            ")  Msk ", fmt_time(m_msk), " (", fmt_ratio(m_msk, m_ipm), ")")

        push!(rows, (size = M, dof = stats.N0, ipm = m_ipm, hsd = m_hsd, cla = m_cla, msk = m_msk))
    end

    dofs = [r.dof for r in rows]
    println("\nFitted log-log slopes (time vs DOF):")
    for (name, get_t) in [("IPM", r -> r.ipm.t), ("HSD", r -> r.hsd.t), ("Clarabel", r -> r.cla.t), ("Mosek", r -> r.msk.t)]
        sl = loglog_slope(dofs, [get_t(r) for r in rows])
        print("  $name: ", isnan(sl) ? "n/a" : @sprintf("DOF^%.2f", sl))
    end
    println()
end

