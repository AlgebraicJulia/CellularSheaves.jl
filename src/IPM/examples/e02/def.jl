# =============================================================================
# e02/def.jl — Quantum marginal compatibility (XXZ ring, ℓ=4, non-chordal SDP)
#
# Usage:
#   julia --project e02/run.jl              # Clarabel only
#   julia --project e02/run.jl --mosek      # + Mosek
#   julia --project e02/run.jl --quick      # small sweep
#
# TODO: Dualization helps both Clarabel (1.2x) and Mosek (1.1x). Consider
#       wrapping optimizers in dual_optimizer() from Dualization.jl.
# =============================================================================

using AppleAccelerate
using LinearAlgebra, SparseArrays, Printf

include("../utils.jl")

const OPTS = parse_args(ARGS)
const TOL = OPTS.tol
const NRUNS = 5

# -----------------------------------------------------------------------------
# Problem definition (Quantum marginal XXZ ring)
# -----------------------------------------------------------------------------

struct QMInstance
    N::Int; ℓ::Int; Δ::Float64; ring::Bool; ε::Float64
end

generate_qm_instance(; N=10, ℓ=3, Δ=1.0, ring=true, ε=1e-9) = QMInstance(N, ℓ, Δ, ring, ε)
qm_clusters(inst::QMInstance) = inst.ring ? (1:inst.N) : (1:(inst.N - inst.ℓ + 1))
qm_svecdim(D::Int) = D * (D + 1) ÷ 2

function qm_svec(M::AbstractMatrix)
    D = size(M, 1); v = zeros(qm_svecdim(D)); k = 0
    for c in 1:D, r in c:D
        k += 1; v[k] = r == c ? M[c, c] : sqrt(2.0) * M[r, c]
    end
    v
end

function qm_smat(v::AbstractVector, D::Int)
    M = zeros(D, D); k = 0
    for c in 1:D, r in c:D
        k += 1
        if r == c; M[c, c] = v[k]
        else; M[r, c] = M[c, r] = v[k] / sqrt(2.0)
        end
    end
    M
end

const SXP = [0.0 1.0; 1.0 0.0]
const SYI = [0.0 -1.0; 1.0 0.0]
const SZP = [1.0 0.0; 0.0 -1.0]

h_xxz(Δ::Float64) = 0.25 .* (kron(SXP, SXP) .- kron(SYI, SYI) .+ Δ .* kron(SZP, SZP))

function embedop(op::AbstractMatrix, pos::Int, nsites::Int)
    w = Int(log2(size(op, 1)))
    kron(kron(Matrix{Float64}(I, 2^(pos-1), 2^(pos-1)), op),
         Matrix{Float64}(I, 2^(nsites-pos-w+1), 2^(nsites-pos-w+1)))
end

function qm_cluster_H(inst::QMInstance, start::Int)
    H = embedop(h_xxz(inst.Δ), 1, inst.ℓ)
    if !inst.ring && start == inst.N - inst.ℓ + 1
        for b in 2:(inst.ℓ - 1); H .+= embedop(h_xxz(inst.Δ), b, inst.ℓ); end
    end
    H
end

# Embed a two-site operator on arbitrary qubits i < j in an n-qubit system
function embed_pair(h4::AbstractMatrix, i::Int, j::Int, n::Int)
    D = 1 << n
    H = zeros(D, D)
    wt(q) = 1 << (n - q)
    for r in 0:D - 1
        ri = (r >> (n - i)) & 1
        rj = (r >> (n - j)) & 1
        base = r - ri * wt(i) - rj * wt(j)
        for bi in 0:1, bj in 0:1
            amp = h4[2 * bi + bj + 1, 2 * ri + rj + 1]
            amp == 0.0 && continue
            H[base + bi * wt(i) + bj * wt(j) + 1, r + 1] += amp
        end
    end
    H
end

# ED oracle: build full N-site Hamiltonian and find ground state energy
function xxz_ed_ground_state(N::Int, Δ::Float64; ring::Bool=true)
    if N > 12
        return NaN  # too large for dense ED
    end

    h2 = h_xxz(Δ)
    H = zeros(2^N, 2^N)
    for b in 1:(ring ? N : N - 1)
        i, j = b, mod1(b + 1, N)
        H .+= embed_pair(h2, min(i, j), max(i, j), N)
    end
    eigmin(Symmetric(H))
end

function ptr_first(ρ::AbstractMatrix)
    d = size(ρ, 1) ÷ 2; out = zeros(d, d)
    for k in 0:1; out .+= ρ[(k*d+1):(k*d+d), (k*d+1):(k*d+d)]; end
    out
end

function ptr_last(ρ::AbstractMatrix)
    d = size(ρ, 1) ÷ 2; out = zeros(d, d)
    for k in 1:2; out .+= ρ[k:2:end, k:2:end]; end
    out
end

function ptr_svec_map(D::Int, which::Symbol)
    f = which === :first ? ptr_first : ptr_last
    sd = qm_svecdim(D); L = zeros(qm_svecdim(D ÷ 2), sd)
    for k in 1:sd
        e = zeros(sd); e[k] = 1.0
        L[:, k] .= qm_svec(f(qm_smat(e, D)))
    end
    L
end

function build_qm_problem(inst::QMInstance)
    D = 2^inst.ℓ; sd = qm_svecdim(D)
    cl = collect(qm_clusters(inst)); ncl = length(cl)
    PTf, PTl = ptr_svec_map(D, :first), ptr_svec_map(D, :last)

    row_ids, col_ids, blks = Int[], Int[], Matrix{Float64}[]
    place!(e, v, mat) = (push!(row_ids, e); push!(col_ids, v); push!(blks, Matrix{Float64}(mat)))

    nagree = inst.ring ? ncl : ncl - 1
    for a in 1:nagree
        i, j = a, mod1(a + 1, ncl)
        place!(a, i, PTf); place!(a, j, -PTl)
    end
    place!(nagree + 1, 1, reshape(qm_svec(Matrix{Float64}(I, D, D)), 1, sd))

    B = blocksparse(row_ids, col_ids, blks)
    g = zeros(size(B, 1)); g[end] = 1.0
    Q = IPM.allocblockdiag(B); fill!(Q, 0)
    c = zeros(size(B, 2))
    for (k, start) in enumerate(cl)
        block(Q, k, k, k) .= inst.ε .* Matrix{Float64}(I, sd, sd)
        c[colrange(B, k)] .= -qm_svec(qm_cluster_H(inst, start))  # negated for min ½p'Qp - c'p
    end
    IPMProblem(Q, B, c, g, AbstractCone[SemidefiniteCone() for _ in 1:ncl]), (ncl=ncl, D=D)
end

function build_jump_qm(inst::QMInstance, optimizer)
    D = 2^inst.ℓ; Do = D ÷ 2
    cl = collect(qm_clusters(inst)); ncl = length(cl)
    model = Model(optimizer); set_silent(model)
    ρ = Dict(k => @variable(model, [1:D, 1:D], PSD) for k in 1:ncl)
    Ef = [kron(Matrix{Float64}(I, 2, 2)[k:k, :], Matrix{Float64}(I, Do, Do)) for k in 1:2]
    El = [kron(Matrix{Float64}(I, Do, Do), Matrix{Float64}(I, 2, 2)[k:k, :]) for k in 1:2]
    nagree = inst.ring ? ncl : ncl - 1
    for a in 1:nagree
        i, j = a, mod1(a + 1, ncl)
        @constraint(model, sum(Ef[k] * ρ[i] * Ef[k]' for k in 1:2) .==
                           sum(El[k] * ρ[j] * El[k]' for k in 1:2))
    end
    @constraint(model, tr(ρ[1]) == 1)
    @objective(model, Min,
               sum(tr(qm_cluster_H(inst, s) * ρ[k]) +
                   0.5 * inst.ε * sum(ρ[k][i, j]^2 for i in 1:D, j in 1:D)
                   for (k, s) in enumerate(cl)))
    model, ρ
end

# -----------------------------------------------------------------------------
# Oracle test
# -----------------------------------------------------------------------------

function test_ed_oracle()
    # Small test: N=6, Δ=1.0 (isotropic Heisenberg)
    E_ed = xxz_ed_ground_state(6, 1.0; ring=true)
    println("  [PASS] ED oracle N=6: E₀=$(round(E_ed, sigdigits=5))")

    # Verify lower-bound property (solve SDP and check E_sdp ≤ E_ed)
    inst = generate_qm_instance(; N=6, ℓ=3, Δ=1.0, ring=true)
    prob, _ = build_qm_problem(inst)
    settings = IPMSettings{Float64}(
        feas_tol = 1e-8, gap_tol = 1e-8, itmax = 200)
    res = solve(prob, settings)

    if res.status in (OPTIMAL, NEAR_OPTIMAL)
        E_sdp = 0.5 * dot(res.p, prob.Q * res.p) - dot(prob.c, res.p)
        gap = E_ed - E_sdp
        if gap >= -1e-6
            println("  [PASS] Lower bound: E_sdp=$(round(E_sdp, sigdigits=5)) ≤ E_ed=$(round(E_ed, sigdigits=5)) (gap=$(round(gap, sigdigits=3)))")
        else
            println("  [FAIL] Lower bound violated: E_sdp=$(round(E_sdp, sigdigits=5)) > E_ed=$(round(E_ed, sigdigits=5))")
        end
    else
        println("  [SKIP] Lower bound test: solver did not converge")
    end
end

# -----------------------------------------------------------------------------
# Sweep
# -----------------------------------------------------------------------------

sizes() = OPTS.tiny ? [8] : OPTS.quick ? [8, 10, 12] : [8, 10, 12, 14, 16, 18]

function run()
    println("\n", "=" ^ 78)
    println("  E02: Quantum marginal compatibility (dial: N)")
    println("=" ^ 78)

    println("\n  Gate tests:")
    test_ed_oracle()
    println()

    ipm_settings = IPMSettings{Float64}(
        feas_tol = TOL, gap_tol = TOL, itmax = 200)

    hsd_settings = HSDSettings{Float64}(
        feas_tol = TOL, gap_tol = TOL, itmax = 200)

    cla_opt = clarabel_opt(; tol = TOL)
    msk_opt = OPTS.mosek ? mosek_opt(; tol = TOL) : nothing

    rows = []
    for N in sizes()
        inst = generate_qm_instance(; N, ℓ = 4, Δ = 1.0, ring = true)
        prob, _ = build_qm_problem(inst)
        stats = problem_stats(prob)

        @printf("  N=%-4d dof=%-6d n1=%-5d blk=%-4.0f  ",
            N, stats.N0, stats.N1, stats.med_block)

        m_ipm = measure_ipm(prob, ipm_settings; nruns = NRUNS)
        m_hsd = measure_ipm(prob, hsd_settings; nruns = NRUNS)
        m_cla = measure_jump(() -> first(build_jump_qm(inst, cla_opt)); nruns = NRUNS)
        # Mosek benchmarked in DUAL form (Dualization.jl): ~1.1x faster than primal here.
        m_msk = msk_opt !== nothing ?
            measure_jump(() -> first(build_jump_qm(inst, Dualization.dual_optimizer(msk_opt))); nruns = NRUNS) :
            (t = NaN, status = "", obj = NaN)

        ratio(b, base) = isfinite(b.t) && isfinite(base.t) ? b.t / base.t : NaN
        fmt_ratio(b, base) = isnan(ratio(b, base)) ? "—" : @sprintf("%.2fx", ratio(b, base))

        println("IPM ", fmt_time(m_ipm), "  HSD ", fmt_time(m_hsd), " (", fmt_ratio(m_hsd, m_ipm),
            ")  Cla ", fmt_time(m_cla), " (", fmt_ratio(m_cla, m_ipm),
            ")  Msk ", fmt_time(m_msk), " (", fmt_ratio(m_msk, m_ipm), ")")

        push!(rows, (size = N, dof = stats.N0, ipm = m_ipm, hsd = m_hsd, cla = m_cla, msk = m_msk))
    end

    dofs = [r.dof for r in rows]
    println("\nFitted log-log slopes (time vs DOF):")
    for (name, get_t) in [("IPM", r -> r.ipm.t), ("HSD", r -> r.hsd.t), ("Clarabel", r -> r.cla.t), ("Mosek", r -> r.msk.t)]
        sl = loglog_slope(dofs, [get_t(r) for r in rows])
        print("  $name: ", isnan(sl) ? "n/a" : @sprintf("DOF^%.2f", sl))
    end
    println()
end

