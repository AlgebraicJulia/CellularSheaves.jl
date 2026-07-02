######################################################################
# quantum_marginal.jl
#
# QUANTUM MARGINAL ENERGY BOUNDS on XXZ spin chains and rings — the
# 2-RDM row of local_consistency_sheaf.md §4, realized. Certified
# lower bounds on the ground-state energy from overlapping cluster
# reduced density matrices, following Barthel & Hübener, PRL 108,
# 200404 (2012) (parallel: Baumgratz & Plenio, NJP 14, 023027):
#
#     min   Σ_C ⟨H_C, ρ_C⟩  (+ ε/2 Σ‖ρ_C‖_F²)
#     s.t.  Tr_first(ρ_C) = Tr_last(ρ_{C+1})    (overlap agreement)
#           tr(ρ_1) = 1                          (one pin; the rest is
#                                                 implied through the
#                                                 agreement edges)
#           ρ_C ∈ SemidefiniteCone               per cluster
#
# Clusters are ℓ consecutive sites; consecutive clusters overlap in
# ℓ−1 sites, and their agreement is PARTIAL TRACE on both ends — the
# corpus's first genuinely multi-Kraus restriction maps (Kraus family
# {I⊗⟨k|}, k = 0,1; everything before this was a compression, i.e.
# single-Kraus). Every true ground state's marginals are feasible, so
# the optimum is a strict lower bound on E₀ for ANY cluster size.
#
# REAL SYMMETRIC WLOG. The XXZ bond Hamiltonian ¼(σˣσˣ + σʸσʸ + Δσᶻσᶻ)
# is real in the computational basis (σʸ⊗σʸ is real), and Re(ρ) of a
# feasible Hermitian family is feasible at the same energy — so the
# library's real SemidefiniteCone is used without loss. (General
# complex Hamiltonians need the standard 2×-size real embedding.)
#
# THE SCIENTIFIC POINT — a theorem failing. The corpus's classical
# story (local_consistency_sheaf.md §5) is "positive local sections
# glue iff the pattern is chordal / a junction tree." A chain is a
# path — trivially chordal — and yet locally consistent quantum
# marginals need NOT extend to any global state: pin both 2-site RDMs
# of a 3-site chain to the Bell state Φ⁺; each has one-site marginal
# I/2, the overlap agreement holds EXACTLY, and no 3-qubit state has
# both as marginals (monogamy of entanglement; deciding consistency is
# QMA-complete, Liu 2006). The energy shadow of the same fact: on the
# 3-site Heisenberg chain the ℓ=2 relaxation puts BOTH bonds in the
# singlet, E_sdp = −3/2 against exact E₀ = −1. Gluing fails on the
# simplest possible topology — the obstruction has left the base graph
# and moved into the cone.
#
# MEASURED FACTS (quantum_marginal_oracle.py, CLARABEL):
#   • ℓ=2 ring bound = −(2+Δ)/4 per bond EXACTLY (every cluster is
#     independently a singlet-like ground state; all one-site marginals
#     I/2, perfectly consistent — 1-site overlaps carry no physics
#     beyond the single bond).
#   • The ladder improves ONLY AT ODD ℓ (both Δ): per-bond Δ=1:
#     ℓ=2: −0.750, ℓ=3,4: −0.500, ℓ=5,6: −0.4671 (thermodynamic exact
#     ¼ − ln2 ≈ −0.4431). Ring optimum = translation-invariant
#     single-cluster optimum (cyclic symmetrization), verified.
#   • H¹: chain 0, ring 1 — on the ring exactly ONE left-null cycle
#     survives, the trace functional: PTᵀ_first(Y) is the functional
#     ⟨I⊗Y,·⟩ and PTᵀ_last(Y) is ⟨Y⊗I,·⟩, and I⊗Y = Y'⊗I forces
#     Y = cI. So the dual is a 1-dimensional coset on rings; Uzawa's
#     CG pins the basepoint, as in the tensor-mesh saga.
#
# Written against the CellularSheaves.IPM PR-67 API; not executed here.
# quantum_marginal_oracle.py is the numerical ground truth.
######################################################################

using CellularSheaves.IPM
using CellularSheaves.BlockSparseArrays: colrange, rowrange, blocksparse, block, nvtxs
using LinearAlgebra
using Printf

# ---- svec (library convention: col-major LOWER, diag first, ×√2) ----

qm_svecdim(D::Int) = D * (D + 1) ÷ 2

function qm_svec(M::AbstractMatrix)
    D = size(M, 1); v = zeros(qm_svecdim(D)); k = 0
    for c in 1:D, r in c:D
        k += 1
        v[k] = r == c ? M[c, c] : sqrt(2.0) * M[r, c]
    end
    return v
end

function qm_smat(v::AbstractVector, D::Int)
    M = zeros(D, D); k = 0
    for c in 1:D, r in c:D
        k += 1
        if r == c
            M[c, c] = v[k]
        else
            M[r, c] = M[c, r] = v[k] / sqrt(2.0)
        end
    end
    return M
end

# ---- spin operators and the XXZ chain --------------------------------

const SXP = [0.0 1.0; 1.0 0.0]
const SYI = [0.0 -1.0; 1.0 0.0]          # σʸ = i·SYI, so σʸ⊗σʸ = −SYI⊗SYI (real)
const SZP = [1.0 0.0; 0.0 -1.0]

"""Bond Hamiltonian ¼(σˣ⊗σˣ + σʸ⊗σʸ + Δ σᶻ⊗σᶻ) — real symmetric."""
h_xxz(Δ::Float64) = 0.25 .* (kron(SXP, SXP) .- kron(SYI, SYI) .+ Δ .* kron(SZP, SZP))

"""Embed `op` (width w) acting on sites pos..pos+w−1 of nsites (site 1 = slowest)."""
function embedop(op::AbstractMatrix, pos::Int, nsites::Int)
    w = Int(log2(size(op, 1)))
    kron(kron(Matrix{Float64}(I, 2^(pos - 1), 2^(pos - 1)), op),
         Matrix{Float64}(I, 2^(nsites - pos - w + 1), 2^(nsites - pos - w + 1)))
end

"""Exact ground energy by dense diagonalization (fine for N ≤ 12)."""
function exact_ground(N::Int, Δ::Float64; ring::Bool = true)
    H = sum(embedop(h_xxz(Δ), i, N) for i in 1:(N - 1))
    if ring && N > 2
        for (A, B, w) in ((SXP, SXP, 0.25), (SYI, SYI, -0.25), (SZP, SZP, 0.25Δ))
            H .+= w .* (embedop(A, N, N) * embedop(B, 1, N))
        end
    end
    return eigvals(Symmetric(H))[1]
end

# ---- partial traces in svec coordinates ------------------------------

function ptr_first(ρ::AbstractMatrix)
    d = size(ρ, 1) ÷ 2
    out = zeros(d, d)
    for k in 0:1
        out .+= ρ[(k * d + 1):(k * d + d), (k * d + 1):(k * d + d)]
    end
    return out
end

function ptr_last(ρ::AbstractMatrix)
    d = size(ρ, 1) ÷ 2
    out = zeros(d, d)
    for k in 1:2
        out .+= ρ[k:2:end, k:2:end]
    end
    return out
end

"""svecdim(D/2) × svecdim(D) matrix of Tr over the first/last site."""
function ptr_svec_map(D::Int, which::Symbol)
    f = which === :first ? ptr_first : ptr_last
    sd = qm_svecdim(D)
    L = zeros(qm_svecdim(D ÷ 2), sd)
    for k in 1:sd
        e = zeros(sd); e[k] = 1.0
        L[:, k] .= qm_svec(f(qm_smat(e, D)))
    end
    return L
end

# ---- instance ---------------------------------------------------------

struct QMInstance
    N::Int          # sites
    ℓ::Int          # cluster size (≥ 2)
    Δ::Float64      # XXZ anisotropy (Δ = 1 is Heisenberg)
    ring::Bool
    ε::Float64
end

generate_qm_instance(; N::Int = 10, ℓ::Int = 3, Δ::Float64 = 1.0,
                     ring::Bool = true, ε::Float64 = 1e-9) =
    (@assert 2 ≤ ℓ ≤ N; QMInstance(N, ℓ, Δ, ring, ε))

qm_clusters(inst::QMInstance) = inst.ring ? (1:inst.N) : (1:(inst.N - inst.ℓ + 1))

"""Bonds owned by the cluster starting at `start`, embedded locally.
Ring: each cluster owns its first bond (covers all N once). Chain: the
last cluster additionally owns its remaining ℓ−2 bonds."""
function qm_cluster_H(inst::QMInstance, start::Int)
    h = h_xxz(inst.Δ)
    H = embedop(h, 1, inst.ℓ)
    if !inst.ring && start == inst.N - inst.ℓ + 1
        for b in 2:(inst.ℓ - 1)
            H .+= embedop(h, b, inst.ℓ)
        end
    end
    return H
end

# ---- native builder ----------------------------------------------------

function build_qm_problem(inst::QMInstance)
    D = 2^inst.ℓ
    sd = qm_svecdim(D)
    cl = collect(qm_clusters(inst))
    ncl = length(cl)
    PTf, PTl = ptr_svec_map(D, :first), ptr_svec_map(D, :last)

    row_ids = Int[]; col_ids = Int[]; blks = Matrix{Float64}[]
    place!(e, v, mat) = (push!(row_ids, e); push!(col_ids, v); push!(blks, Matrix{Float64}(mat)))

    nagree = inst.ring ? ncl : ncl - 1
    for a in 1:nagree
        i, j = a, mod1(a + 1, ncl)
        place!(a, i, PTf); place!(a, j, -PTl)
    end
    place!(nagree + 1, 1, reshape(qm_svec(Matrix{Float64}(I, D, D)), 1, sd))  # trace pin

    B = blocksparse(row_ids, col_ids, blks)
    @assert nvtxs(B) == ncl
    g = zeros(size(B, 1)); g[end] = 1.0

    Q = IPM.allocblockdiag(B); fill!(Q, 0)
    c = zeros(size(B, 2))
    for (k, start) in enumerate(cl)
        block(Q, k, k, k) .= inst.ε .* Matrix{Float64}(I, sd, sd)
        c[colrange(B, k)] .= qm_svec(qm_cluster_H(inst, start))   # tr(Hρ) = ⟨svec H, svec ρ⟩
    end
    prob = IPMProblem(Q, B, c, g, AbstractCone[SemidefiniteCone() for _ in 1:ncl])
    info = (rdmcols = Dict(k => colrange(B, k) for k in 1:ncl), ncl = ncl, D = D,
            h1 = inst.ring ? 1 : 0)     # measured: exactly the trace cycle on rings
    return prob, info
end

"""Objective at the solution, and the ε-corrected certified bound:
the exact ε-optimum v_ε satisfies v_ε ≤ E₀ + ε·ncl/2 (evaluate at the
true ground state's marginals; ‖ρ‖_F ≤ tr ρ = 1), so
E₀ ≥ v_ε − ε·ncl/2, up to solver tolerance."""
function qm_bound(inst::QMInstance, prob, res, info)
    v = prob.c' * res.p + 0.5 * inst.ε * sum(abs2, res.p)
    return v, v - 0.5 * inst.ε * info.ncl
end

# ---- demos --------------------------------------------------------------

"""Energy ladder on a ring, against exact diagonalization. Expect:
ℓ=2 exactly −(2+Δ)/4 per bond; improvement only at odd ℓ; every bound
below E_exact."""
function qm_ladder_demo(; N::Int = 10, Δ::Float64 = 1.0, ℓs = 2:4, kwargs...)
    Eex = exact_ground(N, Δ; ring = true)
    @printf("XXZ Δ=%.2f  N=%d ring   exact E0 = %.8f  [%+.6f /bond]\n", Δ, N, Eex, Eex / N)
    out = []
    for ℓ in ℓs
        inst = generate_qm_instance(; N, ℓ, Δ, ring = true, kwargs...)
        prob, info = build_qm_problem(inst)
        solve(prob, tensor_settings_qm())
        t = @elapsed (res = solve(prob, tensor_settings_qm()))
        v, cert = qm_bound(inst, prob, res, info)
        @printf("  ℓ=%d: E_sdp = %.8f  [%+.6f /bond]  gap %.4f  H1=%d  %.1f ms it=%d %s\n",
                ℓ, cert, cert / N, Eex - cert, info.h1, 1e3t, res.niter, res.status)
        push!(out, (ℓ, cert))
    end
    ℓ2 = findfirst(x -> x[1] == 2, out)
    ℓ2 !== nothing &&
        @printf("  check: ℓ=2 analytic −(2+Δ)/4·N = %.6f\n", -(2 + Δ) / 4 * N)
    return out
end

"""The monogamy gap, in energy form: 3-site Heisenberg chain, ℓ=2. The
relaxation places BOTH bonds in the singlet (their shared one-site
marginals are I/2 — consistent!), certifying −3/2 where the truth is
−1. The feasibility form of the same witness (double Bell, globally
inextendable) is in the oracle."""
function qm_monogamy_demo(; kwargs...)
    inst = generate_qm_instance(; N = 3, ℓ = 2, Δ = 1.0, ring = false, kwargs...)
    prob, info = build_qm_problem(inst)
    res = solve(prob, tensor_settings_qm())
    _, cert = qm_bound(inst, prob, res, info)
    Eex = exact_ground(3, 1.0; ring = false)
    @printf("3-site Heisenberg chain: exact %.6f   ℓ=2 bound %.6f   monogamy gap %.4f\n",
            Eex, cert, Eex - cert)
    return inst, res, info
end

"""Solver settings: same Uzawa recipe as the spline corpus; defined
here so the file stays standalone."""
tensor_settings_qm(; raug = 1e5) =
    IPMSettings{Float64}(kkt = UzawaSettings{Float64}(raug = raug),
                         feas_tol = 1e-8, gap_tol = 1e-8, itmax = 200)

# ---- benchmark -----------------------------------------------------------

function run_qm_benchmark(; optimizer = nothing, dual_optimizer = nothing,
                          nwarmup::Int = 2, nruns::Int = 5)
    optimizer === nothing && error("Pass optimizer, e.g. run_qm_benchmark(optimizer=Mosek.Optimizer)")
    # (N, ℓ, ring, label, raug) — raug tuned per case
    cases = [
        (6,  2, true,  "N6  ℓ2 ring", 1e3),
        (6,  3, true,  "N6  ℓ3 ring", 1e6),
        (10, 2, true,  "N10 ℓ2 ring", 1e1),
        (10, 3, true,  "N10 ℓ3 ring", 1e6),
        (10, 4, true,  "N10 ℓ4 ring", 1e3),
        (12, 3, true,  "N12 ℓ3 ring", 1e5),
        (12, 4, true,  "N12 ℓ4 ring", 1e4),
    ]
    println("="^85)
    println("Quantum Marginal Benchmark (SDP): Sheaf IPM vs Mosek vs Mosek-D")
    println("="^85)
    if dual_optimizer !== nothing
        @printf("%-14s %5s %5s %5s %9s %9s %9s %7s %7s\n",
                "Config", "DOF", "|V|", "H1", "IPM(ms)", "Mosek", "Mosek-D", "P/IPM", "D/IPM")
    else
        @printf("%-14s %5s %5s %5s %9s %9s %8s\n",
                "Config", "DOF", "|V|", "H1", "IPM(ms)", "Mosek", "speedup")
    end
    println("-"^85)
    for (N, ℓ, ring, label, raug) in cases
        inst = generate_qm_instance(; N, ℓ, ring)
        prob, info = build_qm_problem(inst)
        settings = tensor_settings_qm(; raug)
        for _ in 1:nwarmup; solve(prob, settings); end
        t_ipm = minimum([@elapsed solve(prob, settings) for _ in 1:nruns])
        for _ in 1:nwarmup
            m, _ = build_jump_qm(inst, optimizer); optimize!(m)
        end
        t_mosek = minimum([@elapsed begin
            m, _ = build_jump_qm(inst, optimizer); optimize!(m)
        end for _ in 1:nruns])

        if dual_optimizer !== nothing
            for _ in 1:nwarmup
                m, _ = build_jump_qm(inst, dual_optimizer); optimize!(m)
            end
            t_dual = minimum([@elapsed begin
                m, _ = build_jump_qm(inst, dual_optimizer); optimize!(m)
            end for _ in 1:nruns])
            @printf("%-14s %5d %5d %5d %9.1f %9.1f %9.1f %6.2fx %6.2fx\n",
                    label, size(prob.B, 2), info.ncl, info.h1,
                    t_ipm * 1000, t_mosek * 1000, t_dual * 1000,
                    t_mosek / t_ipm, t_dual / t_ipm)
        else
            @printf("%-14s %5d %5d %5d %9.1f %9.1f %7.2fx\n",
                    label, size(prob.B, 2), info.ncl, info.h1,
                    t_ipm * 1000, t_mosek * 1000, t_mosek / t_ipm)
        end
    end
end

# =====================================================================
# JuMP reference (PSD matrix variables; partial traces as Kraus sums)
# =====================================================================

using JuMP

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
    return model, ρ
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

    run_qm_benchmark(; optimizer, dual_optimizer,
                     nwarmup = opts.nwarmup, nruns = opts.nruns)
end
