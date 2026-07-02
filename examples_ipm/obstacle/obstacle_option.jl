######################################################################
# obstacle_option.jl
#
# OBSTACLE PROBLEMS / AMERICAN OPTIONS — the example where Q stops
# being regularization and IS the problem: the objective is the PDE
# energy, the cells are overlapping subdomain patches (Schwarz domain
# decomposition semantics: Lions 1988; Badea 1991/2006; Tai's
# linear-convergence theory), and the orthant cell is the option's
# EARLY-EXERCISE PREMIUM w = u − ψ ≥ 0, vanishing identically on the
# exercise region. No relaxation gap anywhere — this is the corpus's
# pure-Q, pure-glue flagship, and its solver novelty is a MACROSCOPIC
# ACTIVE SET against dense Q (strict complementarity degenerating
# exactly at the free boundary).
#
# Perpetual American put (Black–Scholes, log-price x = log S),
# symmetrized by the integrating factor ρ = e^{2μx/σ²}, μ = r − σ²/2:
#   ρLu = −(pu')' + rρu,  p = ½σ²ρ   ⇒   the VI is the QP
#   min ½wᵀAw + (Aψ)ᵀw,  w ≥ 0,  A = weighted stiffness + lumped mass
# (A is PD natively — ε = 0; the mass term is the floor). The obstacle
# ψ = (K − eˣ)⁺ rides entirely in shifts (c and the boundary pins).
#
# One precision the literature forces (and the md expands): the
# Schwarz papers analyze ITERATIVE subdomain sweeps whose rate depends
# on the overlap; this builder is the MONOLITHIC QP rewritten on
# patches — exact by construction, overlap-invariant (measured:
# ≤ 3.5e-10 across overlaps 1–8). The iteration converges to us.
#
# MEASURED (obstacle_option_oracle.py, CLARABEL; K=1, r=0.05, σ=0.3):
#   • Merton closed form: S* = γK/(1+γ), γ = 2r/σ² = 1.111,
#     S* = 0.526316; BS-ODE residual of the formula 2.1e-5 (FD-limited).
#   • n=801, 4 patches, overlap 4: patched ≡ monolithic (value 2.5e-10,
#     ‖Δu‖∞ 5.0e-7); ‖u − PSOR‖∞ 1.1e-5 (Cryer); ‖u − Merton‖∞ 1.2e-6;
#     free boundary within one grid cell of S* (0.5300 vs 0.5263,
#     h_S ≈ 0.0029). Right boundary pinned to the closed form (the
#     domain-truncation value at S = 4 is 0.0498, NOT 0).
#   • Contact measure (the obstacle reaction, = cone dual on premium
#     cells): ≥ +1.9e-4 on the contact set, ≤ 2.1e-7 off it.
#   • American put T=1: implicit-Euler chain of obstacle QPs (200
#     steps) = 0.098677 vs CRR(2000) 0.098694, |Δ| = 1.7e-5.
#   • Game option (Kifer; Zeng–Zhou two-sided Schwarz): both contact
#     sets live (514 lower / 1 upper at δ = 0.02).
#   • 1D IDENTITY: the minimal-surface obstacle solution equals the
#     Dirichlet one (‖Δ‖∞ = 2.8e-7) — any energy ∫φ(u') is piecewise
#     linear off contact with tangency at contact; the free-boundary
#     difference between energies is a ≥ 2D phenomenon. The SOC mode
#     is therefore an EXACT cross-check of SecondOrderCone leaves
#     against native Q.
#   • H¹ of the patch chain: 0.
#
# Written against the CellularSheaves.IPM PR-67 API; not executed here.
# obstacle_option_oracle.py is the numerical ground truth.
######################################################################

using CellularSheaves.IPM
using CellularSheaves.BlockSparseArrays: colrange, rowrange, blocksparse, block, nvtxs
using LinearAlgebra
using Printf

# ---- Black–Scholes data -------------------------------------------------------

struct ObstacleInstance
    K::Float64; r::Float64; σ::Float64
    xlo::Float64; xhi::Float64; n::Int
    P::Int; olap::Int
    mode::Symbol                     # :dirichlet | :game
    δ::Float64                       # game-option penalty (mode = :game)
end

obstacle_instance(; K = 1.0, r = 0.05, σ = 0.3, xlo = log(0.05), xhi = log(4.0),
                  n = 801, P = 4, olap = 4, mode = :dirichlet, δ = 0.02) =
    ObstacleInstance(K, r, σ, xlo, xhi, n, P, olap, mode, δ)

function merton_perpetual(K, r, σ)
    γ = 2r / σ^2
    Sstar = γ * K / (1 + γ)
    V(S) = S ≤ Sstar ? K - S : (K - Sstar) * (S / Sstar)^(-γ)
    return Sstar, V
end

grid(inst) = collect(range(inst.xlo, inst.xhi, length = inst.n))
payoff(inst, x) = max.(inst.K .- exp.(x), 0.0)

"""P overlapping index windows covering 1..n with ≥ olap shared nodes."""
function obstacle_patches(n, P, olap)
    cuts = round.(Int, range(0, n, length = P + 1))
    [(max(1, cuts[p] + 1 - (p > 1 ? olap : 0)), cuts[p + 1]) for p in 1:P]
end

"""Per-patch Grams: element-split weighted stiffness (each interval by
1/#owning-patches — the tensor_spline splitting) + node-split lumped
mass. Σ_p embed(Q_p) = A exactly. `G` transform: Q_p ↦ mass_p + dt·Q_p
for time stepping."""
function patch_grams(inst; dt = nothing)
    x = grid(inst); n = inst.n; h = x[2] - x[1]
    μ = inst.r - inst.σ^2 / 2
    ρ = exp.(2μ .* x ./ inst.σ^2)
    pmid = [0.5 * inst.σ^2 * exp(2μ * (x[e] + h / 2) / inst.σ^2) for e in 1:(n - 1)]
    pat = obstacle_patches(n, inst.P, inst.olap)
    eown = zeros(n - 1); nown = zeros(n)
    for (a, b) in pat
        eown[a:(b - 1)] .+= 1; nown[a:b] .+= 1
    end
    Qs = Matrix{Float64}[]; Ms = Matrix{Float64}[]
    for (a, b) in pat
        m = b - a + 1
        Qp = zeros(m, m); Mp = zeros(m, m)
        for e in a:(b - 1)
            k = pmid[e] / h / eown[e]; i = e - a + 1
            Qp[i, i] += k; Qp[i + 1, i + 1] += k
            Qp[i, i + 1] -= k; Qp[i + 1, i] -= k
        end
        for i in a:b
            Mp[i - a + 1, i - a + 1] = ρ[i] * h / nown[i]
            Qp[i - a + 1, i - a + 1] += inst.r * ρ[i] * h / nown[i]
        end
        push!(Qs, dt === nothing ? Qp : Mp .+ dt .* Qp)
        push!(Ms, Mp)
    end
    return pat, Qs, Ms
end

# ---- builders -------------------------------------------------------------------

"""Premium cells w_p = u_p − ψ_p ∈ PositiveCone (:dirichlet), or u cells
∈ CofreeCone with two PositiveCone corridor slack cells (:game).
Agreement = shared node values on overlaps; boundary pins in g (right
end pinned to the Merton value — the truncation is NOT zero).
For time stepping pass dt and uprev: Q ↦ M + dt·A, c ↦ Gψ − M·uprev."""
function build_obstacle_problem(inst::ObstacleInstance; dt = nothing,
                                uprev::Union{Nothing, Vector{Float64}} = nothing)
    x = grid(inst); n = inst.n
    ψ = payoff(inst, x)
    Sstar, Vf = merton_perpetual(inst.K, inst.r, inst.σ)
    ubc = (ψ[1], Vf(exp(x[end])))
    pat, Qs, Ms = patch_grams(inst; dt)
    P = length(pat)
    game = inst.mode === :game
    ψhi = ψ .+ inst.δ

    row_ids = Int[]; col_ids = Int[]; blks = Matrix{Float64}[]
    place!(e, v, mat) = (push!(row_ids, e); push!(col_ids, v); push!(blks, Matrix{Float64}(mat)))
    ec = Ref(0); gval = Dict{Int, Vector{Float64}}()
    vslo(p) = P + p; vshi(p) = 2P + p
    for p in 1:(P - 1)                                  # overlap agreement
        a1, b1 = pat[p]; a2, b2 = pat[p + 1]
        m = b1 - a2 + 1
        id = (ec[] += 1)
        L = zeros(m, b1 - a1 + 1); R = zeros(m, b2 - a2 + 1)
        for (t, i) in enumerate(a2:b1)
            L[t, i - a1 + 1] = 1.0; R[t, i - a2 + 1] = -1.0
        end
        place!(id, p, L); place!(id, p + 1, R)
        game || (gval[id] = zeros(m))                   # w agree; u agree same
    end
    for (p, (a, b)) in enumerate(pat)                   # boundary pins
        if a == 1
            id = (ec[] += 1); sel = zeros(1, b - a + 1); sel[1] = 1.0
            place!(id, p, sel); gval[id] = [game ? ubc[1] : ubc[1] - ψ[1]]
        end
        if b == n
            id = (ec[] += 1); sel = zeros(1, b - a + 1); sel[end] = 1.0
            place!(id, p, sel); gval[id] = [game ? ubc[2] : ubc[2] - ψ[n]]
        end
    end
    if game                                             # corridor edges
        for (p, (a, b)) in enumerate(pat)
            m = b - a + 1
            id = (ec[] += 1)
            place!(id, p, Matrix{Float64}(I, m, m)); place!(id, vslo(p), -Matrix{Float64}(I, m, m))
            gval[id] = ψ[a:b]
            id = (ec[] += 1)
            place!(id, p, Matrix{Float64}(I, m, m)); place!(id, vshi(p), Matrix{Float64}(I, m, m))
            gval[id] = ψhi[a:b]
        end
    end
    B = blocksparse(row_ids, col_ids, blks)
    g = zeros(size(B, 1))
    for (e, v) in gval
        g[rowrange(B, e)] .= v
    end
    nv = game ? 3P : P
    Q = IPM.allocblockdiag(B); fill!(Q, 0)
    c = zeros(size(B, 2))
    for (p, (a, b)) in enumerate(pat)
        block(Q, p, p, p) .= Qs[p]                       # native PD Q: ε = 0
        if game
            c[colrange(B, p)] .= 0.0
        else
            cp = Qs[p] * ψ[a:b]
            uprev !== nothing && (cp .-= Ms[p] * uprev[a:b])
            c[colrange(B, p)] .= cp
        end
    end
    if game
        for p in (P + 1):(3P)
            d = length(colrange(B, p))
            block(Q, p, p, p) .= 1e-9 .* Matrix{Float64}(I, d, d)   # slack floor
        end
    end
    cones = AbstractCone[]
    for p in 1:nv
        push!(cones, game ? (p ≤ P ? CofreeCone() : PositiveCone()) : PositiveCone())
    end
    prob = IPMProblem(Q, B, c, g, cones)
    info = (pat = pat, wcols = Dict(p => colrange(B, p) for p in 1:P),
            ψ = ψ, x = x, Sstar = Sstar, Vf = Vf, nv = nv)
    return prob, info
end

"""Assemble the global u from patch premiums (or u cells in :game)."""
function obstacle_solution(inst, res, info)
    u = copy(info.ψ)
    for (p, (a, b)) in enumerate(info.pat)
        w = res.p[info.wcols[p]]
        u[a:b] = inst.mode === :game ? w : info.ψ[a:b] .+ w
    end
    return u
end

# ---- minimal-surface mode (SOC leaves; 1D identity cross-check) -----------------

"""min Σ h·t_e, t_e ≥ ‖(1, Δu/h)‖, u ≥ ψ, u(0)=u(1)=0: premium cells +
one 3-dim SecondOrderCone leaf per element (bound first), wired by
[y₁ = 1; y₂ − Δw/h = Δψ/h]. Must reproduce the Dirichlet obstacle
solution exactly in 1D (oracle: 2.8e-7) — SOC leaves vs native Q."""
function build_minimal_surface_problem(; n = 301, amp = 0.45, wid = 8.0, ε = 1e-9)
    x = collect(range(0, 1, length = n)); h = x[2] - x[1]
    ψ = max.(0.0, amp .- wid .* (x .- 0.5) .^ 2)
    row_ids = Int[]; col_ids = Int[]; blks = Matrix{Float64}[]
    place!(e, v, m) = (push!(row_ids, e); push!(col_ids, v); push!(blks, Matrix{Float64}(m)))
    ec = Ref(0); gval = Dict{Int, Vector{Float64}}()
    # vertex 1: premium w (dim n, PositiveCone); vertices 1+e: leaves
    for e in 1:(n - 1)
        id = (ec[] += 1)
        place!(id, 1 + e, [0.0 1.0 0.0; 0.0 0.0 1.0])
        D = zeros(2, n); D[2, e] = -1 / h; D[2, e + 1] = 1 / h
        place!(id, 1, -D)                                # y₂ − Δw/h = Δψ/h
        gval[id] = [1.0, (ψ[e + 1] - ψ[e]) / h]
    end
    id = (ec[] += 1); sel = zeros(1, n); sel[1] = 1.0; place!(id, 1, sel); gval[id] = [0.0]
    id = (ec[] += 1); sel = zeros(1, n); sel[n] = 1.0; place!(id, 1, sel); gval[id] = [0.0]
    B = blocksparse(row_ids, col_ids, blks)
    g = zeros(size(B, 1))
    for (e, v) in gval
        g[rowrange(B, e)] .= v
    end
    Q = IPM.allocblockdiag(B); fill!(Q, 0)
    c = zeros(size(B, 2))
    for v in 1:nvtxs(B)
        d = length(colrange(B, v))
        block(Q, v, v, v) .= ε .* Matrix{Float64}(I, d, d)
        v > 1 && (c[colrange(B, v)[1]] = h)              # cost on t
    end
    cones = AbstractCone[PositiveCone()]
    for _ in 1:(n - 1); push!(cones, SecondOrderCone()); end
    prob = IPMProblem(Q, B, c, g, cones)
    return prob, (x = x, ψ = ψ, wcols = colrange(B, 1))
end

# ---- demos ------------------------------------------------------------------------

obstacle_settings() = IPMSettings{Float64}(kkt = UzawaSettings{Float64}(raug = 1e5))

"""Perpetual put: value vs Merton, free boundary vs patch seams, and
the contact measure read from the cone dual on premium cells."""
function perpetual_put_demo(; kwargs...)
    inst = obstacle_instance(; kwargs...)
    prob, info = build_obstacle_problem(inst)
    res = solve(prob, obstacle_settings())
    u = obstacle_solution(inst, res, info)
    Vex = info.Vf.(exp.(info.x))
    w = u .- info.ψ
    wtol = 1e-4 * maximum(w)
    bidx = findfirst(>(wtol), w)
    seams = [exp(info.x[a]) for (a, _) in info.pat[2:end]]
    @printf("perpetual put: ‖u − Merton‖∞ = %.2e   boundary S = %.4f (S* = %.4f)   it=%d %s\n",
            maximum(abs.(u .- Vex)), exp(info.x[bidx]), info.Sstar, res.niter, res.status)
    @printf("  patch seams at S ≈ %s — the contact set crosses them freely\n",
            string(round.(seams, digits = 3)))
    return inst, res, info
end

"""Finite-horizon American put: implicit-Euler chain of obstacle QPs.
Q = M + Δt·A is FIXED across steps — only c changes (a factorization-
reuse showcase). Oracle: 0.098677 vs CRR(2000) 0.098694."""
function american_put_demo(; T = 1.0, nt = 200, n = 601, kwargs...)
    inst = obstacle_instance(; n, kwargs...)
    dt = T / nt
    u = payoff(inst, grid(inst))
    local res, info
    for m in 1:nt
        prob, info = build_obstacle_problem(inst; dt, uprev = u)
        res = solve(prob, obstacle_settings())
        u = obstacle_solution(inst, res, info)
    end
    V0 = u[argmin(abs.(grid(inst)))]                     # S = 1 node
    @printf("American put T=%.1f: V(S=1) = %.6f   (oracle: 0.098677; CRR 0.098694)\n",
            T, V0)
    return u
end

function game_option_demo(; kwargs...)
    inst = obstacle_instance(; mode = :game, kwargs...)
    prob, info = build_obstacle_problem(inst)
    res = solve(prob, obstacle_settings())
    u = obstacle_solution(inst, res, info)
    lo = info.ψ; hi = info.ψ .+ inst.δ
    @printf("game option δ=%.2f: contact lower %d / upper %d nodes   it=%d %s\n",
            inst.δ, count(u .≤ lo .+ 1e-6), count(u .≥ hi .- 1e-6),
            res.niter, res.status)
    return u
end

function minimal_surface_demo(; n = 301)
    prob, info = build_minimal_surface_problem(; n)
    res = solve(prob, obstacle_settings())
    u = info.ψ .+ res.p[info.wcols]
    @printf("minimal surface: value %.6f   (must equal Dirichlet in 1D, oracle 2.8e-7)\n",
            prob.c' * res.p)
    return u
end

# =====================================================================
# JuMP reference (:dirichlet)
# =====================================================================

using JuMP

function build_jump_obstacle(inst::ObstacleInstance, optimizer)
    x = grid(inst); n = inst.n
    ψ = payoff(inst, x)
    _, Vf = merton_perpetual(inst.K, inst.r, inst.σ)
    pat, Qs, _ = patch_grams(inst)
    model = Model(optimizer); set_silent(model)
    w = [@variable(model, [1:(b - a + 1)], lower_bound = 0) for (a, b) in pat]
    for p in 1:(length(pat) - 1)
        a1, b1 = pat[p]; a2, b2 = pat[p + 1]
        for i in a2:b1
            @constraint(model, w[p][i - a1 + 1] == w[p + 1][i - a2 + 1])
        end
    end
    @constraint(model, w[1][1] == 0.0)
    @constraint(model, w[end][end] == Vf(exp(x[n])) - ψ[n])
    @objective(model, Min,
               sum(0.5 * w[p]' * Qs[p] * w[p] + (Qs[p] * ψ[a:b])' * w[p]
                   for (p, (a, b)) in enumerate(pat)))
    return model, w
end

# =====================================================================
# Benchmark: Sheaf IPM vs Mosek
# =====================================================================

"""Benchmark IPM vs Mosek on obstacle problems.
Dense Q per patch (PDE stiffness) + macroscopic active set."""
function run_obstacle_benchmark(; optimizer = nothing, dual_optimizer = nothing,
                                 solver_name::String = "Mosek", nwarmup::Int = 2, nruns::Int = 5)
    optimizer === nothing && error("Pass optimizer")

    cases = [
        # (n, P, olap, mode, label, raug) — low raug for dense PDE stiffness
        (401, 4, 4, :dirichlet, "perpetual 4-patch", 1e-2),
        (801, 4, 4, :dirichlet, "perpetual 801", 1e-2),
        (801, 8, 4, :dirichlet, "perpetual 8-patch", 1e-2),
        (1601, 8, 8, :dirichlet, "perpetual 1601", 1e-2),
        (801, 4, 4, :game, "game δ=0.02", 1e-2),
    ]

    sname = rpad(solver_name, 9)
    sname_d = rpad(solver_name * "-D", 9)
    println("="^95)
    println("Obstacle/American-Option Benchmark: Sheaf IPM vs $solver_name")
    println("="^95)
    if dual_optimizer !== nothing
        @printf("%-18s %6s %5s %4s %5s %9s %9s %9s %7s %7s\n",
                "Config", "raug", "DOF", "|V|", "H1", "IPM(ms)", sname, sname_d, "P/IPM", "D/IPM")
    else
        @printf("%-18s %6s %5s %4s %5s %9s %9s %8s\n",
                "Config", "raug", "DOF", "|V|", "H1", "IPM(ms)", sname, "speedup")
    end
    println("-"^95)

    for (n, P, olap, mode, label, raug) in cases
        inst = obstacle_instance(; n, P, olap, mode)
        prob, info = build_obstacle_problem(inst)
        dof = size(prob.B, 2)
        h1 = dof - rank(Matrix(prob.B))

        settings = IPMSettings{Float64}(
            kkt = UzawaSettings{Float64}(raug = raug),
            feas_tol = 1e-8, gap_tol = 1e-8, itmax = 200)

        for _ in 1:nwarmup; solve(prob, settings); end
        t_ipm = minimum([@elapsed solve(prob, settings) for _ in 1:nruns])

        for _ in 1:nwarmup
            m, _ = build_jump_obstacle(inst, optimizer); optimize!(m)
        end
        t_mosek = minimum([@elapsed begin
            m, _ = build_jump_obstacle(inst, optimizer); optimize!(m)
        end for _ in 1:nruns])

        if dual_optimizer !== nothing
            for _ in 1:nwarmup
                m, _ = build_jump_obstacle(inst, dual_optimizer); optimize!(m)
            end
            t_dual = minimum([@elapsed begin
                m, _ = build_jump_obstacle(inst, dual_optimizer); optimize!(m)
            end for _ in 1:nruns])
            @printf("%-18s %6.0e %5d %4d %5d %9.1f %9.1f %9.1f %6.2fx %6.2fx\n",
                    label, raug, dof, info.nv, h1, t_ipm * 1000, t_mosek * 1000, t_dual * 1000,
                    t_mosek / t_ipm, t_dual / t_ipm)
        else
            @printf("%-18s %6.0e %5d %4d %5d %9.1f %9.1f %7.2fx\n",
                    label, raug, dof, info.nv, h1, t_ipm * 1000, t_mosek * 1000,
                    t_mosek / t_ipm)
        end
    end

    # Time-stepping benchmark (factorization reuse)
    println("\n--- Time Stepping (200 obstacle QPs, same Q) ---")
    inst = obstacle_instance(; n = 601, P = 4, olap = 4)
    T, nt = 1.0, 200
    dt = T / nt

    # IPM chain
    u = payoff(inst, grid(inst))
    for _ in 1:nwarmup
        prob, info = build_obstacle_problem(inst; dt, uprev = u)
        solve(prob, obstacle_settings())
    end
    u = payoff(inst, grid(inst))
    t_ipm_chain = @elapsed begin
        for _ in 1:nt
            prob, info = build_obstacle_problem(inst; dt, uprev = u)
            res = solve(prob, obstacle_settings())
            u = obstacle_solution(inst, res, info)
        end
    end

    # Mosek chain
    u = payoff(inst, grid(inst))
    t_mosek_chain = @elapsed begin
        for _ in 1:nt
            prob, info = build_obstacle_problem(inst; dt, uprev = u)
            m, w = build_jump_obstacle(inst, optimizer)
            optimize!(m)
            for (p, (a, b)) in enumerate(info.pat)
                u[a:b] = info.ψ[a:b] .+ value.(w[p])
            end
        end
    end

    @printf("American put chain: IPM %.2fs   %s %.2fs   speedup %.2fx\n",
            t_ipm_chain, solver_name, t_mosek_chain, t_mosek_chain / t_ipm_chain)
end

# =====================================================================
# Command-line interface
# =====================================================================

if abspath(PROGRAM_FILE) == @__FILE__
    include(joinpath(@__DIR__, "..", "benchmark_utils.jl"))
    using Dualization

    opts = parse_benchmark_args(ARGS)
    println("Obstacle/American-Option benchmark")
    println("Solver: $(opts.mosek ? "Mosek" : "OSQP (open-source)")")
    println()

    if opts.mosek
        using MosekTools
        optimizer = Mosek.Optimizer
        dual_optimizer = Dualization.dual_optimizer(Mosek.Optimizer)
        solver_name = "Mosek"
    else
        using OSQP
        optimizer = OSQP.Optimizer
        dual_optimizer = Dualization.dual_optimizer(OSQP.Optimizer)
        solver_name = "OSQP"
    end

    run_obstacle_benchmark(; optimizer, dual_optimizer, solver_name,
                           nwarmup = opts.nwarmup, nruns = opts.nruns)
end
