# =============================================================================
# e13.jl — Minimum-fuel planetary soft landing with tracking
#              (dense-Q SOC chain; the dense-Q habitat where the sheaf solver leads)
#
# Usage:
#   julia --project e13.jl              # Clarabel baseline
#   julia --project e13.jl --mosek      # + Mosek (benchmarked primal)
#   julia --project e13.jl --quick      # quick (smaller) sweep
#
# Problem. 3-DoF powered descent, fixed mass, exact-ZOH discretization — the
# convex (SOC/LCvx) relaxation of Açıkmeşe–Ploen minimum-fuel powered-descent
# guidance (JGCD 2007; the G-FOLD line) — plus a dense quadratic TRACKING term:
#
#   min  Σ_t σ_t Δt  +  (γ/2) Σ_{t=1}^{N+1} (x_t − x̄_t)ᵀ W (x_t − x̄_t)
#   s.t. x_{t+1} = A x_t + Bd(u_t + g),  x_1, x_{N+1} pinned,
#        ‖u_t‖ ≤ σ_t,  ρ1 ≤ σ_t ≤ ρ2   (LCvx thrust: the nonconvex ‖u‖ ≥ ρ1
#                                        relaxed through the slack σ),
#        z_t ≥ tan(γ_gs)‖(x_t, y_t)‖    (glideslope cone).
#
# W = inv(Van Loan dispersion of white accelerometer noise through the double
# integrator over τ_c = 8 s): a dense 6×6 with real r–v cross terms (velocity
# error integrates into position error), normalized to unit spectral norm so γ
# is interpretable. Reference x̄: a straight-line descent corridor (constant
# velocity −r0/t_f), deliberately endpoint-inconsistent with the pins.
#
# WHY THIS EXAMPLE EXISTS. The dense-Q habitat: a per-state dense 6×6 Q block on
# an SOC chain (cf. E03/E09). The quadratic folds cheaply into the sheaf solver's
# block Hessian but is a reformulation tax for the conic solvers — so the sheaf
# method LEADS (HSD < Clarabel < Mosek, gap growing with N; see run.jl),
# inverting the Q ≡ 0 ranking where Mosek wins. (The r–v cross terms need
# positions and velocities in ONE stalk — a genuinely dense-Q construction.)
#
# RELAXATION-TIGHTNESS CAVEAT (measured): tracking pulls ‖u‖ below ρ1 at some
# nodes — the LCvx relaxation is NOT tight under tracking (unlike pure min-fuel,
# whose census is 0). The benchmark object is honestly "the convex relaxation
# with tracking"; the census gate asserts the violation count is nonzero.
#
# Validation (tools/tracking_oracle.py + tools/soft_landing_oracle.py, executed):
#   * INVARIANCE: x̄ := the min-fuel optimum simultaneously minimizes both terms
#     ⇒ x(γ) ≡ x̄ ∀γ (exact; matched at oracle accuracy 2.3e-5 m, Δfuel 1.8e-10).
#   * γ → 0 recovers pure min-fuel at rate γ^1.00; the fuel↑ / tracking↓ frontier
#     is monotone in γ; boundary feasibility is objective-independent (statuses
#     flip at t_f^min unchanged by tracking).
#   * MARGIN LAW: the min-fuel optimum is a NEAR-degenerate face — the strict-
#     complementarity margin over ρ1-active arcs collapses as N^−2.00 (measured
#     cross-solver at 1e-11, reproduced to the digit); γ lifts it ×8–×29.
#
# THE AFFINE-IPM STALL. The plain primal-dual IPM STALLS at N ≥ 240 (a limit
# cycle at the near-degenerate min-fuel face); HSD solves every N. The N^−2
# margin law is real but does NOT govern the stall — lifting it ×29 with γ does
# not cure it (refuted). Mechanism open; every formulation-level candidate
# (degeneracy, corner-smoothing, encoding, augmentation) is falsified by
# experiment. HSD is the solver at scale here.
# METROLOGY (owned): boundary quantities (active sets, margins) are measured at
# TIGHT tolerance (1e-11 Clarabel, tight_margin); benchmark-tol iterates sit
# ~1e-4 off the boundary with unidentifiable active sets.
#
# Dial: grid N at fixed t_f = 60 s, γ = γ_base (count dial). DOF = 15N + 3, med
# block 3–4. Physics dials: γ (frontier), t_f (feasibility boundary).
# =============================================================================

using AppleAccelerate
using LinearAlgebra, SparseArrays, Printf

include("utils.jl")

using CellularSheaves.IPM: SecondOrderCone, CofreeCone, PRIMAL_INFEASIBLE
using CellularSheaves.BlockSparseArrays: rowrange

const OPTS = parse_args(ARGS)
const TOL = OPTS.tol
const NRUNS = 5

# Oracle constant (tools/soft_landing_oracle.py): the minimum landing time at N = 40.
const TFMIN_ORACLE_N40 = 35.5376   # s, base 3-D instance

# -----------------------------------------------------------------------------
# Problem definition
# -----------------------------------------------------------------------------

struct LandingInstance
    r0::Vector{Float64}; v0::Vector{Float64}   # initial position, velocity [m, m/s]
    rho1::Float64; rho2::Float64               # thrust-accel bounds [m/s²]
    gmag::Float64                              # gravity magnitude [m/s²]
    tan_gs::Float64                            # tan(glideslope angle)
end

# Physical constants after Luo–Echigo–Açıkmeşe (2024) §6.1 / Sostaric–Rea
# (2005): m = 1000 kg, T_max = 7500 N, range [0.3, 0.8] T_max, lunar gravity.
# γ_gs = 55°: the unconstrained approach dips to 51.6°, so the cone is live.
landing_instance(; r0 = [600.0, -400.0, 2000.0], v0 = [30.0, -10.0, -60.0],
                 rho1 = 2.25, rho2 = 6.0, gmag = 1.62, gs_deg = 55.0) =
    LandingInstance(r0, v0, rho1, rho2, gmag, tand(gs_deg))

const TF_BASE = 60.0
const N_BASE = 60

"Exact ZOH discretization of the double integrator."
function zoh(tf, N)
    dt = tf / N
    A = [I(3) dt * I(3); zeros(3, 3) I(3)]
    Bd = [dt^2 / 2 * I(3); dt * I(3)]
    return dt, Matrix{Float64}(A), Matrix{Float64}(Bd)
end

gvec(inst) = [0.0, 0.0, -inst.gmag]

# -----------------------------------------------------------------------------
# Build (state chain + thrust SOC + bound NOC + glideslope SOC)
# -----------------------------------------------------------------------------
#
# Stalk order: state 1..N+1 (dim 6) | thrust N+1+k, k = 1..N (dim 4) |
# bounds 2N+1+k, k = 1..N (dim 2) | glideslope 3N+1+(t−1), t = 2..N (dim 3)
# | [mode = :dev only] deviation stalk (dim 7, (t0, e) ∈ SOC(7)) last.
#
# mode = :fuel — the landing problem (terminal state pinned to 0);
# mode = :dev  — the minimum-landing-error / distance-to-feasibility problem
#                (Blackmore–Açıkmeşe–Scharf 2010): terminal pin replaced by
#                e = x_{N+1} into a SOC(7) epigraph, min t0. Always feasible;
#                optimal value 0 iff the :fuel problem is feasible. Note the
#                free terminal node carries no glideslope stalk (ties run
#                t = 2..N); the oracle confirms the deviation value is
#                unchanged by this at the base instance (8.6e-7).

function build_landing(inst::LandingInstance, tf, N; glideslope = true, mode = :fuel)
    dt, A, Bd = zoh(tf, N)
    vstate(t) = t                           # t = 1..N+1
    vthrust(k) = N + 1 + k                  # k = 1..N
    vbnd(k) = 2N + 1 + k                    # k = 1..N
    vgs(t) = 3N + 1 + (t - 1)               # t = 2..N  (only if glideslope)
    vdev = (glideslope ? 4N : 3N + 1) + 1   # mode = :dev only; vertex ids are
                                            # contiguous in either case

    I6 = Matrix{Float64}(I, 6, 6)
    row_ids, col_ids, blocks = Int[], Int[], Matrix{Float64}[]
    rhs_val = Dict{Int, Vector{Float64}}()
    e = 0

    # ---- dynamics: x_{t+1} − A x_t − Bd u_t = Bd g_vec (6 dense rows each)
    Bu = zeros(6, 4); Bu[:, 2:4] .= -Bd
    for t in 1:N
        e += 1
        push!(row_ids, e); push!(col_ids, vstate(t)); push!(blocks, -copy(A))
        push!(row_ids, e); push!(col_ids, vstate(t + 1)); push!(blocks, copy(I6))
        push!(row_ids, e); push!(col_ids, vthrust(t)); push!(blocks, copy(Bu))
        rhs_val[e] = Bd * gvec(inst)
    end

    # ---- initial pin: x_1 = (r0, v0)
    e += 1
    push!(row_ids, e); push!(col_ids, vstate(1)); push!(blocks, copy(I6))
    rhs_val[e] = vcat(inst.r0, inst.v0)

    # ---- terminal: pin (fuel) or epigraph tie (dev)
    e += 1
    if mode == :fuel
        push!(row_ids, e); push!(col_ids, vstate(N + 1)); push!(blocks, copy(I6))
        rhs_val[e] = zeros(6)
    else
        Ad = zeros(6, 7); Ad[:, 2:7] .= I6
        push!(row_ids, e); push!(col_ids, vdev); push!(blocks, Ad)
        push!(row_ids, e); push!(col_ids, vstate(N + 1)); push!(blocks, -copy(I6))
        rhs_val[e] = zeros(6)
    end

    # ---- thrust bounds: σ_k − a_k = ρ1, σ_k + b_k = ρ2 (NOC slacks)
    As = [1.0 0.0 0.0 0.0; 1.0 0.0 0.0 0.0]
    Ab = [-1.0 0.0; 0.0 1.0]
    for k in 1:N
        e += 1
        push!(row_ids, e); push!(col_ids, vthrust(k)); push!(blocks, copy(As))
        push!(row_ids, e); push!(col_ids, vbnd(k)); push!(blocks, copy(Ab))
        rhs_val[e] = [inst.rho1, inst.rho2]
    end

    # ---- glideslope ties: w0 = z_t / tan, (w1, w2) = (x_t, y_t), t = 2..N
    if glideslope
        Ax = zeros(3, 6)
        Ax[1, 3] = -1.0 / inst.tan_gs; Ax[2, 1] = -1.0; Ax[3, 2] = -1.0
        for t in 2:N
            e += 1
            push!(row_ids, e); push!(col_ids, vgs(t)); push!(blocks, Matrix{Float64}(I, 3, 3))
            push!(row_ids, e); push!(col_ids, vstate(t)); push!(blocks, copy(Ax))
            rhs_val[e] = zeros(3)
        end
    end

    B = blocksparse(row_ids, col_ids, blocks)

    g_rhs = zeros(size(B, 1))
    for (edge, val) in rhs_val
        g_rhs[rowrange(B, edge)] .= val
    end

    # ---- Q ≡ 0 (density lives in the maps — the E05–E07 cell); fuel in c
    Q = IPM.allocblockdiag(B); fill!(Q, 0)
    c = zeros(size(B, 2))
    if mode == :fuel
        for k in 1:N
            c[first(colrange(B, vthrust(k)))] = -dt   # min Σ σ_k Δt
        end
    else
        c[first(colrange(B, vdev))] = -1.0            # min t0 (deviation ONLY)
    end

    conekinds = vcat(fill(:free, N + 1), fill(:soc, N), fill(:noc, N),
                     glideslope ? fill(:soc, N - 1) : Symbol[],
                     mode == :dev ? [:soc] : Symbol[])
    K_cones = IPM.AbstractCone[k == :free ? CofreeCone() :
                               k == :noc ? PositiveCone() :
                               SecondOrderCone() for k in conekinds]
    prob = IPMProblem(Q, B, c, g_rhs, K_cones)

    ctx = (; inst, tf, N, dt, A, Bd, mode, glideslope,
           conekinds, vstate, vthrust, vbnd, vgs, vdev)
    return prob, ctx
end

"State trajectory (N+1) × 6 from the solution vector."
function extract_x(prob, ctx, p::AbstractVector)
    X = zeros(ctx.N + 1, 6)
    for t in 1:ctx.N + 1
        X[t, :] .= p[colrange(prob.B, ctx.vstate(t))]
    end
    return X
end

"Thrust vectors N × 3 and slacks σ (length N) from the solution vector."
function extract_u(prob, ctx, p::AbstractVector)
    U = zeros(ctx.N, 3); s = zeros(ctx.N)
    for k in 1:ctx.N
        blk = p[colrange(prob.B, ctx.vthrust(k))]
        s[k] = blk[1]; U[k, :] .= blk[2:4]
    end
    return U, s
end

# -----------------------------------------------------------------------------
# Baseline (same conic program via JuMP)
# -----------------------------------------------------------------------------

function build_jump_landing(prob, ctx, optimizer)
    model = Model(optimizer)
    set_silent(model)
    nv = nvtxs(prob.B)
    xs = Vector{Vector{VariableRef}}(undef, nv)
    for v in 1:nv
        d = length(colrange(prob.B, v))
        xv = @variable(model, [1:d])
        kind = ctx.conekinds[v]
        if kind == :soc
            @constraint(model, xv in JuMP.SecondOrderCone())
        elseif kind == :noc
            @constraint(model, xv .>= 0)
        end
        xs[v] = xv
    end
    x = reduce(vcat, xs)
    Qs = sparse(prob.Q); Bs = sparse(prob.B)
    @objective(model, Min, 0.5 * x' * Qs * x - prob.c' * x)
    @constraint(model, Bs * x .== prob.g)
    return model, x
end

# -----------------------------------------------------------------------------
# Farkas certificate quality (HSD, PRIMAL_INFEASIBLE)
# -----------------------------------------------------------------------------
#
# For (P) min ½pᵀQp − cᵀp s.t. Bp = g, p ∈ K, a primal-infeasibility
# certificate is y with −Bᵀy ∈ K* and gᵀy > 0. The HSD result returns y
# normalized (||y|| = 1) with pobj = gᵀy (result/hsd.jl C5 carve-out). Dual
# cones here: (CofreeCone)* = {0}, NOC and SOC self-dual. Returns
# (gy, free_res, cone_res): gy must be positive, both residuals ~ feas_tol.

sizes() = OPTS.tiny ? [120] :
    OPTS.quick ? [240, 480, 960] :
    [240, 480, 960, 1920]

# -----------------------------------------------------------------------------
# Tracking metric, reference corridor, base weight
# -----------------------------------------------------------------------------

"W = inv(Van Loan P + eps·I), unit spectral norm. Dense by nature (r–v
cross terms from integrating accelerometer noise)."
function dispersion_W(; tau_c = 8.0, q = 1.0, eps_reg = 1e-3)
    Ac = [zeros(3, 3) I(3); zeros(3, 6)]
    L = [zeros(3, 3); I(3)]
    M = [-Ac (q .* (L * L')); zeros(6, 6) Ac']
    F = exp(Matrix{Float64}(M) .* tau_c)
    P = F[7:12, 7:12]' * F[1:6, 7:12]
    P = 0.5 .* (P .+ P') .+ eps_reg .* I(6)
    W = inv(P)
    W = 0.5 .* (W .+ W')
    return W ./ opnorm(W, 2)
end

const W_BASE = dispersion_W()
const GAMMA_BASE = 2e-3

# POLICY (2026-07-25): Mosek is benchmarked in PRIMAL form — measured as its best
# form here (primal ~2x faster than the dualized problem via Dualization.jl; the
# earlier dual-only policy was based on an un-remeasured guess and is retired). If the
# problem class changes materially, re-measure the dual once before trusting this choice.

"Straight-line descent corridor: r linear r0 → 0, v = −r0/t_f constant.
Deliberately endpoint-inconsistent with the pins (v0 ≠ v̄ ≠ 0)."
function corridor(inst::LandingInstance, tf, N)
    xbar = zeros(N + 1, 6)
    vbar = -inst.r0 ./ tf
    for t in 1:N + 1
        xbar[t, 1:3] .= inst.r0 .* (1 - (t - 1) / N)
        xbar[t, 4:6] .= vbar
    end
    return xbar
end

# -----------------------------------------------------------------------------
# Build: the min-fuel constraint set (build_landing) + a post-build dense-Q/c mutation
# -----------------------------------------------------------------------------
#
# min ½pᵀQp − cᵀp convention: tracking ½γ(x−x̄)ᵀW(x−x̄) contributes
# Q_state = γW and c_state = +γWx̄_t (plus a constant ½γ Σ x̄ᵀWx̄ NOT in the
# solver objective — the assembly-identity gate accounts for it exactly).

function build_tracking(inst::LandingInstance, tf, N; gamma = GAMMA_BASE,
                        xbar = nothing, W = W_BASE, glideslope = true)
    xbar === nothing && (xbar = corridor(inst, tf, N))
    prob, ctx = build_landing(inst, tf, N; glideslope)
    for t in 1:N + 1
        v = ctx.vstate(t)
        block(prob.Q, v, v, v) .= gamma .* W
        prob.c[colrange(prob.B, v)] .+= gamma .* (W * xbar[t, :])
    end
    ctx2 = merge(ctx, (; gamma, xbar, W))
    return prob, ctx2
end

"Objective constant dropped by the solver form: ½γ Σ_t x̄_tᵀ W x̄_t."
objective_const(ctx) =
    0.5 * ctx.gamma * sum(dot(ctx.xbar[t, :], ctx.W * ctx.xbar[t, :])
                          for t in 1:ctx.N + 1)

"Direct evaluation of fuel + tracking (assembly check)."
function tracking_objective(ctx, X, s)
    tr = sum(dot(X[t, :] .- ctx.xbar[t, :], ctx.W * (X[t, :] .- ctx.xbar[t, :]))
             for t in 1:ctx.N + 1)
    return ctx.dt * sum(s) + 0.5 * ctx.gamma * tr
end

"Tracking cost alone (frontier metric)."
tracking_cost(ctx, X) =
    sum(dot(X[t, :] .- ctx.xbar[t, :], ctx.W * (X[t, :] .- ctx.xbar[t, :]))
        for t in 1:ctx.N + 1)

"Strict-complementarity margin at a SOLVED ITERATE. WARNING: at benchmark
tolerance the iterate sits ~1e-4 off the boundary and the active set is
unidentifiable — this function is only valid on tightly-solved iterates.
Gates use tight_margin (1e-11 Clarabel) instead."
function compl_margin(prob, ctx, res; tol = 1e-6)
    margin = Inf; both = 0; nact = 0
    for k in 1:ctx.N
        seg = colrange(prob.B, ctx.vbnd(k))
        a = res.p[seg][1]; λ = res.d[seg][1]
        if a < tol * (1 + ctx.inst.rho1)
            nact += 1
            margin = min(margin, a + λ)
            λ < tol && (both += 1)
        end
    end
    return (margin = margin, both = both, nact = nact)
end

"JuMP build that also returns the per-vertex NOC constraint refs (for
tight-tolerance dual extraction)."
function build_jump_with_bnd_refs(prob, ctx, optimizer)
    model = Model(optimizer)
    set_silent(model)
    nv = nvtxs(prob.B)
    xs = Vector{Vector{VariableRef}}(undef, nv)
    bnd_cons = Dict{Int, Any}()
    for v in 1:nv
        d = length(colrange(prob.B, v))
        xv = @variable(model, [1:d])
        kind = ctx.conekinds[v]
        if kind == :soc
            @constraint(model, xv in JuMP.SecondOrderCone())
        elseif kind == :noc
            bnd_cons[v] = @constraint(model, xv .>= 0)
        end
        xs[v] = xv
    end
    x = reduce(vcat, xs)
    Qs = sparse(prob.Q); Bs = sparse(prob.B)
    @objective(model, Min, 0.5 * x' * Qs * x - prob.c' * x)
    @constraint(model, Bs * x .== prob.g)
    return model, xs, bnd_cons
end

"Margin metrology done right: TIGHT Clarabel solve (1e-11), active set and
duals read at the near-exact optimum (adjudicated: at 1e-8 the
active set is empty at threshold 1e-6; at 1e-11 the N^-2.00 law reproduces
to the digit)."
function tight_margin(inst, tf, N; gamma = 0.0, tol = 1e-11)
    prob, ctx = gamma == 0.0 ? build_landing(inst, tf, N) :
                               build_tracking(inst, tf, N; gamma)
    model, xs, bnd = build_jump_with_bnd_refs(prob, ctx, clarabel_opt(; tol))
    optimize!(model)
    margin = Inf; nact = 0; both = 0
    for k in 1:N
        v = ctx.vbnd(k)
        a = value(xs[v][1])
        lam = dual(bnd[v][1])
        if a < 1e-6 * (1 + inst.rho1)
            nact += 1
            margin = min(margin, a + lam)
            lam < 1e-6 && (both += 1)
        end
    end
    return (; margin, nact, both)
end

# -----------------------------------------------------------------------------
# Gate tests (oracle figures in messages; tools/tracking_oracle.py)
# -----------------------------------------------------------------------------

function test_assembly_identity(prob, ctx, res)
    X = extract_x(prob, ctx, res.p)
    _, s = extract_u(prob, ctx, res.p)
    o_direct = tracking_objective(ctx, X, s)
    o_ipm = ipm_objective(prob, res) + objective_const(ctx)
    rel = abs(o_direct - o_ipm) / abs(o_direct)
    @assert rel < 1e-5 "assembly mismatch: $o_direct vs $o_ipm"
    println("  [PASS] objective identity: F(x) = $(round(o_direct, digits = 4)) ",
        "(rel $(round(rel, sigdigits = 2)); constant term accounted)")
    return X
end

function test_invariance(inst, hsd_settings)
    # EXACT structural property: x̄ := the E13 optimum minimizes fuel AND the tracking term
    # (which is 0 at x̄), so x(γ) ≡ x̄. Solved TIGHT with Clarabel (1e-10) — at benchmark HSD
    # tolerance min-fuel's flat optimal face and the large-γ dense Q hide it (cf. tight_margin;
    # metrology: the property is real, the benchmark-tol iterate is the artifact).
    cla = clarabel_opt(; tol = 1e-10)
    prob14, ctx14 = build_landing(inst, TF_BASE, N_BASE)
    m14, xv14 = build_jump_landing(prob14, ctx14, cla)
    optimize!(m14)
    X14 = extract_x(prob14, ctx14, value.(xv14))
    fuel14 = objective_value(m14)
    for gamma in (0.1, 10.0)
        prob, ctx = build_tracking(inst, TF_BASE, N_BASE; gamma, xbar = X14)
        m, xv = build_jump_landing(prob, ctx, cla)
        optimize!(m)
        X = extract_x(prob, ctx, value.(xv))
        _, s = extract_u(prob, ctx, value.(xv))
        dx = maximum(abs.(X .- X14))
        df = abs(ctx.dt * sum(s) - fuel14) / fuel14
        @assert dx < 1e-3 "invariance broken at γ=$gamma: ‖Δx‖ $dx (oracle: 2.3e-5)"
        @assert df < 1e-6 "invariance fuel drift at γ=$gamma: $df (oracle: 1.8e-10)"
    end
    println("  [PASS] invariance (tight Clarabel): x̄ := E13 optimum ⇒ x(γ) ≡ x̄ at γ ∈ {0.1, 10} ",
        "(simultaneous minimizer; exact structural gate)")
    return X14
end

function test_gamma0_limit(inst, X14, hsd_settings)
    # tight Clarabel (X14 is the tight-Clarabel optimum from test_invariance; the γ→0 trajectory
    # limit only converges cleanly at oracle accuracy — the flat fuel face otherwise adds tol noise)
    cla = clarabel_opt(; tol = 1e-10)
    errs = Float64[]
    for gamma in (1e-4, 1e-6)
        prob, ctx = build_tracking(inst, TF_BASE, N_BASE; gamma)
        m, xv = build_jump_landing(prob, ctx, cla)
        optimize!(m)
        push!(errs, maximum(abs.(extract_x(prob, ctx, value.(xv)) .- X14)))
    end
    @assert errs[2] < errs[1] / 50 "γ→0 rate broken: $errs (oracle: rate γ^1.00, ×100)"
    println("  [PASS] γ → 0 recovers E13: ‖Δx‖∞ $(round(errs[1], sigdigits = 2)) ",
        "→ $(round(errs[2], sigdigits = 2)) over γ = 1e-4 → 1e-6 (oracle rate γ^1.00)")
end

function test_frontier(inst, hsd_settings)
    fuels, trs = Float64[], Float64[]
    for gamma in (2e-4, 2e-3, 2e-2)
        prob, ctx = build_tracking(inst, TF_BASE, N_BASE; gamma)
        res = solve(prob, hsd_settings)
        @assert res.status in (OPTIMAL, NEAR_OPTIMAL)
        X = extract_x(prob, ctx, res.p)
        _, s = extract_u(prob, ctx, res.p)
        push!(fuels, ctx.dt * sum(s)); push!(trs, tracking_cost(ctx, X))
    end
    @assert issorted(fuels) "fuel not monotone: $fuels (oracle: 182.6 → 230.8)"
    @assert issorted(trs; rev = true) "tracking not monotone: $trs"
    println("  [PASS] frontier: fuel $(round(fuels[1], digits = 1)) → ",
        "$(round(fuels[3], digits = 1)) up, tracking ",
        "$(round(trs[1], digits = 0)) → $(round(trs[3], digits = 0)) down ",
        "(γ = 2e-4 → 2e-2; oracle: 182.6 → 230.8, 173237 → 61636)")
end

function test_margin_lift(inst)
    # TIGHT metrology only: 1e-11 Clarabel; at benchmark tolerance the active
    # set is empty (iterate ~1e-4 off the boundary) and this gate is invalid.
    m0 = tight_margin(inst, TF_BASE, N_BASE)
    mg = tight_margin(inst, TF_BASE, N_BASE; gamma = GAMMA_BASE)
    @assert m0.nact > 20 "active set not identified — tolerance? nact = $(m0.nact) (oracle: 46)"
    @assert mg.margin > 5 * m0.margin "margin lift broken: $(m0.margin) → $(mg.margin) (oracle: ×8)"
    @assert m0.both == 0 && mg.both == 0 "exact degeneracy appeared: $(m0.both), $(mg.both) (oracle: 0)"
    println("  [PASS] margin lift (tight, 1e-11): min(a + λ) $(round(m0.margin, sigdigits = 3)) → ",
        "$(round(mg.margin, sigdigits = 3)) at γ_base (×$(round(mg.margin / m0.margin, digits = 1)); ",
        "oracle ×8, $(m0.nact) active arcs); doubly-degenerate 0/0. NOTE: the ",
        "margin LAW is real but does NOT govern the affine stall (refuted)")
end

function test_onset_scaling(inst)
    # TIGHT metrology only; cross-solver adjudicated: at 1e-11
    # the law reproduces to the digit; at 1e-8 it is invisible.
    Ns = (60, 120, 240)
    ms = [tight_margin(inst, TF_BASE, N).margin for N in Ns]
    slope = log(ms[end] / ms[1]) / log(Ns[end] / Ns[1])
    @assert -2.4 < slope < -1.6 "onset law broken: N^$slope (oracle: N^-2.00 at 1e-11)"
    println("  [PASS] onset law (tight, 1e-11): γ = 0 margin ~ N^$(round(slope, digits = 2)) ",
        "over N = 60 → 240 (oracle: N^-2.00; ",
        join(round.(ms, sigdigits = 3), ", "), "). Real law of the exact optimum; ",
        "NOT the stall mechanism (refuted)")
end

function test_stall_outcome(inst, ipm_settings, hsd_settings)
    # HISTORY: this gate originally ASSERTED the cure prediction (margin
    # lift ⇒ affine convergence at N = 480, γ_base). First run 2026-07-24:
    # FALSIFIED — STALLED at BOTH γ = 0 and γ_base (adjudicated;
    # the margin law is real but not causal). The gate now PINS the refuted
    # state, ledger-style: it flips loudly if the behavior ever changes.
    N = OPTS.quick ? 240 : 480
    probg, _ = build_tracking(inst, TF_BASE, N)
    res_i = solve(probg, ipm_settings)
    res_h = solve(probg, hsd_settings)
    @assert res_h.status in (OPTIMAL, NEAR_OPTIMAL) "HSD on E13 broke: $(res_h.status)"
    prob0, _ = build_landing(inst, TF_BASE, N)
    res_i0 = solve(prob0, ipm_settings)
    if res_i.status in (OPTIMAL, NEAR_OPTIMAL)
        println("  [FLIP] stall outcome: affine IPM now solves E13 at N = $N ",
            "($(res_i.status), $(res_i.niter) iters; γ = 0: $(res_i0.status)) — ",
            "a fix landed or the mechanism shifted; update this gate")
    else
        println("  [PASS] stall outcome pinned: affine IPM $(res_i.status) at ",
            "γ_base, $(res_i0.status) at γ = 0 (N = $N; refuted cure recorded); ",
            "HSD $(res_h.status) in $(res_h.niter) iters")
    end
end

function test_lcvx_census_tracking(prob, ctx, res)
    # Reads OUR HSD solution (a solver-facing gate). The robust, solver-independent claim is that
    # the LCvx relaxation is NOT tight under tracking (0 < violations < N); the exact count is
    # tol/solver-dependent (cvxpy oracle 11, HSD ~33 near the ‖u‖ = ρ1 boundary) so it is not gated.
    U, s = extract_u(prob, ctx, res.p)
    unorm = [norm(U[k, :]) for k in 1:ctx.N]
    bad = count(unorm .< ctx.inst.rho1 * (1 - 1e-6))
    @assert 0 < bad < ctx.N "LCvx census out of band: $bad of $(ctx.N)"
    println("  [PASS] LCvx census (the tightness caveat, measured on our solution): $bad of ",
        "$(ctx.N) nodes with ‖u‖ < ρ1 at γ_base (cvxpy oracle 11; count is tol-dependent) — ",
        "the relaxation is NOT tight on E13; the benchmark object is the relaxation")
end

function test_tracking_vs_clarabel(prob, ctx, res, X_ipm)
    # Solver-facing gate on OUR HSD solve. Gate the robust invariant (the objective); report the
    # trajectory (tol-limited where the weakly-tracked fuel face is near-flat, cf. test_invariance).
    model, xv = build_jump_landing(prob, ctx, clarabel_opt(; tol = TOL))
    optimize!(model)
    dx = maximum(abs.(X_ipm .- extract_x(prob, ctx, value.(xv))))
    rel = abs(ipm_objective(prob, res) - objective_value(model)) / (1 + abs(objective_value(model)))
    # 5-digit agreement: on the near-flat weakly-tracked face both solvers' objective accuracy is
    # ~1e-5 at TOL = 1e-8 (below tol), so gate at 1e-4 — the correctness check is that HSD lands
    # on the same optimum as Clarabel, not that either hits machine precision here.
    @assert rel < 1e-4 "HSD vs Clarabel objective mismatch: rel $rel (dense-Q tracking)"
    println("  [PASS] HSD vs Clarabel (same program, dense Q): objective agrees to rel ",
        "$(round(rel, sigdigits = 2)); ‖ΔX‖∞ = $(round(dx, sigdigits = 2)) m")
end

function test_boundary_invariance(inst, hsd_settings)
    for (f, want) in ((0.98, PRIMAL_INFEASIBLE), (1.02, OPTIMAL))
        tf = f * TFMIN_ORACLE_N40
        prob, _ = build_tracking(inst, tf, 40; xbar = corridor(inst, tf, 40))
        res = solve(prob, hsd_settings)
        ok = want == OPTIMAL ? res.status in (OPTIMAL, NEAR_OPTIMAL) :
             res.status == PRIMAL_INFEASIBLE
        @assert ok "boundary moved under tracking at $f·tf_min: $(res.status)"
    end
    println("  [PASS] boundary invariance: feasibility is objective-independent; ",
        "statuses flip with E13 at 0.98/1.02·tf_min")
end

# -----------------------------------------------------------------------------
# Sweep
# -----------------------------------------------------------------------------

function run()
    println("\n", "=" ^ 78)
    println("  E13: soft landing with tracking, dense-Q ",
        "(dial: N at tf = $(TF_BASE), γ = $(GAMMA_BASE))")
    println("=" ^ 78)

    ipm_settings = IPMSettings{Float64}(
        feas_tol = TOL, gap_tol = TOL, itmax = 200)
    hsd_settings = HSDSettings{Float64}(
        feas_tol = TOL, gap_tol = TOL, itmax = 200)

    inst = landing_instance()

    println("\n  Gate tests (base instance, N = $(N_BASE), γ = $(GAMMA_BASE)):")
    prob, ctx = build_tracking(inst, TF_BASE, N_BASE)
    res = solve(prob, hsd_settings)
    @assert res.status in (OPTIMAL, NEAR_OPTIMAL) "HSD base: $(res.status)"
    X = test_assembly_identity(prob, ctx, res)
    X14 = test_invariance(inst, hsd_settings)
    test_gamma0_limit(inst, X14, hsd_settings)
    test_frontier(inst, hsd_settings)
    test_margin_lift(inst)
    test_onset_scaling(inst)
    test_lcvx_census_tracking(prob, ctx, res)
    test_tracking_vs_clarabel(prob, ctx, res, X)
    test_boundary_invariance(inst, hsd_settings)
    OPTS.tiny || test_stall_outcome(inst, ipm_settings, hsd_settings)
    println()

    cla_opt = clarabel_opt(; tol = TOL)
    msk_opt = OPTS.mosek ? mosek_opt(; tol = TOL) : nothing

    rows = []
    for N in sizes()
        probN, ctxN = build_tracking(inst, TF_BASE, N)
        stats = problem_stats(probN)

        @printf("  N=%-5d dof=%-6d n1=%-5d blk=%-4.0f  ",
            N, stats.N0, stats.N1, stats.med_block)

        m_ipm = measure_ipm(probN, ipm_settings; nruns = NRUNS)
        m_hsd = measure_ipm(probN, hsd_settings; nruns = NRUNS)
        m_cla = measure_jump(() -> first(build_jump_landing(probN, ctxN, cla_opt)); nruns = NRUNS)
        # Mosek benchmarked in PRIMAL form: ~2x faster than dualized here.
        m_msk = msk_opt !== nothing ?
            measure_jump(() -> first(build_jump_landing(probN, ctxN, mosek_opt(; tol = TOL))); nruns = NRUNS) :
            (t = NaN, status = "", obj = NaN)

        ratio(b, base) = isfinite(b.t) && isfinite(base.t) ? b.t / base.t : NaN
        fmt_ratio(b, base) = isnan(ratio(b, base)) ? "—" : @sprintf("%.2fx", ratio(b, base))

        println("IPM ", fmt_time(m_ipm), "  HSD ", fmt_time(m_hsd), " (", fmt_ratio(m_hsd, m_ipm),
            ")  Cla ", fmt_time(m_cla), " (", fmt_ratio(m_cla, m_ipm),
            ")  Msk ", fmt_time(m_msk), " (", fmt_ratio(m_msk, m_ipm), ")")

        push!(rows, (size = N, dof = stats.N0, ipm = m_ipm, hsd = m_hsd, cla = m_cla,
                     msk = m_msk))
    end

    dofs = [r.dof for r in rows]
    println("\nFitted log-log slopes (time vs DOF):")
    for (name, get_t) in [("IPM", r -> r.ipm.t), ("HSD", r -> r.hsd.t), ("Clarabel", r -> r.cla.t),
                          ("Mosek", r -> r.msk.t)]
        sl = loglog_slope(dofs, [get_t(r) for r in rows])
        print("  $name: ", isnan(sl) ? "n/a" : @sprintf("DOF^%.2f", sl))
    end
    println()
end

if abspath(PROGRAM_FILE) == @__FILE__
    run()
end

# =============================================================================
# Sample run: 2026-07-25 (--quick --mosek). Mosek is benchmarked PRIMAL —
# its best form on the dense-Q tracking objective (primal ~2x faster than the
# dualized problem, and far better-behaved at scale; the earlier dual-only
# policy was an un-remeasured guess and is retired).
# Ratios below are vs HSD (the affine IPM stalls at N >= 240 — see def.jl banner).
# -----------------------------------------------------------------------------
# N=240   dof=3603    IPM —        HSD   27.9ms  Cla   32.8ms (1.18x)  Msk   47.5ms (1.70x)
# N=480   dof=7203    IPM —        HSD   57.5ms  Cla   73.3ms (1.27x)  Msk  100.1ms (1.74x)
# N=960   dof=14403   IPM —        HSD  114.4ms  Cla  157.8ms (1.38x)  Msk  197.1ms (1.72x)
# Slopes: HSD DOF^1.02  Clarabel DOF^1.13  Mosek DOF^1.03
#
# Dense-Q verdict: HSD < Clarabel < Mosek at every N, gap steady (~1.7x for
# Mosek). The quadratic tracking Q folds cheaply into the block Hessian for the
# sheaf solver but is a reformulation tax for the conic solvers. This INVERTS the
# Q ≡ 0 ranking, where Mosek led. (Dualizing Mosek was tried and rejected — it
# was competitive at N=240 but blew up to 454ms by N=960, slope 1.54.) The affine
# IPM stalls past N=240 (a near-degenerate-face limit cycle; mechanism open,
# candidates refuted by experiment — see def.jl banner). HSD is the solver at scale here.
# =============================================================================
