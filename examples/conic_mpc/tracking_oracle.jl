#
# Conic tracking oracle (M2a) — the conic replacement for the learned-MPC QP expert.
#
# The learned_mpc policy today imitates the linear harmonic readout u* = E_u·x_B, which
# minimizes the windowed sheaf energy (tracking + consensus + control cost) with NO input
# constraint — so it can and does command thrust beyond the actuator limit. This oracle
# keeps the same windowed position-tracking objective but adds a hard second-order-cone
# actuator ball ‖u_t‖₂ ≤ ūcap, solved by SheafSDP's interior-point method. The first
# control u_i^0 of the solve is the label the conic policy will learn.
#
# Planar-quadrotor dynamics (nx=6, nu=2), matching train/train_neural_solver.jl. This
# first cut is per-agent (independent tracking); soft sheaf consensus is a later refinement.
#
# Env: examples/conic_mpc (SheafSDP + BlockSparseArrays). No CellularSheaves dependency —
# ZOH discretization is done here via the matrix exponential.
#

module TrackingOracle

using LinearAlgebra
using SheafSDP
using CommonSolve: solve
using BlockSparseArrays: colrange, rowrange, blocksparse, block

export quad_dynamics, QuadTrack, conic_track_control, unconstrained_track_control

const NX = 6
const NU = 2
const PIDX = [1, 2]                      # (y, z) position indices tracked / consensus-coupled

# Planar-quadrotor continuous dynamics → ZOH discrete (matches the learned-MPC setup).
function quad_dynamics(; g=9.81, m=0.5, Iq=0.01, ell=0.25, h=0.05)
    Ac = [0.0 0.0 0.0 1.0 0.0 0.0;
          0.0 0.0 0.0 0.0 1.0 0.0;
          0.0 0.0 0.0 0.0 0.0 1.0;
          0.0 0.0 -g  0.0 0.0 0.0;
          0.0 0.0 0.0 0.0 0.0 0.0;
          0.0 0.0 0.0 0.0 0.0 0.0]
    Bc = [0.0 0.0; 0.0 0.0; 0.0 0.0; 0.0 0.0; 1/m 1/m; ell/(2Iq) -ell/(2Iq)]
    Mexp = exp([Ac Bc; zeros(NU, NX + NU)] .* h)
    Ad = Mexp[1:NX, 1:NX]
    Bd = Mexp[1:NX, NX+1:NX+NU]
    return Ad, Bd
end

"""
    QuadTrack(Ad, Bd, W; ρ, P)

Window tracking-oracle config: horizon `W`, control regularization `ρ`, position
selector `P` (the (y,z) rows). `P` defaults to picking indices `PIDX`.
"""
struct QuadTrack
    Ad::Matrix{Float64}
    Bd::Matrix{Float64}
    W::Int
    ρ::Float64
    P::Matrix{Float64}
end
function QuadTrack(Ad, Bd, W; ρ=1e-3)
    P = zeros(length(PIDX), NX)
    for (r, c) in enumerate(PIDX); P[r, c] = 1.0; end
    return QuadTrack(Ad, Bd, W, ρ, P)
end

# ---------------------------------------------------------------------------
# Build the per-agent windowed conic tracking program.
#   min Σ_t ‖P x^t − P τ^t‖² + ρ‖u^t‖²   s.t.  x^0 = x0,  dyn,  ‖u^t‖₂ ≤ ūcap
#
# Ball formulation (epigraph-slack, raw SOC coords): the SOC block ζ^t = (s^t; u^t)
# gives s^t ≥ ‖u^t‖₂ directly (socdet on the raw block), and a POS slack caps it,
# s^t + sl^t = ūcap ⟹ ‖u^t‖₂ ≤ ūcap. No isometric 1/√2 rescale is needed because
# the objective lives on the states, not on s — so the dynamics read u from ζ's tail
# verbatim and the control regularization stays exactly ρ‖u‖² (matching the ball-off QP).
# With ūcap = Inf the ball is dropped and controls are bare :NOC (plain tracking QP).
# τ :: Vector of length W+1 of target states (nx-vectors); only P τ is used.
# ---------------------------------------------------------------------------
function build_track_problem(cfg::QuadTrack, x0::AbstractVector, τ; ucap::Float64=Inf)
    W = cfg.W; P = cfg.P
    ball = isfinite(ucap)
    dζ = 1 + NU

    # cols: x^t (NOC,NX) t=0..W ; per t: control block(s)
    #   ball:  ζ^t (SOC,1+NU), sl^t (POS,1)   |   no ball: u^t (NOC,NU)
    bctrl = ball ? 2 : 1
    cx(t)  = t + 1
    cζ(t)  = (W + 1) + bctrl * t + 1        # control / ζ block
    csl(t) = (W + 1) + bctrl * t + 2        # ball slack (ball only)

    # rows: init(NX); per t: dyn(NX); ball ⇒ cap(1)
    rctrl = ball ? 2 : 1
    rinit() = 1
    rdyn(t) = 1 + rctrl * t + 1
    rcap(t) = 1 + rctrl * t + 2

    I_nx = Matrix(1.0I, NX, NX)
    I_nu = Matrix(1.0I, NU, NU)
    Bsel = ball ? hcat(zeros(NX, 1), cfg.Bd) : cfg.Bd   # reads u from ζ tail (raw)
    cap_head = reshape([1.0; zeros(NU)], 1, dζ)         # picks s = ζ head
    I_1 = reshape([1.0], 1, 1)

    ri, ci, bl = Int[], Int[], Matrix{Float64}[]
    push!(ri, rinit()); push!(ci, cx(0)); push!(bl, I_nx)
    for t in 0:W-1
        push!(ri, rdyn(t)); push!(ci, cx(t));     push!(bl, -cfg.Ad)
        push!(ri, rdyn(t)); push!(ci, cx(t + 1)); push!(bl, I_nx)
        push!(ri, rdyn(t)); push!(ci, cζ(t));     push!(bl, -Bsel)
        if ball
            push!(ri, rcap(t)); push!(ci, cζ(t));  push!(bl, cap_head)
            push!(ri, rcap(t)); push!(ci, csl(t)); push!(bl, I_1)
        end
    end
    B = blocksparse(ri, ci, bl)

    g = zeros(size(B, 1))
    g[rowrange(B, rinit())] .= x0
    if ball
        for t in 0:W-1; g[rowrange(B, rcap(t))] .= ucap; end
    end

    # Objective: Q = 2PᵀP on states, ρ on controls; c = -2 PᵀP τ on states
    PtP = P' * P
    c = zeros(size(B, 2))
    curv = Tuple{Int,Matrix{Float64}}[]
    for t in 0:W
        push!(curv, (cx(t), 2.0 .* PtP))
        c[colrange(B, cx(t))] .= -2.0 .* (P' * (P * τ[t + 1]))
    end
    for t in 0:W-1
        if ball
            Qc = zeros(dζ, dζ); for k in 2:dζ; Qc[k, k] = cfg.ρ; end   # ρ‖u‖² on ζ tail (raw)
            push!(curv, (cζ(t), Qc))
        else
            push!(curv, (cζ(t), cfg.ρ .* I_nu))
        end
    end

    nv = (W + 1) + bctrl * W
    cones = Vector{Symbol}(undef, nv)
    for t in 0:W; cones[cx(t)] = :NOC; end
    for t in 0:W-1
        if ball; cones[cζ(t)] = :SOC; cones[csl(t)] = :POS
        else;    cones[cζ(t)] = :NOC end
    end

    u0_extract(pvec) = ball ? pvec[colrange(B, cζ(0))[2:end]] : pvec[colrange(B, cζ(0))]
    return c, g, B, cones, curv, u0_extract
end

function assemble_Q(B, curv)
    Q = SheafSDP.allocblockdiag(B); fill!(Q, 0)
    for (b, M) in curv
        blk = block(Q, b, b, b); blk .+= M
    end
    return Q
end

# Gap tolerance kept well off the cone boundary (1e-4): a hard-binding SOC ball
# drives the iterate onto socdet ≈ 0, where socscale!'s sqrt sees a
# roundoff-negative determinant. Stopping short (1e-4) avoids it on ~99.7% of
# binding instances and is plenty accurate for control labels (which are
# standardized for training anyway). The rare remainder falls back to CLIPPED.
const SETTINGS = IPMSettings{Float64}(kkt=UzawaSettings{Float64}(raug=1e4),
                                      feas_tol=1e-6, gap_tol=1e-4, itmax=400)

const OK_STATUS = (OPTIMAL, NEAR_OPTIMAL)

"""
    conic_track_control(cfg, x0, τ; ucap) -> (u0, status)

First optimal control u^0 of the conic (SOC-ball) tracking solve. If the IPM
fails on a hard-binding instance, fall back to the unconstrained control
projected onto the ball (`status = :CLIPPED`) so datagen never crashes.
"""
function conic_track_control(cfg::QuadTrack, x0, τ; ucap::Float64)
    c, g, B, cones, curv, u0 = build_track_problem(cfg, x0, τ; ucap=ucap)
    try
        res = solve(IPMProblem(c, g, B, assemble_Q(B, curv), cones), SETTINGS)
        if res.status in OK_STATUS
            u = u0(res.p)
            nu = norm(u)
            return (nu > ucap ? u .* (ucap / nu) : u), res.status   # guard tiny overshoot
        end
    catch err
        err isa DomainError || rethrow(err)
    end
    u_un, _ = unconstrained_track_control(cfg, x0, τ)
    nu = norm(u_un)
    return (nu > ucap ? u_un .* (ucap / nu) : u_un), :CLIPPED
end

"First optimal control u^0 of the UNCONSTRAINED quadratic tracking solve (ball off)."
function unconstrained_track_control(cfg::QuadTrack, x0, τ)
    c, g, B, cones, curv, u0 = build_track_problem(cfg, x0, τ; ucap=Inf)
    res = solve(IPMProblem(c, g, B, assemble_Q(B, curv), cones), SETTINGS)
    return u0(res.p), res.status
end

end # module TrackingOracle
