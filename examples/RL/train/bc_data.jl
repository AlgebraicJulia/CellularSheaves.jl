#
# Phase A (sheaf-native, NO [A B]/LTI encoding) — collect behavior-cloning data for the layered
# ECC tracking architecture. Homogeneous agents, NO drift: single integrator  q̇ = u.
#
#   Outer (analytic, reused):  spatial sheaf harmonic extension  q*_i = harmonic_extension(targets)
#   Expert (BC label):         the paper's analytic law          u_i = -k · (L z)_i
#   Policy target (trainer):   π_θ(x_i, e_i),  e_i = q*_i - x_i
#
# This is the deliberate "no-drift canary": the expert is the paper's decentralized feedback law,
# and BC tests whether the harmonic error e_i is a sufficient statistic to clone it. The reference
# comes from the spatial sheaf Laplacian (Fairbanks formulation: agents interior, targets boundary)
# — explicitly NOT the time-expanded LTI (Ad,Bd) readout used by train_bc_conic.jl.
#
# Diverse scenarios over team sizes / topologies / target trajectories / seeds. Threaded; JLD2.
# Runs on the CLUSTER — CPU + threads, embarrassingly parallel (each rollout is independent).
#
#   julia -t auto --project=examples/RL examples/RL/train/bc_data.jl
# Env: BCS_ND, BCS_NS, BCS_NSCN, BCS_K, BCS_DT, BCS_KGAIN, BCS_SEED, BCS_EXTRA_EDGES, BCS_OUT
#
using CellularSheaves                       # EuclideanSheaf, add_sheaf_edge!, coboundary_map
import CellularSheaves.NetworkSheaves.EuclideanSheaves: harmonic_extension   # not re-exported
using LinearAlgebra, Random, JLD2, Printf, Statistics
import Base.Threads: @threads, nthreads

const ND     = parse(Int,     get(ENV,"BCS_ND","2"))          # state dim per vertex (2 = planar canary)
const NS     = parse.(Int,    split(get(ENV,"BCS_NS","3,4,5"), ","))   # team sizes
const NSCN   = parse(Int,     get(ENV,"BCS_NSCN","200"))      # scenarios per team size
const K      = parse(Int,     get(ENV,"BCS_K","400"))         # steps per rollout
const DT     = parse(Float64, get(ENV,"BCS_DT","0.02"))
const KGAIN  = parse(Float64, get(ENV,"BCS_KGAIN","8.0"))     # analytic-law proportional gain
const SEED   = parse(Int,     get(ENV,"BCS_SEED","1"))
const XEDGES = parse(Int,     get(ENV,"BCS_EXTRA_EDGES","1")) # random consensus edges beyond the ring
const CWEIGHT= parse(Float64, get(ENV,"BCS_CWEIGHT","1.0"))  # consensus restriction scale (stress test: >1)
const DRIFT  = parse(Float64, get(ENV,"BCS_DRIFT","0.0"))    # nonlinear swirl-drift amplitude (Stage 2: >0)
const OUT    = get(ENV,"BCS_OUT", joinpath(@__DIR__,"..","cache","bc_sheaf_dataset.jld2"))

# Unknown nonlinear agent drift f(q) (swirl; divergence-free). 0 → pure single integrator (canary).
drift(q) = DRIFT == 0.0 ? zero(q) : DRIFT .* vcat(-sin(q[2]), sin(q[1]), zeros(length(q) - 2))

# --- sheaf construction (public API; replicates Joana's build_agent_sheaf) -------------------
# Agent–agent consensus edges carry identity restriction maps; pinning edges (agent i ↔ target j)
# carry √w · I, so the coboundary residual weights tracking against consensus.
eye(nd) = Matrix{Float64}(I, nd, nd)
function build_tracking_sheaf(N::Int, M::Int, nd::Int, agent_edges, pin)
    s = EuclideanSheaf{Float64}(fill(nd, N + M))                 # verts 1..N agents, N+1..N+M targets
    cw = sqrt(CWEIGHT) * eye(nd)                                  # consensus restriction (scaled for stress test)
    for (i, j) in agent_edges
        add_sheaf_edge!(s, i, j, cw, cw)
    end
    for i in 1:N, j in 1:M
        w = pin[i, j]
        if w != 0.0
            R = sqrt(abs(w)) * eye(nd)
            add_sheaf_edge!(s, i, N + j, R, R)
        end
    end
    return s, collect(1:N), collect(N+1:N+M)
end

# --- random homogeneous scenario: N agents, M=N targets (1:1 pinning), ring + random consensus --
function make_topology(N::Int, rng)
    edges = Tuple{Int,Int}[(i, i % N + 1) for i in 1:N]          # ring
    for _ in 1:XEDGES
        i = rand(rng, 1:N); j = rand(rng, 1:N)
        i != j && push!(edges, (min(i, j), max(i, j)))
    end
    pin = zeros(N, N)
    for i in 1:N; pin[i, i] = 1.0; end                           # agent i tracks target i
    return edges, pin
end

# --- smooth random target trajectory (nd × (K+1)) --------------------------------------------
function target_traj(rng, nd::Int)
    c     = (rand(rng, nd) .- 0.5) .* 6.0                        # center
    A     = 0.5 .+ rand(rng, nd) .* 2.0                          # amplitudes
    omega = (0.3 .+ rand(rng, nd) .* 0.7) .* (2π / (K * DT))     # frequencies
    phi   = rand(rng, nd) .* 2π
    drift = (rand(rng, nd) .- 0.5) .* 0.3
    P = Array{Float64}(undef, nd, K + 1)
    for t in 0:K
        tt = t * DT
        @. P[:, t + 1] = c + A * sin(omega * tt + phi) + drift * tt
    end
    return P
end

# --- one rollout: roll the analytic expert, record (state, q*, control) per agent per step ----
function rollout(N::Int, seed::Int)
    rng = MersenneTwister(seed)
    edges, pin = make_topology(N, rng)
    s, _averts, tverts = build_tracking_sheaf(N, N, ND, edges, pin)
    d = coboundary_map(s)

    Tg = Array{Float64}(undef, ND, K + 1, N)                     # target trajectories
    for j in 1:N; Tg[:, :, j] = target_traj(rng, ND); end

    X = Array{Float64}(undef, ND, K, N)                          # agent states
    Q = Array{Float64}(undef, ND, K, N)                          # harmonic references q*_i
    U = Array{Float64}(undef, ND, K, N)                          # expert controls (labels)
    x = [Tg[:, 1, i] .+ 0.5 .* randn(rng, ND) for i in 1:N]      # init near own target

    for t in 1:K
        # cochain z = [agent positions ; current target positions]
        z = zeros(ND * (N + N))
        for i in 1:N; z[(i-1)*ND+1 : i*ND]        = x[i];        end
        for j in 1:N; z[(N+j-1)*ND+1 : (N+j)*ND]  = Tg[:, t, j]; end
        Lz = d' * (d * z)

        # harmonic reference: pin current targets as boundary, solve agents (interior)
        boundary = Dict(tverts[j] => Tg[:, t, j] for j in 1:N)
        qstar, _null = harmonic_extension(s, boundary)
        qvec = Vector(qstar)                                     # agents are verts 1..N

        for i in 1:N
            u = -KGAIN .* Lz[(i-1)*ND+1 : i*ND]
            X[:, t, i] = x[i]
            Q[:, t, i] = qvec[(i-1)*ND+1 : i*ND]
            U[:, t, i] = u
            x[i] = x[i] .+ DT .* (u .+ drift(x[i]))              # q̇ = u + f(q)  (f=0 → pure integrator)
        end
    end
    return (; N, seed, Tg, X, Q, U)
end

# --- one-shot validation (single-threaded, before the sweep): prove the math in the log ------
# At the harmonic extension q*, the interior (agent) Laplacian rows are ~0 by construction, so the
# expert control u=-k(Lz) evaluated AT q* must be ~0. We also report the magnitudes BC will see.
function validate_one()
    println("=== VALIDATION (single scenario, N=$(first(NS))) ===========================")
    r = rollout(first(NS), SEED)
    N = r.N
    # rebuild sheaf to re-check the harmonic residual at t=1
    rng = MersenneTwister(SEED); edges, pin = make_topology(N, rng)
    s, _av, tverts = build_tracking_sheaf(N, N, ND, edges, pin); d = coboundary_map(s)
    boundary = Dict(tverts[j] => r.Tg[:, 1, j] for j in 1:N)
    qstar, _ = harmonic_extension(s, boundary); qvec = Vector(qstar)
    zstar = zeros(ND * (N + N))
    for i in 1:N; zstar[(i-1)*ND+1:i*ND] = qvec[(i-1)*ND+1:i*ND]; end
    for j in 1:N; zstar[(N+j-1)*ND+1:(N+j)*ND] = r.Tg[:, 1, j]; end
    Lzstar = d' * (d * zstar)
    interior_res = maximum(abs, Lzstar[1:ND*N])                      # should be ~0 (harmonic interior)
    track_err = mean(norm(r.X[:, t, i] .- r.Tg[:, t, i]) for t in 1:K, i in 1:N)
    ctrl_mag  = mean(norm(r.U[:, t, i]) for t in 1:K, i in 1:N)
    err_mag   = mean(norm(r.Q[:, t, i] .- r.X[:, t, i]) for t in 1:K, i in 1:N)
    nbad      = count(!isfinite, r.X) + count(!isfinite, r.U) + count(!isfinite, r.Q)
    @printf("  harmonic interior residual ‖(Lz*)_agents‖∞ = %.3e   (expect ~0)\n", interior_res)
    @printf("  mean expert |u|            = %.4f\n", ctrl_mag)
    @printf("  mean |e=q*-x|              = %.4f\n", err_mag)
    @printf("  mean tracking err |x-tgt|  = %.4f   (final-time settles below this)\n", track_err)
    @printf("  non-finite entries         = %d   (must be 0)\n", nbad)
    println("============================================================================")
    flush(stdout)
    @assert interior_res < 1e-6 "harmonic interior residual too large — sheaf/boundary wiring is wrong"
    @assert nbad == 0 "non-finite values in rollout — expert diverged"
    @assert ctrl_mag < 1e4 "analytic law DIVERGED (mean|u|=$ctrl_mag): config is dynamically unstable. " *
        "Discrete stability needs DT·KGAIN·λmax(L) < 2, and λmax(L) grows with CWEIGHT / graph density. " *
        "Lower BCS_KGAIN (or BCS_DT) when raising BCS_CWEIGHT / BCS_EXTRA_EDGES."
end

# --- collect (threaded) ----------------------------------------------------------------------
function main()
    @info "Sheaf BC data collection — config" ND NS NSCN K DT KGAIN seed=SEED threads=nthreads() out=OUT
    println("Julia $(VERSION), nthreads=$(nthreads())"); flush(stdout)
    validate_one()

    jobs = [(N, SEED * 1_000_003 + 7919 * ni + s)
            for (ni, N) in enumerate(NS) for s in 1:NSCN]
    @info "Starting threaded sweep" total=length(jobs)
    res  = Vector{Any}(undef, length(jobs))
    done = Threads.Atomic{Int}(0)
    t0 = time()
    @threads for q in 1:length(jobs)
        (N, sd) = jobs[q]
        res[q]  = rollout(N, sd)
        c = Threads.atomic_add!(done, 1) + 1
        if c % 50 == 0
            @info "  progress" c total=length(jobs) elapsed_s=round(time() - t0, digits=1)
            flush(stdout); flush(stderr)
        end
    end

    # aggregate diagnostics across the whole dataset (so the log proves data health)
    track = mean(mean(norm(r.X[:, t, i] .- r.Tg[:, t, i]) for t in 1:K, i in 1:r.N) for r in res)
    umag  = mean(mean(norm(r.U[:, t, i]) for t in 1:K, i in 1:r.N) for r in res)
    nbad  = sum(count(!isfinite, r.X) + count(!isfinite, r.U) + count(!isfinite, r.Q) for r in res)
    @info "Sweep complete" elapsed_s=round(time() - t0, digits=1) mean_track_err=round(track, digits=4) mean_ctrl=round(umag, digits=4) nonfinite=nbad
    nbad == 0 || @warn "Dataset contains non-finite values — investigate before training" nbad

    mkpath(dirname(OUT))
    jldsave(OUT;
        Ns    = [r.N    for r in res],
        seeds = [r.seed for r in res],
        Tg    = [r.Tg   for r in res],
        X     = [r.X    for r in res],
        Q     = [r.Q    for r in res],
        U     = [r.U    for r in res],
        meta  = Dict("ND"=>ND, "K"=>K, "DT"=>DT, "KGAIN"=>KGAIN, "CWEIGHT"=>CWEIGHT, "DRIFT"=>DRIFT,
                     "XEDGES"=>XEDGES, "NS"=>NS, "NSCN"=>NSCN,
                     "kind"=> DRIFT==0.0 ? "sheaf_bc_nodrift" : "sheaf_bc_drift"))
    npairs = sum(r.N * K for r in res)
    @printf("Saved %d scenarios (%d agent-step pairs) -> %s\n", length(res), npairs, OUT)
    flush(stdout)
end

main()
