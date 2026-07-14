#
# sheaf_rl.jl — shared library for the sheaf-coordinated learned-control demonstrations.
#
# One readable implementation of the pieces every demo shares, so the demo scripts are thin
# configuration files instead of near-identical copies:
#
#   * the fixed 13-agent / 4-target tracking geometry (consensus graph + pinning),
#   * the cellular sheaf (homogeneous identity maps, or heterogeneous rotation/projection maps),
#   * the target trajectories (3-D figure-eights) and the unknown time-varying drift field,
#   * Layer 1   — the harmonic-extension reference q* (a linear sheaf-Laplacian solve),
#   * Layer 2   — the base control law, the learned residual, and the hard actuator-ball projection,
#   * the environment (reset / step / observation) and a generic TD3 trainer,
#   * evaluation (analytic vs drift-oracle vs sheaf+RL) and model save/load.
#
# Everything that VARIES between demonstrations lives in a single `Config` struct; the demo
# scripts set it and call `train_and_save` / `evaluate_policy`.
#
module SheafRL

using CellularSheaves
import CellularSheaves.NetworkSheaves.EuclideanSheaves: harmonic_extension
using Flux, Random, Statistics, LinearAlgebra, Printf, JLD2

# Layer-1 coordination via Richard's sheaf interior-point method (the conic harmonic extension).
# Provides build_ipm_harmonic / ipm_harmonic_problem / update_boundary! / extract_qstar.
include(joinpath(@__DIR__, "ipm_harmonic.jl"))

export Config, SheafEnv, reset!, step!, observation, base_control, drift_field, project_ball,
       push_history!, train_td3, evaluate_policy, train_and_save, save_model,
       build_reference_bank, reference_episode, load_reference_bank, sample_cone_config,
       N_AGENTS, N_TARGETS, STALK_DIM, AGENT_EDGES, PINNING, build_sheaf, make_targets,
       obs_dim, pinned_target, edge_angle

# ----------------------------------------------------------------------------------------------
# Fixed problem geometry — the heterogeneous 13-agent / 4-target tracking structure.
# ----------------------------------------------------------------------------------------------
const STALK_DIM = 3        # each agent/target lives in ℝ³
const N_AGENTS  = 13
const N_TARGETS = 4

# consensus graph (agent–agent edges): four triangles/squares joined by cross-links
const AGENT_EDGES = [(1,2),(2,3),(3,1), (4,5),(5,6),(6,4), (7,8),(8,9),(9,7),
                     (10,11),(11,12),(12,13),(13,10), (1,4),(6,13),(9,11),(7,10),(3,12)]

# pinning: which agent tracks which target. Agents absent here are *interior* nodes, anchored to a
# target only through the consensus graph.
const PINNING = let p = zeros(N_AGENTS, N_TARGETS)
    for (agent, target) in [(2,1), (5,2),(6,2), (7,3),(8,3),(9,3), (10,4),(13,4)]
        p[agent, target] = 1.0
    end
    p
end

"Target index that agent `i` is pinned to, or `0` if it is an interior (consensus-only) agent."
pinned_target(i) = (k = findfirst(!=(0.0), PINNING[i, :]); k === nothing ? 0 : k)

# ----------------------------------------------------------------------------------------------
# Configuration — everything that differs between the demonstrations.
# ----------------------------------------------------------------------------------------------
Base.@kwdef struct Config
    drift_amp::Float64      = 0.0                       # swirl amplitude a (0 ⇒ no drift)
    drift_freq::Float64     = 3.0                       # swirl temporal frequency ω
    history::Int            = drift_amp > 0 ? 6 : 0     # history-window length H (needed iff drift is time-varying)
    actuator_cap::Float64   = Inf                       # hard SOC cap ‖u‖ ≤ ū (Inf ⇒ unconstrained)
    residual_scale::Float64 = 2.0                       # residual scale α
    gain::Float64           = 2.0                       # base-controller gain k
    dt::Float64             = 0.05
    horizon::Int            = 300                       # steps per episode K
    rotate_consensus::Bool  = false                     # heterogeneous: per-edge rotation maps on consensus edges
    planar_agents::Set{Int} = Set{Int}()               # heterogeneous: agents that pin to the x-y projection only
    coordination::Symbol    = :analytic                 # Layer 1: :analytic (linear harmonic ext) | :conic (IPM)
    coord_cones::Any        = nothing                   # per-vertex cones for :conic (nothing ⇒ free ⇒ = analytic)
    coord_ball::Any         = nothing                   # Dict{Int,Float64} vertex→r: SOC reachability cap ‖q*_v‖≤r
    base_law::Symbol        = :laplacian                # Layer 2 base: :laplacian (sheaf flow, equilibrates at the
                                                        #   UNCONSTRAINED q*) | :reference (track the Layer-1 q*,
                                                        #   consistent with a conic reference — no steady-state offset)
    ref_bank::Int           = 0                         # >0: precompute this many reference episodes in parallel
                                                        #   (threads) and sample them in training — takes the IPM
                                                        #   out of the hot loop (q* depends only on targets + cone)
    drift_rand::Bool        = false                     # TRAINING only (via step!): sample a fresh per-episode
                                                        #   drift amplitude in [0, drift_amp] (w/ drift_zero_prob
                                                        #   chance of exactly 0) instead of always using drift_amp.
                                                        #   evaluate_policy/quick_eval always use the fixed drift_amp.
    drift_zero_prob::Float64 = 0.2                      # P(a no-drift episode) when drift_rand — teaches the
                                                        #   residual when NOT to act (fixes the no-drift regression)
    phase_feature::Bool     = false                     # append [sin(ωt), cos(ωt)] (known freq, NOT amplitude) to
                                                        #   the observation — a deterministic clock signal so the
                                                        #   residual can infer drift phase instead of blind-guessing
                                                        #   from a short history window.
end

"Observation dimension: own state, reference error, the (x,u) history window, and (if enabled) the
2-D phase-clock feature."
obs_dim(c::Config) = 2*STALK_DIM + 2*STALK_DIM*c.history + (c.phase_feature ? 2 : 0)

# ----------------------------------------------------------------------------------------------
# Restriction maps and the cellular sheaf.
# ----------------------------------------------------------------------------------------------
identity_map() = Matrix{Float64}(I, STALK_DIM, STALK_DIM)
rotation_z(θ)  = [cos(θ) -sin(θ) 0.0; sin(θ) cos(θ) 0.0; 0.0 0.0 1.0]
const PROJECT_XY = [1.0 0.0 0.0; 0.0 1.0 0.0]           # ℝ³ → ℝ² (drops z)

"Per-edge consensus rotation angle (deterministic in the edge index), spread over [-0.5, 0.5] rad."
edge_angle(edge_index) = 0.25 * (((edge_index - 1) % 5) - 2)

"""
Build the cellular sheaf for `cfg`. Consensus edges use identity maps, or per-edge z-rotations when
`cfg.rotate_consensus`; pinning edges use identity, or the x-y projection for `cfg.planar_agents`.
Returns the sheaf and the list of target vertices.
"""
function build_sheaf(cfg::Config)
    s = EuclideanSheaf{Float64}(fill(STALK_DIM, N_AGENTS + N_TARGETS))
    for (e, (i, j)) in enumerate(AGENT_EDGES)
        right = cfg.rotate_consensus ? rotation_z(edge_angle(e)) : identity_map()
        add_sheaf_edge!(s, i, j, identity_map(), right)
    end
    for i in 1:N_AGENTS, k in 1:N_TARGETS
        PINNING[i, k] == 0 && continue
        if i in cfg.planar_agents
            add_sheaf_edge!(s, i, N_AGENTS + k, PROJECT_XY, PROJECT_XY)   # track x-y only
        else
            add_sheaf_edge!(s, i, N_AGENTS + k, identity_map(), identity_map())
        end
    end
    s, collect(N_AGENTS+1 : N_AGENTS+N_TARGETS)
end

# ----------------------------------------------------------------------------------------------
# Targets and drift.
# ----------------------------------------------------------------------------------------------
"Four well-separated targets, each tracing a 3-D figure-eight (x∼sin, y∼sin2, z a quarter out of phase)."
function make_targets(cfg::Config, rng)
    centers = [[-3.0,-3,0.0], [3.0,-3,0.0], [3.0,3,0.0], [-3.0,3,0.0]]
    K = cfg.horizon
    Tg = Array{Float64}(undef, STALK_DIM, K+1, N_TARGETS)
    for k in 1:N_TARGETS
        center = centers[k] .+ 0.3 .* randn(rng, STALK_DIM)
        ω  = (0.6 + 0.4*rand(rng)) * (2π / (K*cfg.dt))
        φ  = rand(rng) * 2π
        Ax = 1.8 + 0.5*rand(rng); Ay = 1.2 + 0.4*rand(rng); Az = 1.3 + 0.5*rand(rng)
        for t in 0:K
            τ = t * cfg.dt
            Tg[1, t+1, k] = center[1] + Ax*sin(ω*τ + φ)
            Tg[2, t+1, k] = center[2] + Ay*sin(2ω*τ + φ)
            Tg[3, t+1, k] = center[3] + Az*sin(ω*τ + φ + π/2)
        end
    end
    Tg
end

"The unknown, time-varying swirl drift f(q,t) = amp·[-sin(q_y+ωt), sin(q_x+ωt), 0] (horizontal only).
`amp` defaults to `cfg.drift_amp` (unchanged behavior for eval/viz call sites); training's `step!`
passes the env's per-episode-randomized amplitude explicitly when `cfg.drift_rand`."
function drift_field(cfg::Config, q, t, amp::Float64 = cfg.drift_amp)
    amp == 0.0 && return zero(q)
    phase = cfg.drift_freq * t * cfg.dt
    amp .* [-sin(q[2] + phase), sin(q[1] + phase), 0.0]
end

"Closed-form projection onto the actuator ball ‖u‖ ≤ ū (the only place the SOC constraint enters)."
project_ball(cfg::Config, u) = (n = norm(u); n > cfg.actuator_cap ? u .* (cfg.actuator_cap / n) : u)

# ----------------------------------------------------------------------------------------------
# Environment.
# ----------------------------------------------------------------------------------------------
mutable struct SheafEnv
    cfg::Config
    sheaf
    target_vertices::Vector{Int}
    coboundary
    targets::Array{Float64,3}     # STALK_DIM × (K+1) × N_TARGETS
    reference::Array{Float64,3}   # STALK_DIM × K × N_AGENTS — the harmonic q* at each step
    x::Vector{Vector{Float64}}    # current agent positions
    hist_x::Array{Float64,3}      # STALK_DIM × H × N_AGENTS
    hist_u::Array{Float64,3}
    t::Int
    bank::Any                     # precomputed (targets, reference) episodes to sample, or nothing
    drift_amp_ep::Float64         # this episode's drift amplitude (= cfg.drift_amp unless cfg.drift_rand)
end

"Create an environment for `cfg` (builds the fixed sheaf once)."
function SheafEnv(cfg::Config)
    s, tv = build_sheaf(cfg)
    SheafEnv(cfg, s, tv, coboundary_map(s),
             zeros(0,0,0), zeros(0,0,0), Vector{Float64}[], zeros(0,0,0), zeros(0,0,0), 1, nothing,
             cfg.drift_amp)
end

"Layer 1: the harmonic extension q* = H⁻¹Bp for the current targets, at every step of the episode.
`cfg.coordination = :analytic` uses the linear solve; `:conic` solves the same coordination as a
sheaf conic QP through the IPM (identical q* with free cones; feasibility-constrained q* with cones)."
function compute_reference!(env::SheafEnv)
    K = env.cfg.horizon
    env.reference = Array{Float64}(undef, STALK_DIM, K, N_AGENTS)
    bd(t) = Dict(env.target_vertices[k] => env.targets[:, t, k] for k in 1:N_TARGETS)
    if env.cfg.coordination === :conic
        # assemble the sheaf conic QP once (structure fixed); only the boundary g changes per step
        prob, info = ipm_harmonic_problem(env.sheaf, bd(1); vertex_cones = env.cfg.coord_cones, ball = env.cfg.coord_ball)
        for t in 1:K
            update_boundary!(prob, info, bd(t))
            q = extract_qstar(solve(prob, default_ipm_settings()), info)
            for i in 1:N_AGENTS; env.reference[:, t, i] = q[(i-1)*STALK_DIM+1 : i*STALK_DIM]; end
        end
    else
        for t in 1:K
            q = Vector(harmonic_extension(env.sheaf, bd(t))[1])
            for i in 1:N_AGENTS; env.reference[:, t, i] = q[(i-1)*STALK_DIM+1 : i*STALK_DIM]; end
        end
    end
end

"""
Sample one cone from a diverse multi-family set: SOC reachability (varying radius on all agents),
SOC on a random subset, positive-orthant on a random subset, or free. Returns a NamedTuple of the
cone-relevant `Config` fields (splat into `Config`). This is the constraint diversity that trains a
cone-GENERAL residual policy: the policy never sees the cone, only the reference it induces, so
covering the induced-reference distribution is exactly what generalization needs.
"""
function sample_cone_config(rng)
    fam = rand(rng, (:ball, :ball_subset, :orthant, :free))
    if fam === :free
        (coordination = :conic, coord_ball = nothing, coord_cones = nothing, cone_label = "free")
    elseif fam === :ball
        r = 2.0 + 1.5 * rand(rng)
        (coordination = :conic, coord_ball = Dict(i => r for i in 1:N_AGENTS), coord_cones = nothing,
         cone_label = "ball_r$(round(r, digits=2))")
    elseif fam === :ball_subset
        r = 2.0 + rand(rng)
        sub = [i for i in 1:N_AGENTS if rand(rng) < 0.5]; isempty(sub) && push!(sub, rand(rng, 1:N_AGENTS))
        (coordination = :conic, coord_ball = Dict(i => r for i in sub), coord_cones = nothing,
         cone_label = "ballsub_r$(round(r, digits=2))_n$(length(sub))")
    else
        NV = N_AGENTS + N_TARGETS
        cones = AbstractCone[CofreeCone() for _ in 1:NV]; capped = 0
        for i in 1:N_AGENTS; rand(rng) < 0.4 && (cones[i] = PositiveCone(); capped += 1); end
        capped == 0 && (cones[rand(rng, 1:N_AGENTS)] = PositiveCone())
        (coordination = :conic, coord_ball = nothing, coord_cones = cones, cone_label = "orthant")
    end
end

"""
One reference episode `(targets, reference)` under `cfg`'s cone — a single, self-contained call with
no threads. This is the unit of the sharded datagen: run it in independent PROCESSES (a SLURM array),
not threads, because the IPM allocates heavily and 16 threads serialize on Julia's stop-the-world GC
(the threaded `build_reference_bank` measured ~3× SLOWER per solve than serial). Separate processes =
separate heaps = near-linear scaling. `cfg`'s cone (`coord_ball` / `coord_cones`) may vary per call,
so a multi-cone bank is just a stream of episodes each built under a sampled cone.
"""
function reference_episode(cfg::Config, rng)
    env = SheafEnv(cfg)
    env.targets = make_targets(cfg, rng)
    compute_reference!(env)
    (env.targets, env.reference)
end

"Load and concatenate reference-bank shard files (each a jld2 holding key \"bank\")."
function load_reference_bank(paths)
    bank = NTuple{2, Array{Float64,3}}[]
    for p in paths; append!(bank, JLD2.load(p, "bank")); end
    bank
end

"Convenience threaded bank build (single node). Prefer the process-sharded datagen for large banks."
function build_reference_bank(cfg::Config, n::Int; seed=0)
    bank = Vector{NTuple{2, Array{Float64,3}}}(undef, n)
    nblas = BLAS.get_num_threads(); BLAS.set_num_threads(1)
    try
        Threads.@threads for b in 1:n
            bank[b] = reference_episode(cfg, MersenneTwister(seed + b))
        end
    finally
        BLAS.set_num_threads(nblas)
    end
    bank
end

"Reset to a fresh episode: new targets + reference (sampled from the bank when one is attached,
recomputed otherwise), agents started near q*(t=1)."
function reset!(env::SheafEnv, rng)
    c = env.cfg
    if env.bank !== nothing
        env.targets, env.reference = env.bank[rand(rng, 1:length(env.bank))]
    else
        env.targets = make_targets(c, rng)
        compute_reference!(env)
    end
    env.x = [env.reference[:, 1, i] .+ 0.4 .* randn(rng, STALK_DIM) for i in 1:N_AGENTS]
    env.hist_x = zeros(STALK_DIM, c.history, N_AGENTS)
    env.hist_u = zeros(STALK_DIM, c.history, N_AGENTS)
    env.t = 1
    env.drift_amp_ep = !c.drift_rand ? c.drift_amp :
        rand(rng) < c.drift_zero_prob ? 0.0 : c.drift_amp * rand(rng)
    env
end

"Layer 2 base law. `:laplacian` is the sheaf feedback -k (L_F z)_i, whose equilibrium is the
UNCONSTRAINED harmonic extension; `:reference` tracks the Layer-1 q* directly, -k (x_i − q*_i),
so the base is consistent with a conic (cone-projected) reference — no steady-state offset."
function base_control(env::SheafEnv)
    env.cfg.base_law === :reference &&
        return [-env.cfg.gain .* (env.x[i] .- env.reference[:, env.t, i]) for i in 1:N_AGENTS]
    z = zeros(STALK_DIM * (N_AGENTS + N_TARGETS))
    for i in 1:N_AGENTS;  z[(i-1)*STALK_DIM+1 : i*STALK_DIM] = env.x[i];               end
    for k in 1:N_TARGETS; z[(N_AGENTS+k-1)*STALK_DIM+1 : (N_AGENTS+k)*STALK_DIM] = env.targets[:, env.t, k]; end
    Lz = env.coboundary' * (env.coboundary * z)
    [-env.cfg.gain .* Lz[(i-1)*STALK_DIM+1 : i*STALK_DIM] for i in 1:N_AGENTS]
end

"Per-agent observation o_i = [x_i; q*_i − x_i; history; (phase clock if enabled)].
The phase clock is [sin(ωt), cos(ωt)] using the KNOWN frequency `cfg.drift_freq` — a deterministic
timing signal (any real controller has a clock), not the unknown drift amplitude or value itself."
function observation(env::SheafEnv, i)
    base = vcat(env.x[i], env.reference[:, env.t, i] .- env.x[i],
               vec(env.hist_x[:, :, i]), vec(env.hist_u[:, :, i]))
    if env.cfg.phase_feature
        phase = env.cfg.drift_freq * env.t * env.cfg.dt
        base = vcat(base, sin(phase), cos(phase))
    end
    Float32.(base)
end

push_history!(env::SheafEnv, i, u) = begin
    env.cfg.history == 0 && return
    env.hist_x[:, 1:end-1, i] = env.hist_x[:, 2:end, i]; env.hist_x[:, end, i] = env.x[i]
    env.hist_u[:, 1:end-1, i] = env.hist_u[:, 2:end, i]; env.hist_u[:, end, i] = u
end

"Apply per-agent residuals: u_i = Π_ball(base_i + α·Δ_i); advance the plant with drift. Returns (reward, done)."
function step!(env::SheafEnv, residual::Matrix{Float64})
    c = env.cfg
    base = base_control(env)
    reward = zeros(Float32, N_AGENTS)
    for i in 1:N_AGENTS
        u = project_ball(c, base[i] .+ c.residual_scale .* residual[:, i])
        reward[i] = Float32(-norm(env.x[i] .- env.reference[:, env.t, i]))
        push_history!(env, i, u)
        env.x[i] = env.x[i] .+ c.dt .* (u .+ drift_field(c, env.x[i], env.t, env.drift_amp_ep))
    end
    env.t += 1
    reward, env.t > c.horizon
end

# ----------------------------------------------------------------------------------------------
# TD3 — twin-critic, target-smoothed, delayed-actor DDPG with a residual anchor.
# ----------------------------------------------------------------------------------------------
make_actor(c::Config)  = Chain(Dense(obs_dim(c), 256, relu), Dense(256, 256, relu), Dense(256, STALK_DIM, tanh))
make_critic(c::Config) = Chain(Dense(obs_dim(c) + STALK_DIM, 256, relu), Dense(256, 256, relu), Dense(256, 1))

mutable struct Replay
    s::Matrix{Float32}; a::Matrix{Float32}; r::Vector{Float32}
    s2::Matrix{Float32}; done::Vector{Float32}
    capacity::Int; idx::Int; full::Bool
end
Replay(capacity, obs, act) = Replay(zeros(Float32, obs, capacity), zeros(Float32, act, capacity),
    zeros(Float32, capacity), zeros(Float32, obs, capacity), zeros(Float32, capacity), capacity, 0, false)

n_stored(rb::Replay) = rb.full ? rb.capacity : rb.idx
function store!(rb::Replay, s, a, r, s2, done)
    for i in 1:length(r)
        p = (rb.idx % rb.capacity) + 1
        rb.s[:, p] = s[:, i]; rb.a[:, p] = a[:, i]; rb.r[p] = r[i]; rb.s2[:, p] = s2[:, i]; rb.done[p] = done[i]
        rb.idx += 1; rb.idx >= rb.capacity && (rb.full = true)
    end
end
function sample(rb::Replay, batch)
    idx = rand(1:n_stored(rb), batch)
    rb.s[:, idx], rb.a[:, idx], rb.r[idx], rb.s2[:, idx], rb.done[idx]
end

"Mean tracking error ‖x−q*‖ of the *deterministic* `:rl` policy over `n` fixed eval episodes — a
cheap, deterministic score for checkpoint selection (same `rng` ⇒ same episodes ⇒ comparable)."
function quick_eval(actor, cfg::Config, evalenv::SheafEnv, rng; n=5)
    total = 0.0; cnt = 0
    for _ in 1:n
        reset!(evalenv, rng)
        for _ in 1:cfg.horizon
            base = base_control(evalenv)
            for i in 1:N_AGENTS
                u = project_ball(cfg, base[i] .+ cfg.residual_scale .* Float64.(vec(actor(observation(evalenv, i)))))
                total += norm(evalenv.x[i] .- evalenv.reference[:, evalenv.t, i]); cnt += 1
                push_history!(evalenv, i, u)
                evalenv.x[i] = evalenv.x[i] .+ cfg.dt .* (u .+ drift_field(cfg, evalenv.x[i], evalenv.t))
            end
            evalenv.t += 1
        end
    end
    total / cnt
end

"""
Train a residual TD3 policy on `cfg`. Many environments roll in parallel; each agent is an
independent sample of the same weight-shared policy (so the controller is decentralized and
size-agnostic). TD3 wanders, so we eval a deterministic score every `eval_every` steps and return
the BEST-scoring actor, not the last one (final ≠ best; the last policy is often past its peak).
"""
function train_td3(cfg::Config; total_steps, n_envs=32, seed=1,
                   γ=0.99f0, τ=0.005f0, actor_lr=3f-4, critic_lr=1f-3, batch=256,
                   explore=0.15f0, target_noise=0.2f0, target_clip=0.5f0,
                   policy_delay=2, warmup=2000, replay_capacity=300_000,
                   anchor_hi=0.1f0, anchor_lo=0.02f0, anchor_anneal_frac=0.5f0,
                   eval_every=25_000, eval_scenarios=5)
    # anchor anneal: start conservative (stability near the base law), relax over the first
    # `anchor_anneal_frac` of training so late updates get more room to exploit (measured: anchor
    # 0.1→0.633, 0.05→0.482, 0.02→0.439 final score — lower anchor = closer to the oracle, but early
    # instability needs the higher value; annealing gets both).
    anchor_at(step) = anchor_hi + (anchor_lo - anchor_hi) * min(1f0, Float32(step) / (anchor_anneal_frac * total_steps))
    rng = MersenneTwister(seed)
    actor = make_actor(cfg); critic1 = make_critic(cfg); critic2 = make_critic(cfg)
    actor_t = deepcopy(actor); critic1_t = deepcopy(critic1); critic2_t = deepcopy(critic2)
    opt_a  = Flux.setup(Adam(actor_lr), actor)
    opt_c1 = Flux.setup(Adam(critic_lr), critic1); opt_c2 = Flux.setup(Adam(critic_lr), critic2)
    replay = Replay(replay_capacity, obs_dim(cfg), STALK_DIM)
    polyak!(dst, src) = for (p, q) in zip(Flux.params(dst), Flux.params(src)); p .= (1-τ).*p .+ τ.*q; end

    bank = nothing
    if cfg.ref_bank > 0
        bankdir  = get(ENV, "RLJ_BANKDIR", "")     # a directory of datagen shards (multi-cone bank)
        bankfile = get(ENV, "RLJ_BANKFILE", "")    # a single cached bank file
        if !isempty(bankdir) && isdir(bankdir)
            shards = sort(filter(f -> endswith(f, ".jld2"), readdir(bankdir; join=true)))
            bank = load_reference_bank(shards)
            @printf("  reference bank: %d episodes from %d shards in %s\n", length(bank), length(shards), bankdir)
        elseif !isempty(bankfile) && isfile(bankfile)
            bank = JLD2.load(bankfile, "bank")
            @printf("  reference bank: %d episodes loaded from %s\n", length(bank), bankfile)
        else
            tb = time()
            bank = build_reference_bank(cfg, cfg.ref_bank; seed = seed + 31)
            @printf("  reference bank: %d episodes on %d threads in %.0fs\n", cfg.ref_bank, Threads.nthreads(), time()-tb)
            isempty(bankfile) || (jldsave(bankfile; bank); @printf("  reference bank: cached to %s\n", bankfile))
        end
    end
    envs = map(1:n_envs) do _
        e = SheafEnv(cfg); e.bank = bank; reset!(e, rng)
    end
    all_obs() = reduce(hcat, [reduce(hcat, [observation(e, i) for i in 1:N_AGENTS]) for e in envs])

    evalenv = SheafEnv(cfg); evalenv.bank = bank      # deterministic checkpoint-selection env
    best_actor = deepcopy(actor); best_score = Inf; last_eval = 0

    obs = all_obs(); steps = 0; last_loss = NaN32; t0 = time()
    while steps < total_steps
        act = clamp.(actor(obs) .+ explore .* randn(Float32, STALK_DIM, size(obs, 2)), -1f0, 1f0)
        rewards = Float32[]; dones = Float32[]; col = 1
        for e in envs
            Δ = Float64.(act[:, col:col+N_AGENTS-1]); col += N_AGENTS
            r, done = step!(e, Δ)
            append!(rewards, r); append!(dones, fill(done ? 1f0 : 0f0, N_AGENTS))
            done && reset!(e, rng)
        end
        obs2 = all_obs()
        store!(replay, obs, act, rewards, obs2, dones); obs = obs2; steps += size(act, 2)

        if n_stored(replay) >= max(batch, warmup)
            s, a, r, s2, d = sample(replay, batch)
            ε  = clamp.(target_noise .* randn(Float32, STALK_DIM, batch), -target_clip, target_clip)
            a2 = clamp.(actor_t(s2) .+ ε, -1f0, 1f0)
            y  = r .+ γ .* (1f0 .- d) .* min.(vec(critic1_t(vcat(s2, a2))), vec(critic2_t(vcat(s2, a2))))
            _, g1 = Flux.withgradient(m -> mean((vec(m(vcat(s, a))) .- y).^2), critic1); Flux.update!(opt_c1, critic1, g1[1])
            l, g2 = Flux.withgradient(m -> mean((vec(m(vcat(s, a))) .- y).^2), critic2); Flux.update!(opt_c2, critic2, g2[1])
            last_loss = Float32(l)
            if (steps ÷ size(act, 2)) % policy_delay == 0
                anchor = anchor_at(steps)
                _, ga = Flux.withgradient(actor) do m
                    aa = m(s); -mean(critic1(vcat(s, aa))) + anchor * mean(aa.^2)   # maximise Q, stay near the base controller
                end
                Flux.update!(opt_a, actor, ga[1])
                polyak!(actor_t, actor); polyak!(critic1_t, critic1); polyak!(critic2_t, critic2)
            end
        end
        if steps % (n_envs*N_AGENTS*50) < n_envs*N_AGENTS
            @printf("  steps %7d  buffer %7d  critic_loss %.4f  anchor %.3f  %.0fs\n",
                    steps, n_stored(replay), last_loss, anchor_at(steps), time()-t0)
            flush(stdout)
        end
        if n_stored(replay) >= max(batch, warmup) && steps - last_eval >= eval_every
            last_eval = steps
            score = quick_eval(actor, cfg, evalenv, MersenneTwister(seed + 999); n = eval_scenarios)
            score < best_score && (best_score = score; best_actor = deepcopy(actor))
            @printf("  eval @ %7d  rl=%.4f  best=%.4f  %s\n", steps, score, best_score,
                    score ≈ best_score ? "← kept" : "")
            flush(stdout)
        end
    end
    @printf("  best eval score over training: %.4f\n", best_score)
    best_actor
end

# ----------------------------------------------------------------------------------------------
# Evaluation and persistence.
# ----------------------------------------------------------------------------------------------
"""
Mean tracking error ‖x − q*‖ over `n_scenarios` rollouts, for each controller:
`:analytic` (base law), `:oracle` (base law minus the exact drift), `:rl` (base + learned residual).
All three are projected onto the actuator ball.
"""
function evaluate_policy(actor, cfg::Config; rng, n_scenarios=30)
    env = SheafEnv(cfg)
    rollout(mode) = begin
        total = 0.0; count = 0
        for _ in 1:n_scenarios
            reset!(env, rng)
            for _ in 1:cfg.horizon
                base = base_control(env)
                for i in 1:N_AGENTS
                    raw = mode === :analytic ? base[i] :
                          mode === :oracle   ? base[i] .- drift_field(cfg, env.x[i], env.t) :
                                               base[i] .+ cfg.residual_scale .* Float64.(vec(actor(observation(env, i))))
                    u = project_ball(cfg, raw)
                    total += norm(env.x[i] .- env.reference[:, env.t, i]); count += 1
                    push_history!(env, i, u)
                    env.x[i] = env.x[i] .+ cfg.dt .* (u .+ drift_field(cfg, env.x[i], env.t))
                end
                env.t += 1
            end
        end
        total / count
    end
    (analytic = rollout(:analytic), oracle = rollout(:oracle), rl = rollout(:rl))
end

"Save the actor weights plus everything a renderer needs to rebuild the exact env."
function save_model(path, actor, cfg::Config, result)
    net = cpu(actor)
    jldsave(path;
        Ws = [Array(l.weight) for l in net.layers], bs = [Array(l.bias) for l in net.layers],
        N = N_AGENTS, M = N_TARGETS, ND = STALK_DIM, EDGES = AGENT_EDGES, PIN = PINNING,
        DRIFT = cfg.drift_amp, DRIFTW = cfg.drift_freq, HIST = cfg.history, KGAIN = cfg.gain,
        DT = cfg.dt, K = cfg.horizon, RESS = cfg.residual_scale, OBS = obs_dim(cfg), UCAP = cfg.actuator_cap,
        rotate_consensus = cfg.rotate_consensus, planar_agents = collect(cfg.planar_agents),
        COORD = String(cfg.coordination),                                    # Layer-1 solver: analytic | conic
        BASELAW = String(cfg.base_law),                                      # Layer-2 base: laplacian | reference
        PHASEFEAT = cfg.phase_feature,                # obs includes [sin(ωt),cos(ωt)] — MUST match to reconstruct obs_dim
        DRIFTRAND = cfg.drift_rand, DRIFTZEROP = cfg.drift_zero_prob,        # training-time drift randomization (eval is fixed)
        BALLV = cfg.coord_ball === nothing ? Int[]     : collect(keys(cfg.coord_ball)),   # capped vertices …
        BALLR = cfg.coord_ball === nothing ? Float64[] : collect(values(cfg.coord_ball)), # … and their radii
        eval = Dict("analytic"=>result.analytic, "oracle"=>result.oracle, "rl"=>result.rl))
end

"Train on `cfg`, evaluate, save to `out`, and print the result. The standard demo entry point.
Extra kwargs (`anchor`, `eval_every`, `explore`, …) are forwarded to `train_td3`."
function train_and_save(cfg::Config; total_steps, out, seed=1, n_envs=32, n_scenarios=30, label="", train_kw...)
    @info "sheaf-track RL$(isempty(label) ? "" : " — "*label)" drift=cfg.drift_amp cap=cfg.actuator_cap history=cfg.history rotate=cfg.rotate_consensus planar=collect(cfg.planar_agents) total_steps
    actor = train_td3(cfg; total_steps, n_envs, seed, train_kw...)
    result = evaluate_policy(actor, cfg; rng = MersenneTwister(seed + 777), n_scenarios)
    @printf("\n=== mean ‖x−q*‖ over %d scenarios (cap ū=%.3g, drift a=%.3g) ===\n", n_scenarios, cfg.actuator_cap, cfg.drift_amp)
    @printf("  analytic = %.4f\n  drift-oracle = %.4f\n  sheaf+RL = %.4f\n", result.analytic, result.oracle, result.rl)
    save_model(out, actor, cfg, result)
    @info "saved" out
    actor, result
end

end # module SheafRL
