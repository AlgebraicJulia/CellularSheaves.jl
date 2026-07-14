#
# Demonstration 9: CONIC coordination. Layer 1 is solved by the interior-point method (a conic
# program over the sheaf cochain complex) instead of the linear harmonic extension. A second-order
# cone REACHABILITY cap ‖q*_i‖ ≤ r is imposed on every agent: when a target orbits beyond radius r,
# the coordinated reference is projected onto the reachable ball — a feasibility-constrained
# harmonic extension — and the residual policy learns to track that conic q*.
#
# With r = ∞ (or coordination = :analytic) this is byte-identical to the linear solve; the IPM only
# earns its keep once the cap binds. On the fixed 13-agent / 4-target geometry, r = 2.5 binds on the
# far cluster (agents 5,7,8,9, whose free ‖q*‖ ≈ 2.6–3.4) and propagates through the consensus graph.
#
#   julia --project=examples/RL examples/RL/train/conic_reach.jl
# Env: RLJ_REACH (cap r, the canonical EVAL cone — training may see a whole bank via RLJ_BANKDIR),
#      RLJ_DRIFT, RLJ_RESSCALE, RLJ_TOTAL, RLJ_NENVS, RLJ_OUT,
#      RLJ_DRIFTRAND (1 ⇒ per-episode TRAINING drift in [0,RLJ_DRIFT] incl. zero-drift episodes;
#                      eval/report always uses the fixed RLJ_DRIFT — see sheaf_rl.jl drift_rand),
#      RLJ_DRIFTZEROP (P(zero-drift episode) when RLJ_DRIFTRAND), RLJ_PHASEFEAT (1 ⇒ obs gets the
#      [sin ωt, cos ωt] clock feature), RLJ_ANCHOR_HI/LO/ANNEAL (anchor anneals hi→lo over the first
#      ANNEAL-fraction of training — see train_td3's anchor_at).
#
# NOTE: :conic re-solves the IPM at every step of every reset (K solves per episode) UNLESS a
# reference bank is attached (RLJ_BANKDIR/RLJ_BANKFILE) — this is a CLUSTER job either way. Run via sbatch.
#
include(joinpath(@__DIR__, "..", "lib", "sheaf_rl.jl"))
using .SheafRL

reach = parse(Float64, get(ENV, "RLJ_REACH", "2.5"))
cfg = Config(
    drift_amp       = parse(Float64, get(ENV, "RLJ_DRIFT", "2.0")),
    residual_scale  = parse(Float64, get(ENV, "RLJ_RESSCALE", "2.0")),
    coordination    = :conic,
    coord_ball      = Dict(i => reach for i in 1:N_AGENTS),           # canonical EVAL cone (report/compare)
    base_law        = Symbol(get(ENV, "RLJ_BASELAW", "reference")),   # track the CONIC q* (no offset)
    ref_bank        = parse(Int, get(ENV, "RLJ_BANK", "1000")),       # precompute q* episodes on threads
    drift_rand      = get(ENV, "RLJ_DRIFTRAND", "0") == "1",
    drift_zero_prob = parse(Float64, get(ENV, "RLJ_DRIFTZEROP", "0.2")),
    phase_feature   = get(ENV, "RLJ_PHASEFEAT", "0") == "1",
)
train_and_save(cfg;
    total_steps        = parse(Int, get(ENV, "RLJ_TOTAL", "1500000")),
    n_envs             = parse(Int, get(ENV, "RLJ_NENVS", "32")),
    seed               = parse(Int, get(ENV, "RLJ_SEED", "1")),
    anchor_hi          = parse(Float32, get(ENV, "RLJ_ANCHOR_HI", "0.1")),
    anchor_lo          = parse(Float32, get(ENV, "RLJ_ANCHOR_LO", "0.02")),
    anchor_anneal_frac = parse(Float32, get(ENV, "RLJ_ANCHOR_ANNEAL", "0.5")),
    eval_every         = parse(Int, get(ENV, "RLJ_EVALEVERY", "25000")),
    out                = get(ENV, "RLJ_OUT", joinpath(@__DIR__, "..", "cache", "rl_conic_reach.jld2")),
    label              = "conic reachability cap ‖q*‖≤$(reach) (IPM Layer 1)")
