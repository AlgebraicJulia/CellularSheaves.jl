#
# Demonstration 7: the SAME 13-agent / 4-target structure, but a genuinely heterogeneous sheaf —
# per-edge rotation maps on the consensus edges, and x-y projection pinning for three "planar"
# agents (5, 8, 10) that track only their target's x-y.
#
#   julia --project=examples/RL examples/RL/train/heterogeneous.jl
# Env: RLJ_DRIFT, RLJ_RESSCALE, RLJ_TOTAL, RLJ_NENVS, RLJ_OUT
#
include(joinpath(@__DIR__, "..", "lib", "sheaf_rl.jl"))
using .SheafRL

cfg = Config(
    drift_amp        = parse(Float64, get(ENV, "RLJ_DRIFT", "0.0")),
    residual_scale   = parse(Float64, get(ENV, "RLJ_RESSCALE", "3.0")),
    rotate_consensus = true,
    planar_agents    = Set([5, 8, 10]),
)
train_and_save(cfg;
    total_steps = parse(Int, get(ENV, "RLJ_TOTAL", "400000")),
    n_envs      = parse(Int, get(ENV, "RLJ_NENVS", "32")),
    out         = get(ENV, "RLJ_OUT", joinpath(@__DIR__, "..", "cache", "rl_hetero.jld2")),
    label       = "heterogeneous sheaf (rotation + projection maps)")
