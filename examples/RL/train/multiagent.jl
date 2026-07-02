#
# Demonstrations 4–5: the 13-agent / 4-target structure with identity restriction maps.
#   Flavor 1 (no drift):              RLJ_DRIFT=0
#   Flavor 2 (unknown time-varying):  RLJ_DRIFT=2.0
#
#   julia --project=examples/RL examples/RL/train/multiagent.jl
# Env: RLJ_DRIFT, RLJ_RESSCALE, RLJ_TOTAL, RLJ_NENVS, RLJ_OUT
#
include(joinpath(@__DIR__, "..", "lib", "sheaf_rl.jl"))
using .SheafRL

cfg = Config(
    drift_amp      = parse(Float64, get(ENV, "RLJ_DRIFT", "0.0")),
    residual_scale = parse(Float64, get(ENV, "RLJ_RESSCALE", "3.0")),
)
train_and_save(cfg;
    total_steps = parse(Int, get(ENV, "RLJ_TOTAL", "400000")),
    n_envs      = parse(Int, get(ENV, "RLJ_NENVS", "32")),
    out         = get(ENV, "RLJ_OUT", joinpath(@__DIR__, "..", "cache", "rl_multiagent.jld2")),
    label       = "multi-agent (identity restriction maps)")
