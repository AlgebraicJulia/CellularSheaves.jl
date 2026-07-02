#
# Demonstration 8: hard actuator constraint, NO MPC. Every control is projected onto the ball
# ‖u‖ ≤ ū (a second-order cone) — closed form, no solver, no horizon. The cone is the only place
# the constraint enters; Layer 1 stays a linear solve and Layer 2 is feedback + a projection.
#
#   julia --project=examples/RL examples/RL/train/actuator_cap.jl
# Env: RLJ_UCAP (cap ū), RLJ_DRIFT, RLJ_RESSCALE, RLJ_TOTAL, RLJ_NENVS, RLJ_OUT
#
include(joinpath(@__DIR__, "..", "lib", "sheaf_rl.jl"))
using .SheafRL

cfg = Config(
    drift_amp      = parse(Float64, get(ENV, "RLJ_DRIFT", "2.0")),
    actuator_cap   = parse(Float64, get(ENV, "RLJ_UCAP", "1.5")),
    residual_scale = parse(Float64, get(ENV, "RLJ_RESSCALE", "2.0")),
)
train_and_save(cfg;
    total_steps = parse(Int, get(ENV, "RLJ_TOTAL", "1500000")),
    n_envs      = parse(Int, get(ENV, "RLJ_NENVS", "32")),
    out         = get(ENV, "RLJ_OUT", joinpath(@__DIR__, "..", "cache", "rl_ballcap.jld2")),
    label       = "ball cap ‖u‖≤ū (conic, no MPC)")
