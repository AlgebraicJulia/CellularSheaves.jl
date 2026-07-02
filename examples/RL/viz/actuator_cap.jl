#
# Render Demonstration 8 (hard actuator cap ‖u‖≤ū, no MPC) from a saved model. The cap is read
# from the model file (saved at training time).
#
#   julia --project=examples/RL examples/RL/viz/actuator_cap.jl
# Env: AJ_MODEL, AJ_TAG, AJ_OUTDIR
#
get(ENV, "DISPLAY", "") == "" && (ENV["GKSwstype"] = "100")
include(joinpath(@__DIR__, "..", "lib", "sheaf_rl.jl")); using .SheafRL
include(joinpath(@__DIR__, "..", "lib", "render.jl"))

render_all(get(ENV, "AJ_MODEL", joinpath(@__DIR__, "..", "cache", "rl_ballcap_u15.jld2"));
    tag    = get(ENV, "AJ_TAG", "u15"),
    outdir = get(ENV, "AJ_OUTDIR", joinpath(@__DIR__, "..", "cache", "viz_ballcap")),
    header = "hard actuator cap (no MPC)")
