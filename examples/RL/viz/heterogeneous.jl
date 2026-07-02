#
# Render Demonstration 7 (genuinely heterogeneous sheaf) from a saved model. Adds a 3-D animation
# (the only view where the planar/projection agents' off-z behaviour is visible). The rotation and
# projection maps are set here because older model files do not record them.
#
#   julia --project=examples/RL examples/RL/viz/heterogeneous.jl
# Env: AJ_MODEL, AJ_TAG, AJ_OUTDIR
#
get(ENV, "DISPLAY", "") == "" && (ENV["GKSwstype"] = "100")
include(joinpath(@__DIR__, "..", "lib", "sheaf_rl.jl")); using .SheafRL
include(joinpath(@__DIR__, "..", "lib", "render.jl"))

render_all(get(ENV, "AJ_MODEL", joinpath(@__DIR__, "..", "cache", "rl_hetero_f2.jld2"));
    tag    = get(ENV, "AJ_TAG", "f2"),
    outdir = get(ENV, "AJ_OUTDIR", joinpath(@__DIR__, "..", "cache", "viz_hetero")),
    header = "heterogeneous sheaf",
    three_d          = true,
    rotate_consensus = true,
    planar_agents    = Set([5, 8, 10]))
