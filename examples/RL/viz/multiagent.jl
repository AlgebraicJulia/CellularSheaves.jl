#
# Render the multi-agent formation demonstration (identity restriction maps) from a saved model.
#   julia --project=. examples/RL/viz/multiagent.jl
# Env: AJ_MODEL (policy .jld2), AJ_TAG (output label), AJ_OUTDIR, AJ_DRIFT (override drift
#      amplitude — set to 0 for the nominal/no-disturbance rollout).
#
get(ENV, "DISPLAY", "") == "" && (ENV["GKSwstype"] = "100")
include(joinpath(@__DIR__, "..", "lib", "sheaf_rl.jl")); using .SheafRL
include(joinpath(@__DIR__, "..", "lib", "render.jl"))

overrides = haskey(ENV, "AJ_DRIFT") ? (; drift_amp = parse(Float64, ENV["AJ_DRIFT"])) : (;)

render_all(get(ENV, "AJ_MODEL", joinpath(@__DIR__, "..", "cache", "rl_multiagent_f2.jld2"));
    tag    = get(ENV, "AJ_TAG", "f2"),
    outdir = get(ENV, "AJ_OUTDIR", joinpath(@__DIR__, "..", "cache", "viz_multiagent")),
    header = "multi-agent structure (13 agents, 4 targets)",
    overrides...)
