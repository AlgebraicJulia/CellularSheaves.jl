# Run the deployed Tikhonov experiment from the repository's existing
# multi-agent target-tracking environment.

using Pkg

const ROOT = normpath(joinpath(@__DIR__, "..", ".."))
Pkg.activate(joinpath(ROOT, "examples", "multi_agent_target_tracking"))
include(joinpath(ROOT, "examples", "multi_agent_target_tracking", "tikhonov_experiment.jl"))
