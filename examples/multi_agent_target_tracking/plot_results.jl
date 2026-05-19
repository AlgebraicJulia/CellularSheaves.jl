# examples/multi_agent_target_tracking/plot_results.jl

include(joinpath(@__DIR__, "src", "visualization", "plotter.jl"))

using .Plotter

if abspath(PROGRAM_FILE) == abspath((@__FILE__)) || isinteractive()
    Plotter.results(show=true)
end
