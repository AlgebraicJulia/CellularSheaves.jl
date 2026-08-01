include(joinpath(@__DIR__, "src", "CellularSheavesBenchmarks.jl"))

using .CellularSheavesBenchmarks

render_report_from_env!()
