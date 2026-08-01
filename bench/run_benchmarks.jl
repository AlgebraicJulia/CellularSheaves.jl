include(joinpath(@__DIR__, "src", "CellularSheavesBenchmarks.jl"))

using .CellularSheavesBenchmarks

run_benchmarks_from_env!()
