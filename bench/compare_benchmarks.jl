include(joinpath(@__DIR__, "src", "CellularSheavesBenchmarks.jl"))

using .CellularSheavesBenchmarks

compare_benchmarks_from_env!()
