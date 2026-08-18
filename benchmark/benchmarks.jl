include(joinpath(@__DIR__, "..", "bench", "src", "CellularSheavesBenchmarks.jl"))

using .CellularSheavesBenchmarks

const SUITE = filtered_suite_from_env()
