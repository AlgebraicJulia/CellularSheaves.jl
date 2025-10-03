using Test

@testset "Code Quality (Aqua.jl)" begin
  include("aqua.jl")
end

@testset "Network Sheaves" begin
  include("network_sheaves/ADT.jl")
  include("network_sheaves/Parser.jl")
end
