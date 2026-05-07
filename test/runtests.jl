using Test

@testset "Code Quality (Aqua.jl)" begin
  include("aqua.jl")
end

@testset "Network Sheaves" begin
  include("network_sheaves/ADT.jl")
  include("network_sheaves/Parser.jl")
  include("network_sheaves/SheafLaplacian.jl")
  include("network_sheaves/HarmonicExtension.jl")
  include("network_sheaves/Morphisms.jl")
  include("network_sheaves/GraphHomomorphisms.jl")
  include("network_sheaves/Pushforwards.jl")
end
