using Test

@testset "Code Quality (Aqua.jl)" begin
  include("aqua.jl")
end

@testset "Network Sheaves" begin
  include("network_sheaves/Asynch.jl")
  include("network_sheaves/ADT.jl")
  include("network_sheaves/Parser.jl")
  include("network_sheaves/SheafLaplacian.jl")
  include("network_sheaves/HarmonicExtension.jl")
  include("network_sheaves/Morphisms.jl")
  include("network_sheaves/GraphHomomorphisms.jl")
  include("network_sheaves/Morphisms.jl")
  include("network_sheaves/Pushforwards.jl")
  include("network_sheaves/Pushouts.jl")
  include("network_sheaves/TrajectorySheaf.jl")
  include("network_sheaves/ControlledTrajectorySheaf.jl")
  include("network_sheaves/ControlledOptimalControl.jl")
  include("network_sheaves/NullspaceTrajectoryFamily.jl")
  include("network_sheaves/ControlledTrajectoryExamples.jl")
end
