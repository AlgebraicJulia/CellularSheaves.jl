using Test

@testset "Code Quality (Aqua.jl)" begin
  include("aqua.jl")
end

@testset "Docs Consistency" begin
  include("docs_consistency.jl")
end

@testset "Tikhonov Filter" begin
  include("ControlSheaves/Tikhonov.jl")
end

@testset "Block Sparse Arrays" begin
  include("block_sparse_arrays.jl")
end

@testset "Network Sheaves" begin
  include("network_sheaves/Asynch.jl")
  include("network_sheaves/ADT.jl")
  include("network_sheaves/Parser.jl")
  include("network_sheaves/RestrictionMaps.jl")
  include("network_sheaves/SheafLaplacian.jl")
  include("network_sheaves/HarmonicExtension.jl")
  include("network_sheaves/DistributedSolve.jl")
  include("network_sheaves/SchwarzMethods.jl")
  include("network_sheaves/GridSchwarz.jl")
  include("network_sheaves/GridMultigrid.jl")
  include("network_sheaves/MultilevelSchwarz.jl")
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
  include("network_sheaves/MultiAgentTracking.jl")
  include("network_sheaves/MultiAgentTrackingMPC.jl")
  include("network_sheaves/per_agent_and_target_dynamics.jl")
  include("network_sheaves/TrackingDSL.jl")
  include("network_sheaves/QuadraticCosts.jl")
  include("network_sheaves/EscortTracking.jl")
  include("network_sheaves/Formations.jl")
  include("network_sheaves/test_feedforward_control.jl")
  include("network_sheaves/TransferMapLift.jl")
  include("network_sheaves/LayeredEscort.jl")
end

@testset "Control Sheaves" begin
  include("ControlSheaves/AgentControllers.jl")
  include("ControlSheaves/NestedSystems.jl")
  include("ControlSheaves/NestedDSL.jl")
  include("ControlSheaves/CoordinationBenchmarks.jl")
  include("ControlSheaves/DoubleIntegratorHJB.jl")
  include("ControlSheaves/MechanicalHJB.jl")
  include("ControlSheaves/HJBMultigrid.jl")
end
