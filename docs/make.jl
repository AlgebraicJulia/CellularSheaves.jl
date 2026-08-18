# Render GR/Plots figures offscreen during the build so executing the examples
# never pops an on-screen gksqt window that has to be closed by hand. Must be
# set before anything loads Plots/GR (including the CellularSheavesPlots
# extension). "100" = offscreen file rendering; respects an explicit override.
get!(ENV, "GKSwstype", "100")

using Documenter
using Literate

const literate_dir = joinpath(@__DIR__, "literate")
const generated_dir = joinpath(@__DIR__, "src", "generated")

@info "Loading CellularSheaves"
using CellularSheaves

const no_literate = "--no-literate" in ARGS
if !no_literate
  @info "Building Literate.jl docs"

  # Set Literate.jl config if not being compiled on recognized service.
  config = Dict{String,String}()
  if !(haskey(ENV, "GITHUB_ACTIONS") || haskey(ENV, "GITLAB_CI"))
    config["nbviewer_root_url"] = "https://nbviewer.jupyter.org/github/AlgebraicJulia/CellularSheaves.jl/blob/gh-pages/dev"
    config["repo_root_url"] = "https://github.com/AlgebraicJulia/CellularSheaves.jl/blob/main/docs"
  end

  for (root, dirs, files) in walkdir(literate_dir)
    out_dir = joinpath(generated_dir, relpath(root, literate_dir))
    for file in files
      f, l = splitext(file)
      if l == ".jl" && !startswith(f, "_")
        src_path = joinpath(root, file)
        execute = !(contains(src_path, "asynch") || contains(src_path, "distributed") || contains(src_path, "escort.jl") || contains(src_path, "escort_feedforward") || contains(src_path, "scenario5"))
        Literate.markdown(src_path, out_dir;
          execute=execute, config=config, documenter=false, credit=false)
        Literate.notebook(src_path, out_dir;
          execute=false, documenter=false, credit=false)
      end
    end
  end
end

let src = joinpath(@__DIR__, "figures", "asynch"),
    dst = joinpath(@__DIR__, "src", "assets", "figures", "asynch")
  mkpath(dirname(dst))
  cp(src, dst; force=true)
end

let src = joinpath(@__DIR__, "figures", "distributed_solve"),
    dst = joinpath(@__DIR__, "src", "assets", "figures", "distributed_solve")
  if isdir(src)
    mkpath(dirname(dst))
    cp(src, dst; force=true)
  end
end

let src = joinpath(@__DIR__, "figures", "tikhonov"),
    dst = joinpath(@__DIR__, "src", "assets", "figures", "tikhonov")
  if isdir(src)
    mkpath(dirname(dst))
    cp(src, dst; force=true)
  end
end

@info "Building Documenter.jl docs"
makedocs(
  modules=[
    CellularSheaves,
    CellularSheaves.ControlSheaves,
    CellularSheaves.ControlSheaves.Tikhonov,
    CellularSheaves.ControlSheaves.AgentControllers,
    CellularSheaves.ControlSheaves.SafetyFilters,
    CellularSheaves.ControlSheaves.DistributedLayeredControl,
    CellularSheaves.ControlSheaves.Layered,
    CellularSheaves.ControlSheaves.NestedSystems,
    CellularSheaves.ControlSheaves.NestedDSL,
    CellularSheaves.ControlSheaves.NestedDSL.NestedDSLTerm,
    CellularSheaves.ControlSheaves.NestedDSL.NestedDSLParser,
    CellularSheaves.ControlSheaves.NestedDSL.NestedDSLValidator,
    CellularSheaves.ControlSheaves.NestedDSL.NestedDSLLowering,
    CellularSheaves.ControlSheaves.MultiAgentTracking,
    CellularSheaves.ControlSheaves.MultiAgentTracking.QuadraticCosts,
    CellularSheaves.ControlSheaves.CoordinationBenchmarks,
    CellularSheaves.ControlSheaves.TrackingDSL,
    CellularSheaves.ControlSheaves.TrackingDSL.TrackingDSLTerm,
    CellularSheaves.ControlSheaves.TrackingDSL.TrackingDSLParser,
    CellularSheaves.ControlSheaves.TrackingDSL.TrackingDSLValidator,
    CellularSheaves.ControlSheaves.TrackingDSL.TrackingDSLResolver,
    CellularSheaves.ControlSheaves.TrackingDSL.TrackingDSLLowering,
    CellularSheaves.AsynchSheaves,
    CellularSheaves.SheafInterface,
    CellularSheaves.NetworkSheaves.EuclideanSheaves,
    CellularSheaves.NetworkSheaves.GraphHomomorphisms,
    CellularSheaves.NetworkSheaves.SheafMorphisms,
    CellularSheaves.NetworkSheaves.Pushforwards,
    CellularSheaves.NetworkSheaves.Pushouts,
    CellularSheaves.NetworkSheaves.CellularSheafParser,
    CellularSheaves.NetworkSheaves.PotentialSheaves,
    CellularSheaves.NetworkSheaves.TrajectorySheaves,
    CellularSheaves.NetworkSheaves.DistributedSolve,
    CellularSheaves.NetworkSheaves.Formations,
    CellularSheaves.BlockSparseArrays,
    CellularSheaves.IPM,
  ],
  draft=false,
  format=Documenter.HTML(
    assets=["assets/benchtables.css"],
    size_threshold=nothing
  ),
  sitename="CellularSheaves.jl",
  doctest=false,
  checkdocs=:none,
  warnonly=[:cross_references],
  pages=Any[
    "CellularSheaves.jl"=>"index.md",
    "Examples"=>Any[
      "generated/literate_example.md",
      "generated/sheaf_morphisms.md",
      "generated/nullspace.md",
      "generated/nearest_global_section_iterative.md",
      "generated/pushforward.md",
      "generated/trajectory_sheaf.md",
      "Control Examples"=>Any[
        "generated/control/simple_integrator.md",
        "generated/control/controlled_double_integrator.md",
        "generated/control/controlled_vehicle_platoon.md",
        "generated/control/controlled_planar_quadrotor.md",
        "generated/control/controlled_mass_spring_damper_chain.md",
        "control/single_integrator_target_tracking.md",
        "generated/control/multi_quadrotor_target_tracking.md",
        "generated/control/mpc_target_tracking.md"],
      "Layered Control Architecture Examples" => Any[
        "generated/layered/distributed_harmonic_tracking.md",
        "generated/layered/tikhonov_harmonic_tracking.md",
        "generated/layered/layered_scenario5.md",
        "generated/layered/escort.md",
        "generated/layered/escort_feedforward.md",
        "generated/layered/escort_cbf_clf.md",
        "generated/layered/diffusion_vs_direct.md",
        "generated/layered/multilayer_escort.md",
        ],
      "Nested Systems Examples" => Any[
        "generated/nested/centroid_formation_tracking.md",
        "generated/nested/wheel_formation.md",
        "generated/nested/n_ring_formation.md",
        "generated/nested/rescaling_formation.md",
        "generated/nested/hierarchical_fleet.md",
        ],
      "Asynchronous Diffusion"=>Any[
        "generated/asynch/convergence_vs_delay.md",
        "generated/asynch/step_size_comparison.md",
        "generated/asynch/restriction_map_comparison.md",
        "generated/asynch/orthogonal_projection.md",
      ]
    ],
    "Feature Guides"=>Any[
      "features/core_sheaf_workflows.md",
      "features/morphisms_and_pushforwards.md",
      "features/trajectory_and_control.md",
      "Distributed Sheaf Solve"=>Any[
        "features/distributed_sheaf_solve/index.md",
        "features/distributed_sheaf_solve/theory_coordination.md",
        "features/distributed_sheaf_solve/theory_multifrontal.md",
        "features/distributed_sheaf_solve/theory_distributed.md",
        "features/distributed_sheaf_solve/comparison.md",
        "features/distributed_sheaf_solve/benchmarks.md",
        "features/distributed_sheaf_solve/api.md",
      ],
      "Tikhonov Harmonic Tracking"=>Any[
        "features/tikhonov_harmonic_tracking/index.md",
        "features/tikhonov_harmonic_tracking/harmonic_manifold.md",
        "features/tikhonov_harmonic_tracking/stability.md",
        "features/tikhonov_harmonic_tracking/feedforward.md",
        "features/tikhonov_harmonic_tracking/scenarios.md",
        "features/tikhonov_harmonic_tracking/api.md",
      ],
    ],
    "Benchmarks"=>Any[
      "benchmarks.md",
      "benchmark_report.md",
    ],
    "API Reference"=>Any[
      "api/index.md",
      "api/network_sheaves.md",
      "api/sheaf_interface.md",
      "api/euclidean_sheaves.md",
      "api/distributed_solve.md",
      "api/graph_homomorphisms.md",
      "api/sheaf_morphisms.md",
      "api/pushforwards.md",
      "api/pushouts.md",
      "api/dsl.md",
      "api/block_sparse_arrays.md",
      "api/potential_sheaves.md",
      "api/formations.md",
      "api/trajectory_sheaves.md",
      "api/tikhonov.md",
      "api/agent_controllers.md",
      "api/safety_filters.md",
      "api/distributed_layered_control.md",
      "api/layered.md",
      "api/nested_systems.md",
      "api/nested_dsl.md",
      "api/coordination_benchmarks.md",
      "api/herding_platoon.md",
      "api/multi_agent_tracking.md",
      "api/quadratic_costs.md",
      "api/asynch.md",
      "api/tracking_dsl.md",
      "api/ipm.md"
    ],
  ]
)

@info "Deploying docs"
deploydocs(
  target="build",
  repo="github.com/AlgebraicJulia/CellularSheaves.jl.git",
  branch="gh-pages",
  devbranch="main",
  push_preview=true
)
