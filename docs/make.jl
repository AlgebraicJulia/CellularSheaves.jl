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
        Literate.markdown(src_path, out_dir;
          config=config, documenter=true, credit=false)
        execute = !contains(src_path, "asynch")
        Literate.notebook(src_path, out_dir;
          execute=execute, documenter=true, credit=false)
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

@info "Building Documenter.jl docs"
makedocs(
  modules=[CellularSheaves, CellularSheaves.ControlSheaves.MultiAgentTracking, CellularSheaves.ControlSheaves.MultiAgentTracking.QuadraticCosts, CellularSheaves.AsynchSheaves],
  draft=false,
  format=Documenter.HTML(assets=["assets/benchtables.css"]),
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
        "generated/control/mpc_target_tracking.md",
        "generated/control/distributed_harmonic_tracking.md"

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
      "api/trajectory_sheaves.md",
      "api/herding_platoon.md",
      "api/multi_agent_tracking.md",
      "api/quadratic_costs.md",
      "api/asynch.md"
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
