# ACC26 Experiment Runner
#
# Generates the pre-computed CSV data consumed by the literate docs pages in
# docs/literate/asynch/. Run this script from the project root to regenerate
# all experiment data:
#
#   julia --project=docs docs/scripts/acc26_experiments.jl
#
# Output is written to docs/data/asynch/. Subdirectories are created as needed.
# The script uses Random.seed!(8) for reproducibility.
#
# NOTE: These experiments use the new AsynchSheaves API introduced in this PR,
# where compute_trajectory and compute_trajectory_asynch take an EuclideanSheaf
# as their first argument rather than a pre-computed Laplacian matrix.

using CellularSheaves
using Distributions
using Random
using Graphs
using LinearAlgebra
using Plots

#Random.seed!(8)
Random.seed!(42)

const data_path = joinpath(@__DIR__, "../data/asynch/")

function ensure_dir(path)
    isdir(path) || mkpath(path)
end

# ---------------------------------------------------------------------------
# Convergence vs. delay experiment
# Generates: data_path/converge_vs_delay/<topology>/{loss,error,proj_error}<i>.csv
# ---------------------------------------------------------------------------

function generate_convergence_trajectories(file_path::String, s::EuclideanSheaf)
    ensure_dir(file_path)
    n_agents = length(vertex_stalks(s))
    L = sheaf_laplacian_matrix_direct(s)
    K = opnorm(Array(L), 2)
    γ = 1 / K
    x0 = rand(5.0:0.001:15.0, 4 * n_agents)
    T = Int(1e5)
    f = energy_function(L)
    x_orthog = nearest_global_section(s, x0)

    Bs = [1, 10, 50, 100, 200]
    for (i, B) in enumerate(Bs)
        update_d1 = Normal(0.05 * B, 0.005 * B)
        update_d2 = Normal(0.5 * B, 0.05 * B)
        update_mixture = MixtureModelParams([update_d1, update_d2], [0.5, 0.5])

        broadcast_d1 = Normal(0.1 * B, 0.01 * B)
        broadcast_d2 = Normal(0.8 * B, 0.05 * B)
        broadcast_mixture = MixtureModelParams([broadcast_d1, broadcast_d2], [0.5, 0.5])

        traj = compute_trajectory_asynch(s, x0, γ, update_mixture, broadcast_mixture; max_iters=T, B=B)
        loss = f.(traj)
        save_trajectory(joinpath(file_path, "loss$i.csv"), loss)

        x_star = nearest_global_section(s, traj[end])
        rel_error(x) = norm(x - x_star) / norm(x_star)
        save_trajectory(joinpath(file_path, "error$i.csv"), rel_error.(traj))

        proj_error(x) = norm(x - x_orthog) / norm(x_orthog)
        save_trajectory(joinpath(file_path, "proj_error$i.csv"), proj_error.(traj))
    end
end

function generate_all_convergence_trajectories(file_path::String)
    ensure_dir(file_path)
    n_agents = 20
    gs = [random_regular_graph(20, 4), erdos_renyi(20, 0.3), star_graph(20)]
    rm_generator(stalk_dim) = rand(1, stalk_dim)
    sheaves = [sheaf_from_graph(g, 4, rm_generator) for g in gs]
    paths = [joinpath(file_path, "regular_graph/"), joinpath(file_path, "er_graph/"), joinpath(file_path, "complete_graph/")]
    for (s, p) in zip(sheaves, paths)
        generate_convergence_trajectories(p, s)
    end
end

# ---------------------------------------------------------------------------
# Step size comparison experiment
# Generates: data_path/ss_vs_converge/{loss_synch,loss_synch_ss,loss_asynch_ss}.csv
# ---------------------------------------------------------------------------

function generate_stepsize_trajectories(file_path::String)
    ensure_dir(file_path)
    n_agents = 20
    g = random_regular_graph(n_agents, 4)
    s = sheaf_from_graph(g, 4, dim -> random_semi_orthogonal_matrix(2, dim))
    L = sheaf_laplacian_matrix_direct(s)
    K = opnorm(Array(L), 2)
    γ1 = 1 / K
    B = 50
    γ2 = 0.1 * γ1
    x0 = rand(5.0:0.001:15.0, 4 * n_agents)
    T = Int(1e6)
    f = energy_function(L)

    update_d1 = Normal(0.05 * B, 0.005 * B)
    update_d2 = Normal(0.5 * B, 0.05 * B)
    update_mixture = MixtureModelParams([update_d1, update_d2], [0.8, 0.2])

    broadcast_d1 = Normal(0.1 * B, 0.01 * B)
    broadcast_d2 = Normal(0.8 * B, 0.05 * B)
    broadcast_mixture = MixtureModelParams([broadcast_d1, broadcast_d2], [0.5, 0.5])

    trajs = compute_trajectory_asynch(s, x0, [γ1, γ2], update_mixture, broadcast_mixture; max_iters=T, B=B)
    losses = [f.(traj) for traj in trajs]
    save_trajectory(joinpath(file_path, "loss_synch_ss.csv"), losses[1])
    save_trajectory(joinpath(file_path, "loss_asynch_ss.csv"), losses[2])

    traj_synch = compute_trajectory(s, x0, γ1; max_iters=T)
    save_trajectory(joinpath(file_path, "loss_synch.csv"), f.(traj_synch))
end

# ---------------------------------------------------------------------------
# Restriction map type experiment
# Generates: data_path/rm_experiment/{constant_sheaf,vector_bundle,matrix_weighted}/...
# ---------------------------------------------------------------------------

function run_rm_experiment(file_path::String)
    ensure_dir(file_path)
    n_agents = 20
    stalk_dim = 4
    g = random_regular_graph(n_agents, 4)
    rm_identity(d) = Matrix{Float64}(I(d))
    rm_rand(d) = rand(1, d)
    rm_matrix = matrix_weighted_edge_generator(; pd_prob=0.8)

    sheaves = [
        sheaf_from_graph(g, stalk_dim, rm_identity),
        sheaf_from_graph(g, stalk_dim, rm_rand),
        sheaf_from_graph(g, stalk_dim, rm_matrix; symmetric_edges=true),
    ]
    paths = [
        joinpath(file_path, "constant_sheaf/"),
        joinpath(file_path, "vector_bundle/"),
        joinpath(file_path, "matrix_weighted/"),
    ]
    for (s, p) in zip(sheaves, paths)
        generate_convergence_trajectories(p, s)
    end
end

# ---------------------------------------------------------------------------
# Orthogonal projection / convergence bounds experiment
# Generates: data_path/op_experiment/results.csv
# ---------------------------------------------------------------------------

function run_op_experiment(file_path::String)
    ensure_dir(file_path)
    n_agents = 20
    g = random_regular_graph(n_agents, 4)
    s = sheaf_from_graph(g, 4, dim -> random_semi_orthogonal_matrix(1, dim))
    L = sheaf_laplacian_matrix_direct(s)
    K = opnorm(Array(L), 2)
    γ = 1 / K
    T = Int(1e5)

    x0 = rand(Normal(0, 10), 4 * n_agents)
    x_star = nearest_global_section(s, x0)

    results = Float64[]
    for i in 1:15
        B = 2^(i - 1)
        sols = []
        for _ in 1:3
            traj = compute_trajectory_asynch(
                s, x0, γ,
                ones(Int, n_agents), repeat([B], n_agents);
                B=B, max_iters=T,
            )
            push!(sols, traj[end])
        end
        r = mean(norm(sol - x_star) for sol in sols)
        push!(results, r)
    end
    save_trajectory(joinpath(file_path, "results.csv"), results)
end

# ---------------------------------------------------------------------------
# Run all experiments
# ---------------------------------------------------------------------------

@info "Generating convergence vs. delay data..."
generate_all_convergence_trajectories(joinpath(data_path, "converge_vs_delay/"))

@info "Generating step size comparison data..."
generate_stepsize_trajectories(joinpath(data_path, "ss_vs_converge/"))

@info "Generating restriction map type data..."
run_rm_experiment(joinpath(data_path, "rm_experiment/"))

@info "Generating orthogonal projection / convergence bounds data..."
run_op_experiment(joinpath(data_path, "op_experiment/"))

@info "Done. Data written to $data_path"
