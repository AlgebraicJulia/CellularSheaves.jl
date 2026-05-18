"""
Module for simulating asynchronous cellular sheaf diffusion.

Agents communicate over a sheaf's underlying graph using a gossip-like protocol:
each agent maintains a local estimate of the global cochain, periodically
applies a gradient descent step on its own component, and periodically
broadcasts its local state to neighboring agents.

## Scheduling models

Three scheduling models govern *when* each agent updates and broadcasts:

- **Probabilistic** ([`ProbabilisticModelParams`](@ref)): each agent independently
  decides to update or broadcast at each iteration with random probabilities.
- **Deterministic** (period vectors): each agent `i` updates every
  `update_periods[i]` iterations and broadcasts every `broadcast_periods[i]`
  iterations, with a random phase offset.
- **Mixture** ([`MixtureModelParams`](@ref)): periods are sampled from a
  Gaussian mixture at the start of the simulation and held fixed thereafter.

## Key functions

- [`compute_trajectory`](@ref) — synchronous gradient descent baseline.
- [`compute_trajectory_asynch`](@ref) — asynchronous diffusion (all scheduling
  variants; single-run and batch forms).

## Utilities

- [`random_psd`](@ref), [`random_pd`](@ref), [`random_semi_orthogonal_matrix`](@ref)
  — random matrix generators for constructing test sheaves.
- [`matrix_weighted_rm`](@ref), [`matrix_weighted_edge_generator`](@ref)
  — helpers for building matrix-weighted sheaves.
- [`save_trajectory`](@ref), [`save_trajectories`](@ref), [`load_trajectory`](@ref)
  — CSV I/O for trajectory data.
- [`empty_experiment_plot`](@ref), [`plot_loss_curve!`](@ref),
  [`plot_log_loss_curve!`](@ref) — plotting helpers (requires `using Plots`).
"""
module AsynchSheaves

export MixtureModelParams, ProbabilisticModelParams
export compute_trajectory, compute_trajectory_asynch
export random_psd, random_pd, random_semi_orthogonal_matrix
export matrix_weighted_rm, matrix_weighted_edge_generator
export save_trajectory, save_trajectories, load_trajectory
export empty_experiment_plot, plot_loss_curve!, plot_log_loss_curve!

using ArgCheck: @argcheck
using BlockArrays
using CSV
using Distributions: Normal, MixtureModel
using Tables
using LinearAlgebra

using ..EuclideanSheaves: EuclideanSheaf, sheaf_laplacian_matrix, energy_function
using ..SheafInterface: vertex_stalks

# ---------------------------------------------------------------------------
# Parameter structs
# ---------------------------------------------------------------------------

"""
    MixtureModelParams(dists, weights)

Parameters for a Gaussian mixture distribution over agent periods.

Each agent's update or broadcast period is sampled by:
1. Selecting component `k` with probability `weights[k]`.
2. Drawing from `dists[k]` (a `Normal`), rounding up to the nearest integer.
3. Clamping to `[1, B]`.

Construct with:
```julia
using Distributions: Normal
d1 = Normal(5.0, 0.5)
d2 = Normal(50.0, 5.0)
params = MixtureModelParams([d1, d2], [0.5, 0.5])
```
"""
struct MixtureModelParams
    dists::Vector{Normal{Float64}}
    weights::Vector{Float64}
end

"""
    ProbabilisticModelParams(update_prob_upper_bound, broadcast_prob_upper_bound)

Parameters for coin-flip-based update and broadcast scheduling.

At each iteration, agent `i` updates with probability `pᵢ` drawn uniformly
from `[0, update_prob_upper_bound]`, and broadcasts with probability `qᵢ` drawn
from `[0, broadcast_prob_upper_bound]`. These probabilities are sampled once at
the start of the simulation.
"""
struct ProbabilisticModelParams
    update_prob_upper_bound::Float64
    broadcast_prob_upper_bound::Float64
end

# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

function _block_structure(sheaf::EuclideanSheaf)
    stalk_dims = vertex_stalks(sheaf)
    return stalk_dims, length(stalk_dims)
end

function _sample_mixture_periods(model::MixtureModelParams, nblocks::Int, B::Int)
    mixture = MixtureModel(model.dists, model.weights)
    periods = ceil.(Int, rand(mixture, nblocks))
    periods = clamp.(periods, 1, B)
    return periods
end

# ---------------------------------------------------------------------------
# Random matrix generators
# ---------------------------------------------------------------------------

"""
    random_psd(dim) -> Matrix{Float64}

Return a random `dim × dim` positive semidefinite matrix.
"""
function random_psd(dim::Int)
    @argcheck dim >= 1
    A = rand(ceil(Int, dim / 2), dim)
    return A' * A
end

"""
    random_pd(dim) -> Matrix{Float64}

Return a random `dim × dim` positive definite matrix.
"""
function random_pd(dim::Int)
    @argcheck dim >= 1
    A = rand(ceil(Int, dim / 2), dim)
    return A' * A + I(dim)
end

"""
    random_semi_orthogonal_matrix(n, m) -> Matrix{Float64}

Return a random `n × m` semi-orthogonal matrix (columns are orthonormal when
`n ≥ m`, rows are orthonormal when `m > n`).
"""
function random_semi_orthogonal_matrix(n::Int, m::Int)
    @argcheck n >= 1 && m >= 1
    N = max(n, m)
    A = rand(N, N)
    Q, _ = qr(A)
    return Matrix(Q)[1:n, 1:m]
end

"""
    matrix_weighted_rm(A) -> Matrix{Float64}

Extract a restriction map from `A` via QR factorization (returns the R factor).
"""
function matrix_weighted_rm(A::AbstractMatrix)
    _, U = qr(A)
    return Matrix(U)
end

"""
    matrix_weighted_edge_generator(; pd_prob=0.5) -> Function

Return a function `stalk_dim -> Matrix` that generates restriction maps from
either a random PD or PSD matrix with probability `pd_prob` and `1-pd_prob`
respectively.

Useful as an `rm_generator` argument to `sheaf_from_graph`.
"""
function matrix_weighted_edge_generator(; pd_prob::Float64=0.5)
    @argcheck 0.0 <= pd_prob <= 1.0
    return stalk_dim -> begin
        A = rand() < pd_prob ? random_pd(stalk_dim) : random_psd(stalk_dim)
        return matrix_weighted_rm(A)
    end
end

# ---------------------------------------------------------------------------
# I/O utilities
# ---------------------------------------------------------------------------

"""
    save_trajectory(filename, trajectory)

Write a trajectory (vector of scalars or vectors) to a CSV file at `filename`.
"""
function save_trajectory(filename::AbstractString, trajectory)
    CSV.write(filename, Tables.table(trajectory))
end

"""
    save_trajectories(path, experiment_name, trajectories)

Write each trajectory in `trajectories` to a separate CSV file under `path`,
named `<experiment_name>_traj1.csv`, `<experiment_name>_traj2.csv`, etc.
"""
function save_trajectories(path::AbstractString, experiment_name::AbstractString, trajectories)
    for (i, t) in enumerate(trajectories)
        save_trajectory(joinpath(path, "$(experiment_name)_traj$i.csv"), t)
    end
end

"""
    load_trajectory(filename) -> Matrix

Read a trajectory from a CSV file at `filename`. Returns a matrix (rows are
iterations, columns are state components).
"""
function load_trajectory(filename::AbstractString)
    return CSV.File(filename) |> Tables.matrix
end

# ---------------------------------------------------------------------------
# Synchronous baseline
# ---------------------------------------------------------------------------

"""
    compute_trajectory(sheaf, x0, γ; max_iters=1000, tol=1e-8) -> Vector

Run synchronous gradient descent on `sheaf_laplacian_matrix(sheaf)` from
initial cochain `x0` with step size `γ`.

Returns a vector of cochains (one per iteration, including `x0`). Stops early
when energy drops below `tol`.
"""
function compute_trajectory(
    sheaf::EuclideanSheaf,
    x0::AbstractVector,
    γ::Real;
    max_iters::Int=1000,
    tol::Real=1e-8,
)
    @argcheck γ > 0
    @argcheck max_iters >= 1
    @argcheck tol > 0
    L = sheaf_laplacian_matrix(sheaf)
    f = energy_function(L)
    traj = [Vector(x0)]
    x_curr = Vector(x0)
    for _ in 1:max_iters
        x_curr = x_curr - γ * (L * x_curr)
        push!(traj, copy(x_curr))
        f(x_curr) < tol && break
    end
    return traj
end

# ---------------------------------------------------------------------------
# Asynchronous diffusion — internal core
# ---------------------------------------------------------------------------

function _asynch_core!(
    L::AbstractMatrix,
    global_state::BlockVector,
    local_states::BlockMatrix,
    γ::Real,
    update_periods::Vector{Int},
    broadcast_periods::Vector{Int},
    update_phases::Vector{Int},
    broadcast_phases::Vector{Int},
    f,
    stalk_dims::Vector{Int},
    nblocks::Int;
    max_iters::Int,
    tol::Real,
)
    traj = [Vector(global_state)]
    for t in 1:max_iters
        g = BlockArray(L * local_states, stalk_dims, ones(Int, nblocks))
        for i in 1:nblocks
            if t % update_periods[i] == update_phases[i]
                local_states[Block(i), Block(i)] -= γ * g[Block(i), Block(i)]
                global_state[Block(i)] = vec(local_states[Block(i), Block(i)])
                update_phases[i] = rand(0:update_periods[i]-1)
            end
            if t % broadcast_periods[i] == broadcast_phases[i]
                x = copy(local_states[Block(i), Block(i)])
                local_states[Block(i), :] .= x
                broadcast_phases[i] = rand(0:broadcast_periods[i]-1)
            end
        end
        push!(traj, Vector(global_state))
        f(traj[end]) < tol && break
    end
    return traj
end

function _init_local_states(x0::AbstractVector, stalk_dims::Vector{Int}, nblocks::Int)
    global_state = BlockArray(Vector(x0), stalk_dims)
    local_states = BlockArray(
        repeat(Vector(x0), 1, nblocks),
        stalk_dims,
        ones(Int, nblocks),
    )
    return global_state, local_states
end

# ---------------------------------------------------------------------------
# Asynchronous diffusion — public API
# ---------------------------------------------------------------------------

"""
    compute_trajectory_asynch(sheaf, x0, γ, params::ProbabilisticModelParams;
                               max_iters=1000, tol=1e-8, B=50) -> Vector

Asynchronous sheaf diffusion with coin-flip scheduling.

At the start of the simulation each agent `i` is assigned update probability
`pᵢ ∈ [0, params.update_prob_upper_bound]` and broadcast probability
`qᵢ ∈ [0, params.broadcast_prob_upper_bound]`, drawn uniformly at random.
At each iteration, agent `i` independently flips these coins.

Returns a vector of global cochains (one per iteration).
"""
function compute_trajectory_asynch(
    sheaf::EuclideanSheaf,
    x0::AbstractVector,
    γ::Real,
    params::ProbabilisticModelParams;
    max_iters::Int=1000,
    tol::Real=1e-8,
    B::Int=50,
)
    @argcheck γ > 0
    @argcheck max_iters >= 1
    @argcheck 0 < params.update_prob_upper_bound <= 1
    @argcheck 0 < params.broadcast_prob_upper_bound <= 1
    L = sheaf_laplacian_matrix(sheaf)
    stalk_dims, nblocks = _block_structure(sheaf)
    @argcheck length(x0) == sum(stalk_dims)
    f = energy_function(L)

    broadcast_probs = rand(nblocks) .* params.broadcast_prob_upper_bound
    update_probs = rand(nblocks) .* params.update_prob_upper_bound

    global_state = BlockArray(Vector(x0), stalk_dims)
    local_states = BlockArray(repeat(Vector(x0), 1, nblocks), stalk_dims, ones(Int, nblocks))

    traj = [Vector(x0)]
    for t in 1:max_iters
        g = BlockArray(L * local_states, stalk_dims, ones(Int, nblocks))
        for i in 1:nblocks
            if rand() < update_probs[i]
                local_states[Block(i), Block(i)] -= γ * g[Block(i), Block(i)]
                global_state[Block(i)] = vec(local_states[Block(i), Block(i)])
            end
            if rand() < broadcast_probs[i]
                x = copy(local_states[Block(i), Block(i)])
                local_states[Block(i), :] .= x
            end
        end
        push!(traj, Vector(global_state))
        f(traj[end]) < tol && break
    end
    return traj
end

"""
    compute_trajectory_asynch(sheaf, x0, γ,
                               update_periods::Vector{Int},
                               broadcast_periods::Vector{Int};
                               max_iters=1000, tol=1e-8, B=50) -> Vector

Asynchronous sheaf diffusion with deterministic periodic scheduling.

Agent `i` updates every `update_periods[i]` iterations and broadcasts every
`broadcast_periods[i]` iterations, starting from a random phase in `[0, period-1]`.

Returns a vector of global cochains (one per iteration).
"""
function compute_trajectory_asynch(
    sheaf::EuclideanSheaf,
    x0::AbstractVector,
    γ::Real,
    update_periods::AbstractVector{Int},
    broadcast_periods::AbstractVector{Int};
    max_iters::Int=1000,
    tol::Real=1e-8,
    B::Int=50,
)
    @argcheck γ > 0
    @argcheck max_iters >= 1
    L = sheaf_laplacian_matrix(sheaf)
    stalk_dims, nblocks = _block_structure(sheaf)
    @argcheck length(x0) == sum(stalk_dims)
    @argcheck length(update_periods) == nblocks
    @argcheck length(broadcast_periods) == nblocks
    f = energy_function(L)

    update_phases = [rand(0:update_periods[i]-1) for i in 1:nblocks]
    broadcast_phases = [rand(0:broadcast_periods[i]-1) for i in 1:nblocks]

    global_state, local_states = _init_local_states(x0, stalk_dims, nblocks)
    return _asynch_core!(
        L, global_state, local_states, γ,
        Vector(update_periods), Vector(broadcast_periods),
        update_phases, broadcast_phases,
        f, stalk_dims, nblocks;
        max_iters=max_iters, tol=tol,
    )
end

"""
    compute_trajectory_asynch(sheaf, x0, γ,
                               update_model::MixtureModelParams,
                               broadcast_model::MixtureModelParams;
                               max_iters=1000, tol=1e-8, B=50) -> Vector

Asynchronous sheaf diffusion with mixture-model period scheduling.

Agent periods are sampled once at the start from `update_model` and
`broadcast_model` (Gaussian mixture distributions), then clamped to `[1, B]`.
Each agent starts at a random phase within its period.

Returns a vector of global cochains (one per iteration).
"""
function compute_trajectory_asynch(
    sheaf::EuclideanSheaf,
    x0::AbstractVector,
    γ::Real,
    update_model::MixtureModelParams,
    broadcast_model::MixtureModelParams;
    max_iters::Int=1000,
    tol::Real=1e-8,
    B::Int=50,
)
    @argcheck γ > 0
    @argcheck max_iters >= 1
    @argcheck B >= 1
    L = sheaf_laplacian_matrix(sheaf)
    stalk_dims, nblocks = _block_structure(sheaf)
    @argcheck length(x0) == sum(stalk_dims)
    f = energy_function(L)

    update_periods = _sample_mixture_periods(update_model, nblocks, B)
    broadcast_periods = _sample_mixture_periods(broadcast_model, nblocks, B)
    update_phases = [rand(0:update_periods[i]-1) for i in 1:nblocks]
    broadcast_phases = [rand(0:broadcast_periods[i]-1) for i in 1:nblocks]

    global_state, local_states = _init_local_states(x0, stalk_dims, nblocks)
    return _asynch_core!(
        L, global_state, local_states, γ,
        update_periods, broadcast_periods,
        update_phases, broadcast_phases,
        f, stalk_dims, nblocks;
        max_iters=max_iters, tol=tol,
    )
end

"""
    compute_trajectory_asynch(sheaf, x0s::Vector{<:AbstractVector}, γ,
                               update_model, broadcast_model;
                               max_iters=1000, tol=1e-8, B=50) -> Vector{Vector}

Batch asynchronous diffusion over multiple initial cochains.

A single schedule (periods and phases) is sampled once and reused across all
initial conditions in `x0s`. This isolates the effect of initialization from
the effect of the schedule.

Returns a vector of trajectories, one per element of `x0s`.
"""
function compute_trajectory_asynch(
    sheaf::EuclideanSheaf,
    x0s::AbstractVector{<:AbstractVector},
    γ::Real,
    update_model::MixtureModelParams,
    broadcast_model::MixtureModelParams;
    max_iters::Int=1000,
    tol::Real=1e-8,
    B::Int=50,
)
    @argcheck γ > 0
    @argcheck max_iters >= 1
    @argcheck B >= 1
    @argcheck !isempty(x0s)
    L = sheaf_laplacian_matrix(sheaf)
    stalk_dims, nblocks = _block_structure(sheaf)
    f = energy_function(L)

    update_periods = _sample_mixture_periods(update_model, nblocks, B)
    broadcast_periods = _sample_mixture_periods(broadcast_model, nblocks, B)

    return map(x0s) do x0
        @argcheck length(x0) == sum(stalk_dims)
        update_phases = [rand(0:update_periods[i]-1) for i in 1:nblocks]
        broadcast_phases = [rand(0:broadcast_periods[i]-1) for i in 1:nblocks]
        global_state, local_states = _init_local_states(x0, stalk_dims, nblocks)
        _asynch_core!(
            L, global_state, local_states, γ,
            copy(update_periods), copy(broadcast_periods),
            update_phases, broadcast_phases,
            f, stalk_dims, nblocks;
            max_iters=max_iters, tol=tol,
        )
    end
end

"""
    compute_trajectory_asynch(sheaf, x0, γs::AbstractVector{<:Real},
                               update_model, broadcast_model;
                               max_iters=1000, tol=1e-8, B=50) -> Vector{Vector}

Batch asynchronous diffusion over multiple step sizes.

A single schedule (periods and phases) is sampled once and reused across all
step sizes in `γs`. This isolates the effect of the step size from
the effect of the schedule.

Returns a vector of trajectories, one per element of `γs`.
"""
function compute_trajectory_asynch(
    sheaf::EuclideanSheaf,
    x0::AbstractVector,
    γs::AbstractVector{<:Real},
    update_model::MixtureModelParams,
    broadcast_model::MixtureModelParams;
    max_iters::Int=1000,
    tol::Real=1e-8,
    B::Int=50,
)
    @argcheck all(>(0), γs)
    @argcheck max_iters >= 1
    @argcheck B >= 1
    @argcheck !isempty(γs)
    L = sheaf_laplacian_matrix(sheaf)
    stalk_dims, nblocks = _block_structure(sheaf)
    @argcheck length(x0) == sum(stalk_dims)
    f = energy_function(L)

    update_periods = _sample_mixture_periods(update_model, nblocks, B)
    broadcast_periods = _sample_mixture_periods(broadcast_model, nblocks, B)

    return map(γs) do γ
        update_phases = [rand(0:update_periods[i]-1) for i in 1:nblocks]
        broadcast_phases = [rand(0:broadcast_periods[i]-1) for i in 1:nblocks]
        global_state, local_states = _init_local_states(x0, stalk_dims, nblocks)
        _asynch_core!(
            L, global_state, local_states, γ,
            copy(update_periods), copy(broadcast_periods),
            update_phases, broadcast_phases,
            f, stalk_dims, nblocks;
            max_iters=max_iters, tol=tol,
        )
    end
end

# ---------------------------------------------------------------------------
# Plotting (implemented in ext/CellularSheavesPlots.jl when Plots is loaded)
# ---------------------------------------------------------------------------

"""
    empty_experiment_plot(; kwargs...) -> Plot

Create a blank plot with standard experiment styling.
Requires `using Plots` before calling.
"""
function empty_experiment_plot end

"""
    plot_loss_curve!(plt, losses, label; kwargs...)

Add a linear-scale loss curve to `plt`.
Requires `using Plots` before calling.
"""
function plot_loss_curve! end

"""
    plot_log_loss_curve!(plt, losses, label; kwargs...)

Add a log₁₀-scale loss curve to `plt`.
Requires `using Plots` before calling.
"""
function plot_log_loss_curve! end

end # module AsynchSheaves
