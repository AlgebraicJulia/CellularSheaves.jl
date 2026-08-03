module DistributedLayeredControl

using Distributed
using LinearAlgebra
using SparseArrays
using CliqueTrees.Multifrontal
using ...NetworkSheaves: EuclideanSheaf
using ...NetworkSheaves.EuclideanSheaves: _harmonic_extension_restricted_laplacian
using ...NetworkSheaves.DistributedSolve: partition_tree, distributed_tree_solve
using ..Tikhonov: tikhonov_step!
using ..AgentControllers: AbstractAgentState, AgentState, FeedforwardAgentState,
                           AbstractAgentDynamics, QuadrotorDynamics,
                           AbstractAgentController, FeedforwardLQRController, step_agent!

export init_distributed_agents!, run_layered_simulation, LayeredControlProblem, LayeredSimulationResult, animate_layered_escort, animate_scenario5

function animate_layered_escort end
function animate_scenario5 end

# Worker-local state storage. 
# We use a Dict to allow a worker to potentially simulate multiple agents if NA > nworkers().
const LOCAL_AGENTS = Dict{Int, AbstractAgentState}()


"""
    LayeredControlProblem

Problem definition for a distributed layered control simulation.

`pos_dim` controls how many leading dimensions of each agent's stalk reference vector
are passed to the agent's physical controller. For a sheaf using homogeneous/affine stalks
(e.g. SE(3) with D=4) set `pos_dim = D-1` to drop the homogeneous coordinate. For a pure
Euclidean sheaf (e.g. D=2 spatial positions) set `pos_dim = D` to pass all coordinates.
Defaults to `D-1` for backward compatibility with the escort mission.
"""
struct LayeredControlProblem
    sheaf::EuclideanSheaf{Float64}
    target_nodes::Vector{Int}
    target_trajectory_func::Any
    target_velocity_func::Any
    agent_configs::Vector
    dt::Float64
    steps::Int
    r_ring::Float64
    pos_dim::Int
end

# Positional constructor: coerces agent_configs, infers pos_dim as D-1 (escort convention)
function LayeredControlProblem(sheaf, target_nodes, target_trajectory_func, agent_configs, dt, steps, r_ring; target_velocity_func=nothing)
    D = sheaf.vertex_stalks[1]
    return LayeredControlProblem(sheaf, target_nodes, target_trajectory_func, target_velocity_func,
        collect(agent_configs), dt, steps, r_ring, D - 1)
end

# Keyword constructor: accepts optional r_ring, target_velocity_func, and explicit pos_dim
function LayeredControlProblem(; sheaf, target_nodes, target_trajectory_func, agent_configs, dt, steps, r_ring=0.0, target_velocity_func=nothing, pos_dim=nothing)
    D = sheaf.vertex_stalks[1]
    resolved_pos_dim = isnothing(pos_dim) ? D - 1 : pos_dim
    return LayeredControlProblem(sheaf, target_nodes, target_trajectory_func, target_velocity_func,
        collect(agent_configs), dt, steps, r_ring, resolved_pos_dim)
end

struct LayeredSimulationResult
    problem::LayeredControlProblem
    sim_data::Array{Float64, 3}
    qstar_history::Array{Float64, 3}
end

"""
    _init_agent_on_worker!(agent_id::Int, x0::Vector{Float64}, dyn::AbstractAgentDynamics, dt::Float64, K_lqr::Matrix{Float64}, eps::Float64)

Internal function called on worker to initialize standard feedback agent state.
"""
function _init_agent_on_worker!(agent_id::Int, x0::Vector{Float64}, dyn::AbstractAgentDynamics, dt::Float64, K_lqr::Matrix{Float64}, eps::Float64)
    LOCAL_AGENTS[agent_id] = AgentState(x0, dyn, dt, K_lqr, eps)
    return nothing
end

"""
    _init_feedforward_agent_on_worker!(agent_id::Int, x0::Vector{Float64}, dyn::AbstractAgentDynamics, dt::Float64, ctrl::FeedforwardLQRController, eps::Float64)

Internal function called on worker to initialize feedforward agent state.
"""
function _init_feedforward_agent_on_worker!(agent_id::Int, x0::Vector{Float64}, dyn::AbstractAgentDynamics, dt::Float64, ctrl::FeedforwardLQRController, eps::Float64)
    LOCAL_AGENTS[agent_id] = FeedforwardAgentState(x0, dyn, dt, ctrl, eps)
    return nothing
end

"""
    _step_agent_on_worker!(agent_id::Int, qstar_i::Vector{Float64}, dt::Float64)

Step standard feedback agent on worker.
"""
function _step_agent_on_worker!(agent_id::Int, qstar_i::Vector{Float64}, dt::Float64)
    w = LOCAL_AGENTS[agent_id]
    return step_agent!(w, qstar_i, dt)
end

"""
    _step_agent_on_worker!(agent_id::Int, qstar_i::Vector{Float64}, qstar_dot_i::Vector{Float64}, dt::Float64)

Step feedforward agent on worker.
"""
function _step_agent_on_worker!(agent_id::Int, qstar_i::Vector{Float64}, qstar_dot_i::Vector{Float64}, dt::Float64)
    w = LOCAL_AGENTS[agent_id]
    return step_agent!(w, qstar_i, qstar_dot_i, dt)
end

"""
    init_distributed_agents!(pids, agent_configs, dt::Float64, eps::Float64)

Deploys agent states to the specified worker processes. 
Supports both standard feedback `(x0, dyn, K_lqr)` and feedforward `(x0, dyn, FeedforwardLQRController)` configurations.
"""
function init_distributed_agents!(pids::Vector{Int}, agent_configs::Vector, dt::Float64, eps::Float64)
    np = length(pids)
    for (i, cfg) in enumerate(agent_configs)
        pid = pids[(i - 1) % np + 1]
        if cfg[3] isa FeedforwardLQRController
            x0, dyn, ctrl = cfg
            remotecall_fetch(_init_feedforward_agent_on_worker!, pid, i, x0, dyn, dt, ctrl, eps)
        else
            x0, dyn, K_lqr = cfg
            remotecall_fetch(_init_agent_on_worker!, pid, i, x0, dyn, dt, K_lqr, eps)
        end
    end
end

"""
    run_layered_simulation(prob::LayeredControlProblem, pids::Vector{Int}; mode=:distributed, nx::Int=10)

Runs the layered control simulation defined by `prob`.
If `prob.target_velocity_func` is supplied, performs a secondary solve for reference velocity `qstar_dot`
and passes it to feedforward agent state controllers.
"""
function run_layered_simulation(prob::LayeredControlProblem, pids::Vector{Int}; mode=:distributed, nx::Int=10)
    sheaf = prob.sheaf
    target_nodes = prob.target_nodes
    target_trajectory_func = prob.target_trajectory_func
    target_velocity_func = prob.target_velocity_func
    dt = prob.dt
    steps = prob.steps
    
    NA = length(sheaf.vertex_stalks) - length(target_nodes)
    D = sheaf.vertex_stalks[1]
    
    sim_data = zeros(steps, NA, nx)
    qstar_history = zeros(steps, NA, D)

    # Initial target positions
    boundary0 = Dict(node => target_trajectory_func(node, 0.0) for node in target_nodes)
    _, _, Hraw, LIBraw = _harmonic_extension_restricted_laplacian(sheaf, boundary0)
    H = Matrix(Hraw)
    LIB = Matrix(LIBraw)
    
    H_pinv = pinv(H)
    
    # Precompute chordal factorisation
    F = cholesky!(ChordalCholesky(sparse(H)), NoPivot())
    Lfac = F.L
    
    # Partition tree for distributed solve
    nchunk = min(NA, length(pids))
    partition = partition_tree(Lfac, nchunk)
    solve_pids = [pids[(i - 1) % length(pids) + 1] for i in 1:length(partition.chunks)]

    np = length(pids)
    has_velocity = target_velocity_func !== nothing
    
    for t_idx in 1:steps
        t = t_idx * dt
        
        # 1. Primary harmonic extension solve: H * qstar = -LIB * b_t
        b_t = vcat([target_trajectory_func(node, t) for node in target_nodes]...)
        rhs = Vector(-LIB * b_t)
        
        if mode == :centralised
            qstar_full = H_pinv * rhs
        else
            rhs_p = Vector(F.P' \ rhs)
            y_sol = distributed_tree_solve(Lfac, rhs_p, length(partition.chunks); pids = solve_pids)
            qstar_full = F.P \ y_sol
        end
        
        qstar = [qstar_full[(i-1)*D+1:i*D] for i in 1:NA]
        for i in 1:NA
            qstar_history[t_idx, i, :] = qstar[i]
        end
        
        # 2. Secondary solve for harmonic extension velocity if velocity function provided:
        # H * qstar_dot = -LIB * b_dot_t
        qstar_dot = Vector{Vector{Float64}}(undef, NA)
        if has_velocity
            b_dot_t = vcat([target_velocity_func(node, t) for node in target_nodes]...)
            rhs_dot = Vector(-LIB * b_dot_t)
            
            if mode == :centralised
                qstar_dot_full = H_pinv * rhs_dot
            else
                rhs_dot_p = Vector(F.P' \ rhs_dot)
                y_dot_sol = distributed_tree_solve(Lfac, rhs_dot_p, length(partition.chunks); pids = solve_pids)
                qstar_dot_full = F.P \ y_dot_sol
            end
            for i in 1:NA
                qstar_dot[i] = qstar_dot_full[(i-1)*D+1:i*D]
            end
        end
        
        # 3. Step agents on workers
        step_futures = []
        for i in 1:NA
            pid = pids[(i - 1) % np + 1]
            if has_velocity
                push!(step_futures, remotecall(_step_agent_on_worker!, pid, i, qstar[i][1:prob.pos_dim], qstar_dot[i][1:prob.pos_dim], dt))
            else
                push!(step_futures, remotecall(_step_agent_on_worker!, pid, i, qstar[i][1:prob.pos_dim], dt))
            end
        end
        
        step_results = fetch.(step_futures)
        for i in 1:NA
            x_act, x_ref = step_results[i]
            sim_data[t_idx, i, :] = x_act
        end
    end
    
    return LayeredSimulationResult(prob, sim_data, qstar_history)
end

end # module
