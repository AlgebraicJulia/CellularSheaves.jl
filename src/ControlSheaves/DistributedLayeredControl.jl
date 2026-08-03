module DistributedLayeredControl

using Distributed
using LinearAlgebra
using SparseArrays
using CliqueTrees.Multifrontal
using ...NetworkSheaves: EuclideanSheaf
using ...NetworkSheaves.EuclideanSheaves: _harmonic_extension_restricted_laplacian
using ...NetworkSheaves.DistributedSolve: partition_tree, distributed_tree_solve
using ..Tikhonov: tikhonov_step!
using ..AgentControllers: AgentState, AbstractAgentDynamics, QuadrotorDynamics, step_agent!

export init_distributed_agents!, run_layered_simulation, LayeredControlProblem, LayeredSimulationResult, animate_layered_escort

function animate_layered_escort end

# Worker-local state storage. 
# We use a Dict to allow a worker to potentially simulate multiple agents if NA > nworkers().
const LOCAL_AGENTS = Dict{Int, AgentState}()

struct LayeredControlProblem
    sheaf::EuclideanSheaf
    target_nodes::Vector{Int}
    target_trajectory_func::Function
    agent_configs::Vector{Tuple{Vector{Float64}, AbstractAgentDynamics, Matrix{Float64}}}
    dt::Float64
    steps::Int
    r_ring::Float64
end

function LayeredControlProblem(sheaf, target_nodes, target_trajectory_func, agent_configs::Vector{<:Tuple}, dt, steps, r_ring)
    configs = Tuple{Vector{Float64}, AbstractAgentDynamics, Matrix{Float64}}[c for c in agent_configs]
    return LayeredControlProblem(sheaf, target_nodes, target_trajectory_func, configs, dt, steps, r_ring)
end

# Keyword constructor to allow optional r_ring and arbitrary dynamics
function LayeredControlProblem(; sheaf, target_nodes, target_trajectory_func, agent_configs, dt, steps, r_ring=0.0)
    return LayeredControlProblem(sheaf, target_nodes, target_trajectory_func, agent_configs, dt, steps, r_ring)
end

struct LayeredSimulationResult
    problem::LayeredControlProblem
    sim_data::Array{Float64, 3}
    qstar_history::Array{Float64, 3}
end

"""
    _init_agent_on_worker!(agent_id::Int, x0::Vector{Float64}, dyn::QuadrotorDynamics, dt::Float64, K_lqr::Matrix{Float64}, eps::Float64)

Internal function called on the worker process to initialize the agent state.
"""
function _init_agent_on_worker!(agent_id::Int, x0::Vector{Float64}, dyn::AbstractAgentDynamics, dt::Float64, K_lqr::Matrix{Float64}, eps::Float64)
    LOCAL_AGENTS[agent_id] = AgentState(x0, dyn, dt, K_lqr, eps)
    return nothing
end

"""
    _step_agent_on_worker!(agent_id::Int, qstar_4d_i::Vector{Float64}, dt::Float64)

Internal function called on the worker process to step the agent dynamics.
"""
function _step_agent_on_worker!(agent_id::Int, qstar_4d_i::Vector{Float64}, dt::Float64)
    w = LOCAL_AGENTS[agent_id]
    return step_agent!(w, qstar_4d_i, dt)
end


"""
    init_distributed_agents!(pids, agent_configs, dt::Float64, eps::Float64)

Deploys agent states to the specified worker processes. 
`agent_configs` is a Vector of Tuples `(x0, dyn, K_lqr)`.
If `length(pids)` < `length(agent_configs)`, agents are distributed across available pids.
"""
function init_distributed_agents!(pids::Vector{Int}, agent_configs::Vector{<:Tuple{Vector{Float64}, <:AbstractAgentDynamics, Matrix{Float64}}}, dt::Float64, eps::Float64)
    np = length(pids)
    for (i, (x0, dyn, K_lqr)) in enumerate(agent_configs)
        pid = pids[(i - 1) % np + 1]
        remotecall_fetch(_init_agent_on_worker!, pid, i, x0, dyn, dt, K_lqr, eps)
    end
end

"""
    run_layered_simulation(prob::LayeredControlProblem, pids::Vector{Int}; mode=:distributed, nx::Int=10)

Runs the layered control simulation defined by `prob`.
Returns `LayeredSimulationResult`.
"""
function run_layered_simulation(prob::LayeredControlProblem, pids::Vector{Int}; mode=:distributed, nx::Int=10)
    sheaf = prob.sheaf
    target_nodes = prob.target_nodes
    target_trajectory_func = prob.target_trajectory_func
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
    # Use min(NA, length(pids)) to avoid over-partitioning if we have fewer workers
    nchunk = min(NA, length(pids))
    partition = partition_tree(Lfac, nchunk)
    solve_pids = [pids[(i - 1) % length(pids) + 1] for i in 1:length(partition.chunks)]

    np = length(pids)
    
    for t_idx in 1:steps
        t = t_idx * dt
        
        # Build boundary vector
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
        
        # Step agents
        step_futures = []
        for i in 1:NA
            pid = pids[(i - 1) % np + 1]
            push!(step_futures, remotecall(_step_agent_on_worker!, pid, i, qstar[i][1:D-1], dt))
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
