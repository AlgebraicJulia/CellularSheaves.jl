"""
    Layered

Generic Team-of-Teams Layered Escort Sheaves, Pushforward Solvers, and Hierarchical Reference Generators.
"""
module Layered

export RingSpec, SupportSpec, LayeredEscortSpec, LayeredFiberBases,
       HomogeneousDynamics, TeamHomogeneousDynamics, IndividualizedDynamics,
       LayeredEscortProblem, LayeredEscortResult,
       build_layered_escort_sheaf, build_layered_homomorphism, build_layered_fiber_bases,
       world_to_pf_stalk, solve_high_level_harmonic, solve_mid_level_harmonic,
       solve_direct_harmonic, run_layered_escort_simulation, compute_formation_radius_error,
       animate_comprehensive_escort

using LinearAlgebra
using Graphs
using ArgCheck
using ...NetworkSheaves: vertex_stalks, get_restriction_map, EuclideanSheaf
using ...NetworkSheaves.EuclideanSheaves: add_sheaf_edge!, _harmonic_extension_restricted_laplacian
using ...NetworkSheaves.Formations: build_escort_ring, se3_translation_matrix
using ...NetworkSheaves.GraphHomomorphisms: GraphHomomorphism, fiber_vertices
using ...NetworkSheaves.Pushforwards: pushforward_sheaf, all_fiber_bases
using ..AgentControllers: AbstractAgentDynamics, QuadrotorDynamics, AgentState, step_agent!

function animate_comprehensive_escort end

# ---------------------------------------------------------------------------
# Topology Specification Types
# ---------------------------------------------------------------------------

"""
    RingSpec(target_id::Int, n_agents::Int, radius::Float64)

Specification for an escort ring around a target vertex.
"""
struct RingSpec
    target_id::Int
    n_agents::Int
    radius::Float64
end

"""
    SupportSpec(u_ring::Int, v_ring::Int, n_agents::Int)

Specification for a support pool of agents bridging Escort Ring `u_ring` and Escort Ring `v_ring`.
"""
struct SupportSpec
    u_ring::Int
    v_ring::Int
    n_agents::Int
end

"""
    LayeredEscortSpec(rings::Vector{RingSpec}, supports::Vector{SupportSpec})

Full multi-team topology specification for arbitrary escort rings and support edge pools.
"""
struct LayeredEscortSpec
    rings::Vector{RingSpec}
    supports::Vector{SupportSpec}
    n_rings::Int
    n_supports::Int
    n_targets::Int
    n_agents::Int
    total_nodes::Int
    ring_node_ranges::Vector{UnitRange{Int}}
    support_node_ranges::Vector{UnitRange{Int}}
    target_nodes::Vector{Int}
end

function LayeredEscortSpec(rings::Vector{RingSpec}, supports::Vector{SupportSpec})
    n_rings = length(rings)
    n_supports = length(supports)
    n_targets = n_rings

    ring_ranges = UnitRange{Int}[]
    curr_node = 1
    for r in rings
        push!(ring_ranges, curr_node:(curr_node + r.n_agents - 1))
        curr_node += r.n_agents
    end

    support_ranges = UnitRange{Int}[]
    for s in supports
        push!(support_ranges, curr_node:(curr_node + s.n_agents - 1))
        curr_node += s.n_agents
    end

    n_agents = curr_node - 1
    target_nodes = [curr_node + i - 1 for i in 1:n_targets]
    total_nodes = n_agents + n_targets

    return LayeredEscortSpec(rings, supports, n_rings, n_supports, n_targets, n_agents, total_nodes, ring_ranges, support_ranges, target_nodes)
end

# ---------------------------------------------------------------------------
# Structured Fiber Basis Management
# ---------------------------------------------------------------------------

"""
    LayeredFiberBases

Structured data type managing fiber basis matrices and target sub-basis inverses across pushforward sheaf `PfF`.
"""
struct LayeredFiberBases
    hom::GraphHomomorphism
    fiber_bases::Vector{Matrix{Float64}}
    target_subbases::Dict{Int, Matrix{Float64}} # maps ring_idx => target sub-basis matrix
end

"""
    build_layered_fiber_bases(hom::GraphHomomorphism, F::EuclideanSheaf, spec::LayeredEscortSpec) -> LayeredFiberBases

Construct structured fiber bases and extract target sub-basis matrices for world target mapping.
"""
function build_layered_fiber_bases(hom::GraphHomomorphism, F::EuclideanSheaf, spec::LayeredEscortSpec)
    bases = all_fiber_bases(hom, F)
    target_subbases = Dict{Int, Matrix{Float64}}()

    for r_idx in 1:spec.n_rings
        tv = r_idx
        fverts = fiber_vertices(hom, tv)
        t_node = spec.target_nodes[r_idx]
        
        # Row offset of target node within the concatenated fiber stalk space
        row_offset = 0
        for v in fverts
            if v == t_node
                break
            end
            row_offset += vertex_stalks(F)[v]
        end
        
        d_target = vertex_stalks(F)[t_node]
        target_subbases[r_idx] = bases[tv][row_offset+1 : row_offset+d_target, :]
    end

    return LayeredFiberBases(hom, bases, target_subbases)
end

"""
    world_to_pf_stalk(bases::LayeredFiberBases, ring_idx::Int, p_world::Vector{Float64}) -> Vector{Float64}

Convert a 4D world target vector `p_world` into the `PfF` stalk basis coordinates for ring `ring_idx`.
"""
function world_to_pf_stalk(bases::LayeredFiberBases, ring_idx::Int, p_world::Vector{Float64})
    B_target = bases.target_subbases[ring_idx]
    return B_target \ p_world
end

# ---------------------------------------------------------------------------
# Dynamics Specification Types (Memory-Efficient Shared Dynamics)
# ---------------------------------------------------------------------------

abstract type AbstractDynamicsSpec end

"""
    HomogeneousDynamics(dyn::AbstractAgentDynamics)

Single dynamics model shared by all agents in the system (allocates 1 dynamics object).
"""
struct HomogeneousDynamics <: AbstractDynamicsSpec
    dyn::AbstractAgentDynamics
end

"""
    TeamHomogeneousDynamics(team_dyns::Dict{Int, AbstractAgentDynamics})

Dynamics models shared per team/ring.
"""
struct TeamHomogeneousDynamics <: AbstractDynamicsSpec
    team_dyns::Dict{Int, AbstractAgentDynamics}
end

"""
    IndividualizedDynamics(dyns::Vector{<:AbstractAgentDynamics})

Distinct dynamics models for each individual agent.
"""
struct IndividualizedDynamics <: AbstractDynamicsSpec
    dyns::Vector{AbstractAgentDynamics}
end

# ---------------------------------------------------------------------------
# Sheaf Topology & Homomorphism Builders
# ---------------------------------------------------------------------------

"""
    build_layered_escort_sheaf(spec::LayeredEscortSpec) -> EuclideanSheaf

Construct the full cellular sheaf `F` for an arbitrary layered escort specification.
"""
function build_layered_escort_sheaf(spec::LayeredEscortSpec)
    F = EuclideanSheaf{Float64}(fill(4, spec.total_nodes))

    # 1. Escort Rings
    for (r_idx, r) in enumerate(spec.rings)
        r_range = spec.ring_node_ranges[r_idx]
        t_node = spec.target_nodes[r_idx]
        target_local_idx = r.n_agents + 1
        
        ring = build_escort_ring(r.n_agents, target_local_idx, r.radius; observers=[1])
        for e in edges(ring.underlying_graph)
            u_loc, v_loc = src(e), dst(e)
            u = u_loc == target_local_idx ? t_node : r_range[u_loc]
            v = v_loc == target_local_idx ? t_node : r_range[v_loc]
            add_sheaf_edge!(F, u, v, get_restriction_map(ring, u_loc, v_loc), get_restriction_map(ring, v_loc, u_loc))
        end
    end

    # 2. Support Pools
    for (s_idx, s) in enumerate(spec.supports)
        s_range = spec.support_node_ranges[s_idx]
        
        # Support consensus ring / chain
        for i in 1:s.n_agents
            u = s_range[i]
            v = s_range[i % s.n_agents + 1]
            add_sheaf_edge!(F, u, v, Matrix{Float64}(I, 4, 4), Matrix{Float64}(I, 4, 4))
        end

        # Direct tracking edges to target u and target v
        if s.n_agents >= 1
            add_sheaf_edge!(F, s_range[1], spec.target_nodes[s.u_ring], Matrix{Float64}(I, 4, 4), Matrix{Float64}(I, 4, 4))
        end
        if s.n_agents >= 2
            add_sheaf_edge!(F, s_range[2], spec.target_nodes[s.v_ring], Matrix{Float64}(I, 4, 4), Matrix{Float64}(I, 4, 4))
        end
    end

    return F
end

"""
    build_layered_homomorphism(spec::LayeredEscortSpec) -> GraphHomomorphism

Construct graph homomorphism `f : G -> H` mapping ring nodes & target to coarse vertex `v_H`.
"""
function build_layered_homomorphism(spec::LayeredEscortSpec)
    v_map = zeros(Int, spec.total_nodes)

    for (r_idx, r_range) in enumerate(spec.ring_node_ranges)
        v_map[r_range] .= r_idx
        v_map[spec.target_nodes[r_idx]] = r_idx
    end

    for (s_idx, s_range) in enumerate(spec.support_node_ranges)
        v_map[s_range] .= spec.n_rings + s_idx
    end

    n_target_H = spec.n_rings + spec.n_supports
    return GraphHomomorphism(v_map, n_target_H)
end

# ---------------------------------------------------------------------------
# Hierarchical & Direct Solvers
# ---------------------------------------------------------------------------

"""
    solve_high_level_harmonic(PfF, bases::LayeredFiberBases, target_positions::Vector{Vector{Float64}})

Solve harmonic extension on the coarse pushforward sheaf `PfF` given world target positions.
"""
function solve_high_level_harmonic(PfF, bases::LayeredFiberBases, target_positions::Vector{Vector{Float64}})
    boundary = Dict{Int, Vector{Float64}}()
    for (r_idx, p) in enumerate(target_positions)
        boundary[r_idx] = world_to_pf_stalk(bases, r_idx, p)
    end

    _, _, H_mat, B_mat = _harmonic_extension_restricted_laplacian(PfF, boundary)
    
    # Assembly RHS for boundary nodes
    p_boundary = vcat([boundary[r_idx] for r_idx in 1:length(target_positions)]...)
    q_interior = vec(H_mat \ (-Matrix(B_mat) * p_boundary))

    q_H = Vector{Vector{Float64}}(undef, PfF.underlying_graph.ne > 0 ? nv(PfF.underlying_graph) : length(boundary) + length(q_interior)/4)
    for (r_idx, p) in enumerate(target_positions)
        q_H[r_idx] = boundary[r_idx]
    end
    
    int_idx = 1
    for v in 1:length(q_H)
        if !haskey(boundary, v)
            q_H[v] = q_interior[(int_idx-1)*4 + 1 : int_idx*4]
            int_idx += 1
        end
    end

    return q_H
end

"""
    solve_mid_level_harmonic(q_H::Vector{Vector{Float64}}, bases::LayeredFiberBases) -> Matrix{Float64}

Lift coarse pushforward section `q_H` to 4D agent reference states on `G`.
"""
function solve_mid_level_harmonic(q_H::Vector{Vector{Float64}}, bases::LayeredFiberBases)
    n_pf_nodes = length(q_H)
    q_G_components = [bases.fiber_bases[v] * q_H[v] for v in 1:n_pf_nodes]

    # Find total agents from fiber bases
    n_agents = sum(size(bases.fiber_bases[v], 1) ÷ 4 for v in 1:n_pf_nodes) - length(bases.target_subbases)
    q_agents = zeros(n_agents, 4)

    curr_agent = 1
    for v in 1:n_pf_nodes
        B = bases.fiber_bases[v]
        q_fib = q_G_components[v]
        n_fib_nodes = size(B, 1) ÷ 4
        
        # If this fiber has a target, the last node is the target node
        has_target = haskey(bases.target_subbases, v)
        n_fib_agents = has_target ? n_fib_nodes - 1 : n_fib_nodes

        for i in 1:n_fib_agents
            q_agents[curr_agent, :] .= q_fib[(i-1)*4 + 1 : i*4]
            curr_agent += 1
        end
    end

    return q_agents
end

"""
    solve_direct_harmonic(F::EuclideanSheaf, target_nodes::Vector{Int}, target_positions::Vector{Vector{Float64}}) -> Matrix{Float64}

Solve harmonic extension directly on `F` given pinned target nodes and positions.
"""
function solve_direct_harmonic(F::EuclideanSheaf, target_nodes::Vector{Int}, target_positions::Vector{Vector{Float64}})
    boundary = Dict(target_nodes[i] => target_positions[i] for i in 1:length(target_nodes))
    _, _, H_mat, B_mat = _harmonic_extension_restricted_laplacian(F, boundary)
    
    p_boundary = vcat(target_positions...)
    q_interior = H_mat \ (-Matrix(B_mat) * p_boundary)

    n_total = length(F.vertex_stalks)
    n_targets = length(target_nodes)
    n_agents = n_total - n_targets

    q_full = zeros(n_total, 4)
    
    int_idx = 1
    for v in 1:n_total
        if v in target_nodes
            t_idx = findfirst(==(v), target_nodes)
            q_full[v, :] .= target_positions[t_idx]
        else
            q_full[v, :] .= q_interior[(int_idx-1)*4 + 1 : int_idx*4]
            int_idx += 1
        end
    end

    return q_full[1:n_agents, :]
end

# ---------------------------------------------------------------------------
# Problem & Simulation Setup
# ---------------------------------------------------------------------------

struct LayeredEscortProblem
    spec::LayeredEscortSpec
    sheaf::EuclideanSheaf{Float64}
    hom::GraphHomomorphism
    pf_sheaf::EuclideanSheaf{Float64}
    bases::LayeredFiberBases
    dynamics_spec::AbstractDynamicsSpec
    target_trajectories::Vector{Any} # functions t -> [x,y,z,1.0]
    target_velocities::Vector{Any}   # functions t -> [vx,vy,vz,0.0]
    target_accelerations::Vector{Any}# functions t -> [ax,ay,az,0.0]
    dt::Float64
    steps::Int
end

struct LayeredEscortResult
    problem::LayeredEscortProblem
    sim_data::Vector{Vector{Vector{Float64}}}
    qstar_history::Vector{Matrix{Float64}}
    target_history::Vector{Vector{Vector{Float64}}}
end

"""
    run_layered_escort_simulation(prob::LayeredEscortProblem; use_feedforward::Bool=false) -> LayeredEscortResult

Run complete multi-agent layered escort simulation.
"""
function run_layered_escort_simulation(prob::LayeredEscortProblem; use_feedforward::Bool=false)
    spec = prob.spec
    bases = prob.bases
    STEPS = prob.steps
    DT = prob.dt
    time_grid = 0:DT:(STEPS*DT)

    # Initialize agent states
    init_states = [[0.0, 0.0, 1.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0] for _ in 1:spec.n_agents]
    current_states = copy(init_states)

    sim_data = Vector{Vector{Vector{Float64}}}()
    qstar_history = Vector{Matrix{Float64}}()
    target_history = Vector{Vector{Vector{Float64}}}()

    for step in 1:STEPS
        t = time_grid[step]
        
        # 1. Target positions
        p_targets = [traj(t) for traj in prob.target_trajectories]
        
        # 2. High-level harmonic extension
        q_H = solve_high_level_harmonic(prob.pf_sheaf, bases, p_targets)
        
        # Extract target world positions for visualization
        target_world_centers = Vector{Vector{Float64}}()
        for r_idx in 1:spec.n_rings
            push!(target_world_centers, p_targets[r_idx])
        end
        for s_idx in 1:spec.n_supports
            pf_v = spec.n_rings + s_idx
            center = bases.fiber_bases[pf_v][1:4, :] * q_H[pf_v]
            push!(target_world_centers, center)
        end
        push!(target_history, target_world_centers)

        # 3. Mid-level reference lifting q*
        qstar_agents = solve_mid_level_harmonic(q_H, bases)
        push!(qstar_history, qstar_agents)

        # 4. Integrate Agent Dynamics
        step_states = Vector{Vector{Float64}}()
        for i in 1:spec.n_agents
            x_actual = current_states[i][1:3]
            x_ref = qstar_agents[i, 1:3]
            
            # Simple dynamics integration step (or quadrotor dynamics)
            error = x_actual - x_ref
            current_states[i][1:3] .-= 0.3 * error * DT
            push!(step_states, copy(current_states[i]))
        end
        push!(sim_data, step_states)
    end

    return LayeredEscortResult(prob, sim_data, qstar_history, target_history)
end

"""
    compute_formation_radius_error(agent_positions::Matrix{Float64}, centroid::Vector{Float64}, r_expected::Float64) -> Float64

Compute mean absolute deviation from expected radius `r_expected`.
"""
function compute_formation_radius_error(agent_positions::Matrix{Float64}, centroid::Vector{Float64}, r_expected::Float64)
    N = size(agent_positions, 1)
    if N == 0
        return 0.0
    end
    err_sum = sum(abs(norm(agent_positions[i, 1:2] - centroid[1:2]) - r_expected) for i in 1:N)
    return err_sum / N
end

end # module Layered
