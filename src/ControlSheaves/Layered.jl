"""
    Layered

Generic Team-of-Teams Layered Escort Sheaves, Pushforward Solvers, and Hierarchical Reference Generators.
"""
module Layered

export RingSpec, SupportSpec, LayeredEscortSpec, LayeredFiberBases,
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
using ..AgentControllers: AbstractAgentDynamics, QuadrotorDynamics, AgentState, step_agent!,
                           position_indices, initial_state, _default_initial_position
using ..NestedSystems: SystemBinding, AgentBinding, ResolvedAgent, _fold, _fold_system
import ..NestedSystems: resolve_dynamics

function animate_comprehensive_escort end

# ---------------------------------------------------------------------------
# Topology Specification Types
# ---------------------------------------------------------------------------

"""
    RingSpec(target_id::Int, n_agents::Int, radius::Float64; observers::Vector{Int}=[1])

Specification for an escort ring around a target vertex. `observers` names which ring
agents (by local index `1:n_agents`) are directly pinned to the target; any subset of
ring agents may observe the target, defaulting to just the first agent.
"""
struct RingSpec
    target_id::Int
    n_agents::Int
    radius::Float64
    observers::Vector{Int}
end

RingSpec(target_id::Int, n_agents::Int, radius::Float64; observers=[1]) =
    RingSpec(target_id, n_agents, radius, collect(observers))

"""
    SupportSpec(src_ring::Int, tgt_ring::Int, n_agents::Int; src_observers::Vector{Int}=[1], tgt_observers::Vector{Int}=[1])

Specification for a support pool of agents bridging Escort Ring `src_ring` and Escort Ring `tgt_ring`.
`src_observers`/`tgt_observers` name which support-pool agents (by local index `1:n_agents`)
are directly pinned to the `src_ring`/`tgt_ring` target respectively. Defaults to the first
agent in the pool observing both targets.
"""
struct SupportSpec
    src_ring::Int
    tgt_ring::Int
    n_agents::Int
    src_observers::Vector{Int}
    tgt_observers::Vector{Int}
end

SupportSpec(src_ring::Int, tgt_ring::Int, n_agents::Int; src_observers=[1], tgt_observers=[1]) =
    SupportSpec(src_ring, tgt_ring, n_agents, collect(src_observers), collect(tgt_observers))

"""
    LayeredEscortSpec(rings::Vector{RingSpec}, supports::Vector{SupportSpec}; D::Int=4, affine::Bool=true)

Full multi-team topology specification for arbitrary escort rings and support edge pools.

`D` is the stalk dimension (default 4). When `affine=true` (default), stalks carry `D-1`
homogeneous-affine translation coordinates plus one homogeneous row (e.g. `D=4` for SE(3)
3D translation, `D=3` for planar 2D translation). When `affine=false`, stalks are `D` plain
Euclidean coordinates with no translation offsets representable (ring/support radii must be 0).
"""
struct LayeredEscortSpec
    rings::Vector{RingSpec}
    supports::Vector{SupportSpec}
    D::Int # stalk dimension (default 4 for SE(3) homogeneous coordinates)
    affine::Bool # whether stalks carry a homogeneous affine coordinate
    n_rings::Int
    n_supports::Int
    n_targets::Int
    n_agents::Int
    total_nodes::Int
    ring_node_ranges::Vector{UnitRange{Int}}
    support_node_ranges::Vector{UnitRange{Int}}
    target_nodes::Vector{Int}
end

function LayeredEscortSpec(rings::Vector{RingSpec}, supports::Vector{SupportSpec}; D::Int=4, affine::Bool=true)
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

    return LayeredEscortSpec(rings, supports, D, affine, n_rings, n_supports, n_targets, n_agents, total_nodes, ring_ranges, support_ranges, target_nodes)
end

# ---------------------------------------------------------------------------
# Dynamics binding: LayeredEscortSpec as a depth-1 tree (Issue 013)
# ---------------------------------------------------------------------------

"""
    _flat_children(spec::LayeredEscortSpec) -> Vector{Tuple{Symbol,UnitRange{Int},Int}}

Every ring and support pod of `spec`, named `:ring1, :ring2, …` / `:support1, :support2, …`
matching the ordering of `spec.rings` / `spec.supports`, together with its global agent-index
range and its own agent count. Order matches `spec.ring_node_ranges` followed by
`spec.support_node_ranges` — the same order [`run_layered_escort_simulation`](@ref) iterates.
"""
function _flat_children(spec::LayeredEscortSpec)
    out = Tuple{Symbol,UnitRange{Int},Int}[]
    for (i, r) in enumerate(spec.ring_node_ranges)
        push!(out, (Symbol(:ring, i), r, length(r)))
    end
    for (i, r) in enumerate(spec.support_node_ranges)
        push!(out, (Symbol(:support, i), r, length(r)))
    end
    return out
end

"""
    _resolve_flat_bindings(spec, ctx) -> Vector{Tuple{Symbol,Int,Int,AgentBinding}}

Walk `spec`'s implicit depth-1 tree folding `ctx`'s cascade, exactly like
`resolve_dynamics(::LayeredEscortSpec, ::SystemBinding)`, but stop *before* applying any
initial-position default — callers get the raw folded `AgentBinding` (dynamics guaranteed
non-`nothing`, `K_lqr` still possibly `nothing`, `initial_position` still possibly `nothing`)
together with `(child_name, local_idx, global_idx)`.

This is the shared core both `resolve_dynamics` (which defaults straight to the airstrip
position) and `LayeredEscortProblem`'s three-tier fallback (which checks `initial_positions`
*before* the airstrip) build on, so the precedence for `dynamics`/`K_lqr` cannot drift between
the two call sites.
"""
function _resolve_flat_bindings(spec::LayeredEscortSpec, ctx::SystemBinding)
    !isempty(ctx.agents) && throw(ArgumentError(
        "SystemBinding declares `agents` at the root, but a LayeredEscortSpec's root has no " *
        "direct agents — bind them via a :ring<i>/:support<i> child instead"))

    children = _flat_children(spec)
    valid_names = [name for (name, _, _) in children]
    for name in keys(ctx.children)
        name in valid_names || throw(ArgumentError(
            "SystemBinding names child :$name, but the spec has no such ring/support " *
            "(children: $(join(valid_names, ", ")))"))
    end

    root_level = _fold_system(AgentBinding(), ctx)
    out = Tuple{Symbol,Int,Int,AgentBinding}[]
    for (name, range, n_local) in children
        child_ctx = get(ctx.children, name, SystemBinding())
        !isempty(child_ctx.children) && throw(ArgumentError(
            "SystemBinding declares children at :$name, but :$name is a flat ring/support " *
            "pool with no children of its own — use `agents` to bind individual agents by " *
            "local index"))
        for idx in keys(child_ctx.agents)
            @argcheck 1 <= idx <= n_local "SystemBinding names agent $idx at :$name, out of range 1:$n_local"
        end

        team_level = _fold_system(root_level, child_ctx)
        for (local_idx, global_idx) in enumerate(range)
            b = _fold(team_level, get(child_ctx.agents, local_idx, AgentBinding()))
            b.dynamics === nothing && throw(ArgumentError(
                "no dynamics bound for agent $local_idx at :$name (global index $global_idx) " *
                "— declare a `dynamics` at this agent, this ring/pod, or the root"))
            push!(out, (name, local_idx, global_idx, b))
        end
    end
    return out
end

"""
    resolve_dynamics(spec::LayeredEscortSpec, ctx::SystemBinding) -> Vector{ResolvedAgent}

Resolve dynamics bindings for a flat layered escort specification, treating it as a depth-1
tree: each escort ring and each support pod is a child of a single implicit root.

Children are addressed by the conventional names `:ring1, :ring2, …` and `:support1, :support2,
…`, matching the ordering of `spec.rings` and `spec.supports`. Per-agent overrides use the
agent's **local** index within its ring or pod, consistent with `RingSpec.observers`.

Agents come back in global agent-index order, matching `spec.ring_node_ranges` followed by
`spec.support_node_ranges` — the same ordering [`run_layered_escort_simulation`](@ref) iterates.
Precedence and field independence are identical to the tree-shaped
`NestedSystems.resolve_dynamics(::NestedSystemSpec, ::SystemBinding)`: this reuses that same
`_fold`/`_fold_system` cascade rather than re-deriving the precedence rule.
"""
function resolve_dynamics(spec::LayeredEscortSpec, ctx::SystemBinding)
    resolved = ResolvedAgent[]
    for (name, local_idx, global_idx, b) in _resolve_flat_bindings(spec, ctx)
        K_lqr = something(b.K_lqr, Some(zeros(0, 0)))
        initial_position = something(b.initial_position,
                                     Some(_default_initial_position(b.dynamics, global_idx)))
        push!(resolved, ResolvedAgent(global_idx, [name, Symbol(local_idx)],
                                      b.dynamics, K_lqr, initial_position))
    end
    return resolved
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
    target_subbases_inv::Dict{Int, Matrix{Float64}} # pre-computed pseudoinverses for O(1) projection
end

"""
    build_layered_fiber_bases(hom::GraphHomomorphism, F::EuclideanSheaf, spec::LayeredEscortSpec) -> LayeredFiberBases

Construct structured fiber bases and extract target sub-basis matrices for world target mapping.
"""
function build_layered_fiber_bases(hom::GraphHomomorphism, F::EuclideanSheaf, spec::LayeredEscortSpec)
    bases = all_fiber_bases(hom, F)
    target_subbases = Dict{Int, Matrix{Float64}}()
    target_subbases_inv = Dict{Int, Matrix{Float64}}()

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
        B_sub = bases[tv][row_offset+1 : row_offset+d_target, :]
        target_subbases[r_idx] = B_sub
        target_subbases_inv[r_idx] = pinv(B_sub)
    end

    return LayeredFiberBases(hom, bases, target_subbases, target_subbases_inv)
end

"""
    world_to_pf_stalk(bases::LayeredFiberBases, ring_idx::Int, p_world::Vector{Float64}) -> Vector{Float64}

Convert a world target vector `p_world` into the `PfF` stalk basis coordinates for ring `ring_idx`.
"""
function world_to_pf_stalk(bases::LayeredFiberBases, ring_idx::Int, p_world::Vector{Float64})
    return bases.target_subbases_inv[ring_idx] * p_world
end

# ---------------------------------------------------------------------------
# Sheaf Topology & Homomorphism Builders
# ---------------------------------------------------------------------------

"""
    build_layered_escort_sheaf(spec::LayeredEscortSpec) -> EuclideanSheaf

Construct the full cellular sheaf `F` for an arbitrary layered escort specification.
"""
function build_layered_escort_sheaf(spec::LayeredEscortSpec)
    D = spec.D
    F = EuclideanSheaf{Float64}(fill(D, spec.total_nodes))

    # 1. Escort Rings
    for (r_idx, r) in enumerate(spec.rings)
        r_range = spec.ring_node_ranges[r_idx]
        t_node = spec.target_nodes[r_idx]
        target_local_idx = r.n_agents + 1
        
        ring = build_escort_ring(r.n_agents, target_local_idx, r.radius; observers=r.observers, D=D, affine=spec.affine)
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
        I_D = Matrix{Float64}(I, D, D)

        # Support consensus ring / chain
        if s.n_agents > 1
            for i in 1:s.n_agents
                u = s_range[i]
                v = s_range[i % s.n_agents + 1]
                add_sheaf_edge!(F, u, v, I_D, I_D)
            end
        end

        # Direct tracking edges to the src_ring and tgt_ring targets
        @argcheck all(1 .<= o .<= s.n_agents for o in s.src_observers) "src_observers must be within 1:n_agents"
        @argcheck all(1 .<= o .<= s.n_agents for o in s.tgt_observers) "tgt_observers must be within 1:n_agents"
        for i in s.src_observers
            add_sheaf_edge!(F, s_range[i], spec.target_nodes[s.src_ring], I_D, I_D)
        end
        for i in s.tgt_observers
            add_sheaf_edge!(F, s_range[i], spec.target_nodes[s.tgt_ring], I_D, I_D)
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
    
    p_boundary = vcat([boundary[r_idx] for r_idx in 1:length(target_positions)]...)
    q_interior = vec(H_mat \ (-B_mat * p_boundary))

    n_pf_verts = nv(PfF.underlying_graph) > 0 ? nv(PfF.underlying_graph) : length(boundary) + length(q_interior) ÷ vertex_stalks(PfF)[end]
    q_H = Vector{Vector{Float64}}(undef, n_pf_verts)
    
    for (r_idx, p) in enumerate(target_positions)
        q_H[r_idx] = boundary[r_idx]
    end
    
    int_idx = 1
    for v in 1:length(q_H)
        if !haskey(boundary, v)
            d_v = vertex_stalks(PfF)[v]
            q_H[v] = q_interior[(int_idx-1)*d_v + 1 : int_idx*d_v]
            int_idx += 1
        end
    end

    return q_H
end

"""
    solve_mid_level_harmonic(q_H::Vector{Vector{Float64}}, bases::LayeredFiberBases, D::Int=4) -> Matrix{Float64}

Lift coarse pushforward section `q_H` to D-dimensional agent reference states on `G`.
"""
function solve_mid_level_harmonic(q_H::Vector{Vector{Float64}}, bases::LayeredFiberBases, D::Int=4)
    n_pf_nodes = length(q_H)
    q_G_components = [bases.fiber_bases[v] * q_H[v] for v in 1:n_pf_nodes]

    n_agents = sum(size(bases.fiber_bases[v], 1) ÷ D for v in 1:n_pf_nodes) - length(bases.target_subbases)
    q_agents = zeros(n_agents, D)

    curr_agent = 1
    for v in 1:n_pf_nodes
        B = bases.fiber_bases[v]
        q_fib = q_G_components[v]
        n_fib_nodes = size(B, 1) ÷ D
        
        has_target = haskey(bases.target_subbases, v)
        n_fib_agents = has_target ? n_fib_nodes - 1 : n_fib_nodes

        for i in 1:n_fib_agents
            q_agents[curr_agent, :] .= q_fib[(i-1)*D + 1 : i*D]
            curr_agent += 1
        end
    end

    return q_agents
end

"""
    solve_direct_harmonic(F::EuclideanSheaf, target_nodes::Vector{Int}, target_positions::Vector{Vector{Float64}}, D::Int=F.vertex_stalks[1]) -> Matrix{Float64}

Solve harmonic extension directly on `F` given pinned target nodes and positions.
"""
function solve_direct_harmonic(F::EuclideanSheaf, target_nodes::Vector{Int}, target_positions::Vector{Vector{Float64}}, D::Int=F.vertex_stalks[1])
    boundary = Dict(target_nodes[i] => target_positions[i] for i in 1:length(target_nodes))
    _, _, H_mat, B_mat = _harmonic_extension_restricted_laplacian(F, boundary)
    
    p_boundary = vcat(target_positions...)
    q_interior = H_mat \ (-B_mat * p_boundary)

    n_total = length(F.vertex_stalks)
    n_targets = length(target_nodes)
    n_agents = n_total - n_targets

    q_full = zeros(n_total, D)
    
    int_idx = 1
    for v in 1:n_total
        if v in target_nodes
            t_idx = findfirst(==(v), target_nodes)
            q_full[v, :] .= target_positions[t_idx]
        else
            q_full[v, :] .= q_interior[(int_idx-1)*D + 1 : int_idx*D]
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
    bindings::SystemBinding
    target_trajectories::Vector{Any} # functions t -> [x,y,z,1.0]
    target_velocities::Union{Nothing, Vector{Any}}   # functions t -> [vx,vy,vz,0.0]
    target_accelerations::Union{Nothing, Vector{Any}}# functions t -> [ax,ay,az,0.0]
    dt::Float64
    steps::Int
    initial_positions::Union{Nothing, Vector{Vector{Float64}}} # world position per agent at t=0
end

"""
    LayeredEscortProblem(spec, bindings, target_trajectories; dt=0.05, steps=200, target_velocities=nothing, target_accelerations=nothing, initial_positions=nothing)

High-level constructor that automatically builds `F`, `f`, `PfF`, and structured fiber bases from `spec`.

`bindings` is a `NestedSystems.SystemBinding` resolved against `spec` (see
`NestedSystems.resolve_dynamics(::LayeredEscortSpec, ::SystemBinding)`) — a single root-level
`SystemBinding(dynamics=dyn)` binds every agent to the same model, matching the old
`HomogeneousDynamics(dyn)` behaviour exactly.

Each agent's initial position resolves in three tiers, most specific first:
1. `AgentBinding.initial_position` from `bindings` (a per-agent override)
2. `initial_positions[i]` (this constructor's own field, if supplied)
3. the "airstrip" default (`NestedSystems._default_initial_position`) — agents lined up along
   the first position coordinate with fixed spacing, at a common hover altitude, deliberately
   distinct from the target formation so the control law's convergence into formation is
   visible in the simulation rather than starting there already.

Existing callers that pass `initial_positions` and no per-agent bindings see unchanged
behaviour: `resolve_dynamics` only fills in a position when `bindings` supplies one, so tier 2
governs whenever tier 1 is absent.
"""
function LayeredEscortProblem(
    spec::LayeredEscortSpec,
    bindings::SystemBinding,
    target_trajectories::Vector;
    target_velocities=nothing,
    target_accelerations=nothing,
    dt::Float64=0.05,
    steps::Int=200,
    initial_positions=nothing
)
    sheaf = build_layered_escort_sheaf(spec)
    hom = build_layered_homomorphism(spec)
    pf_sheaf = pushforward_sheaf(hom, sheaf)
    bases = build_layered_fiber_bases(hom, sheaf, spec)

    return LayeredEscortProblem(
        spec, sheaf, hom, pf_sheaf, bases,
        bindings, target_trajectories, target_velocities, target_accelerations,
        dt, steps, initial_positions
    )
end

struct LayeredEscortResult
    problem::LayeredEscortProblem
    sim_data::Vector{Vector{Vector{Float64}}}
    qstar_history::Vector{Matrix{Float64}}
    target_history::Vector{Vector{Vector{Float64}}}
end

"""
    run_layered_escort_simulation(prob::LayeredEscortProblem; use_feedforward::Bool=false) -> LayeredEscortResult

Run complete multi-agent layered escort simulation using closed-loop agent dynamics and optional feedforward tracking.
"""
function run_layered_escort_simulation(prob::LayeredEscortProblem; use_feedforward::Bool=false)
    spec = prob.spec
    bases = prob.bases
    STEPS = prob.steps
    DT = prob.dt
    D = spec.D
    time_grid = 0:DT:(STEPS*DT)

    # Initialize agent states with dynamic configurations. Bindings are resolved once, up
    # front, rather than re-walked per agent per step.
    agent_states = Vector{AgentState}()
    for (name, local_idx, i, b) in _resolve_flat_bindings(spec, prob.bindings)
        dyn = b.dynamics
        K_lqr = something(b.K_lqr, Some(zeros(0, 0)))
        # Three tiers, most specific first: an explicit per-agent binding, then this
        # problem's own `initial_positions`, then the airstrip default.
        pos0 = something(b.initial_position,
                         isnothing(prob.initial_positions) ? nothing : Some(prob.initial_positions[i]),
                         Some(_default_initial_position(dyn, i)))
        x0 = initial_state(dyn, pos0)
        push!(agent_states, AgentState(x0, dyn, DT, K_lqr, 0.02; use_velocity=use_feedforward))
    end

    sim_data = Vector{Vector{Vector{Float64}}}()
    qstar_history = Vector{Matrix{Float64}}()
    target_history = Vector{Vector{Vector{Float64}}}()

    for step in 1:STEPS
        t = time_grid[step]
        
        # 1. Target positions
        p_targets = [traj(t) for traj in prob.target_trajectories]
        
        # 2. High-level harmonic extension
        q_H = solve_high_level_harmonic(prob.pf_sheaf, bases, p_targets)
        
        # Target world positions for visualization
        target_world_centers = Vector{Vector{Float64}}()
        for r_idx in 1:spec.n_rings
            push!(target_world_centers, p_targets[r_idx])
        end
        for s_idx in 1:spec.n_supports
            pf_v = spec.n_rings + s_idx
            center = bases.fiber_bases[pf_v][1:D, :] * q_H[pf_v]
            push!(target_world_centers, center)
        end
        push!(target_history, target_world_centers)

        # 3. Mid-level reference lifting q*
        qstar_agents = solve_mid_level_harmonic(q_H, bases, D)
        push!(qstar_history, qstar_agents)

        # Optional Feedforward Velocity and Acceleration Lifting
        qstar_dot_agents = nothing
        qstar_ddot_agents = nothing
        if use_feedforward && !isnothing(prob.target_velocities)
            v_targets = [v_traj(t) for v_traj in prob.target_velocities]
            v_H = solve_high_level_harmonic(prob.pf_sheaf, bases, v_targets)
            qstar_dot_agents = solve_mid_level_harmonic(v_H, bases, D)
        end
        if use_feedforward && !isnothing(prob.target_accelerations)
            a_targets = [a_traj(t) for a_traj in prob.target_accelerations]
            a_H = solve_high_level_harmonic(prob.pf_sheaf, bases, a_targets)
            qstar_ddot_agents = solve_mid_level_harmonic(a_H, bases, D)
        end

        # 4. Step Agent Dynamics with Feedforward Support
        step_states = Vector{Vector{Float64}}()
        for i in 1:spec.n_agents
            pos_dim_i = length(position_indices(agent_states[i].dyn))
            q_ref_i = qstar_agents[i, 1:pos_dim_i]

            if use_feedforward && !isnothing(qstar_dot_agents) && !isnothing(qstar_ddot_agents)
                q_dot_i = qstar_dot_agents[i, 1:pos_dim_i]
                q_ddot_i = qstar_ddot_agents[i, 1:pos_dim_i]
                step_agent!(agent_states[i], q_ref_i, q_dot_i, q_ddot_i, DT)
            elseif use_feedforward && !isnothing(qstar_dot_agents)
                q_dot_i = qstar_dot_agents[i, 1:pos_dim_i]
                step_agent!(agent_states[i], q_ref_i, q_dot_i, DT)
            else
                step_agent!(agent_states[i], q_ref_i, DT)
            end

            push!(step_states, copy(agent_states[i].x))
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
