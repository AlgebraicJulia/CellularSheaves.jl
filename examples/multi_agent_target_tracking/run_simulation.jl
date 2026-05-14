# examples/run_main_sheaf.jl

using LinearAlgebra
using JSON3
using CellularSheaves

abstract type AbstractDynamics end
struct NoDynamics <: AbstractDynamics end
struct TrophicDynamics <: AbstractDynamics end
struct CurrentAgentDynamics <: AbstractDynamics end
struct CurrentTargetDynamics <: AbstractDynamics end

include(joinpath(@__DIR__, "src", "simulation", "integrate.jl"))
include(joinpath(@__DIR__, "src", "simulation", "dynamics.jl"))
include(joinpath(@__DIR__, "src", "core", "entity.jl"))
include(joinpath(@__DIR__, "src", "io", "data_manager.jl"))
include(joinpath(@__DIR__, "src", "core", "cellular_sheaf.jl"))

using .SheafConsensus

# function translating the dynamics_type string in the config to a Julia type
function dynamics_type_from_string(s::AbstractString)
    s = lowercase(String(s))
    if s in ("none", "nodynamics", "no_dynamics")
        return NoDynamics
    elseif s in ("trophic", "trophic_dynamics", "trophicdynamics")
        return TrophicDynamics
    elseif s in ("current_agent", "currentagent", "current_agent_dynamics")
        return CurrentAgentDynamics
    elseif s in ("current_target", "currenttarget", "current_target_dynamics")
        return CurrentTargetDynamics
    else
        error("Unknown dynamics_type=$s.")
    end
end

# this function takes the JSON config and turns the "agents" 
# or "targets" section into a clean list of entity specifications.
function build_entity_specs(base_config, section_key::AbstractString)
    specs = Tuple{Vector{Float64}, Dict{String,Any}}[]
    for item_any in base_config[section_key]
        item = Dict{String,Any}(item_any)
        dyn = String(item["dynamics_type"])
        nstates = Int(get(item, "num_states", get(base_config, "num_states", 3)))
        init_pos = haskey(item, "initial_position") ?
                   Float64.(collect(item["initial_position"])) : zeros(Float64, nstates)
        merged = merge(Dict{String,Any}(base_config), item)
        merged["dynamics_type"] = dyn
        push!(specs, (init_pos, merged))
    end
    return specs
end

# function to build the neighborhood sets for agents and targets from the 
# edge sets in the config
function construct_undirected_neighborhood_set!(entities, edge_set)
    for e in edge_set
        i = Int(e[1]); j = Int(e[2])
        a = entities[i]; b = entities[j]
        b in a.neighbors || push!(a.neighbors, b)
        a in b.neighbors || push!(b.neighbors, a)
    end
end

# Build adjacency lists from the agent edge set.
function agent_neighbor_indices(n_agents::Int, edge_set)
    nbrs = [Int[] for _ in 1:n_agents]
    for e in edge_set
        i = Int(e[1]); j = Int(e[2])
        push!(nbrs[i], j); push!(nbrs[j], i)
    end
    return nbrs
end
# example: if n_agents=5 and edge_set=[(1,2), (2,3), (4,5)], then this returns
# [[2] (neighbors of 1), [1,3] (neighbors of 2), [2] (neighbors of 3), [5] (neighbors of 4), [4] (neighbors of 5)]

function target_neighbor_indices(n_targets::Int, edge_set)
    nbrs = [Int[] for _ in 1:n_targets]
    for e in edge_set
        i = Int(e[1]); j = Int(e[2])
        1 <= i <= n_targets && 1 <= j <= n_targets || continue
        push!(nbrs[i], j); push!(nbrs[j], i)
    end
    return nbrs
end

# Helper to chunk a list of agent indices into groups of 3 for triangle formation.
# example: if idxs=[1,2,3,4,5,6,7], this returns [[1,2,3], [4,5,6]] and ignores the leftover [7].
chunk_triplets(idxs::Vector{Int}) =
    [idxs[i:min(i+2,end)] for i in 1:3:length(idxs) if length(idxs[i:min(i+2,end)]) == 3]

function make_offsets_agents_triangle(agent_specs, target_specs, neighbors,
                                       group::Vector{Int}, target_index::Int)
    nd = Int(agent_specs[1][2]["num_states"])
    N  = length(agent_specs)
    offsets = zeros(Float64, nd, N)
    target_pos = target_specs[target_index][1]
    l = 20.0
    tri_formation = [
        [0.0, l / sqrt(3), 0.0],
        [l / 2, -l * sqrt(3) / 6, 0.0],
        [-l / 2, -l * sqrt(3) / 6, 0.0],
    ]
    g = group[1:3]
    desired_positions = [target_pos .+ v for v in tri_formation]
    for i_local in eachindex(g)
        gi = g[i_local]
        for j in neighbors[gi]
            j in g || continue
            j_local = findfirst(==(j), g)
            delta_ij = desired_positions[j_local] .- desired_positions[i_local]
            offsets[:, gi] .+= delta_ij ./ length(g)
            offsets[:, j]  .-= delta_ij ./ length(g)
        end
    end
    return offsets
end

function make_offsets_agents_square(agent_specs, target_specs, neighbors,
                                     group::Vector{Int}, target_index::Int)
    nd = Int(agent_specs[1][2]["num_states"])
    N  = length(agent_specs)
    offsets = zeros(Float64, nd, N)
    target_pos = target_specs[target_index][1]
    l = 20.0
    square_formation = [
        [-l/2, -l/2, 0.0],
        [ l/2, -l/2, 0.0],
        [ l/2,  l/2, 0.0],
        [-l/2,  l/2, 0.0],
    ]
    g = group[1:4]
    desired_positions = [target_pos .+ v for v in square_formation]
    for i_local in eachindex(g)
        gi = g[i_local]
        for j in neighbors[gi]
            j in g || continue
            j_local = findfirst(==(j), g)
            delta_ij = desired_positions[j_local] .- desired_positions[i_local]
            offsets[:, gi] .+= delta_ij ./ length(g)
            offsets[:, j]  .-= delta_ij ./ length(g)
        end
    end
    return offsets
end

function make_offsets_targets(target_specs, neighbors)
    nd = Int(target_specs[1][2]["num_states"])
    N  = length(target_specs)
    target_offsets = zeros(Float64, nd, N)
    N >= 4 || return target_offsets
    target_pos_1 = target_specs[1][1]
    d = 15.0; h = sqrt(6) * d / 3
    tet_formation = [
        [0.0, 0.0, 0.0],
        [d, 0.0, 0.0],
        [d/2, sqrt(3)*d/2, 0.0],
        [d/2, sqrt(3)*d/6, h],
    ]
    desired_positions = [target_pos_1 .+ v for v in tet_formation]
    for i in 1:4, j in neighbors[i]
        i < j && j <= 4 || continue
        delta_ij = desired_positions[j] .- desired_positions[i]
        target_offsets[:, i] .+= delta_ij ./ 4
        target_offsets[:, j] .-= delta_ij ./ 4
    end
    return target_offsets
end

# =====================================================
# Main sheaf-based simulation
# =====================================================

function run_config_sheaf(path::AbstractString; run_analysis::Bool=false)
    cfg = Dict{String,Any}(JSON3.read(read(path, String), Dict{String,Any}))

    final_time = Float64(cfg["final_time"])
    dt         = Float64(cfg["time_step_delta"])
    steps      = Int(final_time / dt)

    target_specs = build_entity_specs(cfg, "targets")
    agent_specs  = build_entity_specs(cfg, "agents")

    nd        = Int(agent_specs[1][2]["num_states"])
    n_targets = length(target_specs)
    n_agents  = length(agent_specs)

    # ── Formation offsets (identical to baseline) ──────────────────────────
    agent_edge_set  = cfg["agent_edge_set"]
    target_edge_set = get(cfg, "target_edge_set", Any[])
    agent_nbr_idx   = agent_neighbor_indices(n_agents,  agent_edge_set)
    target_nbr_idx  = target_neighbor_indices(n_targets, target_edge_set)

    k_square    = 4
    n_agents < k_square && error("Need at least 4 agents for square formation.")
    square_group = collect(n_agents - k_square + 1 : n_agents)  # ex: if n_agents=10, this is [7,8,9,10]
    tri_pool     = collect(1 : n_agents - k_square)             # ex: if n_agents=10, this is [1,2,3,4,5,6]
    tri_groups   = chunk_triplets(tri_pool)                     # ex: if tri_pool=[1,2,3,4,5,6], this is [[1,2,3], [4,5,6]]

    sq_t_idx  = 1
    tri_t_idx = n_targets >= 2 ? 2 : 1

    offsets_agents = zeros(Float64, nd, n_agents)
    for g in tri_groups
        off = make_offsets_agents_triangle(agent_specs, target_specs,
                                           agent_nbr_idx, g, tri_t_idx)
        offsets_agents[:, g] .+= off[:, g]
    end
    off_sq = make_offsets_agents_square(agent_specs, target_specs,
                                        agent_nbr_idx, square_group, sq_t_idx)
    offsets_agents[:, square_group] .+= off_sq[:, square_group]

    offsets_targets = make_offsets_targets(target_specs, target_nbr_idx)

    # ── Pinning matrix ──────────────────────────────────────────────────────
    pinning_matrix = Float64.(hcat([Float64.(collect(row))
                                    for row in cfg["pinning_matrix"]]...)')
    # pinning_matrix is n_agents × n_targets

    # ── Build SHEAVES ────────────────────────────────────────────────────────
    # Convert edge sets to Vector{Tuple{Int,Int}}
    agent_edges  = [(Int(e[1]), Int(e[2])) for e in agent_edge_set]
    target_edges = [(Int(e[1]), Int(e[2])) for e in target_edge_set]

    println("Building agent–target consensus sheaf...")
    agent_sheaf, agent_verts, target_verts =
        SheafConsensus.build_agent_sheaf(n_agents, n_targets, nd, agent_edges, pinning_matrix)
        # ex: if n_agents=10, n_targets=4, this is a sheaf on a graph with 14 vertices (10 agents + 4 targets)

    println("Building target formation sheaf...")
    target_sheaf = SheafConsensus.build_target_sheaf(n_targets, nd, target_edges)
    # ex: if n_targets=4, this is a sheaf on a graph with 4 vertices (the targets) and edges given by target_edges

    println("Agent–target sheaf: $(length(agent_verts)+length(target_verts)) vertices, ",
            length(agent_edges) + count(!=(0.0), pinning_matrix), " edges")
    println("Target sheaf: $(n_targets) vertices, $(length(target_edges)) edges")

    # ── Instantiate entities ─────────────────────────────────────────────────
    k_agents  = Float64(cfg["agents_proportional_gain"])
    k_targets = Float64(cfg["targets_proportional_gain"])

    targets = AbstractEntity[]
    for (pos, conf) in target_specs
        D = dynamics_type_from_string(String(conf["dynamics_type"]))
        id = String(conf["id"])
        tracking_type = String(get(conf, "tracking_type", ""))
        is_centroid   = Bool(get(conf, "is_centroid", false)) || tracking_type == "f8_dynamics"
        push!(targets, Target{D}(pos, steps, id, nd, dt, k_targets;
                                  is_centroid=is_centroid, tracking_type=tracking_type))
    end
    for j in eachindex(targets)
        targets[j].offsets .= offsets_targets[:, j]
    end

    agents = AbstractEntity[]
    for (pos, conf) in agent_specs
        D = dynamics_type_from_string(String(conf["dynamics_type"]))
        id = String(conf["id"])
        push!(agents, Agent{D}(pos, steps, id, nd, dt, k_agents, targets))
    end
    for i in eachindex(agents)
        agents[i].offsets .= offsets_agents[:, i]
    end

    # Pinning rows
    for i in eachindex(agents)
        agents[i].pin_row .= pinning_matrix[i, :]
    end

    # Neighborhoods (needed for update_dynamics! which reads entity.neighbors
    # only for f8 tracking; not for sheaf control)
    construct_undirected_neighborhood_set!(targets, target_edge_set)
    construct_undirected_neighborhood_set!(agents,  agent_edge_set)

    pos_agents  = [ag.positions  for ag in agents]
    pos_targets = [tg.positions  for tg in targets]

    # ── Data manager ─────────────────────────────────────────────────────────
    outdir = joinpath(@__DIR__, "simulation_data")
    dm = DataManager(outdir)
    init_data_manager!(dm, agents, targets)

    # ── Pre-compute coboundary maps once (they don't change) ─────────────────
    d_agent  = coboundary_map(agent_sheaf)
    d_target = coboundary_map(target_sheaf)
    # ex: if agent_sheaf has 14 vertices and 20 edges, then d_agent is a 20×14 matrix.
    # ex: if target_sheaf has 4 vertices and 4 edges, then d_target is a 4×4 matrix.

    # Helper: assemble agent+target cochain from step `s` columns
    function build_agent_cochain(s::Int)
        z = zeros(Float64, nd * (n_agents + n_targets))
        for i in 1:n_agents
            z[(i-1)*nd+1 : i*nd] = pos_agents[i][:, s] .+ offsets_agents[:, i]
        end
        for j in 1:n_targets
            z[(n_agents+j-1)*nd+1 : (n_agents+j)*nd] = pos_targets[j][:, s]
        end
        return z
    end
    # ex: cochain of length 42 for 10 agents + 4 targets with nd=3, where the first 30 entries are the offset-shifted 
    #agent positions and the last 12 entries are the target positions.
    # z would look like [ag1_x, ag1_y, ag1_z, ag2_x, ag2_y, ag2_z, ..., tg1_x, tg1_y, tg1_z, tg2_x, tg2_y, tg2_z, ...]
    # and every different row would be representing the state of a different entity (agent or target) at step s.

    function build_target_cochain(s::Int)
        z = zeros(Float64, nd * n_targets)
        for j in 1:n_targets
            z[(j-1)*nd+1 : j*nd] = pos_targets[j][:, s] .+ offsets_targets[:, j]
        end
        return z
    end

    # ── Simulation loop ───────────────────────────────────────────────────────
    println("Running sheaf-based simulation for $steps steps...")

    for step in 2:steps
        # Build cochains from previous step positions
        z_agents  = build_agent_cochain(step - 1)
        z_targets = build_target_cochain(step - 1)

        # Sheaf Laplacian action  L z = d^T (d z)
        Lz_agents  = d_agent'  * (d_agent  * z_agents)
        Lz_targets = d_target' * (d_target * z_targets)

        # control outputs
        for i in eachindex(agents)
            err = Lz_agents[(i-1)*nd+1 : i*nd]

            agents[i].synchronization_error .= err

            u = -k_agents * err
            agents[i].control_output .= u

            if nd >= 3
                agents[i].control_output[3] = 0.0
            end
        end


        #    Targets: u_j = -k/100 * (Lz)[target_j_block] + desired_vel
        for j in eachindex(targets)
            time_now = (step - 1) * dt
            desired_vel = zeros(Float64, nd)
            if targets[j].is_centroid && targets[j].tracking_type == "f8_dynamics"
                desired_vel .= f8_dynamics(time_now)
            end

            err = Lz_targets[(j-1)*nd+1 : j*nd]

            targets[j].synchronization_error .= err

            u = -k_targets * err ./ 100.0
            targets[j].control_output .= u .+ desired_vel
        end

        # Integrate dynamics
        for ag in agents;  update_dynamics!(ag, step); end
        for tg in targets; update_dynamics!(tg, step); end

        save_step!(dm, agents, targets, step, dt; python_style=true)

        if step % max(1, Int(floor(steps / 100))) == 0
            pct = round(100 * step / steps, digits=2)
            print("Progress: $pct%\r"); flush(stdout)
        end
    end

    close_all_files!(dm)
    println("\nCSVs saved to $outdir")

    # ── sheaf analysis ───────────────────────────────────────────────
    if run_analysis
        println("\n── Sheaf Analysis ─────────────────────────────────────────")

        # Global sections of the agent–target sheaf
        println("Computing global sections of agent–target sheaf...")
        V = SheafConsensus.global_section_basis(agent_sheaf)
        println("  Global section space dimension: $(size(V, 2))")
        println("  (= number of formation DOF; should be ≥ nd for translation)")

        # Projection of final configuration onto nearest global section
        z_final = build_agent_cochain(steps)
        gs      = SheafConsensus.project_to_global_section(agent_sheaf, z_final; method=:ldl)
        d_mat   = coboundary_map(agent_sheaf)
        residual = norm(d_mat * gs)
        println("  Coboundary residual of projected final state: $residual")
        println("  (→ 0 means agents have reached a global section = ideal formation)")

        # Per-agent tracking error at final step
        println("\n── Per-agent tracking error (final step) ──────────────────")
        for ag in agents
            j = findfirst(!=(0.0), ag.pin_row)
            j === nothing && continue
            tgt     = targets[j]
            dT_xy   = norm(ag.positions[1:2, end] .- tgt.positions[1:2, end])
            z_err   = ag.positions[3, end] - tgt.positions[3, end]
            println("  $(ag.id) → $(tgt.id): XY dist = $(round(dT_xy,digits=3)),  z err = $(round(z_err,digits=3))")
        end
        println("────────────────────────────────────────────────────────────")
    end

    return agents, targets, steps, agent_sheaf, target_sheaf
end

# ── Entry point ───────────────────────────────────────────────────────────────
cfg_path = joinpath(@__DIR__, "configurations", "config_common.json")
run_analysis = "--analysis" in ARGS

agents, targets, steps, agent_sheaf, target_sheaf =
    run_config_sheaf(cfg_path; run_analysis=run_analysis)

println("Done. steps=$steps")
println("Run with `--analysis` flag for sheaf diagnostics (global sections, residuals).")
