# src/sheaf_consensus.jl

module SheafConsensus

export build_agent_sheaf, build_target_sheaf,
       sheaf_control_output, sheaf_target_control_output,
       pack_cochain, unpack_cochain,
       global_section_basis, project_to_global_section

using LinearAlgebra
using CellularSheaves
using BlockArrays: BlockArray


function make_agent_restriction_maps(nd::Int)
    R = Matrix{Float64}(I, nd, nd)
    return R, R
end

"""
    make_pinning_restriction_maps(weight, nd)

For a pinning edge from agent i to target j with scalar weight `w`:
    R_agent = I * sqrt(|w|),   R_target = I * sqrt(|w|)  (same stalk dimension)
so that the coboundary residual carries the weight naturally.
"""
function make_pinning_restriction_maps(weight::Float64, nd::Int)
    sw = sqrt(abs(weight))
    R = sw * Matrix{Float64}(I, nd, nd)
    return R, R        # (R_agent, R_target)
end

# ─────────────────────────────────────────────────────────────────────────────
# Sheaf construction for the agent–(agent+target) graph
# ─────────────────────────────────────────────────────────────────────────────
 
function build_agent_sheaf(n_agents::Int, n_targets::Int, nd::Int,
                            agent_edge_set,   # Vector of (i,j) 1-based pairs
                            pinning_matrix)   # n_agents × n_targets matrix

    n_verts = n_agents + n_targets
    sheaf = EuclideanSheaf{Float64}(fill(nd, n_verts))

    # Agent–agent consensus edges
    for (i, j) in agent_edge_set
        Ri, Rj = make_agent_restriction_maps(nd)
        add_sheaf_edge!(sheaf, i, j, Ri, Rj)
    end

    # Pinning edges: agent i ↔ target j
    for i in 1:n_agents
        for j in 1:n_targets
            w = pinning_matrix[i, j]
            if w != 0.0
                Ragent, Rtarget = make_pinning_restriction_maps(Float64(w), nd)
                agent_v  = i
                target_v = n_agents + j
                add_sheaf_edge!(sheaf, agent_v, target_v, Ragent, Rtarget)
            end
        end
    end

    return sheaf, collect(1:n_agents), collect(n_agents+1:n_agents+n_targets)
end


function build_target_sheaf(n_targets::Int, nd::Int, target_edge_set)
    sheaf = EuclideanSheaf{Float64}(fill(nd, n_targets))
    for (i, j) in target_edge_set
        Ri, Rj = make_agent_restriction_maps(nd)
        add_sheaf_edge!(sheaf, i, j, Ri, Rj)
    end
    return sheaf
end

# ─────────────────────────────────────────────────────────────────────────────
# Cochain packing / unpacking
# ─────────────────────────────────────────────────────────────────────────────

function pack_cochain(positions_agents::Matrix{Float64},
                       offsets_agents::Matrix{Float64},
                       positions_targets::Matrix{Float64},
                       offsets_targets::Matrix{Float64},
                       step::Int,
                       nd::Int) :: Vector{Float64}

    n_agents  = size(positions_agents,  2)
    n_targets = size(positions_targets, 2)
    z = zeros(Float64, nd * (n_agents + n_targets))

    # Agent block: offset-shifted position
    for i in 1:n_agents
        z[(i-1)*nd+1 : i*nd] = positions_agents[:, step] .+ offsets_agents[:, i]
    end

    # Target block: offset-shifted position (targets have their own offsets for
    # the tetrahedral formation, but the pinning term sees raw target position)
    for j in 1:n_targets
        z[(n_agents + j - 1)*nd + 1 : (n_agents + j)*nd] = positions_targets[:, step]
    end

    return z
end

"""
    unpack_cochain(z, n_agents, nd) -> (agent_block, target_block)
"""
function unpack_cochain(z::Vector{Float64}, n_agents::Int, nd::Int)
    agent_block  = reshape(z[1 : n_agents*nd], nd, n_agents)
    target_block = reshape(z[n_agents*nd+1 : end], nd, :)
    return agent_block, target_block
end

# Control law via sheaf Laplacian

function sheaf_control_output(sheaf::EuclideanSheaf,
                               z::Vector{Float64},
                               n_agents::Int,
                               nd::Int,
                               k_agents::Float64) :: Matrix{Float64}

    d   = coboundary_map(sheaf)
    Lz  = d' * (d * z)           # = L * z,  shape: nd*(n_agents+n_targets)

    u = zeros(Float64, nd, n_agents)
    for i in 1:n_agents
        u[:, i] = -k_agents * Lz[(i-1)*nd+1 : i*nd]
    end
    return u
end


function sheaf_target_control_output(sheaf::EuclideanSheaf,
                                      positions_targets::Matrix{Float64},
                                      offsets_targets::Matrix{Float64},
                                      step::Int,
                                      nd::Int,
                                      k_targets::Float64) :: Matrix{Float64}

    n_targets = size(positions_targets, 2)
    z = zeros(Float64, nd * n_targets)
    for j in 1:n_targets
        z[(j-1)*nd+1 : j*nd] = positions_targets[:, step] .+ offsets_targets[:, j]
    end

    d  = coboundary_map(sheaf)
    Lz = d' * (d * z)

    u = zeros(Float64, nd, n_targets)
    for j in 1:n_targets
        u[:, j] = -k_targets * Lz[(j-1)*nd+1 : j*nd]
    end
    return u
end

# ─────────────────────────────────────────────────────────────────────────────
# Analysis helpers

# Return a matrix whose columns span the space of global sections of `sheaf`.
# Each column is a valid equilibrium formation state (up to translation /
# scaling depending on the restriction maps).

function global_section_basis(sheaf::EuclideanSheaf)
    return nullspace_ldlt(sheaf)
end

function project_to_global_section(sheaf::EuclideanSheaf, z::Vector{Float64};
                                   method::Symbol=:ldl)
    gs = nearest_global_section(sheaf, z; method=method)
    return collect(gs)
end

function project_to_global_section(sheaf::EuclideanSheaf, z::BlockArray;
                                   method::Symbol=:ldl)
    return nearest_global_section(sheaf, z; method=method)
end

end # module SheafConsensus
