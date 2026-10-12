# The interaction sheaf of cooperative multi-agent planning, shared by the studies
# `tropical_multiagent_sheaf.jl` (grid agents) and `tropical_multiagent_double_integrator.jl`.
#
#   - Vertex stalk: an agent's own state graph on a shared map (each state has a cell and
#     a velocity), in the gauge of its own cost-to-go d: a move costs its delay
#     c + d(next) - d(current) ≥ 0.
#   - Edge stalk of a pair (i, j): the pushforward of the pair graph (both agents moving at
#     once, collisions removed) along the map that keeps agent i's state and replaces
#     agent j's by
#       * its exact state, when j's block of the map is next to i's (the near region,
#         which contains every collision), or
#       * the pair (Reeb arc of dⱼ, block) otherwise: its progress, and where it is
#         relative to i at the resolution of the blocks.
#     The far coordinate depends on j's state alone and the near/far split on blocks, so
#     the far region is a union of fibres of a product map, and there the pushforward is
#     the tensor product of i's graph with the quotient of j's graph (base change). Only
#     the near region and its boundary are enumerated. The result is the exact
#     pushforward, computed without the pair graph; its cost-to-go Δᵢⱼ bounds the pair's
#     delay.
#   - Gluing: an N-agent plan restricts to every pair and delays are nonnegative, so the
#     sum of Δᵢⱼ over a matching of the interaction graph bounds the total delay.

using CellularSheaves
using SparseArrays

struct GridMap
    cells::Vector{NTuple{2,Int}}
    index::Dict{NTuple{2,Int},Int}
end
GridMap(cells) = GridMap(cells, Dict(c => k for (k, c) in enumerate(cells)))
cell_at(M::GridMap, c, o) = get(M.index, M.cells[c] .+ o, 0)
cheb(o) = max(abs(o[1]), abs(o[2]))
block(M::GridMap, c, b) = cld.(M.cells[c], b)
block_distance(M::GridMap, c, c2, b) = cheb(block(M, c, b) .- block(M, c2, b))

# An agent's own system: `pos[s]` and `vel[s]` of each state, successor lists, a goal state.
struct AgentStalk
    pos::Vector{Int}
    vel::Vector{Int}
    nvel::Int
    state::Dict{Tuple{Int,Int},Int}              # (cell, velocity) -> state
    succ::Vector{Vector{Int}}
    pred::Vector{Vector{Int}}
    goal::Int
    d::Vector{Float64}                           # cost-to-go to the goal
    delays::SparseMatrixCSC{Float64,Int}         # the agent's graph in the gauge of d
    arcs::GraphHomomorphism                      # Reeb arcs of d
end

step_cost(goal, s, s2) = s == s2 == goal ? 0.0 : 1.0
delay(A::AgentStalk, s, s2) = step_cost(A.goal, s, s2) + A.d[s2] - A.d[s]

function AgentStalk(pos, vel, succ, goal)
    S = length(pos)
    pred = [Int[] for _ in 1:S]
    rows, cols, vals = Int[], Int[], Float64[]
    for s in 1:S, s2 in succ[s]
        push!(pred[s2], s)
        push!(rows, s); push!(cols, s2); push!(vals, step_cost(goal, s, s2))
    end
    W = sparse(rows, cols, vals, S, S)
    d = cost_to_go(W, [goal])
    keep = [isfinite(d[r]) && isfinite(d[c]) for (r, c) in zip(rows, cols)]
    delays = sparse(rows[keep], cols[keep], [vals[k] + d[cols[k]] - d[rows[k]] for k in findall(keep)], S, S)
    return AgentStalk(pos, vel, maximum(vel), Dict((pos[s], vel[s]) => s for s in 1:S), succ, pred, goal, d, delays,
        reeb_arcs(W, d))
end

# ---------------------------------------------------------------------------
# The edge stalk of a pair.
# ---------------------------------------------------------------------------

struct PairAbstraction
    b::Int                                       # block size
    nnear::Int                                   # near modes: cells of a 3 × 3 block window × velocities of j
    far::GraphHomomorphism                       # j's state -> (Reeb arc, block)
    far_block::Vector{NTuple{2,Int}}             # the block of each far summary
    Δ::Vector{Float64}                           # coarse delay-to-go over (state of i, mode of j)
end

# j's far summary: its Reeb arc and its block.
function far_summary(M, Aj, b)
    ids = Dict{Tuple{Int,NTuple{2,Int}},Int}()
    vmap = [get!(ids, (Aj.arcs.vertex_map[t], block(M, Aj.pos[t], b)), length(ids) + 1) for t in eachindex(Aj.pos)]
    blocks = Vector{NTuple{2,Int}}(undef, length(ids))
    for ((_, B), k) in ids
        blocks[k] = B
    end
    return GraphHomomorphism(vmap, length(ids)), blocks
end

near(M, Ai, Aj, b, s, t) = block_distance(M, Ai.pos[s], Aj.pos[t], b) <= 1
function near_index(M, Ai, Aj, b, s, t)
    corner = (block(M, Ai.pos[s], b) .- 2) .* b                     # the window's cells are corner .+ (1:3b)
    l = M.cells[Aj.pos[t]] .- corner
    return l[1] + 3b * (l[2] - 1) + (3b)^2 * (Aj.vel[t] - 1)
end

# The mode of agent j (in state t) as seen from agent i (in state s).
pair_mode(M, P::PairAbstraction, Ai, Aj, s, t) =
    near(M, Ai, Aj, P.b, s, t) ? near_index(M, Ai, Aj, P.b, s, t) : P.nnear + P.far.vertex_map[t]
coarse(P::PairAbstraction, Si, s, m) = s + Si * (m - 1)

# The states of agent j near agent i's state s.
function near_states(M, Ai, Aj, b, s)
    corner = (block(M, Ai.pos[s], b) .- 2) .* b
    out = Int[]
    for ly in 1:3b, lx in 1:3b, w in 1:Aj.nvel
        t = get(Aj.state, (get(M.index, corner .+ (lx, ly), 0), w), 0)
        t != 0 && Aj.pos[t] != Ai.pos[s] && isfinite(Aj.d[t]) && push!(out, t)
    end
    return out
end

function PairAbstraction(M::GridMap, Ai::AgentStalk, Aj::AgentStalk, b::Int)
    Si = length(Ai.pos)
    far, far_block = far_summary(M, Aj, b)
    nnear = (3b)^2 * Aj.nvel
    K = nnear + far.n_target
    P = PairAbstraction(b, nnear, far, far_block, Float64[])
    rows, cols, vals = Int[], Int[], Float64[]
    edge!(x, y, c) = (push!(rows, x); push!(cols, y); push!(vals, c))
    # Far region: the tensor product of i's graph with the quotient of j's graph along the
    # far summary (base change), restricted to coarse states whose blocks are apart.
    Q = quotient_graph(Aj.delays, far; loops = true)
    apart(s, c) = cheb(block(M, Ai.pos[s], b) .- far_block[c]) >= 2
    decode(x) = ((x - 1) % Si + 1, (x - 1) ÷ Si + 1)
    T = tensor_product(Ai.delays, Q; keep = (x, y) -> apart(decode(x)...) && apart(decode(y)...))
    r, c, v = findnz(T)
    append!(rows, r .+ Si * nnear); append!(cols, c .+ Si * nnear); append!(vals, v)
    # Near region, exact: from i's state s and j's exact state t near it.
    for s in 1:Si
        isfinite(Ai.d[s]) || continue
        for t in near_states(M, Ai, Aj, b, s)
            x = coarse(P, Si, s, near_index(M, Ai, Aj, b, s, t))
            for s2 in Ai.succ[s], t2 in Aj.succ[t]
                (Ai.pos[s2] == Aj.pos[t2] || (Ai.pos[s2] == Aj.pos[t] && Aj.pos[t2] == Ai.pos[s])) && continue
                (isfinite(Ai.d[s2]) && isfinite(Aj.d[t2])) || continue
                edge!(x, coarse(P, Si, s2, pair_mode(M, P, Ai, Aj, s2, t2)), delay(Ai, s, s2) + delay(Aj, t, t2))
            end
        end
    end
    # Boundary, far to near: j enters the window from a state of its far fibre.
    for s2 in 1:Si
        isfinite(Ai.d[s2]) || continue
        for t2 in near_states(M, Ai, Aj, b, s2)
            y = coarse(P, Si, s2, near_index(M, Ai, Aj, b, s2, t2))
            for s in Ai.pred[s2], t in Aj.pred[t2]
                (isfinite(Ai.d[s]) && isfinite(Aj.d[t]) && !near(M, Ai, Aj, b, s, t)) || continue
                edge!(coarse(P, Si, s, nnear + far.vertex_map[t]), y, delay(Ai, s, s2) + delay(Aj, t, t2))
            end
        end
    end
    Wc = sparse(rows, cols, vals, Si * K, Si * K, min)
    bc = fill(Inf, Si * K)
    bc[coarse(P, Si, Ai.goal, pair_mode(M, P, Ai, Aj, Ai.goal, Aj.goal))] = 0.0
    return PairAbstraction(b, nnear, far, far_block, minplus_closure(Wc, bc))
end

coarse_states(P::PairAbstraction) = count(isfinite, P.Δ)

# The pair's delay bound at states (s, t) of agents i and j, from either side.
pair_bound(M, Pij, Pji, Ai, Aj, s, t) =
    max(Pij.Δ[coarse(Pij, length(Ai.pos), s, pair_mode(M, Pij, Ai, Aj, s, t))],
        Pji.Δ[coarse(Pji, length(Aj.pos), t, pair_mode(M, Pji, Aj, Ai, t, s))])

# ---------------------------------------------------------------------------
# Checks that need the pair graph (small maps only).
# ---------------------------------------------------------------------------

# The pair graph of agents i and j in the delay gauge, with vertex and swap conflicts.
function pair_graph(Ai::AgentStalk, Aj::AgentStalk)
    Si, Sj = length(Ai.pos), length(Aj.pos)
    rows, cols, vals = Int[], Int[], Float64[]
    for t in 1:Sj, s in 1:Si
        Ai.pos[s] == Aj.pos[t] && continue
        (isfinite(Ai.d[s]) && isfinite(Aj.d[t])) || continue
        for s2 in Ai.succ[s], t2 in Aj.succ[t]
            (Ai.pos[s2] == Aj.pos[t2] || (Ai.pos[s2] == Aj.pos[t] && Aj.pos[t2] == Ai.pos[s])) && continue
            (isfinite(Ai.d[s2]) && isfinite(Aj.d[t2])) || continue
            push!(rows, s + Si * (t - 1)); push!(cols, s2 + Si * (t2 - 1))
            push!(vals, delay(Ai, s, s2) + delay(Aj, t, t2))
        end
    end
    return sparse(rows, cols, vals, Si * Sj, Si * Sj)
end
pair_goal(Ai, Aj) = Ai.goal + length(Ai.pos) * (Aj.goal - 1)

# The pushforward of the pair graph along the same map, computed from the pair graph: a
# check that PairAbstraction is exact.
function pushed_pair(M, P::PairAbstraction, Ai, Aj, W = pair_graph(Ai, Aj))
    Si, Sj = length(Ai.pos), length(Aj.pos)
    K = P.nnear + P.far.n_target
    vmap = [Ai.pos[s] == Aj.pos[t] || !isfinite(Aj.d[t]) ? s : coarse(P, Si, s, pair_mode(M, P, Ai, Aj, s, t))
            for t in 1:Sj for s in 1:Si]
    return lower_bound(Abstraction(W, GraphHomomorphism(vmap, Si * K)), [pair_goal(Ai, Aj)])
end
