# The interaction sheaf of cooperative multi-agent planning, shared by the studies
# `tropical_multiagent_sheaf.jl` (grid agents) and `tropical_multiagent_double_integrator.jl`.
#
#   - Vertex stalk: an agent's own state graph on a shared map (each state has a cell and
#     a velocity), in the gauge of its own cost-to-go d: a move costs its delay
#     c + d(next) - d(current) ≥ 0.
#   - Progress map: each state goes to its arc of the Reeb graph of d (its component of the
#     band k ≤ d ≤ k + 1, k = d(state)). The agent's moves pushed forward along it give the
#     arc dynamics; each arc remembers its footprint, the cells its states occupy.
#   - Edge stalk of a pair (i, j): agent i's own states × a mode of agent j: j's exact state
#     (relative cell, velocity) when within distance R of i, where the collision coupling
#     lives, else j's arc. Its coarse graph is built from the two agents' graphs alone and
#     over-approximates the pushforward of the pair graph along that map: far modes exist
#     only where the arc has a cell outside i's window (the far-from-me condition, checked
#     on the footprint), and j entering the window is checked on its actual predecessors.
#     Its cost-to-go Δᵢⱼ bounds the pair's delay.
#   - Gluing: an N-agent plan restricts to every pair, and delays are nonnegative, so the
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
offset(M::GridMap, c, c2) = M.cells[c2] .- M.cells[c]

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
    arc::Vector{Int}                             # Reeb arc of d
    out::Vector{Dict{Int,Float64}}               # pushforward of the delays along the arcs
    stay::Vector{Float64}                        # cheapest delay of a move inside an arc
    footprint::Vector{Vector{Int}}               # cells of each arc's states
end

step_cost(A::AgentStalk, s, s2) = s == s2 == A.goal ? 0.0 : 1.0
delay(A::AgentStalk, s, s2) = step_cost(A, s, s2) + A.d[s2] - A.d[s]

function AgentStalk(pos, vel, succ, goal)
    S = length(pos)
    pred = [Int[] for _ in 1:S]
    rows, cols, vals = Int[], Int[], Float64[]
    for s in 1:S, s2 in succ[s]
        push!(pred[s2], s)
        push!(rows, s); push!(cols, s2); push!(vals, s == s2 == goal ? 0.0 : 1.0)
    end
    d = cost_to_go(sparse(rows, cols, vals, S, S), [goal])
    arc = reeb_arcs(d, [union(succ[s], pred[s]) for s in 1:S])
    C = maximum(arc)
    proto = AgentStalk(pos, vel, maximum(vel), Dict((pos[s], vel[s]) => s for s in 1:S), succ, pred, goal, d, arc,
        Dict{Int,Float64}[], Float64[], Vector{Int}[])
    out = [Dict{Int,Float64}() for _ in 1:C]
    stay = fill(Inf, C)
    footprint = [Int[] for _ in 1:C]
    for s in 1:S
        arc[s] == 0 && continue
        push!(footprint[arc[s]], pos[s])
        for s2 in succ[s]
            arc[s2] == 0 && continue
            c, c2, δ = arc[s], arc[s2], delay(proto, s, s2)
            c == c2 ? (stay[c] = min(stay[c], δ)) : (out[c][c2] = min(get(out[c], c2, Inf), δ))
        end
    end
    return AgentStalk(pos, vel, proto.nvel, proto.state, succ, pred, goal, d, arc, out, stay, unique!.(footprint))
end

# The Reeb graph of d on a graph with neighbour lists `adj`: a level-k state's arc is its
# component of the band {k ≤ d ≤ k + 1}. States that cannot reach the goal get arc 0.
function reeb_arcs(d, adj)
    S = length(d)
    arc, ids = zeros(Int, S), Dict{Tuple{Float64,Int},Int}()
    parent = collect(1:S)
    root(x) = (while parent[x] != x; parent[x] = parent[parent[x]]; x = parent[x]; end; x)
    for k in sort(unique(filter(isfinite, d)))
        band = [s for s in 1:S if k <= d[s] <= k + 1]
        for s in band
            parent[s] = s
        end
        for s in band, u in adj[s]
            if k <= d[u] <= k + 1
                ru, rs = root(u), root(s)
                ru != rs && (parent[ru] = rs)
            end
        end
        for s in band
            d[s] == k && (arc[s] = get!(ids, (k, root(s)), length(ids) + 1))
        end
    end
    return arc
end

# ---------------------------------------------------------------------------
# The edge stalk of a pair.
# ---------------------------------------------------------------------------

struct PairAbstraction
    R::Int
    nnear::Int                                   # (2R + 1)² offsets × velocities of j
    Δ::Vector{Float64}                           # coarse delay-to-go over (state of i, mode of j)
end

near_index(R, o, w) = (o[1] + R) + (2R + 1) * (o[2] + R) + 1 + (2R + 1)^2 * (w - 1)
function near_mode(R, m)
    q = (m - 1) % (2R + 1)^2
    return ((q % (2R + 1)) - R, (q ÷ (2R + 1)) - R), (m - 1) ÷ (2R + 1)^2 + 1
end

# The mode of agent j (in state t) as seen from agent i (in state s).
function pair_mode(M, Ai, Aj, R, s, t)
    o = offset(M, Ai.pos[s], Aj.pos[t])
    return cheb(o) <= R ? near_index(R, o, Aj.vel[t]) : (2R + 1)^2 * Aj.nvel + Aj.arc[t]
end

# Whether arc c has a cell outside the window of radius R around cell p.
far_possible(M, Aj, R, c, p) = any(q -> cheb(offset(M, p, q)) > R, Aj.footprint[c])

function PairAbstraction(M::GridMap, Ai::AgentStalk, Aj::AgentStalk, R::Int)
    Si, nnear = length(Ai.pos), (2R + 1)^2 * Aj.nvel
    K = nnear + length(Aj.stay)
    rows, cols, vals = Int[], Int[], Float64[]
    edge!(x, y, c) = (push!(rows, x); push!(cols, y); push!(vals, c))
    for s in 1:Si, m in 1:K
        isfinite(Ai.d[s]) || continue
        x = s + Si * (m - 1)
        p = Ai.pos[s]
        if m <= nnear                                            # j exact, near i
            o, w = near_mode(R, m)
            t = get(Aj.state, (cell_at(M, p, o), w), 0)
            (o == (0, 0) || t == 0 || Aj.arc[t] == 0) && continue
            for s2 in Ai.succ[s], t2 in Aj.succ[t]
                (Ai.pos[s2] == Aj.pos[t2] || (Ai.pos[s2] == Aj.pos[t] && Aj.pos[t2] == p)) && continue
                Aj.arc[t2] == 0 && continue
                edge!(x, s2 + Si * (pair_mode(M, Ai, Aj, R, s2, t2) - 1), delay(Ai, s, s2) + delay(Aj, t, t2))
            end
        else                                                     # j far, in arc c
            c = m - nnear
            far_possible(M, Aj, R, c, p) || continue             # no state of c is far from i
            for s2 in Ai.succ[s]
                δi, p2 = delay(Ai, s, s2), Ai.pos[s2]
                # stay far: inside c, or along the pushforward of j's moves
                isfinite(Aj.stay[c]) && far_possible(M, Aj, R, c, p2) && edge!(x, s2 + Si * (m - 1), δi + Aj.stay[c])
                for (c2, δj) in Aj.out[c]
                    far_possible(M, Aj, R, c2, p2) && edge!(x, s2 + Si * (nnear + c2 - 1), δi + δj)
                end
                # come near: j's state t2 in i's next window, from a predecessor in c far from i
                for m2 in 1:nnear
                    o2, w2 = near_mode(R, m2)
                    o2 == (0, 0) && continue
                    t2 = get(Aj.state, (cell_at(M, p2, o2), w2), 0)
                    (t2 == 0 || Aj.arc[t2] == 0) && continue
                    δj = minimum((delay(Aj, t, t2) for t in Aj.pred[t2]
                                  if Aj.arc[t] == c && cheb(offset(M, p, Aj.pos[t])) > R); init = Inf)
                    isfinite(δj) && edge!(x, s2 + Si * (m2 - 1), δi + δj)
                end
            end
        end
    end
    Wc = sparse(rows, cols, vals, Si * K, Si * K, min)
    bc = fill(Inf, Si * K)
    bc[Ai.goal + Si * (pair_mode(M, Ai, Aj, R, Ai.goal, Aj.goal) - 1)] = 0.0
    return PairAbstraction(R, nnear, minplus_closure(Wc, bc))
end

coarse_states(P::PairAbstraction) = count(isfinite, P.Δ)

# The pair's delay bound at states (s, t) of agents i and j, from either side.
pair_bound(M, Pij, Pji, Ai, Aj, s, t) =
    max(Pij.Δ[s + length(Ai.pos) * (pair_mode(M, Ai, Aj, Pij.R, s, t) - 1)],
        Pji.Δ[t + length(Aj.pos) * (pair_mode(M, Aj, Ai, Pji.R, t, s) - 1)])

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

# The exact pushforward of the pair graph along the same map: what PairAbstraction
# over-approximates.
function pushed_pair(M, Ai, Aj, R, W = pair_graph(Ai, Aj))
    Si, Sj = length(Ai.pos), length(Aj.pos)
    K = (2R + 1)^2 * Aj.nvel + length(Aj.stay)
    vmap = [Ai.pos[s] == Aj.pos[t] || Aj.arc[t] == 0 ? s : s + Si * (pair_mode(M, Ai, Aj, R, s, t) - 1)
            for t in 1:Sj for s in 1:Si]
    return lower_bound(Abstraction(W, GraphHomomorphism(vmap, Si * K)), [pair_goal(Ai, Aj)])
end
