# N-agent cooperative planning from a sheaf on the interaction graph, never forming the
# joint state graph. Every agent works in the gauge of its own single-agent cost-to-go
# dᵢ (moves cost their *delay* δᵢ = cᵢ + dᵢ(next) - dᵢ(current) ≥ 0), so the joint cost
# is Σᵢ dᵢ plus the total delay, and a heuristic only has to bound the delay.
#
#   - Vertex stalks: each agent's own state graph (the free cells, moves, delays).
#   - Edge stalks: for a pair (i, j), agent i's graph glued to the pushforward of agent j's
#     graph along its *progress map*: each cell goes to its arc of the Reeb graph of dⱼ
#     (its component of the band dⱼ(cell) ≤ dⱼ ≤ dⱼ(cell) + 1; delay-free moves go down
#     the Reeb graph, and equidistant cells on opposite sides of an obstacle stay apart),
#     with agent j kept exact within distance r of agent i, where the
#     collision coupling lives. The pair's coarse graph over-approximates the quotient
#     of the two-agent joint graph, so its cost-to-go Δᵢⱼ is a lower bound on the pair's
#     delay; it is built from the two agents' graphs alone.
#   - Gluing: an N-agent plan restricts to a feasible plan of every pair, and delays are
#     nonnegative, so Σ over a matching of the interaction graph of Δᵢⱼ bounds the total
#     delay. The heuristic is Σᵢ dᵢ + the best matching.
#
# Checked against the exact N-agent cost-to-go and against the exact pairwise delays (the
# best any pairwise gluing can do) on a small map with three agents.
#
#   julia --project=docs docs/scripts/tropical_multiagent_sheaf.jl [n] [queries]

using CellularSheaves
using SparseArrays
using Printf
using Random
using Statistics

const n = length(ARGS) >= 1 ? parse(Int, ARGS[1]) : 10
const NQ = length(ARGS) >= 2 ? parse(Int, ARGS[2]) : 8
const NA = 3                                     # agents
const R = 2                                      # near-field radius
const OUT = joinpath(@__DIR__, "..", "figures")

# A two-cell-thick wall with two one-cell corridors.
const WALL = (n ÷ 2):(n ÷ 2 + 1)
wall(x, y) = x in WALL && !(y in (3, n - 2))
const CELLS = [(x, y) for y in 1:n for x in 1:n if !wall(x, y)]
const INDEX = Dict(c => k for (k, c) in enumerate(CELLS))
const NF = length(CELLS)
const NBRS = [[k; [INDEX[(x + dx, y + dy)] for (dx, dy) in ((1, 0), (-1, 0), (0, 1), (0, -1))
                   if haskey(INDEX, (x + dx, y + dy))]] for (k, (x, y)) in enumerate(CELLS)]

function distances(g)
    d = fill(Inf, NF)
    d[g] = 0
    queue = [g]
    while !isempty(queue)
        u = popfirst!(queue)
        for v in NBRS[u]
            isinf(d[v]) && (d[v] = d[u] + 1; push!(queue, v))
        end
    end
    return d
end

# ---------------------------------------------------------------------------
# Vertex stalks: one agent, its goal and its delays.
# ---------------------------------------------------------------------------

struct AgentStalk
    goal::Int
    d::Vector{Float64}                           # single-agent cost-to-go
    comp::Vector{Int}                            # progress map: cell -> Reeb arc
    out::Vector{Dict{Int,Float64}}               # pushforward of the delays along comp
    stay::Vector{Float64}                        # delay of staying in a component
end
cost(A::AgentStalk, a, a2) = a == a2 == A.goal ? 0.0 : 1.0
delay(A::AgentStalk, a, a2) = cost(A, a, a2) + A.d[a2] - A.d[a]

# The Reeb graph of d: a level-k cell's arc is its component of the band {k ≤ d ≤ k + 1}
# (the discrete level set d = k + 1/2), so equidistant cells joined only through cells much
# nearer to or farther from the goal (opposite sides of an obstacle) are kept apart.
# Returns an arc index per cell.
function level_components(d)
    comp, ids = zeros(Int, NF), Dict{Tuple{Int,Int},Int}()
    parent = collect(1:NF)
    root(x) = (while parent[x] != x; parent[x] = parent[parent[x]]; x = parent[x]; end; x)
    for k in sort(unique(filter(isfinite, d)))
        band = [c for c in 1:NF if k <= d[c] <= k + 1]
        for c in band
            parent[c] = c
        end
        for c in band, u in NBRS[c]
            if k <= d[u] <= k + 1
                ru, rc = root(u), root(c)
                ru != rc && (parent[ru] = rc)
            end
        end
        for c in band
            d[c] == k && (comp[c] = get!(ids, (Int(k), root(c)), length(ids) + 1))
        end
    end
    return comp
end

function AgentStalk(goal)
    proto = AgentStalk(goal, distances(goal), Int[], Dict{Int,Float64}[], Float64[])
    comp = level_components(proto.d)
    C = maximum(comp)
    best = [Dict{Int,Float64}() for _ in 1:C]
    for b in 1:NF, b2 in NBRS[b]
        c, c2 = comp[b], comp[b2]
        c == c2 && continue
        best[c][c2] = min(get(best[c], c2, Inf), delay(proto, b, b2))
    end
    stay = [comp[goal] == c ? 0.0 : 1.0 for c in 1:C]
    return AgentStalk(goal, proto.d, comp, best, stay)
end

# ---------------------------------------------------------------------------
# Edge stalks: the pair abstraction of agent i with a summary of agent j.
# ---------------------------------------------------------------------------

cheb(o) = max(abs(o[1]), abs(o[2]))
const NNEAR = (2R + 1)^2
near_index(o) = (o[1] + R) + (2R + 1) * (o[2] + R) + 1
near_offset(m) = ((m - 1) % (2R + 1) - R, (m - 1) ÷ (2R + 1) - R)
cell_at(a, o) = get(INDEX, CELLS[a] .+ o, 0)

# The mode of agent j as seen from agent i: its offset when near, else its progress component.
mode(Ai, Aj, a, b) = (o = CELLS[b] .- CELLS[a]; cheb(o) <= R ? near_index(o) : NNEAR + Aj.comp[b])

struct PairAbstraction
    K::Int                                       # modes of agent j
    Δ::Vector{Float64}                           # coarse delay-to-go over (cell of i, mode of j)
end
state(P::PairAbstraction, a, m) = a + NF * (m - 1)

function PairAbstraction(Ai::AgentStalk, Aj::AgentStalk)
    K = NNEAR + length(Aj.stay)
    rows, cols, vals = Int[], Int[], Float64[]
    edge!(x, y, c) = (push!(rows, x); push!(cols, y); push!(vals, c))
    for a in 1:NF, m in 1:K
        x = a + NF * (m - 1)
        if m <= NNEAR                                            # j exact, near i
            o = near_offset(m)
            b = cell_at(a, o)
            (o == (0, 0) || b == 0) && continue                  # unoccupied mode
            for a2 in NBRS[a], b2 in NBRS[b]
                (a2 == b2 || (a2 == b && b2 == a)) && continue  # vertex and swap conflicts
                edge!(x, a2 + NF * (mode(Ai, Aj, a2, b2) - 1), delay(Ai, a, a2) + delay(Aj, b, b2))
            end
        else                                                     # j far, in component c
            c = m - NNEAR
            for a2 in NBRS[a]
                δi = delay(Ai, a, a2)
                # stay far: along the pushforward of j's moves, or stay in c
                edge!(x, a2 + NF * (NNEAR + c - 1), δi + Aj.stay[c])
                for (c2, δj) in Aj.out[c]
                    edge!(x, a2 + NF * (NNEAR + c2 - 1), δi + δj)
                end
                # come near: a cell around a2 in c (staying) or in a component c leads to
                for m2 in 1:NNEAR
                    o2 = near_offset(m2)
                    b2 = cell_at(a2, o2)
                    (o2 == (0, 0) || b2 == 0) && continue
                    c2 = Aj.comp[b2]
                    δj = c2 == c ? Aj.stay[c] : get(Aj.out[c], c2, Inf)
                    isfinite(δj) && edge!(x, a2 + NF * (m2 - 1), δi + δj)
                end
            end
        end
    end
    Wc = sparse(rows, cols, vals, NF * K, NF * K, min)
    bc = fill(Inf, NF * K)
    bc[Ai.goal + NF * (mode(Ai, Aj, Ai.goal, Aj.goal) - 1)] = 0.0
    return PairAbstraction(K, minplus_closure(Wc, bc))
end

# The pair's delay bound at the cells (a, b) of agents i and j, from either side.
pair_bound(Pij, Pji, Ai, Aj, a, b) =
    max(Pij.Δ[state(Pij, a, mode(Ai, Aj, a, b))], Pji.Δ[state(Pji, b, mode(Aj, Ai, b, a))])

# ---------------------------------------------------------------------------
# The joint problem (for checking only): NA agents, conflict-free joint moves.
# ---------------------------------------------------------------------------

joint(cells) = 1 + sum((cells[i] - 1) * NF^(i - 1) for i in 1:NA)
cells_of(x) = ntuple(i -> ((x - 1) ÷ NF^(i - 1)) % NF + 1, NA)
conflict(c, c2) = any(c2[i] == c2[j] || (c2[i] == c[j] && c2[j] == c[i]) for i in 1:NA for j in (i + 1):NA)

function joint_graph(agents)
    rows, cols, vals = Int[], Int[], Float64[]
    for x in 1:NF^NA
        c = cells_of(x)
        allunique(c) || continue
        for c2 in Iterators.product((NBRS[ci] for ci in c)...)
            conflict(c, c2) && continue
            push!(rows, x); push!(cols, joint(c2)); push!(vals, sum(cost(agents[i], c[i], c2[i]) for i in 1:NA))
        end
    end
    return sparse(rows, cols, vals, NF^NA, NF^NA)
end

# Exact pairwise delays: the two-agent joint problem of (i, j), others ignored.
function pair_graph(Ai, Aj)
    rows, cols, vals = Int[], Int[], Float64[]
    for b in 1:NF, a in 1:NF
        a == b && continue
        for a2 in NBRS[a], b2 in NBRS[b]
            (a2 == b2 || (a2 == b && b2 == a)) && continue
            push!(rows, a + NF * (b - 1)); push!(cols, a2 + NF * (b2 - 1))
            push!(vals, delay(Ai, a, a2) + delay(Aj, b, b2))
        end
    end
    return sparse(rows, cols, vals, NF^2, NF^2)
end
exact_pair(Ai, Aj) = cost_to_go(pair_graph(Ai, Aj), [Ai.goal + NF * (Aj.goal - 1)])

# The exact pushforward of the pair graph along the same summary (agent i's cell, agent j's
# near offset or Reeb arc): what the cheap construction over-approximates. It needs the pair
# graph, NF² states.
function pushed_pair(Ai, Aj)
    K = NNEAR + length(Aj.stay)
    vmap = [a == b ? a : a + NF * (mode(Ai, Aj, a, b) - 1) for b in 1:NF for a in 1:NF]
    A = Abstraction(pair_graph(Ai, Aj), GraphHomomorphism(vmap, NF * K))
    return lower_bound(A, [Ai.goal + NF * (Aj.goal - 1)])
end

# The best matching of three agents is the best single pair.
const PAIRS = [(i, j) for i in 1:NA for j in (i + 1):NA]
glue_pairs(bounds) = maximum(bounds)

function heuristics(agents, Ps, pushed, exact_pairs)
    N = NF^NA
    dec, sheaf, pushforward, pairwise = fill(Inf, N), fill(Inf, N), fill(Inf, N), fill(Inf, N)
    for x in 1:N
        c = cells_of(x)
        allunique(c) || continue
        base = sum(agents[i].d[c[i]] for i in 1:NA)
        dec[x] = base
        sheaf[x] = base + glue_pairs(pair_bound(Ps[(i, j)], Ps[(j, i)], agents[i], agents[j], c[i], c[j]) for (i, j) in PAIRS)
        pushforward[x] = base + glue_pairs(max(pushed[(i, j)][c[i] + NF * (c[j] - 1)], pushed[(j, i)][c[j] + NF * (c[i] - 1)])
                                           for (i, j) in PAIRS)
        pairwise[x] = base + glue_pairs(exact_pairs[(i, j)][c[i] + NF * (c[j] - 1)] for (i, j) in PAIRS)
    end
    return dec, sheaf, pushforward, pairwise
end

function queries(rng)
    left = [k for (k, (x, y)) in enumerate(CELLS) if x < first(WALL) - 1]
    right = [k for (k, (x, y)) in enumerate(CELLS) if x > last(WALL) + 1]
    cross = map(1:NQ) do _
        # agents 1 and 2 cross left to right, agent 3 right to left
        s12, g3 = left[randperm(rng, length(left))[1:3]][1:2], left[randperm(rng, length(left))[1]]
        g12, s3 = right[randperm(rng, length(right))[1:2]], right[randperm(rng, length(right))[1]]
        starts, goals = (s12..., s3), (g12..., g3)
        (allunique(starts) && allunique(goals)) ? ("cross", starts, goals) : nothing
    end
    randoms = [("random", Tuple(p[1:NA]), Tuple(p[(NA + 1):(2NA)])) for p in (randperm(rng, NF) for _ in 1:NQ)]
    return filter(!isnothing, [cross; randoms])
end

println("grid $(n)², wall with two corridors: $(NF) free cells, $(NA) agents, $(NF^NA) joint states")
rows = NamedTuple[]
for (kind, starts, goals) in queries(MersenneTwister(3))
    agents = [AgentStalk(g) for g in goals]
    t_sheaf = @elapsed Ps = Dict((i, j) => PairAbstraction(agents[i], agents[j]) for i in 1:NA for j in 1:NA if i != j)
    t_pairs = @elapsed exact_pairs = Dict((i, j) => exact_pair(agents[i], agents[j]) for (i, j) in PAIRS)
    t_pushed = @elapsed pushed = Dict((i, j) => pushed_pair(agents[i], agents[j]) for i in 1:NA for j in 1:NA if i != j)
    W = joint_graph(agents)
    goal, start = joint(goals), joint(starts)
    V = cost_to_go(W, [goal])
    isfinite(V[start]) || continue
    # Every pair bound lies below the exact pair delay (the edge stalk is a relaxation).
    for (i, j) in PAIRS, b in 1:NF, a in 1:NF
        a == b && continue
        @assert pair_bound(Ps[(i, j)], Ps[(j, i)], agents[i], agents[j], a, b) <= exact_pairs[(i, j)][a + NF * (b - 1)] + 1e-9
    end
    dec, sheaf, pushforward, pairwise = heuristics(agents, Ps, pushed, exact_pairs)
    reach = isfinite.(V)
    out = sparse(transpose(W))
    coarse = sum(P.K * NF for P in values(Ps))
    @printf("%-6s  optimal %3.0f  decoupled %3.0f  (coordination cost %.0f)\n", kind, V[start], dec[start], V[start] - dec[start])
    for (name, h, build, states) in (("decoupled (sum of distances)", dec, 0.0, 0),
                                     ("sheaf: pair abstractions, glued", sheaf, t_sheaf, coarse),
                                     ("sheaf: exact pair pushforwards, glued", pushforward, t_pushed, coarse),
                                     ("exact pairwise delays, glued", pairwise, t_pairs, length(PAIRS) * NF^2),
                                     ("exact cost-to-go", V, 0.0, NF^NA))
        @assert all(h[reach] .<= V[reach] .+ 1e-9) "$name is not admissible"
        r = shortest_path(W, start, [goal]; heuristic = h, out_edges = out)
        @assert r.cost == V[start]
        push!(rows, (kind = kind, heuristic = name, states = states, optimal = V[start], decoupled = dec[start],
            bound = h[start], expanded = r.expanded, build_seconds = build))
        @printf("    %-40s bound %5.1f  expanded %7d  build %6.3f s\n", name, h[start], r.expanded, build)
    end
end

println("\nsummary (means over queries)")
@printf("%-7s %-40s %13s %10s %11s %9s %9s\n", "queries", "heuristic", "states", "bound/opt", "gap closed", "expanded", "build s")
for kind in ("cross", "random"), name in unique(r.heuristic for r in rows)
    rs = filter(r -> r.kind == kind && r.heuristic == name, rows)
    isempty(rs) && continue
    gap = [(r.bound - r.decoupled) / (r.optimal - r.decoupled) for r in rs if r.optimal > r.decoupled]
    @printf("%-7s %-40s %13d %10.3f %11s %9.0f %9.3f\n", kind, name, round(Int, mean(r.states for r in rs)),
        mean(r.bound / r.optimal for r in rs), isempty(gap) ? "–" : @sprintf("%.2f", mean(gap)),
        mean(r.expanded for r in rs), mean(r.build_seconds for r in rs))
end

mkpath(OUT)
open(joinpath(OUT, "tropical_multiagent_sheaf.csv"), "w") do io
    println(io, join(keys(first(rows)), ","))
    for r in rows
        println(io, join((occursin(',', string(v)) ? "\"$v\"" : v for v in values(r)), ","))
    end
end
