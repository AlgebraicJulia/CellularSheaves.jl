# Cooperative two-agent planning with interaction-aware abstractions. Each agent's
# coarse problem keeps its own cell exactly and a low-dimensional summary of the other
# agent: its relative offset when near (the encounter, kept exactly), its cell when in a
# corridor (where it cannot be passed), and its block of a room otherwise. Pushed
# forward along that summary, the joint graph gives an admissible lower bound on the
# joint cost; how much of the coordination cost (waiting, detours) does it see, against
# the decoupled sum of single-agent distances and the exact joint cost-to-go?
#
# The coarse problems bound the *delay* rather than the cost: each agent is charged
# c + d(next) - d(current) ≥ 0 (d its single-agent distance to its goal), a potential
# shaping under which the joint cost-to-go is exactly V - (d₁ + d₂). Without it, coarsening
# the other agent throws away its travel cost, and the bound never beats the decoupled sum.
#
#   julia --project=docs docs/scripts/tropical_multiagent_study.jl [n] [queries]

using CellularSheaves
using SparseArrays
using Printf
using Random
using Statistics

const n = length(ARGS) >= 1 ? parse(Int, ARGS[1]) : 32
const NQ = length(ARGS) >= 2 ? parse(Int, ARGS[2]) : 8
const OUT = joinpath(@__DIR__, "..", "figures")

# A thick wall (columns 14–19) with two one-cell corridors (rows 16 and 28): agents
# that swap sides must queue through a corridor or detour to the other one.
wall(x, y) = (n ÷ 2 - 2) <= x <= (n ÷ 2 + 3) && !(y in (n ÷ 2, n - 4))
const CELLS = [(x, y) for y in 1:n for x in 1:n if !wall(x, y)]
const INDEX = Dict(c => k for (k, c) in enumerate(CELLS))
const NF = length(CELLS)
const NBRS = [[k; [INDEX[(x + dx, y + dy)] for (dx, dy) in ((1, 0), (-1, 0), (0, 1), (0, -1))
                   if haskey(INDEX, (x + dx, y + dy))]] for (k, (x, y)) in enumerate(CELLS)]

# Joint states (a, b) with a ≠ b, indexed a + NF (b - 1) (states with a = b are unused and
# isolated). A joint move moves both agents one step (or keeps them), with no vertex
# conflict (same cell) and no swap conflict (exchanging cells).
joint(a, b) = a + NF * (b - 1)
agents(x) = ((x - 1) % NF + 1, (x - 1) ÷ NF + 1)

function joint_structure()
    rows, cols = Int[], Int[]
    for b in 1:NF, a in 1:NF
        a == b && continue
        for a2 in NBRS[a], b2 in NBRS[b]
            (a2 == b2 || (a2 == b && b2 == a)) && continue
            push!(rows, joint(a, b)); push!(cols, joint(a2, b2))
        end
    end
    return sparse(rows, cols, ones(length(rows)), NF^2, NF^2)
end

# Each agent pays 1 per step, waiting included, except while it stays at its own goal.
# Returns the joint costs and the shaped joint costs (the agents' delays), on the same
# stored pattern.
function costs(S, g1, g2, d1, d2)
    v, shaped = similar(nonzeros(S)), similar(nonzeros(S))
    rv = rowvals(S)
    for y in 1:size(S, 2), p in nzrange(S, y)
        (a, b), (a2, b2) = agents(rv[p]), agents(y)
        c1 = a == a2 == g1 ? 0.0 : 1.0
        c2 = b == b2 == g2 ? 0.0 : 1.0
        v[p] = c1 + c2
        shaped[p] = (c1 + d1[a2] - d1[a]) + (c2 + d2[b2] - d2[b])
    end
    @assert all(>=(0), shaped)
    mk(w) = SparseMatrixCSC(size(S)..., S.colptr, rv, w)
    return mk(v), mk(shaped)
end

function bfs(g)
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

# Summaries of the other agent, as seen by the agent whose cell is `own`: a mode in 1:K.
cheb(o) = max(abs(o[1]), abs(o[2]))
offset(own, other) = CELLS[other] .- CELLS[own]
near_index(o, r) = (o[1] + r) + (2r + 1) * (o[2] + r) + 1
block_index(c, B) = (CELLS[c][1] - 1) ÷ cld(n, B) + B * ((CELLS[c][2] - 1) ÷ cld(n, B)) + 1

struct Summary{F}
    name::String
    K::Int
    mode::F                 # (own, other, distance of the other to its goal) -> 1:K
end

solo() = Summary("solo (1 mode)", 1, (own, other, d) -> 1)
near(r) = Summary("near r=$r", (2r + 1)^2 + 1, function (own, other, d)
    o = offset(own, other)
    cheb(o) <= r ? near_index(o, r) : (2r + 1)^2 + 1
end)
# Rooms and corridors: the other agent's cell when it is in a corridor (it cannot be
# passed there), else its block of a B × B grid, split by the side of the wall.
const CORRIDOR = [k for (k, (x, y)) in enumerate(CELLS) if (n ÷ 2 - 2) <= x <= (n ÷ 2 + 3)]
const CORRIDOR_INDEX = Dict(k => i for (i, k) in enumerate(CORRIDOR))
room_index(c, B) = haskey(CORRIDOR_INDEX, c) ? CORRIDOR_INDEX[c] :
    length(CORRIDOR) + block_index(c, B) + (CELLS[c][1] > n ÷ 2 ? B^2 : 0)
rooms(B) = Summary("rooms $(B)×$(B) + corridors", length(CORRIDOR) + 2B^2, (own, other, d) -> room_index(other, B))
rooms_near(B, r) = Summary("rooms $(B)×$(B) + corridors + near r=$r", (2r + 1)^2 + length(CORRIDOR) + 2B^2,
    function (own, other, d)
        o = offset(own, other)
        cheb(o) <= r ? near_index(o, r) : (2r + 1)^2 + room_index(other, B)
    end)

# Progress: the other agent's remaining distance to its goal and its side of the wall
# (its cell in a corridor, its offset when near). Along zero-delay moves the distance
# drops by one per step, so the coarse problem keeps the other agent's timing.
const DMAX = 4n
progress_index(c, d) = haskey(CORRIDOR_INDEX, c) ? CORRIDOR_INDEX[c] :
    length(CORRIDOR) + 2 * Int(d[c]) + (CELLS[c][1] > n ÷ 2 ? 2 : 1)
progress_near(r) = Summary("progress + corridors + near r=$r", (2r + 1)^2 + length(CORRIDOR) + 2DMAX + 2,
    function (own, other, d)
        o = offset(own, other)
        cheb(o) <= r ? near_index(o, r) : (2r + 1)^2 + progress_index(other, d)
    end)

# The mode map of agent `i` (1 or 2) on the joint states: (own cell, summary of the other).
function mode_map(s::Summary, i, dists)
    vmap = Vector{Int}(undef, NF^2)
    d_other = dists[3 - i]
    for x in 1:NF^2
        a, b = agents(x)
        own, other = i == 1 ? (a, b) : (b, a)
        vmap[x] = own + NF * (s.mode(own, other, d_other) - 1)
    end
    return GraphHomomorphism(vmap, NF * s.K)
end

# The lower bound of agent i's abstraction along summary s, for costs Wc, and the number
# of coarse states the joint states actually occupy.
function bound(Wc, s, i, goal, dists)
    μ = mode_map(s, i, dists)
    A = Abstraction(Wc, μ)
    occupied = length(unique(μ.vertex_map[x] for x in 1:NF^2 if !=(agents(x)...)))
    return lower_bound(A, [goal]), occupied
end

# Heuristics: name => (W, shaped W, goal, decoupled, (d₁, d₂)) -> (heuristic over the joint states, coarse states).
# Two agents' bounds of the same joint cost combine by max, not sum.
unshaped(summary) = function (W, Ws, goal, dec, dists)
    (h1, c1), (h2, c2) = bound(W, summary, 1, goal, dists), bound(W, summary, 2, goal, dists)
    return max.(dec, h1, h2), c1 + c2
end
shaped(summary) = function (W, Ws, goal, dec, dists)
    (h1, c1), (h2, c2) = bound(Ws, summary, 1, goal, dists), bound(Ws, summary, 2, goal, dists)
    return dec .+ max.(h1, h2), c1 + c2
end
const HEURISTICS = [
    "decoupled (sum of distances)" => (W, Ws, goal, dec, dists) -> (dec, 0),
    "costs, rooms 4×4 + corridors + near 2" => unshaped(rooms_near(4, 2)),
    "delays, near 2" => shaped(near(2)),
    "delays, near 6" => shaped(near(6)),
    "delays, rooms 4×4 + corridors" => shaped(rooms(4)),
    "delays, rooms 4×4 + corridors + near 2" => shaped(rooms_near(4, 2)),
    "delays, rooms 8×8 + corridors + near 2" => shaped(rooms_near(8, 2)),
    "delays, progress + corridors + near 2" => shaped(progress_near(2)),
    "delays, progress + corridors + near 4" => shaped(progress_near(4)),
]

# Queries: agents swapping sides near the corridors (they must coordinate), and random
# start/goal pairs (they mostly need not).
function queries(rng)
    left = [k for (k, (x, y)) in enumerate(CELLS) if n ÷ 2 - 9 <= x <= n ÷ 2 - 4 && abs(y - n ÷ 2) <= 4]
    right = [k for (k, (x, y)) in enumerate(CELLS) if n ÷ 2 + 5 <= x <= n ÷ 2 + 10 && abs(y - n ÷ 2) <= 4]
    swaps = map(1:NQ) do _
        a1, g2 = randperm(rng, length(left))[1:2]
        a2, g1 = randperm(rng, length(right))[1:2]
        ("swap", left[a1], right[a2], right[g1], left[g2])
    end
    randoms = map(1:NQ) do _
        a1, a2, g1, g2 = randperm(rng, NF)[1:4]
        ("random", a1, a2, g1, g2)
    end
    return [swaps; randoms]
end

println("grid $(n)², wall with two 1-cell corridors: $(NF) free cells, $(NF^2) joint states")
t_structure = @elapsed S = joint_structure()
@printf("joint graph: %d moves, built in %.1f s\n\n", nnz(S), t_structure)

rows = NamedTuple[]
rng = MersenneTwister(7)
for (kind, a1, a2, g1, g2) in queries(rng)
    d1, d2 = bfs(g1), bfs(g2)
    W, Ws = costs(S, g1, g2, d1, d2)
    goal = joint(g1, g2)
    start = joint(a1, a2)
    V = cost_to_go(W, [goal])
    dec = [d1[a] + d2[b] for b in 1:NF for a in 1:NF]
    out = sparse(transpose(W))
    exact = shortest_path(W, start, [goal]; heuristic = V, out_edges = out)
    @assert exact.cost == V[start]
    @printf("%-6s  optimal %3.0f  decoupled %3.0f  (coordination cost %.0f)\n", kind, V[start], dec[start], V[start] - dec[start])
    for (name, heuristic) in HEURISTICS
        build = @elapsed h, coarse = heuristic(W, Ws, goal, dec, (d1, d2))
        reachable = isfinite.(V)
        @assert all(h[reachable] .<= V[reachable] .+ 1e-9) "$name is not admissible"
        search = @elapsed r = shortest_path(W, start, [goal]; heuristic = h, out_edges = out)
        @assert r.cost == V[start]
        push!(rows, (kind = kind, heuristic = name, coarse_states = coarse, optimal = V[start],
            decoupled = dec[start], bound = h[start], expanded = r.expanded, exact_expanded = exact.expanded,
            build_seconds = build, search_seconds = search))
        @printf("    %-42s bound %5.1f  expanded %8d  build %5.2f s  search %6.3f s\n",
            name, h[start], r.expanded, build, search)
    end
end

println("\nsummary (means over queries)")
@printf("%-8s %-42s %13s %11s %12s %11s %9s %9s\n", "queries", "heuristic", "coarse states", "bound/opt",
    "gap closed", "expanded", "build s", "search s")
for kind in ("swap", "random"), (name, _) in HEURISTICS
    rs = filter(r -> r.kind == kind && r.heuristic == name, rows)
    gap = [(r.bound - r.decoupled) / (r.optimal - r.decoupled) for r in rs if r.optimal > r.decoupled]
    @printf("%-8s %-42s %13d %11.3f %12s %11.0f %9.2f %9.3f\n", kind, name, round(Int, mean(r.coarse_states for r in rs)),
        mean(r.bound / r.optimal for r in rs), isempty(gap) ? "–" : @sprintf("%.2f", mean(gap)),
        mean(r.expanded for r in rs), mean(r.build_seconds for r in rs), mean(r.search_seconds for r in rs))
end
for kind in ("swap", "random")
    rs = filter(r -> r.kind == kind && r.heuristic == first(HEURISTICS)[1], rows)
    @printf("%-8s %-42s %13s %11.3f %12s %11.0f\n", kind, "exact cost-to-go (reference)", "–", 1.0, "1.00",
        mean(r.exact_expanded for r in rs))
end

mkpath(OUT)
open(joinpath(OUT, "tropical_multiagent.csv"), "w") do io
    println(io, join(keys(first(rows)), ","))
    for r in rows
        println(io, join((occursin(',', string(v)) ? "\"$v\"" : v for v in values(r)), ","))
    end
end
