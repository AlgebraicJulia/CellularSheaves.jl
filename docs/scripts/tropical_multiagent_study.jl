# Cooperative two-agent planning: how much of the coordination cost (waiting, detours)
# can a per-agent abstraction see? Agent i's coarse problem keeps its own cell and a
# summary of the other agent; the joint graph is pushed forward along that map with the
# refinement twisted by the decoupled cost-to-go φ = d₁ + d₂ (each agent's single-agent
# distance to its goal), so the coarse problem bounds only the coordination delay:
# `Abstraction(W, μ; potential = φ)`. The ablation separates the three ingredients:
#
#   - the twist (potential): without it, coarsening the other agent discards its travel cost;
#   - the summary of the other agent: its remaining distance to its goal ("progress"),
#     on which its delay-free moves are constant (level k → k - 1), against its position.
#     A level set alone lets the other agent hop between equidistant cells on opposite
#     sides of a wall; splitting each level into the arcs of the Reeb graph of d (the
#     components of the band k ≤ d ≤ k + 1) keeps it from doing so;
#   - an exact near field: the other agent's offset when within distance r, the support
#     of the collision coupling.
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

# A thick wall with two one-cell corridors: agents that swap sides must queue through a
# corridor or detour to the other one.
const WALL = (n ÷ 2 - 2):(n ÷ 2 + 3)
wall(x, y) = x in WALL && !(y in (n ÷ 2, n - 4))
const CELLS = [(x, y) for y in 1:n for x in 1:n if !wall(x, y)]
const INDEX = Dict(c => k for (k, c) in enumerate(CELLS))
const NF = length(CELLS)
const NBRS = [[k; [INDEX[(x + dx, y + dy)] for (dx, dy) in ((1, 0), (-1, 0), (0, 1), (0, -1))
                   if haskey(INDEX, (x + dx, y + dy))]] for (k, (x, y)) in enumerate(CELLS)]

# Joint states (a, b), a ≠ b, indexed a + NF (b - 1). A joint move moves each agent one
# step or keeps it, with no vertex conflict and no swap conflict.
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
function joint_costs(S, g1, g2)
    v = similar(nonzeros(S))
    rv = rowvals(S)
    for y in 1:size(S, 2), p in nzrange(S, y)
        (a, b), (a2, b2) = agents(rv[p]), agents(y)
        v[p] = (a == a2 == g1 ? 0.0 : 1.0) + (b == b2 == g2 ? 0.0 : 1.0)
    end
    return SparseMatrixCSC(size(S)..., S.colptr, rv, v)
end

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

# Summaries of the other agent seen from `own`: a mode in 1:K, given the other agent's
# distances to its goal.
cheb(o) = max(abs(o[1]), abs(o[2]))
near_offset(own, other) = CELLS[other] .- CELLS[own]
near_index(o, r) = (o[1] + r) + (2r + 1) * (o[2] + r) + 1
block_index(c, B) = (CELLS[c][1] - 1) ÷ cld(n, B) + B * ((CELLS[c][2] - 1) ÷ cld(n, B)) + 1 +
                    (CELLS[c][1] > n ÷ 2 ? B^2 : 0)              # blocks do not straddle the wall
const CORRIDOR = Dict(k => i for (i, k) in enumerate(k for (k, (x, y)) in enumerate(CELLS) if x in WALL))
const LEVELS = 4n                                                # more than any distance

# The single-agent graph (moves and waits on the free cells), for the library's Reeb arcs.
const SOLO = sparse([k for k in 1:NF for _ in NBRS[k]], reduce(vcat, NBRS), ones(sum(length, NBRS)), NF, NF)
level_components(d) = reeb_arcs(SOLO, d).vertex_map

struct Summary{F}
    name::String
    K::Int
    mode::F                 # (own, other, (d, comp) of the other) -> 1:K
end

# A far summary, optionally with an exact near field of radius r and exact corridor cells.
function summary(name, K, far; near = 0, corridors = false)
    nn = near > 0 ? (2near + 1)^2 : 0
    nc = corridors ? length(CORRIDOR) : 0
    label = name * (near > 0 ? " + near $near" : "") * (corridors ? " + corridors" : "")
    return Summary(label, nn + nc + K, function (own, other, info)
        if near > 0
            o = near_offset(own, other)
            cheb(o) <= near && return near_index(o, near)
        end
        corridors && haskey(CORRIDOR, other) && return nn + CORRIDOR[other]
        return nn + nc + far(other, info)
    end)
end
progress(; kw...) = summary("progress levels", LEVELS + 1, (other, info) -> Int(info.d[other]) + 1; kw...)
components(; kw...) = summary("progress (Reeb arcs)", NF, (other, info) -> info.comp[other]; kw...)
position_blocks(B; kw...) = summary("position blocks $(B)×$(B)", 2B^2, (other, info) -> block_index(other, B); kw...)

# Agent i's mode map on the joint states: (own cell, summary of the other).
function mode_map(s::Summary, i, dists)
    d_other = dists[3 - i]
    vmap = Vector{Int}(undef, NF^2)
    for x in 1:NF^2
        a, b = agents(x)
        own, other = i == 1 ? (a, b) : (b, a)
        vmap[x] = own + NF * (s.mode(own, other, d_other) - 1)
    end
    return GraphHomomorphism(vmap, NF * s.K)
end

occupied(μ) = length(unique(μ.vertex_map[x] for x in 1:NF^2 if !=(agents(x)...)))

# The heuristic of a summary: the max of the two agents' bounds (both bound the same
# joint cost) and the decoupled sum, which a twisted bound already dominates.
function heuristic(W, s, goal, dists, φ; twisted = true)
    h, coarse = copy(φ), 0
    for i in 1:2
        μ = mode_map(s, i, dists)
        A = twisted ? Abstraction(W, μ; potential = φ) : Abstraction(W, μ)
        h .= max.(h, lower_bound(A, [goal]))
        coarse += occupied(μ)
    end
    return h, coarse
end

const VARIANTS = [
    ("decoupled (sum of distances)", nothing, true),
    ("untwisted: " * components(near = 2).name, components(near = 2), false),
    ("twisted: " * position_blocks(4; near = 2).name, position_blocks(4; near = 2), true),
    ("twisted: " * progress(near = 2).name, progress(near = 2), true),
    ("twisted: " * progress(near = 2, corridors = true).name, progress(near = 2, corridors = true), true),
    ("twisted: " * components().name, components(), true),
    ("twisted: " * components(near = 2).name, components(near = 2), true),
]

function queries(rng)
    side(lo, hi) = [k for (k, (x, y)) in enumerate(CELLS) if lo <= x <= hi && abs(y - n ÷ 2) <= 4]
    left, right = side(n ÷ 2 - 9, n ÷ 2 - 4), side(n ÷ 2 + 5, n ÷ 2 + 10)
    swaps = map(1:NQ) do _
        a1, g2 = left[randperm(rng, length(left))[1:2]]
        a2, g1 = right[randperm(rng, length(right))[1:2]]
        ("swap", a1, a2, g1, g2)
    end
    randoms = [("random", randperm(rng, NF)[1:4]...) for _ in 1:NQ]
    return [swaps; randoms]
end

println("grid $(n)², wall with two 1-cell corridors: $(NF) free cells, $(NF^2) joint states")
t_structure = @elapsed S = joint_structure()
@printf("joint graph: %d moves, built in %.1f s\n\n", nnz(S), t_structure)

rows = NamedTuple[]
for (kind, a1, a2, g1, g2) in queries(MersenneTwister(7))
    W = joint_costs(S, g1, g2)
    d1, d2 = distances(g1), distances(g2)
    φ = [d1[a] + d2[b] for b in 1:NF for a in 1:NF]
    goal, start = joint(g1, g2), joint(a1, a2)
    V = cost_to_go(W, [goal])
    reachable = isfinite.(V)
    out = sparse(transpose(W))
    exact = shortest_path(W, start, [goal]; heuristic = V, out_edges = out)
    @printf("%-6s  optimal %3.0f  decoupled %3.0f  (coordination cost %.0f)\n", kind, V[start], φ[start], V[start] - φ[start])
    for (name, s, twisted) in VARIANTS
        info = ((d = d1, comp = level_components(d1)), (d = d2, comp = level_components(d2)))
        build = @elapsed h, coarse = s === nothing ? (φ, 0) : heuristic(W, s, goal, info, φ; twisted)
        @assert all(h[reachable] .<= V[reachable] .+ 1e-9) "$name is not admissible"
        search = @elapsed r = shortest_path(W, start, [goal]; heuristic = h, out_edges = out)
        @assert r.cost == V[start]
        push!(rows, (kind = kind, heuristic = name, coarse_states = coarse, optimal = V[start], decoupled = φ[start],
            bound = h[start], expanded = r.expanded, exact_expanded = exact.expanded,
            build_seconds = build, search_seconds = search))
        @printf("    %-52s bound %5.1f  expanded %6d  build %5.2f s\n", name, h[start], r.expanded, build)
    end
end

println("\nsummary (means over queries; coarse states: both agents, occupied)")
@printf("%-7s %-52s %13s %10s %11s %9s\n", "queries", "heuristic", "coarse states", "bound/opt", "gap closed", "expanded")
for kind in ("swap", "random")
    for (name, _, _) in VARIANTS
        rs = filter(r -> r.kind == kind && r.heuristic == name, rows)
        gap = [(r.bound - r.decoupled) / (r.optimal - r.decoupled) for r in rs if r.optimal > r.decoupled]
        @printf("%-7s %-52s %13d %10.3f %11s %9.0f\n", kind, name, round(Int, mean(r.coarse_states for r in rs)),
            mean(r.bound / r.optimal for r in rs), isempty(gap) ? "–" : @sprintf("%.2f", mean(gap)),
            mean(r.expanded for r in rs))
    end
    rs = filter(r -> r.kind == kind && r.heuristic == first(VARIANTS)[1], rows)
    @printf("%-7s %-52s %13d %10.3f %11s %9.0f\n", kind, "exact cost-to-go (reference)", NF^2 - NF, 1.0, "1.00",
        mean(r.exact_expanded for r in rs))
end

mkpath(OUT)
open(joinpath(OUT, "tropical_multiagent.csv"), "w") do io
    println(io, join(keys(first(rows)), ","))
    for r in rows
        println(io, join((occursin(',', string(v)) ? "\"$v\"" : v for v in values(r)), ","))
    end
end
