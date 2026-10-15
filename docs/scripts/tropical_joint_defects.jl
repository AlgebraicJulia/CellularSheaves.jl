# Lifting defects of the joint maps of tropical_multiagent_study.jl: agent i exact,
# agent j summarized, on the *pair* graph (collisions removed), in the gauge of the
# decoupled cost-to-go φ = d₁ + d₂, so that arc costs are delays. This is the map
# whose quotient gives the delay bound in that study, so these are the defects the
# path-lifting theorem is about. (On agent j's own graph the same summaries nearly
# all lift; their failures come from the coupling, near agent i.)
#
# In the delay gauge moves toward the goal are free, so a reverse search from a block
# can flood much of the pair graph. We therefore estimate the shares from random
# samples of coarse arcs (all arcs, and zero-delay arcs separately), each defect
# computed exactly up to the horizon (beyond it, counted as a leak).
#
#   julia --project=docs docs/scripts/tropical_joint_defects.jl [n] [queries] [horizon] [samples]

using CellularSheaves
using CellularSheaves.NetworkSheaves.TropicalAbstraction: _DefectSearch, _defect_search!
using SparseArrays
using Printf
using Random
using Statistics

const n = length(ARGS) >= 1 ? parse(Int, ARGS[1]) : 32
const NQ = length(ARGS) >= 2 ? parse(Int, ARGS[2]) : 4
const HORIZON = length(ARGS) >= 3 ? parse(Float64, ARGS[3]) : 4.0
const SAMPLES = length(ARGS) >= 4 ? parse(Int, ARGS[4]) : 2000
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


# the summaries compared in the study (all twisted by φ)
const SUMMARIES = [position_blocks(4; near = 2), progress(near = 2), progress(near = 2, corridors = true),
                   components(), components(near = 2)]

function gauged(W, φ)
    r, c, v = findnz(W)
    keep = [isfinite(φ[x]) && isfinite(φ[y]) for (x, y) in zip(r, c)]
    return sparse(r[keep], c[keep], max.(v[keep] .+ φ[c[keep]] .- φ[r[keep]], 0.0), size(W)...)
end

println("grid $(n)², wall with two 1-cell corridors: $(NF) free cells, $(NF^2) joint states; horizon $(HORIZON); $(SAMPLES) sampled arcs per kind")
S = joint_structure()
acc = Dict{String,Vector{NTuple{3,Float64}}}()
for (kind, a1, a2, g1, g2) in filter(q -> q[1] == "swap", queries(MersenneTwister(7)))
    W = joint_costs(S, g1, g2)
    d1, d2 = distances(g1), distances(g2)
    φ = [d1[a] + d2[b] for b in 1:NF for a in 1:NF]
    Wg = gauged(W, φ)
    info = ((d = d1, comp = level_components(d1)), (d = d2, comp = level_components(d2)))
    for s in SUMMARIES
        μ = mode_map(s, 1, info)
        I, J, c = findnz(quotient_graph(Wg, μ))
        members = [Int[] for _ in 1:μ.n_target]
        for (x, B) in enumerate(μ.vertex_map)
            push!(members[B], x)
        end
        rng = MersenneTwister(11)
        all_arcs = rand(rng, eachindex(I), min(SAMPLES, length(I)))
        zeros_ = findall(iszero, c)
        zero_arcs = isempty(zeros_) ? Int[] : rand(rng, zeros_, min(SAMPLES, length(zeros_)))
        search = _DefectSearch(size(Wg, 1))
        defect(k) = (_defect_search!(search, Wg, members[J[k]], x -> true, c[k] + HORIZON);
                     max(0.0, maximum(search.dist[x] for x in members[I[k]]) - c[k]))
        lift = mean(iszero(defect(k)) for k in all_arcs)
        ez = [defect(k) for k in zero_arcs]
        leak = isempty(ez) ? 0.0 : count(>(0), ez) / length(ez)
        mleak = all(iszero, ez) ? 0.0 : mean(min.(filter(>(0), ez), HORIZON + 1))
        push!(get!(acc, s.name, NTuple{3,Float64}[]), (lift, leak, mleak))
        @printf("  %-48s lift %5.1f%%   zero-delay arcs that leak %5.1f%%\n", s.name, 100lift, 100leak)
    end
end
println("\nmeans over swap queries")
println("summary of the other agent                         arcs that lift   zero-delay arcs that leak   mean leak (capped at horizon+1)")
rows = []
for s in SUMMARIES
    a = acc[s.name]
    m = [mean(x[i] for x in a) for i in 1:3]
    @printf("%-50s %8.1f%%         %8.1f%%              %5.2f\n", s.name, 100m[1], 100m[2], m[3])
    push!(rows, (s.name, m...))
end
open(joinpath(@__DIR__, "..", "figures", "tropical_joint_defects_n$(n).csv"), "w") do io
    println(io, "summary,lift_fraction,zero_delay_leak_fraction,mean_leak_capped")
    for r in rows
        println(io, join(r, ","))
    end
end
