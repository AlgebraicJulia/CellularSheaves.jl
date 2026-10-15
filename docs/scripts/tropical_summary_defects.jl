# Lifting defects of the summaries of another agent used by the multi-agent studies
# (tropical_multiagent_study.jl, tropical_cover_triangles.jl).
#
# A summary π_j abstracts agent j's single-agent graph. In the gauge of j's distance
# d_j to its goal, a move toward the goal costs 0, waiting costs 1 and a step away
# costs 2: arc costs are j's *delays*. A coarse arc of the summary that does not
# lift lets a summarized agent make, for free, a move its actual state cannot: the
# leaks the studies describe (slipping through a corridor early, hopping across a
# wall, sliding along a level set onto another route). We report, per summary:
#   - coarse states;
#   - the share of coarse arcs that lift (zero defect);
#   - among the *zero-delay* coarse arcs (the ones a delay bound uses for free), the
#     share with positive defect and their mean defect.
#
#   julia --project=docs docs/scripts/tropical_summary_defects.jl [n] [goals]

using CellularSheaves
using SparseArrays
using Printf
using Random
using Statistics

n = length(ARGS) >= 1 ? parse(Int, ARGS[1]) : 32
ngoals = length(ARGS) >= 2 ? parse(Int, ARGS[2]) : 12

# The map of tropical_multiagent_study.jl: a thick wall with two one-cell corridors.
WALL = (n ÷ 2 - 2):(n ÷ 2 + 3)
wall(x, y) = x in WALL && !(y in (n ÷ 2, n - 4))
CELLS = [(x, y) for y in 1:n for x in 1:n if !wall(x, y)]
INDEX = Dict(c => k for (k, c) in enumerate(CELLS))
NF = length(CELLS)
NBRS = [[k; [INDEX[(x + dx, y + dy)] for (dx, dy) in ((1, 0), (-1, 0), (0, 1), (0, -1))
             if haskey(INDEX, (x + dx, y + dy))]] for (k, (x, y)) in enumerate(CELLS)]
SOLO = sparse([k for k in 1:NF for _ in NBRS[k]], reduce(vcat, NBRS), ones(sum(length, NBRS)), NF, NF)

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

# Agent j's delays: arc u → v costs 1 + d(v) − d(u), except waiting at the goal (0).
function delay_graph(d, g)
    r, c, v = findnz(SOLO)
    w = [u == x == g ? 0.0 : 1.0 + d[x] - d[u] for (u, x) in zip(r, c)]
    return sparse(r, c, w, NF, NF)
end

compact(labels) = (ids = Dict{Any,Int}(); GraphHomomorphism([get!(ids, l, length(ids) + 1) for l in labels], length(ids)))
block(k, B) = ((CELLS[k][1] - 1) ÷ B, (CELLS[k][2] - 1) ÷ B, CELLS[k][1] > n ÷ 2)   # blocks do not straddle the wall

function summaries(d, Wd)
    reeb = reeb_arcs(SOLO, d)
    arcblock = compact([(reeb.vertex_map[k], block(k, 2)) for k in 1:NF])
    zero = (s, s2, c) -> c == 0
    return [("position blocks 4×4", compact([block(k, 4) for k in 1:NF])),
            ("level sets of d_j", compact([d[k] for k in 1:NF])),
            ("Reeb arcs of d_j", compact(reeb.vertex_map)),
            ("(Reeb arc, block 2)", arcblock),
            ("(Reeb arc, block 2), future-refined", future_refinement(Wd, arcblock, reeb; arcs = zero)),
            ("(Reeb arc, block 2), stable", stable_refinement(Wd, arcblock; arcs = zero))]
end

rng = MersenneTwister(7)
goals = randperm(rng, NF)[1:ngoals]
names = String[]
acc = Dict{String,Vector{NTuple{4,Float64}}}()
for g in goals
    d = distances(g)
    Wd = delay_graph(d, g)
    for (name, π) in summaries(d, Wd)
        ε = lifting_defect(Wd, π; scope = :global)
        Wq = quotient_graph(Wd, π)
        I, J, c = findnz(Wq)
        e = nonzeros(ε)
        z = findall(iszero, c)                                     # zero-delay coarse arcs
        leak = isempty(z) ? 0.0 : count(>(0), e[z]) / length(z)
        mean_leak = isempty(z) || all(iszero, e[z]) ? 0.0 : mean(filter(>(0), e[z]))
        name in names || push!(names, name)
        push!(get!(acc, name, NTuple{4,Float64}[]), (π.n_target, count(iszero, e) / length(e), leak, mean_leak))
    end
end

@printf("grid %d², wall with two 1-cell corridors, %d free cells; means over %d goals\n\n", n, NF, ngoals)
println("summary of agent j                        coarse states   arcs that lift   zero-delay arcs that leak   mean leak")
rows = []
for name in names
    a = acc[name]
    m = [mean(x[i] for x in a) for i in 1:4]
    @printf("%-40s %10.0f      %8.1f%%         %8.1f%%              %5.2f\n", name, m[1], 100m[2], 100m[3], m[4])
    push!(rows, (name, m...))
end
open(joinpath(@__DIR__, "..", "figures", "tropical_summary_defects_n$(n).csv"), "w") do io
    println(io, "summary,coarse_states,lift_fraction,zero_delay_leak_fraction,mean_leak")
    for r in rows
        println(io, join(r, ","))
    end
end
