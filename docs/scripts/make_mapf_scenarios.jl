# MAPF scenarios from MovingAI single-agent scenarios, in the format of the MAPF
# benchmark (Stern et al. 2019), for maps the MAPF suite does not cover (the
# Warcraft III set, https://movingai.com/benchmarks/wc3maps512/).
#
# The single-agent file lists start/goal pairs, bucketed by octile path length. A
# MAPF scenario is an ordered list of agents; an N-agent instance takes the first N
# entries. Each scenario k = 1:count is built as follows:
#   - keep the pairs whose start and goal both lie in the largest 4-connected
#     component of free cells ('.', 'G', 'S'), since MAPF agents move 4-connected;
#   - shuffle them with seed k;
#   - take them in that order, skipping a pair whose start repeats an earlier
#     start or whose goal repeats an earlier goal (the MAPF suite's distinctness).
# Lines keep the source's octile distance and use its bucket rule, ⌊distance / 4⌋,
# as the MAPF files do. Coordinates stay 0-based, column then row.
#
#   julia docs/scripts/make_mapf_scenarios.jl MAP_FILE SCEN_FILE OUT_DIR [count]
# writes OUT_DIR/<map>-random-<k>.scen, the names tropical_multiagent_cover.jl reads.

using Random
using Printf

function largest_component(path)
    lines = readlines(path)
    h, w = parse(Int, split(lines[2])[2]), parse(Int, split(lines[3])[2])
    free = Set((x, y) for y in 0:(h - 1) for x in 0:(w - 1) if lines[5 + y][x + 1] in ('.', 'G', 'S'))
    best, seen = Set{Tuple{Int,Int}}(), Set{Tuple{Int,Int}}()
    for c in free
        c in seen && continue
        comp, stack = Set([c]), [c]
        push!(seen, c)
        while !isempty(stack)
            (x, y) = pop!(stack)
            for (dx, dy) in ((1, 0), (-1, 0), (0, 1), (0, -1))
                c2 = (x + dx, y + dy)
                (c2 in free && !(c2 in seen)) && (push!(seen, c2); push!(comp, c2); push!(stack, c2))
            end
        end
        length(comp) > length(best) && (best = comp)
    end
    return best, w, h
end

function main(mapfile, scenfile, outdir, count)
    comp, w, h = largest_component(mapfile)
    mapname = basename(mapfile)
    pairs = Tuple{NTuple{4,Int},Float64}[]
    for l in readlines(scenfile)[2:end]
        f = split(l)
        length(f) >= 9 || continue
        sx, sy, gx, gy = parse.(Int, f[5:8])
        ((sx, sy) in comp && (gx, gy) in comp && (sx, sy) != (gx, gy)) || continue
        push!(pairs, ((sx, sy, gx, gy), parse(Float64, f[9])))
    end
    mkpath(outdir)
    base = splitext(mapname)[1]
    sizes = Int[]
    for k in 1:count
        starts, goals = Set{Tuple{Int,Int}}(), Set{Tuple{Int,Int}}()
        open(joinpath(outdir, "$(base)-random-$(k).scen"), "w") do io
            println(io, "version 1")
            for i in randperm(MersenneTwister(k), length(pairs))
                (sx, sy, gx, gy), d = pairs[i]
                ((sx, sy) in starts || (gx, gy) in goals) && continue
                push!(starts, (sx, sy)); push!(goals, (gx, gy))
                @printf(io, "%d\t%s\t%d\t%d\t%d\t%d\t%d\t%d\t%.8f\n", floor(Int, d / 4), mapname, w, h, sx, sy, gx, gy, d)
            end
        end
        push!(sizes, length(starts))
    end
    @printf("%s: %d free cells in the largest 4-connected component; %d of the single-agent pairs usable; %d scenarios of %d–%d agents in %s\n",
        mapname, length(comp), length(pairs), count, minimum(sizes), maximum(sizes), outdir)
end

length(ARGS) >= 3 || error("usage: make_mapf_scenarios.jl MAP_FILE SCEN_FILE OUT_DIR [count]")
main(ARGS[1], ARGS[2], ARGS[3], length(ARGS) >= 4 ? parse(Int, ARGS[4]) : 25)
