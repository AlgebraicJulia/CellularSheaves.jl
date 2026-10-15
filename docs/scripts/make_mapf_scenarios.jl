# MAPF scenarios for maps outside the MAPF benchmark suite, made the way the suite's
# "random" scenarios were (Stern et al. 2019, §4.2–4.3, https://movingai.com/benchmarks/mapf.html):
# "All points in the largest reachable region of each map were randomly paired, and
# then the first 1000 problems were put into the scenario", 25 scenarios per map.
#
# Read against the published files, that is: shuffle the cells of the largest
# 4-connected region of free cells ('.', 'G', 'S'), pair them consecutively into
# (start, goal), and keep the first min(1000, ⌊cells/2⌋) pairs. No cell occurs twice in
# a scenario. (Every suite file checked has ⌊cells/2⌋ agents, capped at 1000, and no
# repeated cell: room-32-32-4 has 682 cells and 341 agents.)
#
# Lines use the suite's format: bucket, map, width, height, start x, start y, goal x,
# goal y (0-based, column then row), and the octile distance (diagonal moves cost √2
# and may not cut corners), with bucket ⌊distance / 4⌋.
#
#   julia -t 8 docs/scripts/make_mapf_scenarios.jl MAP_FILE OUT_DIR [count] [agents]
# writes OUT_DIR/<map>-random-<k>.scen, k = 1:count, the names
# tropical_multiagent_cover.jl and tropical_cover_lazy.jl read.
#   julia docs/scripts/make_mapf_scenarios.jl --check MAP_FILE SCEN_FILE
# recomputes the distances of a published scenario file and compares them.

using Random
using Printf

function read_free(path)
    lines = readlines(path)
    h, w = parse(Int, split(lines[2])[2]), parse(Int, split(lines[3])[2])
    free = falses(w, h)                                              # free[x + 1, y + 1]
    for y in 0:(h - 1), x in 0:(w - 1)
        free[x + 1, y + 1] = lines[5 + y][x + 1] in ('.', 'G', 'S')
    end
    return free, w, h
end

function largest_region(free)
    w, h = size(free)
    label = zeros(Int, w, h)
    best, bestsize, next = 0, 0, 0
    for x in 1:w, y in 1:h
        (free[x, y] && label[x, y] == 0) || continue
        next += 1
        label[x, y] = next
        stack, count = [(x, y)], 0
        while !isempty(stack)
            (a, b) = pop!(stack)
            count += 1
            for (dx, dy) in ((1, 0), (-1, 0), (0, 1), (0, -1))
                a2, b2 = a + dx, b + dy
                (1 <= a2 <= w && 1 <= b2 <= h && free[a2, b2] && label[a2, b2] == 0) || continue
                label[a2, b2] = next
                push!(stack, (a2, b2))
            end
        end
        count > bestsize && ((best, bestsize) = (next, count))
    end
    return [(x - 1, y - 1) for y in 1:h for x in 1:w if label[x, y] == best]   # 0-based, row-major
end

# Octile distance from (sx, sy) to (gx, gy), 0-based: Dijkstra on the 8-neighbour grid,
# a diagonal step allowed only when both orthogonal cells it passes are free.
function octile(free, sx, sy, gx, gy, dist)
    w, h = size(free)
    fill!(dist, Inf)
    dist[sx + 1, sy + 1] = 0.0
    heap = [(0.0, sx + 1, sy + 1)]
    while !isempty(heap)
        d, x, y = popmin!(heap)
        d > dist[x, y] && continue
        (x, y) == (gx + 1, gy + 1) && return d
        for dx in -1:1, dy in -1:1
            (dx, dy) == (0, 0) && continue
            x2, y2 = x + dx, y + dy
            (1 <= x2 <= w && 1 <= y2 <= h && free[x2, y2]) || continue
            dx != 0 && dy != 0 && !(free[x + dx, y] && free[x, y + dy]) && continue
            d2 = d + (dx != 0 && dy != 0 ? sqrt(2.0) : 1.0)
            if d2 < dist[x2, y2]
                dist[x2, y2] = d2
                pushheap!(heap, (d2, x2, y2))
            end
        end
    end
    return Inf
end

function pushheap!(h, x)
    push!(h, x)
    i = length(h)
    while i > 1 && h[i >> 1] > h[i]
        h[i >> 1], h[i] = h[i], h[i >> 1]
        i >>= 1
    end
end

function popmin!(h)
    top = h[1]
    last = pop!(h)
    if !isempty(h)
        h[1] = last
        i, n = 1, length(h)
        while true
            l, r, m = 2i, 2i + 1, i
            l <= n && h[l] < h[m] && (m = l)
            r <= n && h[r] < h[m] && (m = r)
            m == i && break
            h[i], h[m] = h[m], h[i]
            i = m
        end
    end
    return top
end

function generate(mapfile, outdir, count, agents)
    free, w, h = read_free(mapfile)
    region = largest_region(free)
    mapname = basename(mapfile)
    base = splitext(mapname)[1]
    n = min(agents, length(region) ÷ 2)
    mkpath(outdir)
    Threads.@threads for k in 1:count
        cells = shuffle(MersenneTwister(k), region)
        dist = fill(Inf, w, h)
        lines = map(1:n) do i
            (sx, sy), (gx, gy) = cells[2i - 1], cells[2i]
            d = octile(free, sx, sy, gx, gy, dist)
            @sprintf("%d\t%s\t%d\t%d\t%d\t%d\t%d\t%d\t%.8f", floor(Int, d / 4), mapname, w, h, sx, sy, gx, gy, d)
        end
        open(joinpath(outdir, "$(base)-random-$(k).scen"), "w") do io
            println(io, "version 1")
            foreach(l -> println(io, l), lines)
        end
    end
    @printf("%s: largest region %d cells; %d scenarios of %d agents in %s\n", mapname, length(region), count, n, outdir)
end

function check(mapfile, scenfile)
    free, w, h = read_free(mapfile)
    dist = fill(Inf, w, h)
    worst, rows = 0.0, 0
    for l in readlines(scenfile)[2:end]
        f = split(l)
        sx, sy, gx, gy = parse.(Int, f[5:8])
        worst = max(worst, abs(octile(free, sx, sy, gx, gy, dist) - parse(Float64, f[9])))
        rows += 1
    end
    @printf("%s: %d rows, largest distance difference %.2e\n", basename(scenfile), rows, worst)
    worst < 1e-6 || error("distances differ from the published file")
end

if length(ARGS) >= 3 && ARGS[1] == "--check"
    check(ARGS[2], ARGS[3])
elseif length(ARGS) >= 2
    generate(ARGS[1], ARGS[2], length(ARGS) >= 3 ? parse(Int, ARGS[3]) : 25, length(ARGS) >= 4 ? parse(Int, ARGS[4]) : 1000)
else
    error("usage: make_mapf_scenarios.jl MAP_FILE OUT_DIR [count] [agents] | --check MAP_FILE SCEN_FILE")
end
