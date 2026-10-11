# Export the data of the layered-control demo (docs/viz) as JSON (machine-written).
#
#   julia --project=docs docs/scripts/tropical_viz_export.jl [n] [vmax]
#
# The minimum-time double integrator lattice of the studies: value maps at rest, the
# optimistic bounds of the velocity-blind and velocity-mode planners, a dozen queries
# with their optimal and certified layered trajectories (routes, corridors, values
# along the closed loop), the states each search expands, the one-axis phase
# portrait, and the study tables.

using CellularSheaves
using SparseArrays
using Random
using Base64
using Printf

function double_integrator_lattice(n, vmax, obstacle)
    nv = 2vmax + 1
    points = (n, n, nv, nv)
    L = LinearIndices(points)
    rows, cols = Int[], Int[]
    for I in CartesianIndices(points)
        qx, qy, ix, iy = Tuple(I)
        vx, vy = ix - vmax - 1, iy - vmax - 1
        obstacle(qx, qy) && continue
        for ax in -1:1, ay in -1:1
            wx, wy = vx + ax, vy + ay
            (abs(wx) <= vmax && abs(wy) <= vmax) || continue
            px, py = qx + wx, qy + wy
            (1 <= px <= n && 1 <= py <= n) || continue
            (obstacle(px, py) || obstacle(qx + wx ÷ 2, qy + wy ÷ 2)) && continue
            push!(rows, L[I]); push!(cols, L[px, py, wx + vmax + 1, wy + vmax + 1])
        end
    end
    N = prod(points)
    return sparse(rows, cols, ones(length(rows)), N, N), points
end

function heading_speed_modes(vmax, headings, speed_edges)
    M = zeros(Int, 2vmax + 1, 2vmax + 1)
    for ix in 1:(2vmax + 1), iy in 1:(2vmax + 1)
        vx, vy = ix - vmax - 1, iy - vmax - 1
        if vx == 0 && vy == 0
            M[ix, iy] = 1
            continue
        end
        band = findfirst(>=(max(abs(vx), abs(vy))), speed_edges)
        h = mod(round(Int, atan(vy, vx) / (2π / headings)), headings)
        M[ix, iy] = 1 + (band - 1) * headings + h + 1
    end
    return M
end

# ---- a tiny JSON writer (numbers, strings, bools, nothing, vectors, named tuples) ----
js(x::Nothing) = "null"
js(x::Bool) = x ? "true" : "false"
js(x::Integer) = string(x)
js(x::AbstractFloat) = isfinite(x) ? (isinteger(x) && abs(x) < 1e15 ? string(Int(x)) : @sprintf("%.6g", x)) : "null"
js(x::AbstractString) = "\"" * escape_string(x) * "\""
js(x::Symbol) = js(string(x))
js(x::AbstractVector) = "[" * join(js.(x), ",") * "]"
js(x::Tuple) = js(collect(x))
js(x::NamedTuple) = "{" * join(("\"$(k)\":" * js(v) for (k, v) in pairs(x)), ",") * "}"
js(x::AbstractDict) = "{" * join(("\"$(k)\":" * js(v) for (k, v) in x), ",") * "}"
bytes64(v) = base64encode(UInt8.(clamp.(v, 0, 255)))              # one byte per cell
bits64(b) = base64encode(UInt8[sum(UInt8(b[8k + j]) << (j - 1) for j in 1:8 if 8k + j <= length(b); init = 0x00)
                                for k in 0:(cld(length(b), 8) - 1)])

n = length(ARGS) >= 1 ? parse(Int, ARGS[1]) : 128
vmax = length(ARGS) >= 2 ? parse(Int, ARGS[2]) : 6
nv = 2vmax + 1
m = n ÷ 3
walls(qx, qy) = (qx in m:(m + 1) && !(qy in 2:(n ÷ 6))) || (qx in (2m):(2m + 1) && !(qy in (n - n ÷ 6):(n - 1)))
W, points = double_integrator_lattice(n, vmax, walls)
Wout = sparse(transpose(W))
L = LinearIndices(points)
C = CartesianIndices(points)
stop = vmax + 1
goal = [L[qx, qy, stop, stop] for qx in (n - 3):(n - 1) for qy in 2:4]
V = cost_to_go(W, goal)
state(x) = (C[x][1], C[x][2], C[x][3] - vmax - 1, C[x][4] - vmax - 1)

ψ = projection_homomorphism(points, (1, 2))
H = 16
modes = heading_speed_modes(vmax, H, [2, 4, vmax])
T = TruncatedPushforward(W, points, (1, 2), modes)
planners = (blind = LayeredPlanner(W, ψ, (n, n)), modes = LayeredPlanner(T))
bounds = (blind = lower_bound(planners.blind.abstraction, goal), modes = lower_bound(planners.modes.abstraction, goal))
at_rest(v) = [v[L[qx, qy, stop, stop]] for qy in 1:n for qx in 1:n]          # row-major: y outer, x inner
mode_info(mode) = mode == 1 ? (0.0, 0) : (2π * ((mode - 2) % H) / H, (mode - 2) ÷ H + 1)

rng = Random.MersenneTwister(4)
candidates = [x for x in 1:size(W, 1) if isfinite(V[x]) && V[x] > 0]
sources = candidates[randperm(rng, length(candidates))[1:30]]
# A dozen queries: prefer starts left of the first wall, which must find the gaps.
picked = sort([s for s in sources if state(s)[1] < m]; by = s -> -V[s])
picked = vcat(picked, [s for s in sources if !(s in picked)])[1:12]

function explored(f)
    count = zeros(Int, n * n)
    f(x -> (c = state(x); count[(c[2] - 1) * n + c[1]] += 1))
    return count
end

queries = []
for (qi, s) in enumerate(picked)
    nxt = greedy_policy(W, V; out_edges = Wout)
    opt = [s]
    while V[opt[end]] > 0
        push!(opt, nxt[opt[end]])
    end
    exp_d = explored(f -> shortest_path(W, s, goal; out_edges = Wout, on_expand = f))
    exp_b = explored(f -> shortest_path(W, s, goal; heuristic = bounds.blind, out_edges = Wout, on_expand = f))
    exp_m = explored(f -> shortest_path(W, s, goal; heuristic = bounds.modes, out_edges = Wout, on_expand = f))
    runs = []
    for (pname, P) in pairs(planners), (w, mw) in ((4, 1), (16, 1))
        c = certified_path(P, s, goal; width = w, mode_width = mw, lower = bounds[pname])
        U = corridor(P, c.route, c.width, c.mode_width)
        cells = falses(n * n)
        for cs in findall(U)
            q = CartesianIndices((n, n))[P.π.vertex_map[cs]]
            cells[(q[2] - 1) * n + q[1]] = true
        end
        route = map(c.route) do r
            q = CartesianIndices((n, n))[P.π.vertex_map[r]]
            mode = (r - 1) % P.K + 1
            ang, band = P.K == 1 ? (0.0, -1) : mode_info(mode)
            (q[1], q[2], ang, band)
        end
        push!(runs, (planner = string(pname), width = w, mode_width = mw, used_width = c.width, used_mode_width = c.mode_width,
            first_try = c.first_try, cost = c.cost, bound = c.suboptimality_bound, lower = c.lower,
            certified = c.certificate.positive_definite && c.certificate.strict && c.certificate.closed,
            route = route, corridor = bits64(cells), traj = state.(c.path),
            V_U = [c.cost - sum(W[c.path[i], c.path[i + 1]] for i in 1:(k - 1); init = 0.0) for k in eachindex(c.path)],
            V_opt = V[c.path], lower_along = bounds[pname][c.path]))
    end
    push!(queries, (id = qi, src = state(s), cost = V[s], optimal = state.(opt), V_opt = V[opt],
        lower_blind = bounds.blind[opt], lower_modes = bounds.modes[opt],
        explored = (dijkstra = bytes64(exp_d), blind = bytes64(exp_b), modes = bytes64(exp_m)),
        expanded = (dijkstra = sum(exp_d), blind = sum(exp_b), modes = sum(exp_m)), runs = runs))
    println("query $(qi): cost $(V[s]), runs ", [(r.planner, r.width, r.cost, round(r.bound; digits = 2)) for r in runs])
end

# ---- the one-axis phase portrait (classic minimum-time double integrator) ------------
nq, gq = 96, 48
Lp = LinearIndices((nq, nv))
pr, pc = Int[], Int[]
for q in 1:nq, iv in 1:nv, a in -1:1
    v2 = iv - vmax - 1 + a
    abs(v2) <= vmax || continue
    q2 = q + v2
    1 <= q2 <= nq || continue
    push!(pr, Lp[q, iv]); push!(pc, Lp[q2, v2 + vmax + 1])
end
W1 = sparse(pr, pc, ones(length(pr)), nq * nv, nq * nv)
V1 = cost_to_go(W1, [Lp[gq, stop]])
n1 = greedy_policy(W1, V1)
accel = [n1[Lp[q, iv]] == 0 ? 0 : (CartesianIndices((nq, nv))[n1[Lp[q, iv]]][2] - iv) for iv in 1:nv for q in 1:nq]
phase_starts = [(8, 0), (90, 0), (20, 6), (76, -6), (40, -5), (60, 4)]
phase_traj = map(phase_starts) do (q, v)
    x = Lp[q, v + vmax + 1]
    tr = [x]
    while isfinite(V1[tr[end]]) && V1[tr[end]] > 0
        push!(tr, n1[tr[end]])
    end
    [(CartesianIndices((nq, nv))[y][1], CartesianIndices((nq, nv))[y][2] - vmax - 1) for y in tr]
end

# ---- study tables --------------------------------------------------------------------
readcsv(path) = (lines = readlines(path); head = split(lines[1], ","); [Dict(zip(head, split(l, ","))) for l in lines[2:end]])
fig = joinpath(@__DIR__, "..", "figures")
layered = readcsv(joinpath(fig, "tropical_layered.csv"))
certified = readcsv(joinpath(fig, "tropical_certified.csv"))
num(s) = (x = tryparse(Float64, s); isnothing(x) ? s : x)

data = (n = n, vmax = vmax, headings = H, goal = [(state(g)[1], state(g)[2]) for g in goal],
    walls = bits64([walls(qx, qy) for qy in 1:n for qx in 1:n]),
    V_rest = bytes64([isfinite(v) ? Int(v) : 255 for v in at_rest(V)]),
    lower_rest = (blind = bytes64([isfinite(v) ? Int(v) : 255 for v in at_rest(bounds.blind)]),
                  modes = bytes64([isfinite(v) ? Int(v) : 255 for v in at_rest(bounds.modes)])),
    states = size(W, 1), edges = nnz(W), reachable = count(isfinite, V),
    queries = queries,
    phase = (nq = nq, goal = gq, value = [isfinite(v) ? Int(v) : -1 for v in [V1[Lp[q, iv]] for iv in 1:nv for q in 1:nq]],
             accel = accel, trajectories = phase_traj),
    layered = [Dict(k => num(v) for (k, v) in r) for r in layered],
    certified = [Dict(k => num(v) for (k, v) in r) for r in certified])

out = joinpath(@__DIR__, "..", "viz", "data.json")
mkpath(dirname(out))
write(out, js(data))
println("wrote $(out) ($(round(filesize(out) / 1e6; digits = 2)) MB)")
