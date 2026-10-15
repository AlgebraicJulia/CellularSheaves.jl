# The lifting defect of velocity-mode abstractions of the double integrator lattice,
# and the two certificates it gives (docs/papers/path-lifting, Theorems 4.2 and 5.2).
#
# The study works on the *viable* subgraph: the states that can reach the goal.
# Dead states (doomed to leave the grid or hit a wall) would otherwise give every arc
# out of their block an infinite defect, which measures dead ends, not lifting.
#
# For each truncation of the pushforward onto positions (velocity modes per cell):
#   - the defect of every arc of the optimistic quotient (lifting_defect, global
#     scope with a horizon: defects beyond it count as ∞, which keeps the bounds sound);
#   - the lower bound V_ψ ∘ ψ and the defect-inflated upper bound V^{+ε}_ψ ∘ ψ against
#     the exact cost-to-go V;
#   - the pessimistic (∀∃) upper bound for comparison;
#   - the macro-step closed loop of Theorem 5.2(b), run from sampled states: each step
#     follows one arc of an optimal walk of the inflated quotient by a fine path; we
#     check that the inflated value decreases by at least the path's cost and that the
#     loop reaches the goal within its certified cost.
#
#   julia -t 32 --project=docs docs/scripts/tropical_defect_study.jl [n] [vmax] [horizon] [samples]

using CellularSheaves
using CellularSheaves.NetworkSheaves.TropicalAbstraction: _DefectSearch, _defect_search!
using SparseArrays
using Printf
using Statistics
using Random

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

# Arc list of a sparse matrix with an extra sink: the coarse graph with arc costs
# Wq + ε and terminal arcs B → ⋆ of cost bq + εterm (stored zeros kept).
function with_sink(Wq, ε, bq, εterm)
    nc = size(Wq, 1)
    I, J, V = findnz(Wq)
    E = nonzeros(ε)
    t = findall(isfinite, bq)
    return sparse([I; t], [J; fill(nc + 1, length(t))], [V .+ E; bq[t] .+ εterm[t]], nc + 1, nc + 1, min)
end

n = length(ARGS) >= 1 ? parse(Int, ARGS[1]) : 64
vmax = length(ARGS) >= 2 ? parse(Int, ARGS[2]) : 6
horizon = length(ARGS) >= 3 ? parse(Float64, ARGS[3]) : 2.0 * vmax + 2
samples = length(ARGS) >= 4 ? parse(Int, ARGS[4]) : 300

m = n ÷ 3
walls(qx, qy) = (qx in m:(m + 1) && !(qy in 2:(n ÷ 6))) || (qx in (2m):(2m + 1) && !(qy in (n - n ÷ 6):(n - 1)))
Wfull, points = double_integrator_lattice(n, vmax, walls)
L = LinearIndices(points)
stop = vmax + 1
goal_full = [L[qx, qy, stop, stop] for qx in (n - 3):(n - 1) for qy in 2:4]
Vfull = cost_to_go(Wfull, goal_full)
keep = findall(isfinite, Vfull)                                  # the viable states
newindex = zeros(Int, size(Wfull, 1))
newindex[keep] = 1:length(keep)
W = Wfull[keep, keep]
Wt = sparse(transpose(W))
goal = newindex[goal_full]
V = Vfull[keep]
reach = findall(>(0), V)
@printf("lattice %d² × %d², %d states, %d viable (reach the goal); horizon %.0f; %d threads\n",
    n, 2vmax + 1, size(Wfull, 1), length(keep), horizon, Threads.nthreads())
flush(stdout)

# An abstraction of the full lattice, restricted to the viable states (empty blocks dropped).
function restrict_map(μ, keep)
    labels = μ.vertex_map[keep]
    ids = Dict{Int,Int}()
    compact = [get!(ids, l, length(ids) + 1) for l in labels]
    return GraphHomomorphism(compact, length(ids))
end

nv = 2vmax + 1
truncations = [("1 mode (velocity-blind)", ones(Int, nv, nv)),
               ("8 headings × 3 speeds", heading_speed_modes(vmax, 8, [2, 4, vmax])),
               ("16 headings × 3 speeds", heading_speed_modes(vmax, 16, [2, 4, vmax])),
               ("8 headings × $(vmax) speeds", heading_speed_modes(vmax, 8, collect(1:vmax)))]

rng = MersenneTwister(2026)
rows = []
println("\ntruncation                    modes   arcs   lift  ε>h(∞)  mean ε   L/V    certified: inflated / must   U/V: inflated / must   corr(U−L, V−L)  loop ok   build")
for (name, modes) in truncations
    T = TruncatedPushforward(Wfull, points, (1, 2), modes)
    μ = restrict_map(T.μ, keep)
    nc = μ.n_target
    t = @elapsed begin
        Wq = SparseMatrixCSC{Float64,Int}(quotient_graph(W, μ))
        ε = lifting_defect(W, μ; scope = :global, horizon = horizon)
    end
    e = nonzeros(ε)
    lift = count(iszero, e) / length(e)
    over = count(isinf, e) / length(e)
    fin = filter(isfinite, e)

    # lower bound: the optimistic quotient with terminal cost min over the block
    Lb = lower_bound(Abstraction(W, μ), goal)
    # upper bound: inflated quotient; terminal defect from the exact V on the goal blocks
    bq = fill(Inf, nc)
    for x in goal
        bq[μ.vertex_map[x]] = 0.0
    end
    εterm = zeros(nc)
    for x in eachindex(V)
        B = μ.vertex_map[x]
        isfinite(bq[B]) && (εterm[B] = max(εterm[B], V[x] - bq[B]))
    end
    G = with_sink(Wq, ε, bq, εterm)
    Uq = cost_to_go(G, [nc + 1])[1:nc]
    U = Uq[μ.vertex_map]
    # pessimistic (∀∃) quotient on the same viable graph, terminal cost max over the block
    Wm = SparseMatrixCSC{Float64,Int}(robust_quotient(W, μ))
    bm = zeros(nc)
    isgoal = falses(length(V)); isgoal[goal] .= true
    for x in eachindex(V)
        B = μ.vertex_map[x]
        isgoal[x] || (bm[B] = Inf)
    end
    Gm = with_sink(Wm, SparseMatrixCSC{Float64,Int}(sparse(findnz(Wm)[1], findnz(Wm)[2], zeros(nnz(Wm)), nc, nc)), bm, zeros(nc))
    Um = cost_to_go(Gm, [nc + 1])[1:nc][μ.vertex_map]

    @assert all(Lb[x] <= V[x] + 1e-9 for x in eachindex(V))                 # Theorem 4.2, lower side
    @assert all(V[x] <= U[x] + 1e-9 for x in eachindex(V))                  # Theorem 4.2, upper side
    @assert all(U[x] <= Um[x] + 1e-9 for x in eachindex(V))                 # Theorem 4.2, against must

    lr = mean(Lb[reach] ./ V[reach])
    cu = filter(x -> isfinite(U[x]), reach)
    cm = filter(x -> isfinite(Um[x]), reach)
    ur = isempty(cu) ? NaN : mean(U[cu] ./ V[cu])
    mr = isempty(cm) ? NaN : mean(Um[cm] ./ V[cm])
    gap = U[cu] .- Lb[cu]
    loss = V[cu] .- Lb[cu]
    corr = length(cu) > 2 && std(gap) > 0 && std(loss) > 0 ? cor(gap, loss) : NaN

    # The macro-step closed loop of Theorem 5.2(b) from sampled certified states.
    S = _DefectSearch(size(W, 1))
    Gt = sparse(transpose(G))
    members = [Int[] for _ in 1:nc]
    for (x, B) in enumerate(μ.vertex_map)
        push!(members[B], x)
    end
    goalset = BitSet(goal)
    ok = 0
    tried = 0
    for x0 in (isempty(cu) ? Int[] : rand(rng, cu, min(samples, length(cu))))
        tried += 1
        x, spent, good = x0, 0.0, true
        while !(x in goalset)
            B = μ.vertex_map[x]
            # next block on an optimal walk of the inflated quotient (or the sink)
            best, nxt, c = Inf, 0, Inf
            for k in nzrange(Gt, B)
                C = rowvals(Gt)[k]
                w = nonzeros(Gt)[k]
                v = C == nc + 1 ? 0.0 : Uq[C]
                if w + v < best
                    best, nxt, c = w + v, C, w
                end
            end
            # follow it by a cheapest fine path (forward search on the transpose)
            sources = [x]
            _defect_search!(S, Wt, sources, y -> true, c + 1e-9)
            targets = nxt == nc + 1 ? goal : members[nxt]
            y = argmin(z -> S.dist[z], targets)
            d = S.dist[y]
            value_after = nxt == nc + 1 ? 0.0 : Uq[nxt]
            if !(d <= c + 1e-9 && value_after + d <= Uq[B] + 1e-9)
                good = false
                break
            end
            spent += d
            x = y
        end
        ok += good && spent <= U[x0] + 1e-9 && spent >= V[x0] - 1e-9
    end

    @printf("%-28s %5d %7d %5.1f%% %6.1f%% %6.2f  %5.3f   %8.1f%% / %5.1f%%        %6.3f / %6.3f        %6.3f       %d/%d   %6.1f s\n",
        name, T.K, length(e), 100lift, 100over, isempty(fin) ? NaN : mean(fin), lr,
        100length(cu) / length(reach), 100length(cm) / length(reach), ur, mr, corr, ok, tried, t)
    flush(stdout)
    push!(rows, (name, T.K, length(e), lift, over, isempty(fin) ? NaN : mean(fin), isempty(fin) ? NaN : median(fin),
        lr, length(cu) / length(reach), length(cm) / length(reach), ur, mr, corr, ok, tried, t))
end

open(joinpath(@__DIR__, "..", "figures", "tropical_defect_n$(n)_v$(vmax).csv"), "w") do io
    println(io, "truncation,modes,arcs,lift_fraction,beyond_horizon_fraction,mean_finite_defect,median_finite_defect,lower_over_optimal,certified_inflated,certified_must,inflated_over_optimal,must_over_optimal,corr_gap_loss,loop_ok,loop_tried,seconds")
    for r in rows
        println(io, join(r, ","))
    end
end
