# Lyapunov certificates from lax morphisms of min-plus sheaves, on the double
# integrator lattice. For each velocity-mode truncation of the pushforward onto
# positions, the pessimistic (∀∃) quotient makes refinement a lax morphism, and the
# refined coarse value is a control Lyapunov function of the 4-D system: the greedy
# closed loop on it is certified to reach the goal. How much of the state space is
# certified, and how far is the certified cost from the optimum?
#
#   julia --project=docs docs/scripts/tropical_lyapunov_study.jl [n] [vmax]

using CellularSheaves
using SparseArrays
using Printf
using Statistics

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

n = length(ARGS) >= 1 ? parse(Int, ARGS[1]) : 128
vmax = length(ARGS) >= 2 ? parse(Int, ARGS[2]) : 6
m = n ÷ 3
walls(qx, qy) = (qx in m:(m + 1) && !(qy in 2:(n ÷ 6))) || (qx in (2m):(2m + 1) && !(qy in (n - n ÷ 6):(n - 1)))
W, points = double_integrator_lattice(n, vmax, walls)
L = LinearIndices(points)
stop = vmax + 1
goal = [L[qx, qy, stop, stop] for qx in (n - 3):(n - 1) for qy in 2:4]
b = [x in goal ? 0.0 : Inf for x in 1:size(W, 1)]
V = cost_to_go(W, goal)
reach = findall(isfinite, V)
println("lattice $(n)² × $(2vmax + 1)², $(size(W, 1)) states, $(length(reach)) can reach the goal")
@printf("optimal value: control Lyapunov %s; greedy closed loop %s\n", is_control_lyapunov(W, V, b),
    is_lyapunov(greedy_policy(W, V), V, goal))

nv = 2vmax + 1
truncations = [("1 mode (velocity-blind)", ones(Int, nv, nv)),
               ("8 headings × 3 speeds", heading_speed_modes(vmax, 8, [2, 4, vmax])),
               ("16 headings × 3 speeds", heading_speed_modes(vmax, 16, [2, 4, vmax])),
               ("16 headings × 6 speeds", heading_speed_modes(vmax, 16, collect(1:vmax))),
               ("every velocity (isomorphism)", reshape(1:(nv * nv), nv, nv))]
println("\ntruncation                     modes  lax E  certified states (of reachable)  E·Vc / V (mean, worst)  closed loop  build+solve")
rows = []
for (name, modes) in truncations
    T = TruncatedPushforward(W, points, (1, 2), modes)
    t = @elapsed begin
        Wr = robust_quotient(W, T.μ)
        Vup = upper_bound(T, goal)
    end
    E = refinement_matrix(T)
    lax = intertwines(Wr, W, E).lax
    cov = findall(isfinite, Vup)
    nxt = greedy_policy(W, Vup)
    cert = is_lyapunov(nxt, Vup, [x for x in goal if Vup[x] == 0])
    clf = is_control_lyapunov(W, Vup, b)
    ratio = Vup[cov] ./ V[cov]
    ok = clf && cert.positive_definite && cert.strict && cert.closed
    @printf("%-30s %5d  %-5s %10d (%5.1f%%) %20.3f %8.3f %12s %9.1f s\n", name, T.K, lax, length(cov),
        100 * length(cov) / length(reach), mean(ratio), maximum(ratio), ok ? "certified" : "FAILED", t)
    push!(rows, (name, T.K, lax, length(cov), length(cov) / length(reach), mean(ratio), maximum(ratio), ok, t))
end

open(joinpath(@__DIR__, "..", "figures", "tropical_lyapunov.csv"), "w") do io
    println(io, "truncation,modes,refinement_lax,certified_states,certified_fraction,mean_upper_over_optimal,max_upper_over_optimal,certified,seconds")
    for r in rows
        println(io, join(r, ","))
    end
end
