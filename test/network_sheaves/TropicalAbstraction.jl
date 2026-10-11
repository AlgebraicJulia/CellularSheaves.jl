using Test
using CellularSheaves
using SparseArrays
using Random

# The minimum-time double integrator on a lattice: states (qx, qy, vx, vy) with
# positions 1:n and velocities -vmax:vmax; a step picks an acceleration in
# {-1, 0, 1}², v' = v + a, q' = q + v', and costs 1. Steps that leave the grid or
# end in (or pass over the midpoint of) an obstacle cell are not allowed.
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

# Bellman–Ford: the min-plus fixed point by value iteration.
function bellman_ford(W, targets)
    V = fill(Inf, size(W, 1))
    V[targets] .= 0
    rows, cols, vals = findnz(W)
    for _ in 1:size(W, 1)
        changed = false
        for (x, y, w) in zip(rows, cols, vals)
            if V[y] + w < V[x]
                V[x] = V[y] + w
                changed = true
            end
        end
        changed || break
    end
    return V
end

@testset "TropicalAbstraction" begin
    rng = Random.MersenneTwister(3)

    # Min-plus algebra of the pushforward on a random graph: Eᵀ ⊗ E = I, the
    # quotient is Eᵀ ⊗ W ⊗ E, and Dijkstra solves the Bellman equation.
    n = (6, 5)
    N = prod(n)
    W = sprand(rng, N, N, 0.15) .* 3
    ψ = aggregation_homomorphism(n, (2, 2))
    vc = rand(rng, ψ.n_target)
    @test tropical_pushforward(ψ, tropical_pullback(ψ, vc)) == vc
    Wd = fill(Inf, N, N)
    for (x, y, w) in zip(findnz(W)...)
        Wd[x, y] = w
    end
    Q = quotient_graph(W, ψ)
    Qd = fill(Inf, ψ.n_target, ψ.n_target)
    for (b, c, w) in zip(findnz(Q)...)
        Qd[b, c] = w
    end
    fibre(b) = findall(==(b), ψ.vertex_map)
    @test all(Qd[b, c] == minimum(Wd[fibre(b), fibre(c)]) for b in 1:ψ.n_target, c in 1:ψ.n_target if b != c)
    @test all(isinf, Qd[b, b] for b in 1:ψ.n_target)
    targets = [1, 17]
    @test cost_to_go(W, targets).values == bellman_ford(W, targets)

    # Minimum-time double integrator around a wall with a gap.
    wall(qx, qy) = qx in 11:12 && !(qy in 3:5)
    W, points = double_integrator_lattice(24, 2, wall)
    L = LinearIndices(points)
    goal = [L[qx, qy, 3, 3] for qx in 20:22 for qy in 20:22]          # stopped in the goal box
    sol = cost_to_go(W, goal)
    V = sol.values
    @test sol.expanded == count(isfinite, V)
    @test V ≈ bellman_ford(W, goal)

    # A hierarchy of position aggregations, ψ₁ (blocks of 2 × 2), ψ₂ = ψ₁' ∘ ψ₁ (4 × 4), ….
    ψs = GraphHomomorphism[]
    sizes = points
    ψ = GraphHomomorphism(collect(1:prod(points)), prod(points))
    for _ in 1:3
        step = aggregation_homomorphism(sizes, (2, 2, 1, 1))
        ψ = compose(ψ, step)
        push!(ψs, ψ)
        sizes = coarse_points(sizes, (2, 2, 1, 1))
    end

    sources = [L[qx, qy, 3, 3] for (qx, qy) in ((2, 2), (3, 22), (6, 12), (2, 14))]
    @test all(isfinite, V[sources])
    for ψ in ψs
        # The coarse cost-to-go, pulled back, is a lower bound on the fine one.
        Vc = cost_to_go(quotient_graph(W, ψ), unique(ψ.vertex_map[goal])).values
        @test all(tropical_pullback(ψ, Vc) .<= V .+ 1e-12)
        for s in sources
            h = abstraction_heuristic(W, ψ, s)
            @test h[s] == 0
            # Consistent: h(y) ≤ h(x) + W[x, y] on every edge out of a reachable state.
            @test all(isinf(h[x]) || h[y] <= h[x] + w for (x, y, w) in zip(findnz(W)...))
            # A* is exact and expands fewer states than Dijkstra.
            dijkstra = cost_to_go(W, goal; source = s)
            astar = cost_to_go(W, goal; source = s, heuristic = h)
            @test dijkstra.cost == astar.cost == V[s]
            @test astar.expanded < dijkstra.expanded
        end
    end

    @test_throws ArgumentError cost_to_go(-W, goal)
    @test_throws ArgumentError cost_to_go(W, goal; heuristic = zeros(size(W, 1)))
    @test_throws ArgumentError quotient_graph(W, aggregation_homomorphism((6, 5), (2, 2)))
end
