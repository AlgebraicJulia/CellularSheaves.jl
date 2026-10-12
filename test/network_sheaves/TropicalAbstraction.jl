using Test
using CellularSheaves
using SparseArrays
using Random
using PATHPACK

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

# Bellman–Ford: the min-plus closure by value iteration.
function bellman_ford(W, b)
    V = float.(copy(b))
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

function dense_weights(W, zero)
    D = fill(zero, size(W))
    for (x, y, w) in zip(findnz(W)...)
        D[x, y] = w
    end
    return D
end
# Velocity modes: rest, then `headings` direction sectors for each speed band
# (Chebyshev speed ≤ speed_edges[k]).
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

# The min-plus closure by Floyd–Warshall (costs ≥ 0).
function minplus_star(D)
    C = copy(D)
    n = size(C, 1)
    for k in 1:n, i in 1:n, j in 1:n
        C[i, j] = min(C[i, j], C[i, k] + C[k, j])
    end
    for i in 1:n
        C[i, i] = min(C[i, i], 0.0)
    end
    return C
end

minplus(A, B) =[minimum(A[i, k] + B[k, j] for k in axes(A, 2)) for i in axes(A, 1), j in axes(B, 2)]

@testset "TropicalAbstraction" begin
    rng = Random.MersenneTwister(3)

    # ---- the semiring algebra on a random graph --------------------------------
    n = (6, 5)
    N = prod(n)
    W = sprand(rng, N, N, 0.15) .* 3
    ψ = aggregation_homomorphism(n, (2, 2))
    vc = rand(rng, ψ.n_target)
    @test tropical_pushforward(ψ, tropical_pullback(ψ, vc)) == vc                  # Eᵀ ⊗ E = I
    @test tropical_pushforward(ψ, tropical_pullback(ψ, vc); combine = max) == vc
    fibre(b) = findall(==(b), ψ.vertex_map)
    Wd = dense_weights(W, Inf)
    for steps in (1, 2, 3)
        Q = dense_weights(quotient_graph(W, ψ; steps), Inf)
        Wk = reduce(minplus, fill(Wd, steps))
        @test all(Q[b, c] ≈ minimum(Wk[fibre(b), fibre(c)]) for b in 1:ψ.n_target, c in 1:ψ.n_target if b != c)
        @test all(isinf, Q[b, b] for b in 1:ψ.n_target)
    end
    Qmax = dense_weights(quotient_graph(W, ψ; combine = max), -Inf)
    Wdmax = dense_weights(W, -Inf)
    @test all(Qmax[b, c] == maximum(Wdmax[fibre(b), fibre(c)]) for b in 1:ψ.n_target, c in 1:ψ.n_target if b != c)
    b = fill(Inf, N); b[[1, 17]] .= 0; b[9] = 0.5
    @test minplus_closure(W, b) ≈ bellman_ford(W, b)
    @test cost_to_go(W, [1, 17]) == minplus_closure(W, [x in (1, 17) ? 0.0 : Inf for x in 1:N])

    # ---- elimination (Schur complement) and its commutation with aggregation ----
    n = (6, 6)
    N = prod(n)
    W = sprand(rng, N, N, 0.12) .* 3
    D = dense_weights(W, Inf)
    Dstar = minplus_star(D)
    ψ = aggregation_homomorphism(n, (2, 2))
    Sc = sort(randperm(rng, ψ.n_target)[1:4])                  # separator blocks
    S = findall(in(Sc), ψ.vertex_map)                            # saturated: a union of blocks
    I = setdiff(1:N, S)
    WS = schur_complement(W, S)
    @test dense_weights(WS, Inf) ≈ min.(D[S, S], minplus(minplus(D[S, I], minplus_star(D[I, I])), D[I, S]))
    @test minplus_star(dense_weights(WS, Inf)) ≈ Dstar[S, S]     # elimination is exact on S

    # ψ restricted to S, onto the separator blocks.
    position = Dict(b => k for (k, b) in enumerate(Sc))
    ψS = GraphHomomorphism([position[b] for b in ψ.vertex_map[S]], length(Sc))
    offdiag(M) = [M[a, b] for a in axes(M, 1), b in axes(M, 2) if a != b]
    Wq = quotient_graph(W, ψ)
    # Restriction (dropping I) commutes with aggregation exactly: base change.
    @test offdiag(dense_weights(quotient_graph(W[S, S], ψS), Inf)) == offdiag(dense_weights(Wq, Inf)[Sc, Sc])
    # Elimination commutes laxly: aggregate-then-eliminate ≤ eliminate-then-aggregate,
    # both as operators and as closures, and both bound the aggregated exact closure.
    AE = dense_weights(schur_complement(Wq, Sc), Inf)            # aggregate, then eliminate
    EA = dense_weights(quotient_graph(WS, ψS), Inf)              # eliminate, then aggregate
    @test all(offdiag(AE) .<= offdiag(EA) .+ 1e-12)
    exact = [minimum(Dstar[S, S][findall(==(a), ψS.vertex_map), findall(==(b), ψS.vertex_map)])
             for a in 1:length(Sc), b in 1:length(Sc)]
    @test all(minplus_star(AE) .<= minplus_star(EA) .+ 1e-12)
    @test all(minplus_star(EA) .<= exact .+ 1e-12)
    @test any(minplus_star(AE) .< minplus_star(EA) .- 1e-9)     # strictly lax here
    # … and strictly commutes when the aggregation is one-to-one on the eliminated interior.
    ids = [b in Sc ? position[b] : 0 for b in 1:ψ.n_target]
    vmap = [ids[ψ.vertex_map[x]] for x in 1:N]
    vmap[I] .= length(Sc) .+ (1:length(I))
    ψ1 = GraphHomomorphism(vmap, length(Sc) + length(I))
    AE1 = dense_weights(schur_complement(quotient_graph(W, ψ1), collect(1:length(Sc))), Inf)
    @test offdiag(AE1) ≈ offdiag(EA)

    # ---- minimum-time double integrator around a wall with a gap ---------------
    wall(qx, qy) = qx in 11:12 && !(qy in 3:5)
    W, points = double_integrator_lattice(24, 2, wall)
    L = LinearIndices(points)
    goal = [L[qx, qy, 3, 3] for qx in 20:22 for qy in 20:22]          # stopped in the goal box
    V = cost_to_go(W, goal)
    @test V == bellman_ford(W, [x in goal ? 0.0 : Inf for x in 1:size(W, 1)])
    sources = [L[qx, qy, 3, 3] for (qx, qy) in ((2, 2), (3, 22), (6, 12), (2, 14))]
    @test all(isfinite, V[sources])
    Wout = sparse(transpose(W))
    dijkstra = [shortest_path(W, s, goal; out_edges = Wout) for s in sources]
    @test [r.cost for r in dijkstra] == V[sources]

    # A path is a sequence of moves from the source into the goal whose costs add up.
    valid(r, s) = first(r.path) == s && last(r.path) in goal &&
        sum(W[r.path[i], r.path[i + 1]] for i in 1:(length(r.path) - 1); init = 0.0) == r.cost

    # Abstractions: positions only (blocks 2^ℓ, one step), and space–time (velocity
    # blocks k, position blocks k², k steps).
    id = GraphHomomorphism(collect(1:prod(points)), prod(points))
    function coarsen_grid(ψ, sizes, factors)
        step = aggregation_homomorphism(sizes, factors)
        return compose(ψ, step), coarse_points(sizes, factors)
    end
    space = Tuple{GraphHomomorphism,Int}[]
    ψ, sizes = id, points
    for _ in 1:3
        ψ, sizes = coarsen_grid(ψ, sizes, (2, 2, 1, 1))
        push!(space, (ψ, 1))
    end
    spacetime = Tuple{GraphHomomorphism,Int}[]
    ψ, sizes = id, points
    for k in (2, 4)
        ψ, sizes = coarsen_grid(ψ, sizes, (2, 2, 2, 2))
        ψ, sizes = coarsen_grid(ψ, sizes, (2, 2, 1, 1))
        push!(spacetime, (ψ, k))
    end

    for (ψ, steps) in vcat(space, spacetime)
        A = Abstraction(W, ψ; steps)
        h = lower_bound(A, goal)
        @test all(h .<= V)                                            # a lower bound everywhere
        @test all(h[goal] .== 0)
        if steps == 1                                                 # consistent
            @test all(isinf(h[y]) || h[x] <= w + h[y] for (x, y, w) in zip(findnz(W)...))
        end
        for (s, d) in zip(sources, dijkstra)
            r = shortest_path(W, s, goal; heuristic = h, out_edges = Wout)
            @test r.cost == V[s] && valid(r, s)                       # exact
            steps == 1 && ψ === space[1][1] && @test r.expanded < d.expanded
        end
        # PathPack's chordal min-plus LU gives the same coarse closure.
        @test lower_bound(Abstraction(W, ψ; steps, closure = PathPackClosure()), goal) == h
    end
    # Coarsening time with space keeps the bound tight where space alone loses it.
    hs = lower_bound(Abstraction(W, space[2][1]), goal)                # positions 4 × 4, one step
    hst = lower_bound(Abstraction(W, spacetime[1][1]; steps = spacetime[1][2]), goal)          # positions 4 × 4, velocities 2, two steps
    @test sum(hst[sources]) > sum(hs[sources])

    # ---- reach-avoid in the bottleneck semiring --------------------------------
    walls = [(qx, qy) for qx in 1:24, qy in 1:24 if wall(qx, qy)]
    clearance(qx, qy) = minimum(hypot(qx - a, qy - b) for (a, b) in walls) - 0.5
    l = [wall(I[1], I[2]) ? -1.0 : clearance(I[1], I[2]) for I in CartesianIndices(points)][:]
    g = fill(-Inf, length(l)); g[goal] .= Inf
    R = reach_avoid(W, l, g)
    # The Bellman equation, by value iteration in (max, min).
    Rvi = min.(l, g)
    rows, cols = findnz(W)[1:2]
    for _ in 1:size(W, 1)
        changed = false
        for (x, y) in zip(rows, cols)
            d = min(l[x], Rvi[y])
            if d > Rvi[x]
                Rvi[x] = d; changed = true
            end
        end
        changed || break
    end
    @test R == Rvi
    @test all(R .<= l)
    # Reachable at all ⇔ finite minimum time.
    @test (R .> -Inf) == isfinite.(V)
    # R(x) ≥ c ⇔ x reaches the goal through states of clearance ≥ c: a min-plus
    # question on the subgraph of those states.
    for c in (0.5, 1.5, 2.5)
        keep = [l[x] >= c for x in rows]
        Wc = sparse(rows[keep], cols[keep], ones(count(keep)), size(W)...)
        safe_goal = [t for t in goal if l[t] >= c]
        reachable = isempty(safe_goal) ? falses(length(l)) : isfinite.(cost_to_go(Wc, safe_goal)) .& (l .>= c)
        @test (R .>= c) == reachable
    end
    # The max pushforward over-approximates: an upper bound on the reach-avoid value.
    Amax = sparse(rows, cols, l[rows], size(W)...)
    ψ = space[1][1]
    Rc = maxmin_closure(quotient_graph(Amax, ψ; combine = max), tropical_pushforward(ψ, min.(l, g); combine = max))
    @test all(tropical_pullback(ψ, Rc) .>= R)

    # ---- option 1: the coarse all-pairs table (pattern database) ---------------
    for (ψ, steps) in (space[1], spacetime[1])
        h = lower_bound(Abstraction(W, ψ; steps), goal)
        @test lower_bound(Abstraction(W, ψ; steps, table = true), goal) == h
        @test lower_bound(Abstraction(W, ψ; steps, closure = PathPackClosure(), table = true), goal) == h
    end

    # ---- option 2: exact elimination of cell interiors (overlay) ---------------
    cell = [((I[1] - 1) ÷ 8) * 3 + (I[2] - 1) ÷ 8 + 1 for I in CartesianIndices(points)][:]
    O = Overlay(W, cell)
    @test length(O.portals) < size(W, 1)
    @test all(isconcretetype, fieldtypes(typeof(O)))
    hfine = lower_bound(Abstraction(W, space[1][1]), goal)
    for s in sources
        r = overlay_path(O, s, goal)
        @test r.cost == V[s] && first(r.path) == s && last(r.path) in goal
        ra = overlay_path(O, s, goal; heuristic = hfine)
        @test ra.cost == V[s]
    end

    # ---- layered control: a 2-D route, a 4-D search in its corridor ------------
    ψp = projection_homomorphism(points, (1, 2))
    @test ψp.vertex_map == [LinearIndices((24, 24))[I[1], I[2]] for I in CartesianIndices(points)][:]
    P = LayeredPlanner(W, ψp, (24, 24))
    hp = lower_bound(P.abstraction, goal)
    @test all(hp .<= V)

    # A twisted refinement (a potential): zero changes nothing; the exact cost-to-go is
    # reproduced exactly; a consistent bound as potential gives a bound between it and V.
    reach = isfinite.(V)
    @test lower_bound(Abstraction(W, ψp; potential = zeros(size(W, 1))), goal) == hp
    hexact = lower_bound(Abstraction(W, ψp; potential = V), goal)
    @test hexact[reach] ≈ V[reach] && all(isinf, hexact[.!reach])
    htw = lower_bound(Abstraction(W, space[1][1]; potential = hp), goal)
    @test all(htw[reach] .<= V[reach] .+ 1e-9) && all(htw[reach] .>= hp[reach] .- 1e-9)
    @test count(htw[reach] .> hp[reach] .+ 1e-9) > 0
    @test shortest_path(W, sources[1], goal; heuristic = htw).cost == V[sources[1]]
    @test_throws ArgumentError Abstraction(W, ψp; potential = 2 .* V)                 # not consistent

    # Base change: the quotient of a tensor product (two agents moving at once) along a
    # product map is the tensor product of the quotients with their loops.
    rngt = MersenneTwister(5)
    edgeset(M) = Set(zip(findnz(M)...))
    for _ in 1:5
        A, B = sprand(rngt, 7, 7, 0.4) .* 3, sprand(rngt, 6, 6, 0.4) .* 2
        ψA, ψB = GraphHomomorphism(rand(rngt, 1:3, 7), 3), GraphHomomorphism(rand(rngt, 1:2, 6), 2)
        AB = tensor_product(A, B)
        @test nnz(AB) == nnz(A) * nnz(B) && AB[2 + 7 * (3 - 1), 4 + 7 * (5 - 1)] == (A[2, 4] == 0 || B[3, 5] == 0 ? 0 : A[2, 4] + B[3, 5])
        prodmap = GraphHomomorphism([ψA.vertex_map[a] + 3 * (ψB.vertex_map[b] - 1) for b in 1:6 for a in 1:7], 6)
        @test edgeset(quotient_graph(AB, prodmap; loops = true)) ==
              edgeset(tensor_product(quotient_graph(A, ψA; loops = true), quotient_graph(B, ψB; loops = true)))
    end
    @test nnz(tensor_product(sparse([1.0;;]), sparse([2.0;;]); keep = (x, y) -> false)) == 0
    @test_throws ArgumentError quotient_graph(W, ψp; steps = 2, loops = true)

    # Reeb arcs: on a path with the goal in the middle, equidistant states on the two sides
    # are different arcs; around a cycle they meet at the far side; unreachable states
    # form their own arc. Zero-cost edges count as adjacency.
    path = sparse([1:6; 2:7], [2:7; 1:6], ones(12), 7, 7)
    arcs = reeb_arcs(path, abs.((1:7) .- 4.0)).vertex_map
    @test arcs[3] != arcs[5] && arcs[2] != arcs[6] && arcs[1] != arcs[7]
    @test length(unique(arcs)) == 7
    cyc = sparse([1:8; [2:8; 1]], [[2:8; 1]; 1:8], ones(16), 8, 8)
    carcs = reeb_arcs(cyc, Float64[0, 1, 2, 3, 4, 3, 2, 1]).vertex_map
    @test carcs[4] == carcs[6] && carcs[2] != carcs[8] && carcs[3] != carcs[7]
    zarcs = reeb_arcs(sparse([1, 2], [2, 3], [0.0, 1.0], 4, 4), [1.0, 1.0, 0.0, Inf])
    @test zarcs.vertex_map[1] == zarcs.vertex_map[2] && zarcs.vertex_map[4] == zarcs.n_target
    @test_throws ArgumentError lower_bound(Abstraction(W, ψp; potential = fill(Inf, size(W, 1))), goal)
    Wq = P.abstraction.quotient
    for s in sources
        route = plan_route(P, s, goal)
        @test first(route) == ψp.vertex_map[s] && last(route) in ψp.vertex_map[goal]
        @test all(Wq[a, b] > 0 for (a, b) in zip(route[1:(end - 1)], route[2:end]))   # moves of the quotient
        @test all(corridor((24, 24), route, 0)[route])
        @test all(corridor((24, 24), route, 1) .>= corridor((24, 24), route, 0))
        for width in (0, 2, 4)
            r = layered_path(P, s, goal; width, heuristic = hp)
            cells = corridor((24, 24), route, r.width)
            @test isfinite(r.cost) && r.cost >= V[s]                                 # a feasible path, never better than optimal
            @test valid(r, s) && all(cells[ψp.vertex_map[r.path]])                   # made of moves, inside its corridor
            @test r.width == width || !r.first_try
        end
        full = layered_path(P, s, goal; width = 24, heuristic = hp)
        @test full.cost == V[s] && full.first_try                                    # the whole grid: exact
    end
    @test !isfinite(layered_path(P, L[11, 10, 3, 3], goal).cost)                     # a state inside the wall

    # ---- truncated pushforward: velocity modes in the stalks over positions ----
    modes = heading_speed_modes(2, 8, [1, 2])                                         # rest, 8 headings × 2 speeds
    T = TruncatedPushforward(W, points, (1, 2), modes)
    @test T.K == 17
    @test T isa TruncatedPushforward{2} && all(isconcretetype, fieldtypes(typeof(T)))
    # The inner constructor keeps the modes a refinement of the cells (π ∘ μ = ψ).
    πbad = GraphHomomorphism(circshift(T.π.vertex_map, 1), T.π.n_target)
    @test_throws ArgumentError TruncatedPushforward(T.W, T.ψ, T.μ, πbad, T.K, T.quotient, T.cell_points, T.mode_neighbors)
    @test_throws ArgumentError TruncatedPushforward(T.W, T.ψ, T.μ, T.π, T.K, T.quotient, T.cell_points, T.mode_neighbors[1:2])
    @test T.π.vertex_map[T.μ.vertex_map] == ψp.vertex_map                           # stalkwise: over the identity of the cell graph
    E, Pm = refinement_matrix(T), coarsening_matrix(T)
    vc = rand(rng, T.μ.n_target)
    @test coarsen(T, refine(T, vc)) == vc                                            # P E = I
    f = rand(rng, size(W, 1))
    @test all(refine(T, coarsen(T, f)) .<= f)                                         # E P ≤ I
    for _ in 1:20                                                                     # Galois connection E g ≤ f ⟺ g ≤ P f
        g = coarsen(T, f) .+ (rand(rng, T.μ.n_target) .- 0.9) .* 0.1
        @test all(refine(T, g) .<= f) == all(g .<= coarsen(T, f))
    end
    @test intertwines(W, T.quotient, Pm) == (lax = true, oplax = false, strict = false)   # P: lax morphism F → G
    @test intertwines(T.quotient, W, E) == (lax = false, oplax = true, strict = false)    # E: oplax morphism G → F
    # One mode per state: an isomorphism, strict both ways. One mode per stalk: the
    # velocity-blind quotient.
    Tiso = TruncatedPushforward(W, points, (1, 2), reshape(1:25, 5, 5))
    @test intertwines(W, Tiso.quotient, coarsening_matrix(Tiso)).strict
    @test intertwines(Tiso.quotient, W, refinement_matrix(Tiso)).strict
    T1 = TruncatedPushforward(W, points, (1, 2), ones(Int, 5, 5))
    @test offdiag(dense_weights(T1.quotient, Inf)) == offdiag(dense_weights(quotient_graph(W, ψp), Inf))
    # The oplax square makes the refined coarse solution a lower bound, tighter than
    # the velocity-blind one.
    bterm = [x in goal ? 0.0 : Inf for x in 1:size(W, 1)]
    hm = refine(T, minplus_closure(T.quotient, coarsen(T, bterm)))
    @test all(hm .<= V)
    @test all(hm .>= hp)
    @test lower_bound(Abstraction(W, T.μ), goal) == hm
    # The moded planner: routes carry modes; the corridor restricts positions and modes.
    PM = LayeredPlanner(T)
    @test all(isconcretetype, fieldtypes(typeof(PM)))
    @test PM.abstraction.W === T.W
    @test_throws ArgumentError LayeredPlanner(W, ψp, (24, 23))                         # cells ≠ blocks of ψ
    for s in sources
        route = plan_route(PM, s, goal)
        @test first(route) == T.μ.vertex_map[s] && last(route) in T.μ.vertex_map[goal]
        for (width, mode_width) in ((0, 0), (2, 1), (4, 2))
            r = layered_path(PM, s, goal; width, mode_width, heuristic = hm)
            U = corridor(PM, route, r.width, r.mode_width)
            @test isfinite(r.cost) && r.cost >= V[s] && valid(r, s) && all(U[T.μ.vertex_map[r.path]])
        end
        @test layered_path(PM, s, goal; width = 24, mode_width = 17, heuristic = hm).cost == V[s]
    end

    # ---- Lyapunov certificates from lax morphisms of sheaves ---------------------
    # The optimal value is a control Lyapunov function, and the greedy closed loop
    # on it is a Lyapunov morphism into (ℝ≥0, r ↦ {r}) with strict decrease.
    @test is_control_lyapunov(W, V, bterm)
    next = greedy_policy(W, V)
    @test is_lyapunov(next, V, goal) == (positive_definite = true, decrescent = true, strict = true, closed = true,
        domain = count(isfinite, V))
    # The optimistic pushforward's refinement is only oplax: its value is a lower
    # bound, not a control Lyapunov function.
    @test !is_control_lyapunov(W, hm, bterm)
    # The pessimistic (∀∃) pushforward makes refinement a lax morphism, which carries
    # the coarse value to a control Lyapunov function of the fine system.
    Wr = robust_quotient(W, T.μ)
    @test intertwines(Wr, W, E).lax
    Vup = upper_bound(T, goal)
    covered = findall(isfinite, Vup)
    @test length(covered) > length(goal)                                             # certifies more than the goal
    @test all(V[covered] .<= Vup[covered])                                           # an upper bound on the cost
    @test is_control_lyapunov(W, Vup, bterm)
    nup = greedy_policy(W, Vup)
    cert = is_lyapunov(nup, Vup, [x for x in goal if Vup[x] == 0])
    @test cert.positive_definite && cert.decrescent && cert.strict && cert.closed
    for x in covered[1:max(1, length(covered) ÷ 50):end]                             # the closed loop reaches the goal
        cost, y, steps = 0.0, x, 0
        while Vup[y] > 0
            cost += W[y, nup[y]]; y = nup[y]; steps += 1
        end
        @test y in goal && cost <= Vup[x] && steps <= Vup[x]
    end
    # One mode per state: the pessimistic quotient is exact, and so is the certificate.
    @test upper_bound(Tiso, goal) == V
    # The layered controller's low level is certified by the corridor value.
    for s in sources[1:2]
        r = layered_path(PM, s, goal; width = 4, mode_width = 1, widen = false, heuristic = hm)
        isfinite(r.cost) || continue
        c = corridor_certificate(PM, s, goal; width = 4, mode_width = 1)
        @test c.V[s] == r.cost
        inside_goal = [t for t in goal if isfinite(c.V[t])]
        cc = is_lyapunov(c.next, c.V, inside_goal)
        @test cc.positive_definite && cc.strict && cc.closed
    end

    # The certified layered controller: stability from the corridor value (lax side),
    # a suboptimality bound from the optimistic refinement (oplax side).
    for s in sources, (width, mode_width) in ((0, 0), (2, 1), (4, 1))
        c = certified_path(PM, s, goal; width, mode_width, lower = hm)
        l = layered_path(PM, s, goal; width = c.width, mode_width = c.mode_width, widen = false, heuristic = hm)
        @test c.cost == l.cost >= V[s]                                                # the same controller cost
        @test c.certificate.positive_definite && c.certificate.strict && c.certificate.closed
        @test first(c.path) == s && last(c.path) in goal && valid(c, s)              # its closed-loop trajectory
        @test c.suboptimality_bound >= c.cost / V[s] - 1e-12                          # a guaranteed bound on cost / optimum
    end
    @test !isfinite(certified_path(PM, L[11, 10, 3, 3], goal).cost)

    # ---- edge cases: zero-cost moves, stored ∞ costs, unreachable sources ------
    # A zero-cost cycle that ties with the way out, and a zero-cost last move: the
    # certified closed loop still ends in the target.
    Z = sparse([1, 2, 1], [2, 1, 3], [0.0, 0.0, 1.0], 3, 3)
    PZ = LayeredPlanner(Z, GraphHomomorphism([1, 2, 3], 3), (3,))
    cz = certified_path(PZ, 1, [3]; width = 3)
    @test cz.cost == 1 && cz.path[1] == 1 && last(cz.path) == 3
    Z2 = sparse([1, 2], [2, 3], [1.0, 0.0], 3, 3)
    cz2 = certified_path(LayeredPlanner(Z2, GraphHomomorphism([1, 2, 3], 3), (3,)), 1, [3]; width = 3)
    @test cz2.path == [1, 2, 3] && cz2.cost == 1
    # No path: the bound is ∞, not NaN.
    cu = certified_path(PZ, 3, [1]; width = 3, lower = fill(Inf, 3), widen = false)
    @test !isfinite(cu.cost) && cu.suboptimality_bound == Inf
    # A stored ∞ cost (a forbidden move) leaves the k-step quotient a valid sparse matrix.
    WI = sparse([1, 1, 2, 3], [2, 3, 4, 4], [Inf, 1.0, 1.0, 1.0], 4, 4)
    ψI = GraphHomomorphism([1, 2, 2, 3], 3)
    QI = quotient_graph(WI, ψI; steps = 2)
    @test all(allunique(rowvals(QI)[nzrange(QI, j)]) for j in 1:size(QI, 2))
    @test QI[1, 3] == 2.0

    @test_throws ArgumentError cost_to_go(-W, goal)
    @test_throws ArgumentError quotient_graph(W, aggregation_homomorphism((6, 5), (2, 2)))
    @test_throws ArgumentError quotient_graph(W, space[1][1]; steps = 2, combine = max)
    @test_throws TypeError Abstraction(W, space[1][1]; closure = :other)
    @test_throws ArgumentError Abstraction(W, space[1][1]; steps = 0)
    @test_throws ArgumentError Abstraction(W, GraphHomomorphism([1, 2], 2))         # ψ must map the states of W
    @test_throws ArgumentError shortest_path(W, sources[1], goal; heuristic = zeros(3))
end
