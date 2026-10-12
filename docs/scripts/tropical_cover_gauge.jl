# The gauge on the cover of the interaction sheaf (`tropical_pair_sheaf.jl`).
#
# The total delay of a joint plan is Σᵢ delayᵢ. Split each agent's delay among the edges of
# the interaction graph at it, state by state: weights wᵢₑ(s) ≥ 0 with Σₑ wᵢₑ(s) ≤ 1. The
# edge e = {i, j} then carries the pair problem with arc costs wᵢₑ·delayᵢ + wⱼₑ·delayⱼ.
#
#   - A joint plan restricts to a collision-free pair plan on every edge, and the costs of
#     its restrictions add up to at most its total delay. So for every split the sum over
#     edges of the weighted pair optima is a lower bound, at every joint state, and it is
#     consistent (a joint move restricts to a move of every pair).
#   - The split changes how cost is assigned to the cover without changing the cost of any
#     global section: it is a gauge transformation on the cover, as the potential φ is one on
#     an agent's own graph (which adds a coboundary and fixes the cost of every path up to
#     its ends). The agent gauge comes first: in the delay gauge every arc cost is
#     nonnegative and an agent alone costs 0, so the whole bound comes from the edges.
#   - A matching is the split with wᵢₑ ≡ 1 on matched edges. Maximizing over splits, by
#     projected supergradient ascent from the matching, gives a bound at least as good.
#
# The weighted edge stalk is the same pushforward as before (agent i exact; agent j by its
# near state or far summary); the arc cost of j on a far coarse arc is the least weighted
# delay over the fibre, so the stalk stays the exact pushforward of the weighted pair graph.

import CellularSheaves.NetworkSheaves.TropicalAbstraction as TA

struct GaugeFactor
    i::Int                       # the exact side
    j::Int                       # the summarized side
    e::Int                       # edge of the interaction graph
    P::PairAbstraction           # the stalk's coordinates (its Δ is unused)
    n::Int                       # coarse states
    offs::Vector{Int}            # arcs from x are offs[x]:offs[x + 1] - 1
    src::Vector{Int32}
    dst::Vector{Int32}
    s::Vector{Int32}             # i's state at the arc's source
    ci::Vector{Float64}          # i's delay on the arc
    jref::Vector{Int32}          # j's state at the source (> 0), or -q for the far coarse arc q
    cj::Vector{Float64}          # j's delay on the arc (near and boundary arcs)
    fibre_offs::Vector{Int}      # the fibre of far coarse arc q: j's arcs fibre_offs[q]:fibre_offs[q + 1] - 1
    fibre_t::Vector{Int32}
    fibre_c::Vector{Float64}
    target::Int
end

function GaugeFactor(M::GridMap, Ai::AgentStalk, Aj::AgentStalk, i, j, e, b)
    Si = length(Ai.pos)
    far, far_block = far_summary(M, Aj, b)
    velmode = collect(1:Aj.nvel)
    nnear = (3b)^2 * Aj.nvel
    P = PairAbstraction(b, velmode, nnear, far, far_block, Float64[])
    n = Si * (nnear + far.n_target)
    # j's arcs grouped by far coarse arc.
    qid = Dict{Tuple{Int,Int},Int}()
    fibres = Vector{Vector{Tuple{Int,Float64}}}()
    ends = Tuple{Int,Int}[]
    for (t, t2, c) in zip(findnz(Aj.delays)...)
        k = (far.vertex_map[t], far.vertex_map[t2])
        q = get!(qid, k) do
            push!(fibres, Tuple{Int,Float64}[]); push!(ends, k); length(fibres)
        end
        push!(fibres[q], (t, c))
    end
    xs, ys, ss, cis, jr, cjs = Int32[], Int32[], Int32[], Float64[], Int32[], Float64[]
    arc!(x, y, s, ci, r, cj) = (push!(xs, x); push!(ys, y); push!(ss, s); push!(cis, ci); push!(jr, r); push!(cjs, cj))
    # Far region: i's arcs times j's far coarse arcs, both ends apart.
    apart(s, B) = cheb(block(M, Ai.pos[s], b) .- far_block[B]) >= 2
    for (s, s2, c) in zip(findnz(Ai.delays)...), (q, (B, B2)) in enumerate(ends)
        (apart(s, B) && apart(s2, B2)) || continue
        arc!(coarse(P, Si, s, nnear + B), coarse(P, Si, s2, nnear + B2), s, c, -q, 0.0)
    end
    # Near region, exact.
    for s in 1:Si
        isfinite(Ai.d[s]) || continue
        for t in near_states(M, Ai, Aj, b, s)
            x = coarse(P, Si, s, near_index(M, Ai, Aj, b, velmode, s, t))
            for s2 in Ai.succ[s], t2 in Aj.succ[t]
                (Ai.pos[s2] == Aj.pos[t2] || (Ai.pos[s2] == Aj.pos[t] && Aj.pos[t2] == Ai.pos[s])) && continue
                (isfinite(Ai.d[s2]) && isfinite(Aj.d[t2])) || continue
                arc!(x, coarse(P, Si, s2, pair_mode(M, P, Ai, Aj, s2, t2)), s, delay(Ai, s, s2), t, delay(Aj, t, t2))
            end
        end
    end
    # Boundary, far to near: j enters the window from a state of its far fibre.
    for s2 in 1:Si
        isfinite(Ai.d[s2]) || continue
        for t2 in near_states(M, Ai, Aj, b, s2)
            y = coarse(P, Si, s2, near_index(M, Ai, Aj, b, velmode, s2, t2))
            for s in Ai.pred[s2], t in Aj.pred[t2]
                (isfinite(Ai.d[s]) && isfinite(Aj.d[t]) && !near(M, Ai, Aj, b, s, t)) || continue
                arc!(coarse(P, Si, s, nnear + far.vertex_map[t]), y, s, delay(Ai, s, s2), t, delay(Aj, t, t2))
            end
        end
    end
    p = sortperm(xs)
    offs = ones(Int, n + 1)
    for x in xs
        offs[x + 1] += 1
    end
    offs = cumsum(offs)
    fibre_offs = cumsum([1; length.(fibres)])
    target = coarse(P, Si, Ai.goal, pair_mode(M, P, Ai, Aj, Ai.goal, Aj.goal))
    return GaugeFactor(i, j, e, P, n, offs, xs[p], ys[p], ss[p], cis[p], jr[p], cjs[p], fibre_offs,
        Int32[t for f in fibres for (t, _) in f], Float64[c for f in fibres for (_, c) in f], target)
end

# The least weighted delay of j over each far fibre, and the arc attaining it.
function fibre_min(F::GaugeFactor, wj)
    Q = length(F.fibre_offs) - 1
    qw, qa = fill(Inf, Q), zeros(Int, Q)
    for q in 1:Q, k in F.fibre_offs[q]:(F.fibre_offs[q + 1] - 1)
        v = wj[F.fibre_t[k]] * F.fibre_c[k]
        v < qw[q] && ((qw[q], qa[q]) = (v, k))
    end
    return qw, qa
end

arc_cost(F::GaugeFactor, a, wi, wj, qw) = wi[F.s[a]] * F.ci[a] + (F.jref[a] > 0 ? wj[F.jref[a]] * F.cj[a] : qw[-F.jref[a]])

# The weighted pair optimum from coarse state x0, and its supergradient: the delay each
# agent spends at each of its states along the optimal path (j's through the fibre minima).
function factor_solve(F::GaugeFactor, wi, wj, x0)
    qw, qa = fibre_min(F, wj)
    dist = fill(Inf, F.n)
    parc = zeros(Int, F.n)
    dist[x0] = 0.0
    heap = Tuple{Float64,Int}[]
    TA._push!(heap, (0.0, x0))
    while !isempty(heap)
        dx, x = TA._pop!(heap)
        dx > dist[x] && continue
        x == F.target && break
        for a in F.offs[x]:(F.offs[x + 1] - 1)
            y = F.dst[a]
            v = dx + arc_cost(F, a, wi, wj, qw)
            if v < dist[y]
                dist[y] = v
                parc[y] = a
                TA._push!(heap, (v, y))
            end
        end
    end
    gi, gj = Dict{Int,Float64}(), Dict{Int,Float64}()
    y = F.target
    isfinite(dist[y]) || return Inf, gi, gj
    while y != x0
        a = parc[y]
        gi[F.s[a]] = get(gi, F.s[a], 0.0) + F.ci[a]
        if F.jref[a] > 0
            gj[F.jref[a]] = get(gj, F.jref[a], 0.0) + F.cj[a]
        else
            k = qa[-F.jref[a]]
            gj[F.fibre_t[k]] = get(gj, F.fibre_t[k], 0.0) + F.fibre_c[k]
        end
        y = F.src[a]
    end
    return dist[F.target], gi, gj
end

# The weighted edge stalk as a PairAbstraction: its delay-to-go at every coarse state.
function weighted_stalk(F::GaugeFactor, wi, wj)
    qw, _ = fibre_min(F, wj)
    costs = [arc_cost(F, a, wi, wj, qw) for a in eachindex(F.dst)]
    keep = isfinite.(costs)
    Wc = sparse(F.src[keep], F.dst[keep], costs[keep], F.n, F.n, min)
    bc = fill(Inf, F.n)
    bc[F.target] = 0.0
    P = F.P
    return PairAbstraction(P.b, P.velmode, P.nnear, P.far, P.far_block, minplus_closure(Wc, bc))
end

# The heaviest matching of agents 1:N with pair weights w, as a list of pairs.
function matching_pairs(w::AbstractMatrix)
    N = size(w, 1)
    best = zeros(1 << N)
    choice = zeros(Int, 1 << N)
    for mask in 1:((1 << N) - 1)
        i = trailing_zeros(mask) + 1
        rest = mask & ~(1 << (i - 1))
        v, r = best[rest + 1], rest
        while r != 0
            j = trailing_zeros(r) + 1
            r &= r - 1
            u = w[i, j] + best[(rest & ~(1 << (j - 1))) + 1]
            u > v && ((v, choice[mask + 1]) = (u, j))
        end
        best[mask + 1] = v
    end
    out, mask = Tuple{Int,Int}[], (1 << N) - 1
    while mask != 0
        i = trailing_zeros(mask) + 1
        j = choice[mask + 1]
        mask &= ~(1 << (i - 1))
        if j != 0
            push!(out, (min(i, j), max(i, j)))
            mask &= ~(1 << (j - 1))
        end
    end
    return out
end

# The smallest Σᵢ xᵢ over integers xᵢ ≥ 0 with xᵢ + xⱼ ≥ w[i, j]: the weighted vertex cover
# heuristic of Li, Felner et al. (2019). The delays of a joint plan are such an x.
function min_vertex_cover(w::AbstractMatrix)
    N = size(w, 1)
    W = ceil.(Int, w .- 1e-9)
    x = zeros(Int, N)
    best = Ref(sum(maximum(W; dims = 2)))
    function go(k, total)
        total >= best[] && return
        if k > N
            best[] = total
            return
        end
        lo = maximum((W[k, j] - x[j] for j in 1:(k - 1)); init = 0)
        hi = max(lo, maximum(W[k, :]))
        for v in lo:hi
            x[k] = v
            go(k + 1, total + v)
        end
    end
    go(1, 0)
    return float(best[])
end

# Euclidean projection onto {v ≥ 0, Σ v ≤ 1}.
function project_budget!(v)
    v .= max.(v, 0.0)
    sum(v) <= 1 && return v
    u = sort(v; rev = true)
    c, τ = 0.0, 0.0
    for (k, uk) in enumerate(u)
        c += uk
        t = (c - 1) / k
        uk > t && (τ = t)
    end
    v .= max.(v .- τ, 0.0)
    return v
end

struct CoverGauge
    sheaf::InteractionSheaf      # the weighted edge stalks
    weights::Vector{Dict{Int,Vector{Float64}}}
    matching::Float64            # the bound at the start before ascent
    bound::Float64               # the bound at the start after ascent
    history::Vector{Float64}
end

# Maximize the split bound at `starts` by projected supergradient ascent from the
# matching, with Polyak steps toward `target` (an upper bound on the total delay).
function CoverGauge(IS::InteractionSheaf, starts; target = Inf, iters = 60)
    agents, M = IS.agents, IS.map
    N = length(agents)
    weights = [Dict{Int,Vector{Float64}}() for _ in 1:N]
    isempty(IS.pairs) && return CoverGauge(IS, weights, 0.0, 0.0, Float64[])
    b = first(values(IS.stalks)).b
    w = zeros(N, N)
    for (i, j) in IS.pairs
        w[i, j] = w[j, i] = pair_bound(M, IS.stalks[(i, j)], IS.stalks[(j, i)], agents[i], agents[j], starts[i], starts[j])
    end
    matched = Set(matching_pairs(w))
    factors = Tuple{GaugeFactor,GaugeFactor}[]
    for (e, (i, j)) in enumerate(IS.pairs)
        push!(factors, (GaugeFactor(M, agents[i], agents[j], i, j, e, b), GaugeFactor(M, agents[j], agents[i], j, i, e, b)))
        on = (i, j) in matched ? 1.0 : 0.0
        weights[i][e] = fill(on, length(agents[i].pos))
        weights[j][e] = fill(on, length(agents[j].pos))
    end
    x0(F) = coarse(F.P, length(agents[F.i].pos), starts[F.i],
        pair_mode(M, F.P, agents[F.i], agents[F.j], starts[F.i], starts[F.j]))
    function evaluate(wts)
        L = 0.0
        grads = [Dict{Int,Dict{Int,Float64}}() for _ in 1:N]
        for (F, G) in factors
            v, gi, gj = factor_solve(F, wts[F.i][F.e], wts[F.j][F.e], x0(F))
            u, hi, hj = factor_solve(G, wts[G.i][G.e], wts[G.j][G.e], x0(G))
            if u > v
                v, gi, gj, F = u, hi, hj, G
            end
            L += v
            grads[F.i][F.e] = gi
            grads[F.j][F.e] = gj
        end
        return L, grads
    end
    L0, grads = evaluate(weights)
    best, bestw, history = L0, deepcopy(weights), [L0]
    θ, stall = 1.0, 0
    for _ in 1:iters
        target - best < 1e-6 && break
        L, grads = evaluate(weights)
        push!(history, L)
        if L > best + 1e-9
            best, bestw, stall = L, deepcopy(weights), 0
        elseif (stall += 1) >= 5
            θ, stall = θ / 2, 0
        end
        g2 = sum(v^2 for gs in grads for g in values(gs) for v in values(g); init = 0.0)
        g2 == 0 && break
        η = θ * (min(target, best + 2) - L) / g2
        for i in 1:N
            for (e, g) in grads[i], (s, v) in g
                weights[i][e][s] += η * v
            end
            es = collect(keys(weights[i]))
            buf = zeros(length(es))
            for s in eachindex(agents[i].pos)
                for (k, e) in enumerate(es)
                    buf[k] = weights[i][e][s]
                end
                project_budget!(buf)
                for (k, e) in enumerate(es)
                    weights[i][e][s] = buf[k]
                end
            end
        end
    end
    stalks = Dict{Tuple{Int,Int},PairAbstraction}()
    for (F, G) in factors, H in (F, G)
        stalks[(H.i, H.j)] = weighted_stalk(H, bestw[H.i][H.e], bestw[H.j][H.e])
    end
    return CoverGauge(InteractionSheaf(M, agents, IS.pairs, stalks), bestw, L0, best, history)
end

# The split bound at any joint state: the sum over all edges of the weighted pair bounds.
gauge_bound(C::CoverGauge, states) =
    sum((pair_bound(C.sheaf.map, C.sheaf.stalks[(i, j)], C.sheaf.stalks[(j, i)], C.sheaf.agents[i], C.sheaf.agents[j],
        states[i], states[j]) for (i, j) in C.sheaf.pairs); init = 0.0)
