# Numerical check of the path-lifting theorem (docs/issues/019, docs/papers/path-lifting).
#
# Plain Julia, no packages: dense Floyd–Warshall on small random cost graphs.
# For a surjection ψ of the states onto blocks it checks
#
#   (1) V_ψ ∘ ψ ≤ V                 (the optimistic quotient is a lower bound)
#   (2) V ≤ V_ψ^{+ε} ∘ ψ            (lifting bound: quotient with arc costs W_ψ + ε)
#   (3) V_ψ^{+ε} ≤ V^must           (it is never worse than the pessimistic quotient)
#   (4) ε = 0 on every arc between distinct blocks ⇔ V_ψ ∘ ψ = V for every target
#       that is one block (and then for every terminal cost constant on blocks)
#   (5) the local defect, from paths inside B ∪ C only, bounds ε and keeps (2)
#   (6) two-sided: the worst loss over one-block targets at a block B is at least
#       max_C ε(B, C) and at most max_C (V^{+ε}_C(B) − V_ψ,C(B))
#   (7) optimality: h(B, C) = max_{x ∈ B} d(x, C) (all pairs) is a sound upper
#       weighting (V ≤ V_h ∘ ψ for one-block targets), h ≤ W_ψ + ε on the arcs of
#       W_ψ, and lowering any finite entry of h makes it unsound
#
# Run: julia docs/scripts/path_lifting_check.jl
using Random

const INF = Inf

function floyd(W)
    n = size(W, 1)
    d = copy(W)
    for i in 1:n
        d[i, i] = min(d[i, i], 0.0)
    end
    for k in 1:n, i in 1:n, j in 1:n
        d[i, j] = min(d[i, j], d[i, k] + d[k, j])
    end
    return d
end

# Cost-to-go with terminal cost b, from all-pairs distances d.
ctg(d, b) = [minimum(d[x, y] + b[y] for y in eachindex(b)) for x in axes(d, 1)]

blocks_of(ψ, m) = [findall(==(B), ψ) for B in 1:m]

function optimistic(W, ψ, m)
    Wq = fill(INF, m, m)
    for x in axes(W, 1), y in axes(W, 2)
        ψ[x] != ψ[y] && (Wq[ψ[x], ψ[y]] = min(Wq[ψ[x], ψ[y]], W[x, y]))
    end
    return Wq
end

# The pessimistic (∀∃) quotient, as robust_quotient in TropicalAbstraction.jl: an arc
# B → C only when every state of B has a one-step move into C, costing the worst
# such cheapest move.
function pessimistic(W, ψ, m)
    bl = blocks_of(ψ, m)
    Wq = fill(INF, m, m)
    for B in 1:m, C in 1:m
        B == C && continue
        Wq[B, C] = maximum(minimum(W[x, y] for y in bl[C]) for x in bl[B])
    end
    return Wq
end

# Lifting defect ε(B, C) = max_{x ∈ B} (dist(x, C) − W_ψ(B, C))⁺ for the arcs of W_ψ.
function defect(dist, Wq, ψ, m)
    bl = blocks_of(ψ, m)
    ε = fill(INF, m, m)
    for B in 1:m, C in 1:m
        (B == C || isinf(Wq[B, C])) && continue
        ε[B, C] = maximum(max(minimum(dist[x, y] for y in bl[C]) - Wq[B, C], 0.0) for x in bl[B])
    end
    return ε
end

function random_case(rng; zero_internal::Bool)
    n = rand(rng, 6:24)
    m = rand(rng, 2:max(2, n ÷ 2))
    ψ = [rand(rng, 1:m) for _ in 1:n]
    ψ[1:m] .= 1:m                                    # surjective
    W = fill(INF, n, n)
    p = rand(rng, [0.15, 0.3, 0.5])
    for x in 1:n, y in 1:n
        x != y && rand(rng) < p && (W[x, y] = Float64(rand(rng, 0:4)))
    end
    if zero_internal                                 # each block a zero-cost cycle: every arc lifts
        for B in 1:m
            s = findall(==(B), ψ)
            for (i, x) in enumerate(s)
                W[x, s[mod1(i + 1, length(s))]] = length(s) > 1 ? 0.0 : W[x, x]
            end
        end
    end
    return W, ψ, m
end

function check(rng, trials)
    counts = Dict(:lower => 0, :lift => 0, :must => 0, :exact_iff => 0, :local => 0, :zero_defect_cases => 0,
                  :two_sided => 0, :optimal => 0)
    tol = 1e-9
    for t in 1:trials
        W, ψ, m = random_case(rng; zero_internal = t % 3 == 0)
        n = size(W, 1)
        d = floyd(W)
        b = [rand(rng) < 0.3 ? Float64(rand(rng, 0:3)) : INF for _ in 1:n]
        b[rand(rng, 1:n)] = 0.0
        V = ctg(d, b)

        Wq = optimistic(W, ψ, m)
        bq = [minimum(b[x] for x in 1:n if ψ[x] == B) for B in 1:m]
        Vq = ctg(floyd(Wq), bq)

        ε = defect(d, Wq, ψ, m)
        εterm = [maximum(V[x] - bq[B] for x in 1:n if ψ[x] == B; init = -INF) for B in 1:m]
        εterm = [isinf(bq[B]) ? 0.0 : max(εterm[B], 0.0) for B in 1:m]
        Vε = ctg(floyd(Wq .+ ε), bq .+ εterm)

        Wm = pessimistic(W, ψ, m)
        bm = [maximum(b[x] for x in 1:n if ψ[x] == B) for B in 1:m]
        Vm = ctg(floyd(Wm), bm)

        counts[:lower] += all(Vq[ψ[x]] <= V[x] + tol for x in 1:n)
        counts[:lift] += all(V[x] <= Vε[ψ[x]] + tol for x in 1:n)
        counts[:must] += all(Vε[B] <= Vm[B] + tol for B in 1:m)

        # (4) exactness for one-block targets ⇔ zero defect on every inter-block arc
        zero = all(ε[B, C] <= tol for B in 1:m, C in 1:m if B != C && isfinite(Wq[B, C]))  # an unreachable C is an infinite defect
        exact = true
        for C in 1:m
            bC = [ψ[x] == C ? 0.0 : INF for x in 1:n]
            VC = ctg(d, bC)
            VqC = ctg(floyd(Wq), [B == C ? 0.0 : INF for B in 1:m])
            exact &= all(abs(VqC[ψ[x]] - VC[x]) <= tol || (isinf(VqC[ψ[x]]) && isinf(VC[x])) for x in 1:n)
        end
        if zero                                      # then also exact for terminal costs constant on blocks
            bB = [rand(rng) < 0.5 ? Float64(rand(rng, 0:5)) : INF for _ in 1:m]
            bB[rand(rng, 1:m)] = 0.0
            Vb = ctg(d, bB[ψ])
            Vqb = ctg(floyd(Wq), bB)
            exact &= all(abs(Vqb[ψ[x]] - Vb[x]) <= tol || (isinf(Vqb[ψ[x]]) && isinf(Vb[x])) for x in 1:n)
            counts[:zero_defect_cases] += 1
        end
        counts[:exact_iff] += (zero == exact)

        # (5) local defect from paths inside B ∪ C
        bl = blocks_of(ψ, m)
        εloc = fill(INF, m, m)
        for B in 1:m, C in 1:m
            (B == C || isinf(Wq[B, C])) && continue
            S = vcat(bl[B], bl[C])
            dl = floyd(W[S, S])
            nb = length(bl[B])
            εloc[B, C] = maximum(max(minimum(dl[i, j] for j in nb+1:length(S)) - Wq[B, C], 0.0) for i in 1:nb)
        end
        Vloc = ctg(floyd(Wq .+ εloc), bq .+ εterm)
        counts[:local] += all(εloc .+ tol .>= ε) && all(V[x] <= Vloc[ψ[x]] + tol for x in 1:n)

        # (6) two-sided bound on the worst loss over one-block targets
        dq = floyd(Wq)
        dqε = floyd(Wq .+ ε)
        ok6 = true
        for B in 1:m
            loss = -INF
            upper = -INF
            for C in 1:m
                C == B && continue
                for x in bl[B]
                    dxC = minimum(d[x, y] for y in bl[C])
                    (isinf(dxC) && isinf(dq[B, C])) && continue
                    loss = max(loss, dxC - dq[B, C])
                end
                isinf(dq[B, C]) || (upper = max(upper, dqε[B, C] - dq[B, C]))
            end
            lower = maximum((ε[B, C] for C in 1:m if C != B && isfinite(Wq[B, C])); init = -INF)
            ok6 &= (loss >= lower - tol) && (loss <= upper + tol || isinf(upper) || loss == -INF)
        end
        counts[:two_sided] += ok6

        # (7) the Hausdorff weighting h is the least sound upper weighting
        h = fill(INF, m, m)
        for B in 1:m, C in 1:m
            B == C && continue
            h[B, C] = maximum(minimum(d[x, y] for y in bl[C]) for x in bl[B])
        end
        dh = floyd(h)
        ok7 = all(all(minimum(d[x, y] for y in bl[C]) <= dh[ψ[x], C] + tol for x in 1:n) for C in 1:m)
        ok7 &= all(h[B, C] <= Wq[B, C] + ε[B, C] + tol for B in 1:m, C in 1:m if B != C && isfinite(Wq[B, C]))
        finite = [(B, C) for B in 1:m, C in 1:m if B != C && isfinite(h[B, C])]
        if !isempty(finite)
            (B, C) = rand(rng, finite)
            h2 = copy(h); h2[B, C] -= 0.5
            dh2 = floyd(h2)
            ok7 &= any(minimum(d[x, y] for y in bl[C]) > dh2[B, C] + tol for x in bl[B])
        end
        counts[:optimal] += ok7
    end
    return counts
end

counts = check(MersenneTwister(2026), 3000)
for k in (:lower, :lift, :must, :exact_iff, :local, :two_sided, :optimal)
    println(rpad(string(k), 12), counts[k], " / 3000")
end
println("cases with zero defect on every arc: ", counts[:zero_defect_cases])
