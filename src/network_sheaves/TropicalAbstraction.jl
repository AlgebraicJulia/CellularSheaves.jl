"""
    TropicalAbstraction

Coarse grids for deterministic optimal control in the min-plus (tropical)
semiring ``(ℝ ∪ \\{∞\\}, \\min, +)``, where they give exact solutions through
bounds rather than through residual corrections.

**Deterministic optimal control as a min-plus fixed point.** On a finite state
graph whose edge ``x → y`` costs ``W_{xy} ≥ 0`` (a control taking state ``x`` to
``y``), the cost-to-go to a target set ``Γ`` solves the Bellman equation

```math
V(x) = \\min_y \\big( W_{xy} + V(y) \\big), \\qquad V|_Γ = 0,
```

which is linear in min-plus: ``V = W ⊗ V ⊕ b`` with ``b = 0`` on ``Γ`` and ``∞``
elsewhere, so ``V = W^* ⊗ b`` ([`cost_to_go`](@ref), by Dijkstra's algorithm).

**Pushforward and pullback without subtraction.** For a map ``ψ`` sending each
state to a block (an aggregation, e.g.
[`aggregation_homomorphism`](@ref CellularSheaves.NetworkSheaves.GridMultigrid.aggregation_homomorphism)),
the pullback ``E = ψ^*`` copies a block's value onto its states
([`tropical_pullback`](@ref)) and its transpose ``E^T`` takes the min-plus sum,
the minimum, over each fibre ([`tropical_pushforward`](@ref)). Because ``\\min``
is idempotent, ``E^T ⊗ E = I`` with no normalization: the block average of the
real case becomes the block minimum. The Galerkin operator ``W_c = E^T ⊗ W ⊗ E``
is the *quotient graph* ([`quotient_graph`](@ref)): one vertex per block, the
cheapest edge between two blocks, edges inside a block contracted.

**Bounds instead of corrections.** Min-plus has no subtraction, so there is no
residual to correct. Instead the coarse problem is a *relaxation*: every path of
the state graph maps to a walk of the quotient graph that costs no more, so
``E ⊗ V_c ≤ V`` for the coarse cost-to-go ``V_c`` to ``ψ(Γ)``. The pulled-back
coarse solution is moreover *consistent* (``h(x) ≤ W_{xy} + h(y)`` on every
edge), so it is an admissible A* heuristic, and A* guided by it returns the exact
fine cost while expanding fewer states than Dijkstra's algorithm
([`abstraction_heuristic`](@ref), [`cost_to_go`](@ref) with a `heuristic`).
A tighter abstraction (smaller blocks) gives a tighter bound, fewer fine
expansions and a larger coarse problem. This is Hierarchical A* (Holte, Perez,
Zimmer & MacDonald 1996) and the abstraction heuristics of classical planning
(pattern databases, Culberson & Schaeffer 1998; merge-and-shrink, Helmert,
Haslum & Hoffmann 2007): in sheaf terms, the pushforward along ``ψ``.

The same pushforward is the Galerkin coarse operator of the multigrid solvers in
the real semiring (see
[`GridMultigrid`](@ref CellularSheaves.NetworkSheaves.GridMultigrid)), where the
coarse solution is a correction rather than a bound.
"""
module TropicalAbstraction

export tropical_pullback, tropical_pushforward, quotient_graph, cost_to_go, abstraction_heuristic

using ArgCheck: @argcheck
using SparseArrays
using ..GraphHomomorphisms: GraphHomomorphism

"""
    tropical_pullback(ψ, vc) -> Vector

The pullback ``E ⊗ v_c`` along `ψ`: each state takes the value of its block.
"""
function tropical_pullback(ψ::GraphHomomorphism, vc::AbstractVector)
    @argcheck length(vc) == ψ.n_target "one value per block"
    return vc[ψ.vertex_map]
end

"""
    tropical_pushforward(ψ, v) -> Vector

The min-plus pushforward ``E^T ⊗ v`` along `ψ`: each block takes the minimum of
`v` over its states (``∞`` for an empty block). It is a left inverse of
[`tropical_pullback`](@ref) for a surjective ``ψ``.
"""
function tropical_pushforward(ψ::GraphHomomorphism, v::AbstractVector{T}) where {T<:Real}
    @argcheck length(v) == length(ψ.vertex_map) "one value per state"
    vc = fill(typemax(T), ψ.n_target)
    for (x, b) in enumerate(ψ.vertex_map)
        vc[b] = min(vc[b], v[x])
    end
    return vc
end

"""
    quotient_graph(W, ψ) -> SparseMatrixCSC

The min-plus Galerkin operator ``E^T ⊗ W ⊗ E`` of the edge costs `W` (`W[x, y]`
the cost of the edge ``x → y``, absent entries ``∞``) along `ψ`: an edge between
two blocks costs the cheapest edge between their states. Edges inside a block
(self-loops of nonnegative cost, which never shorten a path) are dropped.
"""
function quotient_graph(W::SparseMatrixCSC{T}, ψ::GraphHomomorphism) where {T<:Real}
    _check_costs(W)
    @argcheck size(W, 1) == length(ψ.vertex_map) "ψ must map the states of W"
    rows, cols, vals = findnz(W)
    keep = [ψ.vertex_map[x] != ψ.vertex_map[y] for (x, y) in zip(rows, cols)]
    nc = ψ.n_target
    return sparse(ψ.vertex_map[rows[keep]], ψ.vertex_map[cols[keep]], vals[keep], nc, nc, min)
end

function _check_costs(W)
    @argcheck size(W, 1) == size(W, 2) "W must be square"
    @argcheck all(>=(0), nonzeros(W)) "edge costs must be nonnegative"
end

# A binary min-heap of (key, state) pairs with lazy deletion.
function _push!(heap, item)
    push!(heap, item)
    i = length(heap)
    while i > 1
        p = i >> 1
        heap[p] <= heap[i] && break
        heap[p], heap[i] = heap[i], heap[p]
        i = p
    end
end

function _pop!(heap)
    top = heap[1]
    last = pop!(heap)
    if !isempty(heap)
        heap[1] = last
        i, n = 1, length(heap)
        while true
            l = 2i
            l > n && break
            c = (l < n && heap[l + 1] < heap[l]) ? l + 1 : l
            heap[i] <= heap[c] && break
            heap[i], heap[c] = heap[c], heap[i]
            i = c
        end
    end
    return top
end

"""
    cost_to_go(W, targets; source = nothing, heuristic = nothing)
        -> (values, cost, expanded)

The least solution ``V = W^* ⊗ b`` of the min-plus Bellman equation
``V(x) = \\min_y (W_{xy} + V(y))``, ``V = 0`` on `targets`, by Dijkstra's algorithm
run backwards from the targets over the in-edges of each state (the columns of
`W`). `values[x]` is the cost from `x` to the targets (``∞`` if unreachable),
exact for every state the search settled.

With a `source`, the search stops once the source is settled and `cost` is its
exact cost-to-go. A `heuristic` (with a `source`) turns the search into A*: it
must be a consistent lower bound on the cost *from the source* to each state,
``h(y) ≤ h(x) + W_{xy}`` on every edge and ``h(\\text{source}) = 0``, such as
[`abstraction_heuristic`](@ref); `cost` is still exact. `expanded` counts the
states settled.
"""
function cost_to_go(W::SparseMatrixCSC{T}, targets; source::Union{Nothing,Integer} = nothing,
        heuristic::Union{Nothing,AbstractVector} = nothing) where {T<:Real}
    _check_costs(W)
    N = size(W, 1)
    @argcheck !isempty(targets) && all(t -> 1 <= t <= N, targets) "targets must be states of W"
    @argcheck isnothing(source) || 1 <= source <= N "source must be a state of W"
    @argcheck isnothing(heuristic) || (!isnothing(source) && length(heuristic) == N) "a heuristic needs a source and one value per state"
    S = float(T)
    h = isnothing(heuristic) ? nothing : S.(heuristic)
    g = fill(typemax(S), N)
    settled = falses(N)
    heap = Tuple{S,Int}[]
    for t in targets
        g[t] = zero(S)
        _push!(heap, (isnothing(h) ? zero(S) : h[t], t))
    end
    rows, vals = rowvals(W), nonzeros(W)
    expanded = 0
    while !isempty(heap)
        _, y = _pop!(heap)
        settled[y] && continue
        settled[y] = true
        expanded += 1
        y == source && break
        for k in nzrange(W, y)
            x = rows[k]
            settled[x] && continue
            d = g[y] + vals[k]
            if d < g[x]
                g[x] = d
                key = isnothing(h) ? d : d + h[x]
                isinf(key) || _push!(heap, (key, x))   # h = ∞: the source cannot reach x
            end
        end
    end
    cost = isnothing(source) ? nothing : g[source]
    return (values = g, cost = cost, expanded = expanded)
end

"""
    abstraction_heuristic(W, ψ, source) -> Vector

The A* heuristic of the abstraction `ψ` for searches from `source`: the cost in
the quotient graph ([`quotient_graph`](@ref)) from ``ψ(\\text{source})`` to each
block, pulled back to the states. It is a consistent lower bound on the cost from
`source` to each state, so [`cost_to_go`](@ref) with it is exact; states whose
block the source cannot reach get ``∞`` and are never expanded.
"""
function abstraction_heuristic(W::SparseMatrixCSC, ψ::GraphHomomorphism, source::Integer)
    Wc = quotient_graph(W, ψ)
    from_source = cost_to_go(sparse(transpose(Wc)), [ψ.vertex_map[source]]).values
    return tropical_pullback(ψ, from_source)
end

end # module
