# A tree-parallel triangular solve for CliqueTrees.Multifrontal chordal
# factorizations, built up in three layers — see dear_itay.md for the design
# history.
#
# 1. A per-node-addressed re-implementation of the serial multifrontal
#    `ldiv!` (`tree_forward_ldiv!`/`tree_backward_ldiv!`). CliqueTrees' own
#    forward/backward tree walk (divide.jl) passes the "update matrix"
#    contribution between clique-tree nodes through a single shared LIFO
#    stack (`Mptr`/`Mval`, indexed by a running position `ns`), which cannot
#    be split across independent tasks or processes. These reproduce the
#    same math, storing each node's contribution in a slot keyed by the
#    node's own id instead.
# 2. Tree partitioning (`partition_tree`) and per-worker data extraction
#    (`worker_factorization`) — splitting the clique tree into connected,
#    minimally-cut chunks and giving each one its own copied slice of the
#    factor data, ready to ship to a separate process.
# 3. An actual `Distributed.jl` backend (`distributed_tree_solve`) that runs
#    each chunk's forward/backward pass on its own worker process, with
#    `RemoteChannel`s carrying exactly the boundary corrections across cut
#    edges — real inter-process message passing, not shared memory.
module DistributedSolve

export tree_forward_ldiv!, tree_backward_ldiv!, TreeWorkspace,
    TreePartition, partition_tree, WorkerFactorization, worker_factorization,
    distributed_tree_solve

using ArgCheck: @argcheck
using LinearAlgebra
using LinearAlgebra: ldiv!, mul!
using Distributed: RemoteChannel, remotecall, workers
using Graphs: neighbors, nv
using CliqueTrees.Multifrontal: ChordalTriangular

"""    tree_forward_ldiv!(L, b::AbstractVector) -> b

Solve ``L y = b`` in place, where `L` is the lower-triangular factor of a
`CliqueTrees.Multifrontal` chordal factorization (`ChordalCholesky` or
`ChordalLDLt`). This is the *forward* half of a symmetric solve — the
bottom-up pass over the clique tree that condenses each supernode's local
unknowns and pushes a correction up to its parent (see
[`tree_backward_ldiv!`](@ref) for the top-down half).

Mathematically this reproduces `CliqueTrees.Multifrontal.ldiv!(L, b)`
exactly; the only difference is *how* the intermediate per-node corrections
are stored while walking the tree (indexed by node id here, versus a shared
stack internally in CliqueTrees).
"""
function tree_forward_ldiv!(L::ChordalTriangular, b::AbstractVector{T}) where {T}
    @argcheck length(b) == size(L, 1)
    S = L.S
    store = Vector{Union{Nothing, Dict{Int, T}}}(undef, nv(S.res))
    fill!(store, nothing)

    for j in _roots(S)
        _forward_visit!(b, store, S, L.Dval, L.Lval, j)
    end

    return b
end

function _forward_visit!(b, store, S, Dval, Lval, j)
    resj = collect(neighbors(S.res, j))
    sepj = collect(neighbors(S.sep, j))
    nn = length(resj)
    na = length(sepj)

    c1 = b[resj]
    f2 = zeros(eltype(c1), na)

    for i in neighbors(S.chd, j)
        _forward_visit!(b, store, S, Dval, Lval, i)
        _scatter_contribution!(c1, f2, resj, sepj, store[i])
    end

    Dp = S.Dptr[j]
    D11 = reshape(view(Dval, Dp:(Dp + nn * nn - 1)), nn, nn)
    y1 = LowerTriangular(D11) \ c1
    b[resj] = y1

    if na > 0
        Lp = S.Lptr[j]
        L21 = reshape(view(Lval, Lp:(Lp + na * nn - 1)), na, nn)
        m2 = f2 .- L21 * y1
        store[j] = Dict(gid => m2[k] for (k, gid) in enumerate(sepj))
    else
        store[j] = Dict{Int, eltype(c1)}()
    end

    return nothing
end

function _scatter_contribution!(c1, f2, resj, sepj, contribution)
    for (gid, val) in contribution
        k = findfirst(==(gid), resj)
        if k === nothing
            m = findfirst(==(gid), sepj)
            @argcheck m !== nothing "clique-tree contribution targets neither the residual nor the separator of its parent"
            f2[m] += val
        else
            c1[k] += val
        end
    end
    return nothing
end

"""    tree_backward_ldiv!(L, b::AbstractVector) -> b

Solve ``L^\\mathsf{T} y = b`` in place — the top-down half of a symmetric
solve. Unlike [`tree_forward_ldiv!`](@ref), this pass needs no per-node
store: each node only ever reads the already-finalized values at its own
separator (written by an ancestor earlier in the root-to-leaves traversal)
and writes only to its own residual, so the shared solution vector itself
carries every message.

Calling `tree_forward_ldiv!` then `tree_backward_ldiv!` on the same `b`
reproduces `M \\ b` for the original factored matrix `M`, matching
`CliqueTrees.Multifrontal.ldiv!(L, ldiv!(L', b))` bit for bit.
"""
function tree_backward_ldiv!(L::ChordalTriangular, b::AbstractVector)
    @argcheck length(b) == size(L, 1)
    S = L.S
    b0 = copy(b)

    for j in _roots(S)
        _backward_visit!(b, b0, S, L.Dval, L.Lval, j)
    end

    return b
end

function _backward_visit!(b, b0, S, Dval, Lval, j)
    resj = collect(neighbors(S.res, j))
    sepj = collect(neighbors(S.sep, j))
    nn = length(resj)
    na = length(sepj)

    c1 = b0[resj]

    if na > 0
        Lp = S.Lptr[j]
        L21 = reshape(view(Lval, Lp:(Lp + na * nn - 1)), na, nn)
        c1 = c1 .- L21' * b[sepj]
    end

    Dp = S.Dptr[j]
    D11 = reshape(view(Dval, Dp:(Dp + nn * nn - 1)), nn, nn)
    b[resj] = LowerTriangular(D11)' \ c1

    for i in neighbors(S.chd, j)
        _backward_visit!(b, b0, S, Dval, Lval, i)
    end

    return nothing
end

function _roots(S)
    return (j for j in 1:nv(S.res) if iszero(S.pnt[j]))
end

# ===== Preallocated workspace for repeated solves =====
#
# The recursive `tree_forward_ldiv!`/`tree_backward_ldiv!` above allocate on
# every call — a `Dict` per node for the forward store, `collect`ed neighbor
# lists at each visit, and a linear `findfirst` per scattered contribution.
# That is fine for a one-off solve, but the harmonic-extension use case
# re-solves against the *same* factorization every control step as the
# targets move (only the right-hand side `Bp` changes; the sheaf Laplacian,
# and hence its factor and clique tree, do not). A `TreeWorkspace` hoists all
# of that symbolic, allocation-heavy work out of the per-step hot loop: it is
# built once per factorization and then reused by the three-argument
# `tree_forward_ldiv!`/`tree_backward_ldiv!` methods, which sweep the tree
# flat (postorder / reverse-postorder) with in-place `ldiv!`/`mul!` and no
# per-solve allocation. On grid sheaf Laplacians this closes the gap to
# CliqueTrees' own monolithic `ldiv!` from ~6.7x to ~1.1x.

"""    TreeWorkspace(L, T=Float64) -> TreeWorkspace

Preallocated scratch for repeated flat solves against a fixed
`CliqueTrees.Multifrontal` chordal factor `L` with element type `T`. Build it
**once** per factorization and pass it to the three-argument
[`tree_forward_ldiv!`](@ref) / [`tree_backward_ldiv!`](@ref); those methods
then perform no per-solve allocation, which is what makes re-solving against a
cached factor every control step (as moving targets change only the
right-hand side) cheap.

The workspace caches the symbolic tree — each node's residual, separator and
child lists, a postorder in which every node follows all its descendants, an
offset table locating each node's separator contribution in a single flat
`store` buffer, and a precomputed scatter map `tgt` that replaces the
`findfirst` search of the recursive version — together with the reusable
`store`, `c1` and `f2` numeric buffers.
"""
struct TreeWorkspace{T}
    res::Vector{Vector{Int}}
    sep::Vector{Vector{Int}}
    chd::Vector{Vector{Int}}
    order::Vector{Int}
    off::Vector{Int}
    tgt::Vector{Vector{Int}}
    store::Vector{T}
    c1::Vector{T}
    f2::Vector{T}
end

# Postorder over the (possibly multi-rooted) clique forest: every node is
# emitted only after all of its descendants, so children precede parents and
# the flat forward sweep can process nodes in a single pass. Iterative to
# avoid stack overflow on deep (path-like) trees.
function _postorder(S)
    order = Int[]
    sizehint!(order, nv(S.res))
    stack = Tuple{Int, Bool}[]
    for r in _roots(S)
        push!(stack, (r, false))
        while !isempty(stack)
            j, done = pop!(stack)
            if done
                push!(order, j)
            else
                push!(stack, (j, true))
                for c in neighbors(S.chd, j)
                    push!(stack, (c, false))
                end
            end
        end
    end
    return order
end

function TreeWorkspace(L::ChordalTriangular, ::Type{T}=Float64) where {T}
    S = L.S
    n = nv(S.res)
    res = [collect(neighbors(S.res, j)) for j in 1:n]
    sep = [collect(neighbors(S.sep, j)) for j in 1:n]
    chd = [collect(neighbors(S.chd, j)) for j in 1:n]
    order = _postorder(S)

    off = zeros(Int, n + 1)
    for j in 1:n
        off[j + 1] = off[j] + length(sep[j])
    end

    # tgt[c][k]: where child c's k-th separator entry lands in its parent j —
    # +m into the parent's residual slot m, -m into the parent's separator
    # slot m. This is the scatter that the recursive version recomputes with
    # `findfirst` on every contribution.
    tgt = [Int[] for _ in 1:n]
    for j in 1:n
        pr = Dict(g => k for (k, g) in enumerate(res[j]))
        ps = Dict(g => k for (k, g) in enumerate(sep[j]))
        for c in chd[j]
            tgt[c] = [haskey(pr, g) ? pr[g] : -ps[g] for g in sep[c]]
        end
    end

    mx = isempty(res) ? 0 : maximum(length, res)
    ma = isempty(sep) ? 0 : maximum(length, sep)
    return TreeWorkspace{T}(res, sep, chd, order, off, tgt,
        zeros(T, off[end]), zeros(T, mx), zeros(T, ma))
end

"""    tree_forward_ldiv!(L, b, ws::TreeWorkspace) -> b

Allocation-free forward solve `L y = b` in place, reusing the preallocated
`ws` (see [`TreeWorkspace`](@ref)). Mathematically identical to the
two-argument [`tree_forward_ldiv!`](@ref); it sweeps the clique tree flat in
postorder rather than by recursion, and scatters each node's contribution
through the precomputed `ws.tgt` map instead of searching for it.
"""
function tree_forward_ldiv!(L::ChordalTriangular, b::AbstractVector{T}, ws::TreeWorkspace{T}) where {T}
    @argcheck length(b) == size(L, 1)
    S = L.S
    Dval = L.Dval
    Lval = L.Lval
    @inbounds for j in ws.order
        resj = ws.res[j]
        nn = length(resj)
        na = length(ws.sep[j])
        c1 = view(ws.c1, 1:nn)
        f2 = view(ws.f2, 1:na)
        for k in 1:nn
            c1[k] = b[resj[k]]
        end
        for k in 1:na
            f2[k] = zero(T)
        end
        for c in ws.chd[j]
            m = view(ws.store, (ws.off[c] + 1):ws.off[c + 1])
            tg = ws.tgt[c]
            for k in eachindex(tg)
                t = tg[k]
                t > 0 ? (c1[t] += m[k]) : (f2[-t] += m[k])
            end
        end
        Dp = S.Dptr[j]
        D11 = reshape(view(Dval, Dp:(Dp + nn * nn - 1)), nn, nn)
        ldiv!(LowerTriangular(D11), c1)
        for k in 1:nn
            b[resj[k]] = c1[k]
        end
        if na > 0
            Lp = S.Lptr[j]
            L21 = reshape(view(Lval, Lp:(Lp + na * nn - 1)), na, nn)
            m = view(ws.store, (ws.off[j] + 1):ws.off[j + 1])
            copyto!(m, f2)
            mul!(m, L21, c1, -one(T), one(T))
        end
    end
    return b
end

"""    tree_backward_ldiv!(L, b, ws::TreeWorkspace) -> b

Allocation-free backward solve `Lᵀ y = b` in place, reusing `ws` (see
[`TreeWorkspace`](@ref)). Mathematically identical to the two-argument
[`tree_backward_ldiv!`](@ref). It sweeps the tree in reverse postorder
(every node before its children), so a node's separator entries — which
reference only already-finalized ancestor residuals — are read straight from
`b`; no per-node store is needed, as in the recursive version.
"""
function tree_backward_ldiv!(L::ChordalTriangular, b::AbstractVector{T}, ws::TreeWorkspace{T}) where {T}
    @argcheck length(b) == size(L, 1)
    S = L.S
    Dval = L.Dval
    Lval = L.Lval
    @inbounds for idx in length(ws.order):-1:1
        j = ws.order[idx]
        resj = ws.res[j]
        sepj = ws.sep[j]
        nn = length(resj)
        na = length(sepj)
        c1 = view(ws.c1, 1:nn)
        for k in 1:nn
            c1[k] = b[resj[k]]
        end
        if na > 0
            f2 = view(ws.f2, 1:na)
            for k in 1:na
                f2[k] = b[sepj[k]]
            end
            Lp = S.Lptr[j]
            L21 = reshape(view(Lval, Lp:(Lp + na * nn - 1)), na, nn)
            mul!(c1, L21', f2, -one(T), one(T))
        end
        Dp = S.Dptr[j]
        D11 = reshape(view(Dval, Dp:(Dp + nn * nn - 1)), nn, nn)
        ldiv!(LowerTriangular(D11)', c1)
        for k in 1:nn
            b[resj[k]] = c1[k]
        end
    end
    return b
end

# ===== Tree partitioning (milestone 4) =====
#
# Splits the clique tree into disjoint, connected, rooted subtrees ("chunks")
# to be assigned one per worker/process — the piece that turns the serial
# per-node kernel above into something that can actually be distributed.
# See dear_itay.md §5-6 for the design rationale.

"""    TreePartition

The result of [`partition_tree`](@ref): an assignment of every clique-tree
node to a worker, chosen so that every worker's nodes form one connected
rooted subtree and the number of parent/child edges crossing between workers
is the minimum possible for the resulting worker count (`length(chunks) - 1`
for a tree with a single root).

# Fields
- `owner::Vector{Int}`: `owner[j]` is the worker id (1-based) that node `j`
  is assigned to.
- `chunks::Vector{Vector{Int}}`: `chunks[w]` lists the node ids owned by
  worker `w`. The actual worker count is `length(chunks)`, which may exceed
  the requested `target_workers` — see [`partition_tree`](@ref).
- `cut_edges::Vector{Pair{Int,Int}}`: `child => parent` pairs whose two
  endpoints are owned by different workers. Only these edges require actual
  cross-worker communication; every other edge in the tree is a same-worker
  memory read.
"""
struct TreePartition
    owner::Vector{Int}
    chunks::Vector{Vector{Int}}
    cut_edges::Vector{Pair{Int, Int}}
end

_bag_size(S, j) = length(neighbors(S.res, j)) + length(neighbors(S.sep, j))

"""    partition_tree(L, target_workers::Integer; weight=nothing) -> TreePartition

Partition the clique tree of `L` (a `CliqueTrees.Multifrontal` chordal
factorization's triangular factor) into roughly `target_workers` connected,
rooted subtrees, by a single bottom-up sweep in the tree's own postorder:
each node accumulates the weight of itself plus every child not yet claimed
by an earlier chunk, and is sealed off (together with everything it just
accumulated) into a new chunk as soon as that running total reaches
`total_weight / target_workers`. Every root is always sealed regardless of
weight, since there is nothing above it to merge into.

`weight` defaults to each node's bag size (`|res(j)| + |sep(j)|`), a
standard proxy for triangular-solve cost. Pass a custom `(S, j) -> Int`
function to use a different cost model.

Because sealing only ever happens along an existing parent/child edge, every
chunk is guaranteed connected and the resulting number of cross-worker edges
(`cut_edges`) is provably minimal for the worker count obtained — but that
worker count is only approximately `target_workers` (a next-fit-style
approximation), and can be far short of it: a wide, shallow tree (e.g. a
star sheaf) has nowhere to place an interior cut, so it collapses toward one
large chunk regardless of `target_workers`, only fragmenting into many tiny
single-node chunks once the threshold drops below a single node's weight.
This is not a bug — it is the tree's own separator width limiting available
parallelism, made explicit (see dear_itay.md).
"""
function partition_tree(L::ChordalTriangular, target_workers::Integer; weight=nothing)
    @argcheck target_workers >= 1
    S = L.S
    n = nv(S.res)
    w(j) = weight === nothing ? _bag_size(S, j) : weight(S, j)
    total = sum(w, 1:n)
    threshold = max(1, ceil(Int, total / target_workers))

    owner = zeros(Int, n)
    pending_weight = zeros(Int, n)
    pending_members = Vector{Vector{Int}}(undef, n)
    chunks = Vector{Int}[]

    for j in 1:n
        members = Int[j]
        pw = w(j)
        for c in neighbors(S.chd, j)
            if owner[c] == 0
                pw += pending_weight[c]
                append!(members, pending_members[c])
            end
        end

        if pw >= threshold || iszero(S.pnt[j])
            push!(chunks, members)
            wid = length(chunks)
            for m in members
                owner[m] = wid
            end
        else
            pending_weight[j] = pw
            pending_members[j] = members
        end
    end

    cuts = Pair{Int, Int}[]
    for j in 1:n
        p = S.pnt[j]
        if !iszero(p) && owner[j] != owner[p]
            push!(cuts, j => p)
        end
    end

    return TreePartition(owner, chunks, cuts)
end

"""    WorkerFactorization

A self-contained slice of a chordal factorization, holding only the data one
worker needs to run its portion of [`tree_forward_ldiv!`](@ref) /
[`tree_backward_ldiv!`](@ref): exactly the `Dval`/`Lval` entries for the
nodes it owns (copied out, not aliased into the shared array — this is the
data a separate process, or drone, would actually be shipped), reindexed by
a fresh local `Dptr`/`Lptr` starting at 1.

# Fields
- `ids::Vector{Int}`: the (global) node ids this worker owns, ascending.
- `Dval::Vector`, `Lval::Vector`: this worker's own copy of the diagonal and
  off-diagonal factor blocks, concatenated in `ids` order.
- `Dptr::Vector{Int}`, `Lptr::Vector{Int}`: local offsets into the arrays
  above, `Dptr[k]:Dptr[k+1]-1` giving node `ids[k]`'s diagonal block.
- `res::Dict{Int,Vector{Int}}`, `sep::Dict{Int,Vector{Int}}`,
  `chd::Dict{Int,Vector{Int}}`, `pnt::Dict{Int,Int}`: the symbolic tree
  structure restricted to owned nodes (global vertex/node ids throughout,
  matching `tree_forward_ldiv!`/`tree_backward_ldiv!`'s own indexing).
- `boundary_up::Vector{Int}`: owned nodes whose parent belongs to a
  *different* worker — the forward pass must send that node's contribution
  out, and the backward pass must receive a value for it, over the network.
- `boundary_down::Vector{Int}`: owned nodes with at least one child
  belonging to a different worker — the mirror image of `boundary_up`.
"""
struct WorkerFactorization{TD, TL}
    ids::Vector{Int}
    Dval::Vector{TD}
    Lval::Vector{TL}
    Dptr::Vector{Int}
    Lptr::Vector{Int}
    res::Dict{Int, Vector{Int}}
    sep::Dict{Int, Vector{Int}}
    chd::Dict{Int, Vector{Int}}
    pnt::Dict{Int, Int}
    boundary_up::Vector{Int}
    boundary_down::Vector{Int}
end

"""    worker_factorization(L, partition::TreePartition, worker::Integer) -> WorkerFactorization

Extract the [`WorkerFactorization`](@ref) for one worker of a
[`partition_tree`](@ref) result — the exact packet that worker (or drone)
needs, and nothing else: its own copied slice of `L`'s factor data plus the
symbolic metadata for its owned nodes.
"""
function worker_factorization(L::ChordalTriangular, partition::TreePartition, worker::Integer)
    @argcheck 1 <= worker <= length(partition.chunks)
    S = L.S
    ids = sort(partition.chunks[worker])
    idset = Set(ids)

    Dval = eltype(L.Dval)[]
    Lval = eltype(L.Lval)[]
    Dptr = Vector{Int}(undef, length(ids) + 1)
    Lptr = Vector{Int}(undef, length(ids) + 1)
    Dptr[1] = 1
    Lptr[1] = 1

    res = Dict{Int, Vector{Int}}()
    sep = Dict{Int, Vector{Int}}()
    chd = Dict{Int, Vector{Int}}()
    pnt = Dict{Int, Int}()

    for (k, j) in enumerate(ids)
        Dp, Dp2 = S.Dptr[j], S.Dptr[j + 1] - 1
        append!(Dval, view(L.Dval, Dp:Dp2))
        Dptr[k + 1] = Dptr[k] + (Dp2 - Dp + 1)

        Lp, Lp2 = S.Lptr[j], S.Lptr[j + 1] - 1
        append!(Lval, view(L.Lval, Lp:Lp2))
        Lptr[k + 1] = Lptr[k] + (Lp2 - Lp + 1)

        res[j] = collect(neighbors(S.res, j))
        sep[j] = collect(neighbors(S.sep, j))
        chd[j] = collect(neighbors(S.chd, j))
        pnt[j] = S.pnt[j]
    end

    boundary_up = [j for j in ids if !iszero(pnt[j]) && !(pnt[j] in idset)]
    boundary_down = [j for j in ids if any(c -> !(c in idset), chd[j])]

    return WorkerFactorization(ids, Dval, Lval, Dptr, Lptr, res, sep, chd, pnt, boundary_up, boundary_down)
end

# ===== Distributed backend (milestone 5) =====
#
# Each worker runs its own owned subtree with the exact same math as
# tree_forward_ldiv!/tree_backward_ldiv!, branching only on whether a given
# child/parent is owned locally (a plain recursive call) or by a different
# worker (a RemoteChannel message). One channel per cut edge, keyed by the
# cut edge's child-side node id — the same key is used for both the forward
# message (child → parent's worker) and the backward message (parent's
# worker → child), since each channel only ever needs to carry one value at
# a time and the two phases never overlap.
#
# Two correctness subtleties surfaced only by testing on general (non-star,
# non-path) graphs, where separators can be wider than 1 and reach past the
# immediate parent — star/path never exercise this, since their separators
# never exceed one element:
#   1. A worker's own `solved` dict must absorb every boundary value it
#      receives, not just use it locally — a deeper node in the same chunk
#      can depend on that same pass-through vertex.
#   2. A backward `put!` to an external child must send the sender's whole
#      bag (residual ∪ separator), not just its own residual — the child's
#      separator can reach past the sender into the sender's own separator.
# See test/network_sheaves/DistributedSolve.jl's general-topology tests.

_local_root(wd, idset) = wd.ids[findfirst(j -> iszero(wd.pnt[j]) || !(wd.pnt[j] in idset), wd.ids)]
_idpos(wd) = Dict(j => k for (k, j) in enumerate(wd.ids))

"""    worker_forward_pass(wd::WorkerFactorization, b_owned::Dict{Int,T}, fwd_channels) -> Dict{Int,T}

Run the forward pass (see [`tree_forward_ldiv!`](@ref)) restricted to one
worker's owned subtree. `b_owned` supplies the right-hand-side entries for
this worker's own residual vertices only. `fwd_channels` is a
`node_id => RemoteChannel` map covering every cut edge in the tree (workers
only ever touch the entries for their own boundary nodes); a node whose
child belongs to another worker blocks on `take!` from that child's
channel, and a node whose own parent belongs to another worker `put!`s its
correction there once solved. Returns this worker's solved values, keyed by
global vertex id.
"""
function worker_forward_pass(wd::WorkerFactorization, b_owned::Dict{Int, T}, fwd_channels) where {T}
    idset = Set(wd.ids)
    idpos = _idpos(wd)
    solved = Dict{Int, T}()
    store = Dict{Int, Dict{Int, T}}()

    function visit(j)
        resj = wd.res[j]
        sepj = wd.sep[j]
        nn = length(resj)
        na = length(sepj)

        c1 = [b_owned[v] for v in resj]
        f2 = zeros(T, na)

        for c in wd.chd[j]
            contribution = c in idset ? (visit(c); store[c]) : take!(fwd_channels[c])
            _scatter_contribution!(c1, f2, resj, sepj, contribution)
        end

        k = idpos[j]
        Dp = wd.Dptr[k]
        D11 = reshape(view(wd.Dval, Dp:(Dp + nn * nn - 1)), nn, nn)
        y1 = LowerTriangular(D11) \ c1
        for (idx, gid) in enumerate(resj)
            solved[gid] = y1[idx]
        end

        if na > 0
            Lp = wd.Lptr[k]
            L21 = reshape(view(wd.Lval, Lp:(Lp + na * nn - 1)), na, nn)
            contribution = Dict(gid => (f2 .- L21 * y1)[idx] for (idx, gid) in enumerate(sepj))
            store[j] = contribution
            wd.pnt[j] != 0 && !(wd.pnt[j] in idset) && put!(fwd_channels[j], contribution)
        else
            store[j] = Dict{Int, T}()
        end

        return nothing
    end

    visit(_local_root(wd, idset))
    return solved
end

"""    worker_backward_pass(wd::WorkerFactorization, y1_owned::Dict{Int,T}, bwd_channels) -> Dict{Int,T}

The backward-pass counterpart of [`worker_forward_pass`](@ref) (see
[`tree_backward_ldiv!`](@ref)). `y1_owned` is this worker's own result from
the forward pass. A node whose parent belongs to another worker blocks on
`take!` from `bwd_channels[node]`; a node with a child owned elsewhere
`put!`s its own *whole bag* (residual and separator both — a child's
separator can reach past its immediate parent into the parent's own
separator, so the residual alone is not always enough) to that child's
channel once finished.
"""
function worker_backward_pass(wd::WorkerFactorization, y1_owned::Dict{Int, T}, bwd_channels) where {T}
    idset = Set(wd.ids)
    idpos = _idpos(wd)
    solved = Dict{Int, T}()

    function visit(j)
        resj = wd.res[j]
        sepj = wd.sep[j]
        nn = length(resj)
        na = length(sepj)

        c1 = [y1_owned[v] for v in resj]

        if na > 0
            # A separator can reach past the immediate parent (fill-in can
            # make sep(j) reference a further ancestor), so a value received
            # here from another worker must be merged into `solved`, not
            # just used locally — deeper nodes in this same chunk may need
            # that same pass-through vertex too, and this is their only
            # route to it (the serial version gets this for free by reading
            # off one shared global vector instead of a per-worker dict).
            sepvals_dict = wd.pnt[j] in idset ? Dict(gid => solved[gid] for gid in sepj) : take!(bwd_channels[j])
            merge!(solved, sepvals_dict)
            sepvals = [sepvals_dict[gid] for gid in sepj]
            k = idpos[j]
            Lp = wd.Lptr[k]
            L21 = reshape(view(wd.Lval, Lp:(Lp + na * nn - 1)), na, nn)
            c1 = c1 .- L21' * sepvals
        end

        k = idpos[j]
        Dp = wd.Dptr[k]
        D11 = reshape(view(wd.Dval, Dp:(Dp + nn * nn - 1)), nn, nn)
        y2 = LowerTriangular(D11)' \ c1
        for (idx, gid) in enumerate(resj)
            solved[gid] = y2[idx]
        end

        for c in wd.chd[j]
            if c in idset
                visit(c)
            else
                # An external child's own separator can reach past j into
                # j's separator too (the same pass-through as above), so it
                # is j's *whole bag* — res(j) ∪ sep(j), both already in
                # `solved` at this point — that's guaranteed to cover
                # whatever the child needs, not just j's own residual.
                bag = vcat(resj, sepj)
                put!(bwd_channels[c], Dict(gid => solved[gid] for gid in bag))
            end
        end

        return nothing
    end

    visit(_local_root(wd, idset))
    return solved
end

"""    distributed_tree_solve(L, b::AbstractVector, target_workers::Integer; pids=workers()) -> Vector

Solve against `L` (a `CliqueTrees.Multifrontal` chordal factorization's
triangular factor) genuinely distributed across separate worker processes:
[`partition_tree`](@ref) splits the clique tree into `target_workers`-ish
chunks (see its docstring — the actual chunk count is an approximation, not
exact), each worker gets its own [`worker_factorization`](@ref) slice via
`remotecall`, and [`worker_forward_pass`](@ref)/[`worker_backward_pass`](@ref)
run on each worker with `RemoteChannel` messages carrying exactly the
cross-worker corrections. Errors if `pids` doesn't have enough entries for
the resulting chunk count.

This function does not call `addprocs` itself — set up and provision
worker processes first (`addprocs(N; exeflags="--project=...")`, then
`@everywhere using CellularSheaves` so the worker-side functions above are
loadable on every process), matching how process/cluster provisioning is
the caller's decision, not this package's.
"""
function distributed_tree_solve(L::ChordalTriangular, b::AbstractVector{T}, target_workers::Integer; pids = workers()) where {T}
    @argcheck length(b) == size(L, 1)
    partition = partition_tree(L, target_workers)
    nw = length(partition.chunks)
    @argcheck nw <= length(pids) "need at least $nw worker processes for $nw chunks, got $(length(pids))"

    cut_ids = [child for (child, parent) in partition.cut_edges]
    fwd_channels = Dict(id => RemoteChannel(() -> Channel{Dict{Int, T}}(1)) for id in cut_ids)
    bwd_channels = Dict(id => RemoteChannel(() -> Channel{Dict{Int, T}}(1)) for id in cut_ids)

    wds = [worker_factorization(L, partition, w) for w in 1:nw]
    b_owned = [Dict(v => b[v] for (j, verts) in wd.res for v in verts) for wd in wds]

    fwd_results = fetch.([remotecall(worker_forward_pass, pids[w], wds[w], b_owned[w], fwd_channels) for w in 1:nw])
    bwd_results = fetch.([remotecall(worker_backward_pass, pids[w], wds[w], fwd_results[w], bwd_channels) for w in 1:nw])

    y = Vector{T}(undef, size(L, 1))
    for d in bwd_results, (gid, val) in d
        y[gid] = val
    end
    return y
end

end
