# Data behind the hand-drawn structural figures in the Distributed Sheaf Solve
# feature guide (docs/src/features/distributed_sheaf_solve/). The figures are
# authored as SVG, but every supernode, tree edge, partition assignment, and
# message slot they depict is computed here from the same 12-agent two-ring
# formation as the distributed harmonic tracking worked example. Run this to
# re-derive (or audit) everything the figures show:
#
#   julia --project=docs docs/scripts/distributed_solve_figdata.jl

get!(ENV, "GKSwstype", "100")

using CellularSheaves
using CellularSheaves.NetworkSheaves.DistributedSolve
using CliqueTrees.Multifrontal
using LinearAlgebra
using SparseArrays
using Printf
import CellularSheaves.NetworkSheaves.EuclideanSheaves: _harmonic_extension_restricted_laplacian
import Graphs: neighbors, nv

# ---- the formation (identical to the worked example) ----------------------

const NA, NT = 12, 2
const TV1, TV2 = NA + 1, NA + 2
const D = 3
const I3 = Matrix{Float64}(I, D, D)

consensus_edges = [(1,2),(2,3),(3,4),(4,5),(5,6),(6,1),
                   (7,8),(8,9),(9,10),(10,11),(11,12),(12,7),
                   (1,7)]

sheaf = EuclideanSheaf{Float64}(fill(D, NA + NT))
for (i, j) in consensus_edges
    add_sheaf_edge!(sheaf, i, j, I3, I3)
end
for i in 1:6
    add_sheaf_edge!(sheaf, i, TV1, I3, I3)
end
for i in 7:12
    add_sheaf_edge!(sheaf, i, TV2, I3, I3)
end

boundary0 = Dict(TV1 => zeros(D), TV2 => zeros(D))
_, _, H, B = _harmonic_extension_restricted_laplacian(sheaf, boundary0)

F = cholesky!(ChordalCholesky(sparse(H)), NoPivot())
L = F.L
S = L.S

# ---- supernodes, mapped back to agents -------------------------------------

# S.res/S.sep index PERMUTED rows of H; F.P maps back to original DOFs, and
# original DOF d belongs to agent fld1(d, D).
perm = invperm(F.P.perm)   # permuted row -> original row
agent_of(prow) = fld1(perm[prow], D)

ns = nv(S.res)
println("=== supernodes (j: res_agents | sep_agents | nn, na) ===")
for j in 1:ns
    resj = collect(neighbors(S.res, j))
    sepj = collect(neighbors(S.sep, j))
    ra = sort(unique(agent_of.(resj)))
    sa = sort(unique(agent_of.(sepj)))
    @printf("%2d: res %-12s sep %-12s  nn=%2d na=%2d  parent=%d\n",
        j, string(ra), string(sa), length(resj), length(sepj), S.pnt[j])
end

# ---- partition --------------------------------------------------------------

p = partition_tree(L, 3)
println("\n=== partition into $(length(p.chunks)) workers ===")
for (w, chunk) in enumerate(p.chunks)
    println("worker $w owns supernodes $(sort(chunk))")
end
println("cut edges (child => parent): ", p.cut_edges)

# ---- forward message schedule (greedy half-duplex, sibling-serialized) -----

println("\n=== forward gather schedule (slot: child -> parent) ===")
ready = Dict{Int, Int}()   # slot at which node j's message is ready to send
busy_until = zeros(Int, ns)
order = Int[]              # postorder
function po(j)
    for c in neighbors(S.chd, j); po(c); end
    push!(order, j)
end
for j in 1:ns
    iszero(S.pnt[j]) && po(j)
end
sched = Tuple{Int, Int, Int}[]   # (slot, child, parent)
for j in order
    ch = collect(neighbors(S.chd, j))
    t = 0
    for c in sort(ch; by = c -> ready[c])
        slot = max(ready[c], t) + 1
        push!(sched, (slot, c, j))
        t = slot
    end
    ready[j] = t
end
for (slot, c, par) in sort(sched)
    @printf("slot %2d: node %2d -> node %2d\n", slot, c, par)
end
fwd_makespan = isempty(sched) ? 0 : maximum(first, sched)
@printf("forward makespan = %d slots (round trip = %d)\n", fwd_makespan, 2fwd_makespan)

# ---- factor block sizes for the anatomy figure ------------------------------

println("\n=== factor block bytes ===")
@printf("total Dval %d entries, Lval %d entries, factor %.2f KB\n",
    length(L.Dval), length(L.Lval), (length(L.Dval) + length(L.Lval)) * 8 / 1e3)
