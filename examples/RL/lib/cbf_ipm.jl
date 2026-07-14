# cbf_ipm.jl — decentralized control-barrier-function (CBF) safety filter, solved with the sheaf IPM.
#
# Safety is encoded by the pairwise barrier
#       h_ij(x) = ‖x_i − x_j‖² − d_safe²  ≥ 0        (agents at least d_safe apart).
# For single-integrator agents ẋ_i = u_i, the safe set is forward-invariant if ḣ_ij ≥ −γ·h_ij
# along trajectories (Ames, Coogan, Egerstedt, Notomista, Sreenath, Tabuada, "Control Barrier
# Functions: Theory and Applications", ECC 2019). Splitting the responsibility of each pair evenly
# between its two agents (Wang, Ames, Egerstedt, "Safety Barrier Certificates for Collisions-Free
# Multirobot Systems", T-RO 2017), agent i enforces, for every neighbor j within a sensing radius,
#       2 (x_i − x_j)ᵀ u_i  ≥  −(γ/2) h_ij ,
# and applies the minimally-invasive correction to its nominal command:
#       u* = argmin_u ½‖u − u_nom‖²   s.t. the linear constraints above.
#
# Each per-agent filter is a small conic QP in the IPM's primal form (P):
#       min ½ pᵀQp + cᵀp   s.t.  Bp = g,  p ∈ K,
# with cells p = (u ∈ ℝⁿ cofree, one slack s_j ∈ ℝ₊ per constraint) and one row per neighbor,
#       a_jᵀ u − s_j = b_j ,      a_j = 2(x_i − x_j),  b_j = −(γ/2) h_ij .
# Q = I on the u cell (zero on slacks), c = −u_nom on the u cell. The filter returns u_nom
# untouched whenever it already satisfies every constraint — it only acts near collisions.
#
#   include("examples/RL/lib/cbf_ipm.jl")   # (project must have CellularSheaves available)
using CellularSheaves
using CellularSheaves.IPM: IPMProblem, IPMSettings, UzawaSettings, solve,
                           CofreeCone, PositiveCone, AbstractCone, OPTIMAL, NEAR_OPTIMAL
using CellularSheaves.BlockSparseArrays: blocksparse, block, colrange
using LinearAlgebra

const IPMmod = CellularSheaves.IPM

cbf_ipm_settings() = IPMSettings{Float64}(kkt = UzawaSettings{Float64}(raug = 1e5))

"""
    cbf_filter(u_nom, x, others; d_safe=0.6, γ=6.0, sense=3.0, settings=cbf_ipm_settings())

Minimally-invasive CBF correction of the nominal command `u_nom` for an agent at `x`, given the
positions `others` of the remaining agents. Only neighbors within the sensing radius `sense`
contribute constraints; if `u_nom` already satisfies all of them it is returned unchanged.
`γ` is the class-K gain in ḣ ≥ −γh (larger = more permissive, later intervention). When the
filter must act, `bias` adds a small counterclockwise tangential component that breaks the
symmetric head-on gridlock (set `bias=0` to disable). Solved as a conic QP through
`CellularSheaves.IPM`; falls back to cyclic halfspace projection if the IPM reports a numerical
failure.
"""
function cbf_filter(u_nom::AbstractVector{Float64}, x::AbstractVector{Float64},
                    others::AbstractVector{<:AbstractVector{Float64}};
                    d_safe::Float64 = 0.6, γ::Float64 = 6.0, sense::Float64 = 3.0,
                    bias::Float64 = 0.35, settings = cbf_ipm_settings())
    n = length(u_nom)
    A = Vector{Vector{Float64}}(); b = Float64[]
    for xj in others
        dx = x .- xj
        r2 = sum(abs2, dx)
        r2 > sense^2 && continue
        h = r2 - d_safe^2
        push!(A, 2.0 .* dx)
        push!(b, -(γ / 2) * h)
    end
    isempty(A) && return copy(u_nom)
    all(dot(A[k], u_nom) ≥ b[k] - 1e-12 for k in eachindex(A)) && return copy(u_nom)

    # deadlock breaker: symmetric conflicts (agents pushing head-on into the barrier) can wedge
    # into a gridlock on the safe-set boundary. A small tangential bias on the nominal command
    # makes every agent pass its most-violated conflict counterclockwise, breaking the symmetry
    # (cf. Wang, Ames, Egerstedt, T-RO 2017, §V-B).
    if bias > 0
        k = argmax([b[k] - dot(A[k], u_nom) for k in eachindex(A)])
        dx = A[k]
        txy = hypot(dx[1], dx[2])
        if txy > 1e-9
            tang = [-dx[2] / txy, dx[1] / txy, 0.0]
            u_nom = u_nom .+ bias .* norm(u_nom) .* tang[1:n]
        end
    end

    m = length(A)
    row_ids = Int[]; col_ids = Int[]; blks = Matrix{Float64}[]
    for k in 1:m
        push!(row_ids, k); push!(col_ids, 1);     push!(blks, reshape(A[k], 1, n))
        push!(row_ids, k); push!(col_ids, 1 + k); push!(blks, -ones(1, 1))
    end
    B = blocksparse(row_ids, col_ids, blks)
    Q = IPMmod.allocblockdiag(B); fill!(Q, 0)
    block(Q, 1, 1, 1) .= Matrix{Float64}(I, n, n)
    c = zeros(n + m); c[1:n] .= -u_nom
    g = collect(b)
    K = AbstractCone[CofreeCone()]
    for _ in 1:m; push!(K, PositiveCone()); end

    u = try
        res = solve(IPMProblem(Q, B, c, g, K), settings)
        res.status in (OPTIMAL, NEAR_OPTIMAL) ? res.p[colrange(B, 1)] : nothing
    catch
        nothing
    end
    u === nothing ? _projection_fallback(u_nom, A, b) : u
end

# Cyclic (Kaczmarz-style) projection onto the violated halfspaces {u : aᵀu ≥ b} — not exact for
# several simultaneously-active constraints, but a safe last resort if the IPM solve fails.
function _projection_fallback(u_nom, A, b; passes::Int = 12)
    u = copy(u_nom)
    for _ in 1:passes
        done = true
        for k in eachindex(A)
            v = dot(A[k], u) - b[k]
            if v < 0
                u .-= v .* A[k] ./ sum(abs2, A[k])
                done = false
            end
        end
        done && break
    end
    u
end

"""
    cbf_filter_all(u_noms, xs; kwargs...)

Apply [`cbf_filter`](@ref) to every agent (each sees all the others as neighbors). Returns the
vector of corrected commands. Decentralized: each agent's QP uses only relative positions within
its sensing radius, matching what an onboard implementation would see.
"""
cbf_filter_all(u_noms, xs; kwargs...) =
    [cbf_filter(u_noms[i], xs[i], [xs[j] for j in eachindex(xs) if j != i]; kwargs...)
     for i in eachindex(xs)]
