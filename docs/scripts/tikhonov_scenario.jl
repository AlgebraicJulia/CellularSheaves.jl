module TikhonovGuideScenario

using CellularSheaves
using CellularSheaves.ControlSheaves.Tikhonov
using CellularSheaves.NetworkSheaves.DistributedSolve
using CliqueTrees.Multifrontal
using LinearAlgebra
using SparseArrays
import CellularSheaves.NetworkSheaves.EuclideanSheaves:
    _harmonic_extension_restricted_laplacian

export NA, NT, D, TV1, TV2, AGENT_EDGES, TARGET_PERIOD, H, B,
    target1, target2, target1_rate, target2_rate, harmonic_reference,
    harmonic_rate, planner_rollout, unpack_agents

const NA, NT, D = 12, 2, 2
const TV1, TV2 = NA + 1, NA + 2
const I2 = Matrix{Float64}(I, D, D)
const TARGET_PERIOD = 12.0
const OMEGA = 2pi / TARGET_PERIOD
const AGENT_EDGES = [(1, 2), (2, 3), (3, 4), (4, 5), (5, 6), (6, 1),
                     (7, 8), (8, 9), (9, 10), (10, 11), (11, 12), (12, 7),
                     (1, 7)]

function _build_sheaf()
    sheaf = EuclideanSheaf{Float64}(fill(D, NA + NT))
    for (i, j) in AGENT_EDGES
        add_sheaf_edge!(sheaf, i, j, I2, I2)
    end
    for i in 1:6
        add_sheaf_edge!(sheaf, i, TV1, I2, I2)
    end
    for i in 7:12
        add_sheaf_edge!(sheaf, i, TV2, I2, I2)
    end
    return sheaf
end

target1(t) = [-2.2 + 1.45cos(OMEGA * t), 1.15sin(OMEGA * t)]
target2(t) = [ 2.2 + 1.20cos(OMEGA * t + pi), 0.95sin(2OMEGA * t + pi / 4)]
target1_rate(t) = [-1.45OMEGA * sin(OMEGA * t), 1.15OMEGA * cos(OMEGA * t)]
target2_rate(t) = [-1.20OMEGA * sin(OMEGA * t + pi), 1.90OMEGA * cos(2OMEGA * t + pi / 4)]

const SHEAF = _build_sheaf()
const BOUNDARY0 = Dict(TV1 => target1(0.0), TV2 => target2(0.0))
const _REDUCTION = _harmonic_extension_restricted_laplacian(SHEAF, BOUNDARY0)
const H = sparse(_REDUCTION[3])
const B = sparse(_REDUCTION[4])
const FACTOR = cholesky!(ChordalCholesky(H), NoPivot())
const LFAC = FACTOR.L
const WORKSPACE = TreeWorkspace(LFAC, Float64)

function _cached_solve(rhs)
    state = Vector(FACTOR.P' \ rhs)
    tree_forward_ldiv!(LFAC, state, WORKSPACE)
    tree_backward_ldiv!(LFAC, state, WORKSPACE)
    return Vector(FACTOR.P \ state)
end

harmonic_reference(t) = _cached_solve(-B * vcat(target1(t), target2(t)))
harmonic_rate(t) = _cached_solve(-B * vcat(target1_rate(t), target2_rate(t)))

unpack_agents(q::AbstractVector) = reshape(q, D, NA)
unpack_agents(q::AbstractMatrix) = permutedims(reshape(q, D, NA, size(q, 2)), (1, 3, 2))

"""Simulate the normalized planner on the shared two-target harmonic reference."""
function planner_rollout(times; epsilon, feedforward=false, displaced=false)
    K = length(times)
    K >= 2 || throw(ArgumentError("times must contain at least two samples"))
    ideal = hcat(harmonic_reference.(times)...)
    rates = hcat(harmonic_rate.(times)...)
    x0 = copy(ideal[:, 1])
    if displaced
        x0 .+= repeat([0.65, -0.35], NA)
    end
    filter = TikhonovFilter(x0; epsilon)
    state = similar(ideal)
    state[:, 1] = filter.x
    for k in 1:K-1
        dt = times[k + 1] - times[k]
        if feedforward
            u0 = tikhonov_feedforward_reference(ideal[:, k], rates[:, k], epsilon)
            u1 = tikhonov_feedforward_reference(ideal[:, k + 1], rates[:, k + 1], epsilon)
            tikhonov_step!(filter, u0, u1, dt)
        else
            tikhonov_step!(filter, ideal[:, k], ideal[:, k + 1], dt)
        end
        state[:, k + 1] = filter.x
    end
    error = vec(sqrt.(sum(abs2, state - ideal; dims=1) ./ NA))
    return (; times=collect(times), ideal, rates, state, error)
end

end
