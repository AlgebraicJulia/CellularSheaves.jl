# One axis of the lattice double integrator against the continuous one.
# Lattice: offset D (target minus position) and velocity v; a unit-cost step picks
# v' ∈ {v-1, v, v+1} and sets D' = D - v'. Continuous: ẍ = u, |u| ≤ 1, minimum time
# to rest at the target (bang-bang with one switch). Plain Julia.
#
#   julia docs/scripts/lattice_vs_continuous.jl [Dmax] [vmax]

using Printf

Dmax = length(ARGS) >= 1 ? parse(Int, ARGS[1]) : 60
vq = length(ARGS) >= 2 ? parse(Int, ARGS[2]) : 8
vbig, R = 4vq + 8, 4Dmax + 8vq^2                     # bounds large enough not to bind

T = Dict((0, 0) => 0)
queue = [(0, 0)]
while !isempty(queue)
    (D, v) = popfirst!(queue)
    for vp in (v - 1):(v + 1)                         # predecessors: from (D + v, vp) the step v' = v lands here
        abs(vp) <= vbig || continue
        Dp = D + v
        abs(Dp) <= R || continue
        haskey(T, (Dp, vp)) && continue
        T[(Dp, vp)] = T[(D, v)] + 1
        push!(queue, (Dp, vp))
    end
end

function continuous_time(D, v)
    x = -D                                            # position relative to the target
    s = x + v * abs(v) / 2                            # switching curve x = -v|v|/2
    abs(s) < 1e-12 && return abs(v)
    return s > 0 ? v + 2sqrt(x + v^2 / 2) : -v + 2sqrt(-x + v^2 / 2)
end

diffs, ratios = Float64[], Float64[]
for D in -Dmax:Dmax, v in -vq:vq
    (D, v) == (0, 0) && continue
    t, c = T[(D, v)], continuous_time(D, v)
    push!(diffs, t - c); push!(ratios, t / c)
end
@printf("|D| ≤ %d, |v| ≤ %d: lattice − continuous time in [%.3f, %.3f], mean %.3f; ratio in [%.3f, %.3f]\n",
    Dmax, vq, minimum(diffs), maximum(diffs), sum(diffs) / length(diffs), minimum(ratios), maximum(ratios))
