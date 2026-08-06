# # N-Ring Cyclic Escort Formation
#
# `n` targets orbit a large circle. Each target `i` has its own escort ring of `m_i` agents at a
# small radius (scaled so rings never overlap as `n` grows). The `n` rings are additionally
# wired together at the coarse level through `n` small support pods, one between each
# consecutive pair of rings, forming a **cycle** -- and every edge in that cycle uses
# [`centroid`](@ref) (Issue 011), since neither a ring nor a pod has a single natural
# representative to expose on a cross-ring edge.
#
# The whole spec is written as a function of `n`, so the topology, the escort radius, and the
# dynamics cascade all generalize; this page renders the `n = 3` instance.
#
# ## A note on what the cyclic coupling actually does
#
# Unlike the star-shaped coupling in `centroid_formation_tracking.jl`, a *cycle* has no ring that
# is "close" to its own pin and "far" from the coupling: every ring sits between exactly two
# consensus edges (to its neighboring pods), competing against just one target-pin edge of its
# own. Left as a single `project(1)` pin, that 1-vote-vs-2-votes imbalance drags the whole
# formation toward the group centroid rather than letting each ring sit on its own target.
#
# Each ring's `Observation` here instead uses a **redundant pin**: every other agent around the
# ring observes the target, not just one (only the first uses plain `project`; the rest use
# `translation_pin`, which pins only the translation components -- see the longer note in
# `centroid_formation_tracking.jl` for why several *full* pins to the same target would otherwise
# rescale the whole rigid body). This directly rebalances the vote count in the ring's favor
# without touching the pod coupling at all, and introduces no numerical risk the way naively
# down-weighting the coupling edges did when we tried that (a `RawRestriction` scaled toward zero
# made the harmonic solve's automatic nullity threshold misclassify the weakened edges as null
# directions, blowing up the solution). The result is *not* perfect independence -- most of the
# cohesion cost remains, and that's expected: this cycle has two consensus edges pulling on every
# ring, so one extra vote narrows the gap without closing it, and that's the honest picture, not a
# rounding error to be tuned away. The tracking-error panel shows the (modestly reduced) residual,
# and the topology panel confirms the coupling itself still closes into the declared cycle.

using CellularSheaves
using CellularSheaves.ControlSheaves.NestedSystems
using CellularSheaves.ControlSheaves.AgentControllers
using LinearAlgebra
using Statistics
using Plots
using Printf

## `@__DIR__` resolves relative to Literate.jl's *output* location while a page is being
## executed, not this file's own source directory -- so the shared driver is located from the
## package root instead, which is stable regardless of how this script gets run.
include(joinpath(pkgdir(CellularSheaves), "docs", "literate", "nested", "_nested_simulation.jl"))

# ## Escort radius scaling
#
# Targets sit at angle `2π(i-1)/n` around a circle of radius `R_big`, so the chord distance
# between adjacent targets is `2 R_big sin(π/n)` -- shrinking as `n` grows, as it must for more
# targets to fit around the same circle. Ring `i`'s radius is a `safety` fraction of that chord,
# scaled by its own share of the largest ring's agent count, leaving headroom for the support pod
# between it and its neighbor.

escort_radius(n::Int, R_big::Real, m_i::Int, m_max::Int; safety::Real=0.35) =
    safety * R_big * sin(pi / n) * (m_i / m_max)

# ## Spec construction, parameterized by `n`

"""
    build_n_ring_spec(n; m=fill(5, n), R_big=3.0, support_m=2, safety=0.35, D=4)

Build the `n`-ring cyclic escort `NestedSystemSpec`, its compiled `SheafTower`, and the agent
index ranges for each ring/pod (into `tower.agent_vertices`, in the order the rings/pods were
added -- ring 1, pod 1, ring 2, pod 2, ..., matching `NestedSystems`'s depth-first agent
assignment).
"""
function build_n_ring_spec(n::Int; m::Vector{Int}=fill(5, n), R_big::Real=3.0,
                           support_m::Int=2, safety::Real=0.35, D::Int=4)
    @assert length(m) == n
    m_max = maximum(m)
    rings = [LeafTeam(Symbol(:ring, i), :ring, m[i], escort_radius(n, R_big, m[i], m_max; safety=safety))
             for i in 1:n]
    ## 2-agent pods use :path (one edge) rather than :ring, to avoid the degenerate parallel-edge
    ## 2-cycle a :ring topology gives for exactly 2 agents.
    pods = [LeafTeam(Symbol(:pod, i), :path, support_m, 0.3 * escort_radius(n, R_big, m_max, m_max; safety=safety))
            for i in 1:n]

    children = AbstractSystemNode[]
    for i in 1:n
        push!(children, rings[i])
        push!(children, pods[i])
    end

    cyc_edges = SystemEdge[]
    for i in 1:n
        push!(cyc_edges, SystemEdge(2i - 1, 2i; src_map=centroid(), dst_map=centroid()))          # ring_i -- pod_i
    end
    for i in 1:n
        push!(cyc_edges, SystemEdge(2i, (2i % 2n) + 1; src_map=centroid(), dst_map=centroid()))   # pod_i -- ring_{i+1}
    end

    root = RefinedSystem(:root, children, cyc_edges)
    targets = [TargetSpec(Symbol(:t, i)) for i in 1:n]
    ## Redundant pin: every other agent around each ring observes its own target, rebalancing
    ## the vote count against the pod's two competing consensus edges -- see the note above. Only
    ## the first pin per ring uses plain `project`; the rest use `translation_pin`, which leaves
    ## the homogeneous coordinate alone (see its docstring in `_nested_simulation.jl` -- letting
    ## several full D×D pins fight over the same target drags the homogeneous row away from 1 and
    ## rescales the whole rigid body).
    observations = Observation[]
    for i in 1:n, k in 1:2:m[i]
        push!(observations, Observation([2i - 1], i; system_map=redundant_pin(m[i], D, k)))
    end
    spec = NestedSystemSpec(root, targets, observations, D, true)
    tower = build_sheaf_tower(spec)

    team_sizes = Int[]
    for i in 1:n
        push!(team_sizes, m[i]); push!(team_sizes, support_m)
    end
    ranges = agent_index_ranges(team_sizes)
    ring_ranges = [ranges[2i - 1] for i in 1:n]
    pod_ranges = [ranges[2i] for i in 1:n]

    return spec, tower, ring_ranges, pod_ranges
end

const N = 3
const M = fill(5, N)
const R_BIG = 3.0
const SUPPORT_M = 2
const D = 4

spec, tower, ring_ranges, pod_ranges = build_n_ring_spec(N; m=M, R_big=R_BIG, support_m=SUPPORT_M)
println("Ring radius: ", round(escort_radius(N, R_BIG, M[1], maximum(M)), digits=3), " m ",
       "  (chord = ", round(2R_BIG * sin(pi / N), digits=3), " m)")

# ## Target motion: a rigid n-gon orbiting the big circle

const ω_big, h_alt = 0.3, 1.5

target_traj(i::Int) = t -> [R_BIG * cos(2π * (i - 1) / N + ω_big * t),
                            R_BIG * sin(2π * (i - 1) / N + ω_big * t), h_alt, 1.0]
target_vel(i::Int) = t -> [-R_BIG * ω_big * sin(2π * (i - 1) / N + ω_big * t),
                           R_BIG * ω_big * cos(2π * (i - 1) / N + ω_big * t), 0.0, 0.0]
target_acc(i::Int) = t -> [-R_BIG * ω_big^2 * cos(2π * (i - 1) / N + ω_big * t),
                           -R_BIG * ω_big^2 * sin(2π * (i - 1) / N + ω_big * t), 0.0, 0.0]

target_trajectories = [target_traj(i) for i in 1:N]
target_velocities = [target_vel(i) for i in 1:N]
target_accelerations = [target_acc(i) for i in 1:N]

# ## Dynamics cascade
#
# Every agent gets `QuadrotorDynamics` with a shared gain by default -- but the pods, which have
# no target of their own and exist purely for structural coordination, get a softer gain via a
# per-child override. This is `SystemBinding`'s nested-`Dict` cascade doing real work, not the
# uniform (root-only) case `centroid_formation_tracking.jl` already covers.

DT, STEPS = 0.05, 200
dyn = QuadrotorDynamics()
Ad, Bd = CellularSheaves.AgentControllers.discrete_matrices(dyn, DT)
Q_lqr = Matrix(Diagonal([500.0, 500.0, 500.0, 150.0, 150.0, 100.0, 100.0, 100.0, 5.0, 5.0]))
R_lqr = Matrix(Diagonal([0.005, 0.005, 0.005]))
K = CellularSheaves.AgentControllers.solve_dare(Ad, Bd, 10 * Q_lqr, R_lqr)
K_soft = CellularSheaves.AgentControllers.solve_dare(Ad, Bd, 2 * Q_lqr, R_lqr)

pod_overrides = Dict(Symbol(:pod, i) => SystemBinding(K_lqr=K_soft) for i in 1:N)
bindings = SystemBinding(dynamics=dyn, K_lqr=K, children=pod_overrides)

# ## Run the closed-loop simulation with full feedforward

prob = NestedEscortProblem(tower, bindings, target_trajectories;
                           target_velocities=target_velocities,
                           target_accelerations=target_accelerations,
                           dt=DT, steps=STEPS)
res = run_nested_escort_simulation(prob; use_feedforward=true)
time_grid = 0:DT:(STEPS*DT)

# ## Diagnostics: per-ring formation and tracking error over time

form_err = [zeros(STEPS) for _ in 1:N]
track_err = [zeros(STEPS) for _ in 1:N]
ring_radius = [escort_radius(N, R_BIG, M[i], maximum(M)) for i in 1:N]

for step in 1:STEPS
    t = time_grid[step]
    for i in 1:N
        positions = [res.sim_data[step][a][1:3] for a in ring_ranges[i]]
        c = sum(positions) / length(positions)
        form_err[i][step] = sum(abs(norm(p .- c) - ring_radius[i]) for p in positions) / length(positions)
        track_err[i][step] = norm(c[1:2] .- target_trajectories[i](t)[1:2])
    end
end

# ## Plots

cs = palette(:tab10)
ring_color(i) = cs[mod1(i, 10)]

## Panel 1: top-down view. `top_down_frame` renders one step -- current ring and pod agent
## positions with dashed ring outlines, and each target's trail-and-marker up to that step.
## Reused below for the animated version, so both stay visually consistent.
all_x = Float64[]; all_y = Float64[]
for rngs in (ring_ranges, pod_ranges), rng in rngs, a in rng, step in 1:STEPS
    push!(all_x, res.sim_data[step][a][1]); push!(all_y, res.sim_data[step][a][2])
end
for traj in target_trajectories, t in time_grid[1:STEPS]
    p = traj(t)
    push!(all_x, p[1]); push!(all_y, p[2])
end
const TOPDOWN_XLIM = (minimum(all_x) - 0.5, maximum(all_x) + 0.5)
const TOPDOWN_YLIM = (minimum(all_y) - 0.5, maximum(all_y) + 0.5)

function top_down_frame(step::Int)
    t = time_grid[step]
    p = plot(title=@sprintf("Top-Down View (t = %.2f s)", t), aspect_ratio=1,
            xlabel="x (m)", ylabel="y (m)", xlims=TOPDOWN_XLIM, ylims=TOPDOWN_YLIM,
            legend=:outertopright)
    for i in 1:N
        fx = [res.sim_data[step][a][1] for a in ring_ranges[i]]
        fy = [res.sim_data[step][a][2] for a in ring_ranges[i]]
        scatter!(p, fx, fy, color=ring_color(i), markersize=4, label="ring$i")
        plot!(p, [fx; fx[1]], [fy; fy[1]], color=ring_color(i), linestyle=:dash, label="")

        px = [res.sim_data[step][a][1] for a in pod_ranges[i]]
        py = [res.sim_data[step][a][2] for a in pod_ranges[i]]
        scatter!(p, px, py, color=:gray40, markersize=3, marker=:diamond, label=(i == 1 ? "pods" : ""))

        tp = [target_trajectories[i](tt) for tt in time_grid[1:step]]
        plot!(p, [q[1] for q in tp], [q[2] for q in tp], color=ring_color(i), linestyle=:dot, alpha=0.5, label="")
        pt = target_trajectories[i](t)
        scatter!(p, [pt[1]], [pt[2]], marker=:star5, markersize=9, color=ring_color(i), label="target$i")
    end
    return p
end

p1 = top_down_frame(STEPS)

## Panel 2: per-ring formation radius error.
p2 = plot(title="Formation Radius Error", xlabel="time (s)", ylabel="error (m)")
for i in 1:N
    plot!(p2, time_grid[1:STEPS], form_err[i], color=ring_color(i), label="ring$i", linewidth=2)
end

## Panel 3: per-ring centroid tracking error -- the residual cost of the cyclic coupling that the
## redundant pin doesn't fully cancel, per the note above.
p3 = plot(title="Centroid Tracking Error", xlabel="time (s)", ylabel="error (m)")
for i in 1:N
    plot!(p3, time_grid[1:STEPS], track_err[i], color=ring_color(i), label="ring$i", linewidth=2)
end

## Panel 4: high-level topology snapshot -- ring/pod centroids connected in cycle order, to
## directly validate that the centroid()-wired cycle closes the way it was declared.
p4 = plot(title="High-Level Topology (Ring/Pod Centroids)", aspect_ratio=1, xlabel="x (m)", ylabel="y (m)")
node_xy = Vector{Tuple{Float64,Float64}}()
for i in 1:N
    rc = sum(res.sim_data[end][a][1:2] for a in ring_ranges[i]) / length(ring_ranges[i])
    push!(node_xy, (rc[1], rc[2]))
    scatter!(p4, [rc[1]], [rc[2]], color=ring_color(i), markersize=8, label=(i == 1 ? "ring centroid" : ""))
    pc = sum(res.sim_data[end][a][1:2] for a in pod_ranges[i]) / length(pod_ranges[i])
    push!(node_xy, (pc[1], pc[2]))
    scatter!(p4, [pc[1]], [pc[2]], color=:gray40, marker=:diamond, markersize=6, label=(i == 1 ? "pod centroid" : ""))
end
cyc_x = [p[1] for p in node_xy]; push!(cyc_x, node_xy[1][1])
cyc_y = [p[2] for p in node_xy]; push!(cyc_y, node_xy[1][2])
plot!(p4, cyc_x, cyc_y, color=:black, linewidth=1.5, alpha=0.6, label="cycle")

plot(p1, p2, p3, p4, layout=(2, 2), size=(1100, 900),
    plot_title="$N-Ring Cyclic Escort Formation (centroid-wired)")

# The tracking-error panel should settle to a much smaller residual than a single-pin tower would
# give (compare against `centroid_formation_tracking.jl`'s much larger, unmitigated offset), while
# the topology panel confirms the coupling itself still closes into the declared cycle.

@printf("Mean steady-state tracking error (last 20%% of run): %.3f m (ring radius = %.3f m)\n",
       sum(mean(track_err[i][round(Int, 0.8STEPS):end]) for i in 1:N) / N, ring_radius[1])

# ## Animated top-down view
#
# Same `top_down_frame` function as panel 1 above, stepped through the run -- watch each ring
# hold reasonably close to its own orbiting target while the pods visibly bridge between them.

anim = @animate for step in 1:4:STEPS
    top_down_frame(step)
end
gif(anim, "n_ring_formation_top_down.gif", fps=10)

# ![N-ring formation top-down animation](n_ring_formation_top_down.gif)

println("N-ring cyclic escort formation example complete.")
