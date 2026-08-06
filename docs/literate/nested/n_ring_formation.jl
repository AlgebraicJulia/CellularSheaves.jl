# # N-Ring Cyclic Escort Formation
#
# `n` targets orbit a large circle. Each target `i` has its own escort ring of `m_i` agents at a
# small radius (scaled so rings never overlap as `n` grows). The `n` rings are additionally wired
# together at the coarse level through `n` small support pods, one between each consecutive pair
# of rings, forming a **cycle** -- and every edge in that cycle uses [`centroid`](@ref)
# (Issue 011), since neither a ring nor a pod has a single natural representative to expose on a
# cross-ring edge. The whole spec is written as a function of `n`, so the topology, the escort
# radius, and the dynamics cascade all generalize; this page renders the `n = 3` instance.
#
# Watching the result, each ring sits *near* its own target but not exactly on it -- the next
# section explains precisely why, with a direct measurement rather than a hand-wave.
#
# ## Why a cycle, and not a tree
#
# `wheel_formation.jl`, the previous page in this section, bridges `n` rings with `n` support pods
# too -- but arranged as a **tree** (spokes into a shared, unanchored hub) rather than a cycle, and
# gets noticeably tighter tracking as a result (own-target weight around `0.74`, against `0.6`
# here). That's not an oversight in this page's design; it's a real structural limit.
# `NestedSystemSpec`'s tree of `RefinedSystem`/`LeafTeam` nodes gives every node exactly one
# parent, and `SystemEdge` only ever connects two direct children of the *same* parent. A bridge
# between ring `i` and ring `i+1` needs both of them exposed to whatever bridges them -- fine for
# one bridge, but ring `i` also needs to be exposed to the *previous* bridge (with ring `i-1`), and
# a single node cannot be a child of two different parents at once. There is no way to make a
# cyclic bridge structure a child of one shared parent the way the wheel's spokes can all share the
# hub. Wiring the cycle's edges directly into `root`'s own `internal_edges`, as this page does, is
# the only way to express "ring `i` bridges into both of its neighbors" at all -- and it's exactly
# that direct, unmediated coupling between every adjacent pair that produces the worse cross-talk
# measured below.

using CellularSheaves
using CellularSheaves.ControlSheaves.NestedSystems
using CellularSheaves.ControlSheaves.AgentControllers
using LinearAlgebra
using Statistics
using Plots
using Printf

## `@__DIR__` resolves relative to Literate.jl's *output* location while a page is being
## executed, not this file's own source directory -- so the plotting helpers are located from the
## package root instead, which is stable regardless of how this script gets run. The simulation
## driver and multi-pin helpers themselves live in `NestedSystems` (already `using`-ed above).
include(joinpath(pkgdir(CellularSheaves), "docs", "literate", "nested", "_plot_helpers.jl"))

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
    build_n_ring_spec(n; m=fill(5, n), R_big=3.0, support_m=2, safety=0.35, D=4, pin_scheme=:redundant)

Build the `n`-ring cyclic escort `NestedSystemSpec`, its compiled `SheafTower`, and the agent
index ranges for each ring/pod (into `tower.agent_vertices`, in the order the rings/pods were
added -- ring 1, pod 1, ring 2, pod 2, ..., matching `NestedSystems`'s depth-first agent
assignment).

`pin_scheme` selects how each ring observes its own target: `:redundant` (the default) pins every
other agent around the ring, via [`NestedSystems.redundant_pin`](@ref); `:single` pins only the
first agent, the same default every other example in this section uses. Both are used below --
`:redundant` for the actual simulation, `:single` only as a comparison point for measuring how
much the redundant pin actually buys.
"""
function build_n_ring_spec(n::Int; m::Vector{Int}=fill(5, n), R_big::Real=3.0,
                           support_m::Int=2, safety::Real=0.35, D::Int=4,
                           pin_scheme::Symbol=:redundant)
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
    observations = Observation[]
    for i in 1:n
        ks = pin_scheme == :redundant ? (1:2:m[i]) : (1:1)
        for k in ks
            system_map = pin_scheme == :redundant ? redundant_pin(m[i], D, k) : project(1)
            push!(observations, Observation([2i - 1], i; system_map=system_map))
        end
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

# ## Why rings don't settle on their own targets
#
# Every ring sits between exactly two `centroid()` consensus edges (to its neighboring pods),
# competing against its own target-pin edge -- and the pods have **no target of their own**: at
# the joint least-squares optimum, an unanchored pod's position is *exactly* the average of its
# two neighboring rings' own solved positions (a direct consequence of it appearing in only those
# two quadratic terms). That makes every pod a perfect, undamped relay -- it doesn't just couple a
# ring to its immediate neighbors, it couples every ring to *every other ring's target*, carried
# around the whole cycle with no attenuation from the pods themselves.
#
# Because the whole system is linear in the target positions, this is directly measurable: nudge
# one target at a time and watch how much a ring's own solved position moves.

"""
    target_response_weights(tower, ring_ranges, ring_idx, n) -> Vector{Float64}

The linear response of ring `ring_idx`'s own centroid to a unit displacement of each of the `n`
targets in turn, holding the others fixed at a common base point. Sums to `1` for any ring in a
translation-invariant tower (translating every target by the same vector must translate every
ring's solved position by that same vector) -- so this is a set of blending weights, not just a
sensitivity measure.
"""
function target_response_weights(tower::SheafTower, ring_ranges, ring_idx::Int, n::Int)
    base = fill([0.0, 0.0, 1.5, 1.0], n)
    ## `ring_ranges` holds *local* agent-position ranges (see `agent_index_ranges`) -- the
    ## convention `res.sim_data` uses elsewhere on this page. `solve_hierarchical`'s own output is
    ## indexed by actual sheaf vertex number (targets first, then agents), a different convention
    ## -- so this maps through `tower.agent_vertices` to convert before indexing into it.
    verts = tower.agent_vertices[ring_ranges[ring_idx]]
    ring_centroid(tv) = begin
        q = solve_hierarchical(tower, tv)[end]
        sum(q[v] for v in verts) / length(verts)
    end
    c0 = ring_centroid(base)
    h = 1.0
    return [begin
                perturbed = copy(base)
                perturbed[j] = base[j] .+ [h, 0.0, 0.0, 0.0]
                (ring_centroid(perturbed)[1] - c0[1]) / h
            end
            for j in 1:n]
end

w = target_response_weights(tower, ring_ranges, 1, N)
@printf("ring1's target-response weights: own=%.3f, neighbors=%.3f/%.3f (sum=%.3f)\n", w[1], w[2], w[3], sum(w))

# With the redundant pin used below, ring 1's position is (to the precision printed above) `0.6`
# of its own target plus `0.2` of each neighbor's -- not `1.0`/`0.0`/`0.0` the way a truly
# independent ring would respond. That 60/20/20 split, not any single bad edge or a bug in the
# pin scheme, is *why* the tracking-error panel further down never reaches zero: it's measuring
# exactly this blend. See the "Design notes" section at the end of this page for how that split
# compares to a plain single-pin ring, and why the redundant pin narrows it only modestly rather
# than eliminating it -- and see `wheel_formation.jl` for how much of this blend is specifically a
# cost of the cycle, versus a cost of having *any* bridge at all.

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

## Panel 1: top-down view, and the animated version below both need consistent axis limits --
## computed once, from every agent and target position across the whole run.
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

## Panel 1 proper: a composite of the whole run. Every agent's full path is drawn thin and faint
## (dt is far finer than a static plot needs); ring/pod outlines are drawn at sparse, evenly-spaced
## snapshots (not one per step) with alpha fading from faint (early) to solid (final).
function top_down_composite(; n_snapshots::Int=8)
    snaps = snapshot_steps(STEPS, n_snapshots)
    p = plot(title="Top-Down View (full trajectory)", aspect_ratio=1,
            xlabel="x (m)", ylabel="y (m)", xlims=TOPDOWN_XLIM, ylims=TOPDOWN_YLIM,
            legend=:outertopright)
    for i in 1:N
        for a in ring_ranges[i]
            px = [res.sim_data[s][a][1] for s in 1:STEPS]
            py = [res.sim_data[s][a][2] for s in 1:STEPS]
            plot!(p, px, py, color=ring_color(i), alpha=0.1, linewidth=1, label="")
        end
        for (si, s) in enumerate(snaps)
            a = fade_alpha(si, length(snaps))
            fx = [res.sim_data[s][v][1] for v in ring_ranges[i]]
            fy = [res.sim_data[s][v][2] for v in ring_ranges[i]]
            plot!(p, [fx; fx[1]], [fy; fy[1]], color=ring_color(i), linestyle=:dash, alpha=a,
                 label=(s == snaps[end] ? "ring$i" : ""))
            scatter!(p, fx, fy, color=ring_color(i), markersize=3, alpha=a, label="")

            px = [res.sim_data[s][v][1] for v in pod_ranges[i]]
            py = [res.sim_data[s][v][2] for v in pod_ranges[i]]
            scatter!(p, px, py, color=:gray40, markersize=3, marker=:diamond, alpha=a,
                    label=(i == 1 && s == snaps[end] ? "pods" : ""))
        end

        tp = [target_trajectories[i](t) for t in time_grid[1:STEPS]]
        plot!(p, [q[1] for q in tp], [q[2] for q in tp], color=ring_color(i), linestyle=:dot, alpha=0.5, label="")
        for (si, s) in enumerate(snaps)
            pt = target_trajectories[i](time_grid[s])
            a = fade_alpha(si, length(snaps))
            scatter!(p, [pt[1]], [pt[2]], marker=:star5, markersize=(s == snaps[end] ? 9 : 5),
                    color=ring_color(i), alpha=a, label=(s == snaps[end] ? "target$i" : ""))
        end
    end
    return p
end

p1 = top_down_composite()

## The animated version reuses the same per-step rendering in miniature -- current positions and
## outline only, no accumulated trail (the composite above already shows the whole run at once).
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

## Panel 2: per-ring formation radius error.
p2 = plot(title="Formation Radius Error", xlabel="time (s)", ylabel="error (m)")
for i in 1:N
    plot!(p2, time_grid[1:STEPS], form_err[i], color=ring_color(i), label="ring$i", linewidth=2)
end

## Panel 3: per-ring centroid tracking error -- the 40% (own weight 0.6, not 1.0) that the
## redundant pin leaves on the table, per the "Why rings don't settle on their own targets"
## section above.
p3 = plot(title="Centroid Tracking Error", xlabel="time (s)", ylabel="error (m)")
for i in 1:N
    plot!(p3, time_grid[1:STEPS], track_err[i], color=ring_color(i), label="ring$i", linewidth=2)
end

## Panel 4: high-level topology over time -- ring/pod centroids connected in cycle order, at the
## same sparse snapshots as panel 1 (fading faint-to-solid), so the coarse topology's own arc
## through time is visible directly, not just its final shape. Directly validates that the
## centroid()-wired cycle closes the way it was declared, at every snapshot, not only the last.
p4 = plot(title="High-Level Topology Over Time (Ring/Pod Centroids)", aspect_ratio=1,
         xlabel="x (m)", ylabel="y (m)")
## Skip the very start of the run: agents begin from the airstrip layout (see
## `_default_initial_position`), and that convergence transient is already covered by panels 2/3
## -- including it here would just stretch this panel's axes around one uninformative outlier.
topo_snaps = snapshot_steps(STEPS, 9)[2:end]
for (si, s) in enumerate(topo_snaps)
    a = fade_alpha(si, length(topo_snaps))
    node_xy = Vector{Tuple{Float64,Float64}}()
    for i in 1:N
        rc = sum(res.sim_data[s][v][1:2] for v in ring_ranges[i]) / length(ring_ranges[i])
        push!(node_xy, (rc[1], rc[2]))
        scatter!(p4, [rc[1]], [rc[2]], color=ring_color(i), markersize=(s == topo_snaps[end] ? 8 : 4),
                alpha=a, label=(i == 1 && s == topo_snaps[end] ? "ring centroid" : ""))
        pc = sum(res.sim_data[s][v][1:2] for v in pod_ranges[i]) / length(pod_ranges[i])
        push!(node_xy, (pc[1], pc[2]))
        scatter!(p4, [pc[1]], [pc[2]], color=:gray40, marker=:diamond,
                markersize=(s == topo_snaps[end] ? 6 : 3), alpha=a,
                label=(i == 1 && s == topo_snaps[end] ? "pod centroid" : ""))
    end
    cyc_x = [pt[1] for pt in node_xy]; push!(cyc_x, node_xy[1][1])
    cyc_y = [pt[2] for pt in node_xy]; push!(cyc_y, node_xy[1][2])
    plot!(p4, cyc_x, cyc_y, color=:black, linewidth=(s == topo_snaps[end] ? 1.5 : 1.0), alpha=a * 0.7,
         label=(s == topo_snaps[end] ? "cycle" : ""))
end

plot(p1, p2, p3, p4, layout=(2, 2), size=(1100, 900),
    plot_title="$N-Ring Cyclic Escort Formation (centroid-wired)")

@printf("Mean steady-state tracking error (last 20%% of run): %.3f m (ring radius = %.3f m)\n",
       sum(mean(track_err[i][round(Int, 0.8STEPS):end]) for i in 1:N) / N, ring_radius[1])

# Note too, in the topology panel, that the pods sit closer to the *group centroid* than a naive
# "midpoint between the two targets each pod bridges" would predict -- a pod has no target of its
# own, so it sits exactly at the midpoint of its two neighboring rings' *actual* (already-blended)
# positions, compounding the same effect one step further.

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

# ## Design notes: alternatives considered
#
# Two things were tried and rejected while building this example, kept here for anyone tempted to
# revisit either:
#
# **Weakening the pod coupling with a scaled `RawRestriction`.** Scaling both ends of each cyclic
# edge by a factor `α < 1` (`RawRestriction(α .* materialize_restriction(centroid(), node, D))`)
# does reduce the pull as `α → 0` -- but only down to a floor set by the single-pin boundary
# offset (see `centroid_formation_tracking.jl`), and pushing `α` much smaller than that floor
# needs made the harmonic solve's automatic nullity threshold misclassify the now-very-weak edges
# as null directions, blowing up the solution entirely. Not a usable knob.
#
# **Just adding more redundant pins.** The measurement above explains why this has limited
# payoff: the blend isn't caused by an edge-*count* imbalance that more edges can out-vote, it's
# caused by the pods providing an **undamped** path to every other ring's target, and no amount of
# strengthening one ring's own pin removes that path. To see the size of the effect the redundant
# pin *does* buy, compare against a plain single pin on the same topology:

_, tower_single, ring_ranges_single, _ = build_n_ring_spec(N; m=M, R_big=R_BIG, support_m=SUPPORT_M,
                                                            pin_scheme=:single)
w_single = target_response_weights(tower_single, ring_ranges_single, 1, N)
@printf("single pin:    own=%.3f, neighbors=%.3f/%.3f\n", w_single[1], w_single[2], w_single[3])
@printf("redundant pin: own=%.3f, neighbors=%.3f/%.3f\n", w[1], w[2], w[3])

# The redundant pin moves ring 1's own weight from about `0.56` to `0.6` -- a real, measurable
# improvement, and worth having, but nowhere near the `1.0` a truly independent ring would show.
# That gap is the honest cost of a cyclic topology whose bridges have no anchor of their own.

println("Design notes complete.")
