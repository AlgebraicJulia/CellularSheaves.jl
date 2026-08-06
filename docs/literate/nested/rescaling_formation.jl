# # Emergent Rescaling: A Rigid Team Absorbing Divergent Targets Through Scale
#
# The other two examples in this section go out of their way to keep the affine homogeneous
# coordinate pinned at exactly `1` -- `translation_pin` exists specifically to stop several
# redundant pins from dragging it away and rescaling a rigid body as an unwanted side effect. This
# example does the opposite: it uses **unmodified, fully redundant `project(k)` pins** on purpose,
# and shows that the "bug" is actually a legitimate, useful degree of freedom.
#
# ## The mechanism
#
# A rigid team's space of exact global sections is `D`-dimensional pure translation -- normally
# `D-1` translation coordinates plus a homogeneous row that stays locked at `1`. But when several
# agents are *each* pinned to a *different* target (as opposed to redundantly pinning several
# agents to the *same* target, as in `centroid_formation_tracking.jl`), those pins are jointly
# satisfiable only if the targets happen to sit at exactly the team's nominal relative geometry.
# The moment they don't, the least-squares compromise the harmonic solve finds doesn't just
# translate the team -- it also finds the homogeneous coordinate `w` that best explains the
# mismatch, and because these affine restriction maps scale every translation offset by `w`, that
# compromise **rescales the whole rigid body uniformly**, without ever bending its angles.
#
# ## Setup
#
# Six targets start at a common origin and diverge outward along a rigid hexagon's own six
# nominal directions -- but at six *different* rates, so their divergence is not a uniform
# scaling of the hexagon. A 6-agent rigid ring, one agent pinned to each target via a full
# (untruncated) `project(k)`, has only one shared translation-plus-scale to offer in response.
# We compare it against an **unconstrained** placement -- each of the six targets tracked by its
# own independent, unconnected point, with no shape constraint at all -- which is free to trace
# out whatever irregular shape the six different divergence rates actually produce.

using CellularSheaves
using CellularSheaves.ControlSheaves.NestedSystems
using LinearAlgebra
using Statistics
using Plots
using Printf

## `@__DIR__` resolves relative to Literate.jl's *output* location while a page is being
## executed, not this file's own source directory -- so the shared driver is located from the
## package root instead, which is stable regardless of how this script gets run.
include(joinpath(pkgdir(CellularSheaves), "docs", "literate", "nested", "_nested_simulation.jl"))

# ## Topology: one rigid ring, every agent individually pinned

const N = 6
const D = 4
const h_alt = 1.5

ring = LeafTeam(:formation, :ring, N, 1.0)
root = RefinedSystem(:root, AbstractSystemNode[ring])
targets = [TargetSpec(Symbol(:t, k)) for k in 1:N]

## Deliberately full, untruncated `project(k)` pins -- no `redundant_pin`/`translation_pin` here,
## since the whole point is to let the homogeneous coordinate drift.
observations = [Observation([1], k; system_map=project(k)) for k in 1:N]
spec = NestedSystemSpec(root, targets, observations, D, true)
tower = build_sheaf_tower(spec)

# ## Target motion: a common origin, diverging at six different rates
#
# Target `k` sits at the ring's own nominal angle `2π(k-1)/N`, growing outward at rate
# `base_rate * scale_factor(k)` -- `scale_factor` ranges linearly from `0.6` to `1.4` across the
# six targets, so their divergence is deliberately *not* a uniform scaling of the hexagon.

const base_rate = 0.4
angle(k) = 2π * (k - 1) / N
scale_factor(k) = 0.6 + 0.8 * (k - 1) / (N - 1)

target_traj(k) = t -> [base_rate * scale_factor(k) * t * cos(angle(k)),
                       base_rate * scale_factor(k) * t * sin(angle(k)), h_alt, 1.0]
target_trajectories = [target_traj(k) for k in 1:N]

const T_FINAL, N_STEPS = 10.0, 200
time_grid = range(0.0, T_FINAL, length=N_STEPS)

# ## Sweep the reference solve over time
#
# This is deliberately reference-only -- no closed-loop agent dynamics -- since the point is
# purely about what the harmonic solve geometrically produces, and closed-loop tracking lag would
# only obscure that with an unrelated confound.

rigid = [solve_hierarchical(tower, [traj(t) for traj in target_trajectories])[end] for t in time_grid]
rigid_agents = [[q[v] for v in tower.agent_vertices] for q in rigid]
unconstrained_agents = [[traj(t) for traj in target_trajectories] for t in time_grid]

## `w`, the homogeneous coordinate every agent shares exactly (verified below) -- this is the
## mechanism itself, made visible.
w_history = [rigid_agents[s][1][4] for s in 1:N_STEPS]
@assert all(s -> all(a -> isapprox(a[4], w_history[s]; atol=1e-8), rigid_agents[s]), 1:N_STEPS) "expected every agent to share exactly one homogeneous coordinate at each step"

## Shape regularity: the standard deviation of each agent's distance from the formation's own
## centroid. Exactly zero for a perfectly regular (possibly rescaled) polygon; grows as soon as a
## formation deforms away from its nominal regular shape.
function regularity(agents::Vector{Vector{Float64}})
    positions = [a[1:3] for a in agents]
    c = sum(positions) / length(positions)
    dists = [norm(p .- c) for p in positions]
    return std(dists)
end

rigid_irregularity = [regularity(rigid_agents[s]) for s in 1:N_STEPS]
unconstrained_irregularity = [regularity(unconstrained_agents[s]) for s in 1:N_STEPS]

@printf("Final: w = %.3f, rigid shape irregularity = %.2e, unconstrained shape irregularity = %.3f\n",
       w_history[end], rigid_irregularity[end], unconstrained_irregularity[end])

# ## Plots

colors = (rigid=:steelblue, unconstrained=:darkorange)

all_x = Float64[]; all_y = Float64[]
for agents in (rigid_agents, unconstrained_agents), step_agents in agents, a in step_agents
    push!(all_x, a[1]); push!(all_y, a[2])
end
const TOPDOWN_XLIM = (minimum(all_x) - 0.5, maximum(all_x) + 0.5)
const TOPDOWN_YLIM = (minimum(all_y) - 0.5, maximum(all_y) + 0.5)

## Panel 1: top-down composite -- the rigid formation's regular polygon growing in place across
## sparse snapshots (fading faint-to-solid) against the unconstrained swarm's increasingly
## irregular spread, at the same snapshots.
function formation_snapshot!(p, agents_at_step::Vector{Vector{Float64}}, color, alpha, label)
    fx = [a[1] for a in agents_at_step]; fy = [a[2] for a in agents_at_step]
    plot!(p, [fx; fx[1]], [fy; fy[1]], color=color, linestyle=:dash, alpha=alpha, label=label)
    scatter!(p, fx, fy, color=color, markersize=3, alpha=alpha, label="")
end

p1 = plot(title="Top-Down View (full trajectory)", aspect_ratio=1, xlabel="x (m)", ylabel="y (m)",
         xlims=TOPDOWN_XLIM, ylims=TOPDOWN_YLIM, legend=:outertopright)
for k in 1:N
    rpx = [rigid_agents[s][k][1] for s in 1:N_STEPS]; rpy = [rigid_agents[s][k][2] for s in 1:N_STEPS]
    plot!(p1, rpx, rpy, color=colors.rigid, alpha=0.12, linewidth=1, label="")
    upx = [unconstrained_agents[s][k][1] for s in 1:N_STEPS]; upy = [unconstrained_agents[s][k][2] for s in 1:N_STEPS]
    plot!(p1, upx, upy, color=colors.unconstrained, alpha=0.12, linewidth=1, label="")
end
snaps = snapshot_steps(N_STEPS, 8)
for (si, s) in enumerate(snaps)
    a = fade_alpha(si, length(snaps))
    formation_snapshot!(p1, rigid_agents[s], colors.rigid, a, s == snaps[end] ? "rigid (rescaling)" : "")
    formation_snapshot!(p1, unconstrained_agents[s], colors.unconstrained, a,
                        s == snaps[end] ? "unconstrained" : "")
end

## Panel 2: the homogeneous coordinate `w` over time -- the mechanism made visible.
p2 = plot(title="Homogeneous Coordinate (Scale)", xlabel="time (s)", ylabel="w", legend=false)
plot!(p2, collect(time_grid), w_history, color=colors.rigid, linewidth=2)
hline!(p2, [1.0], color=:gray, linestyle=:dash, alpha=0.6)

## Panel 3: shape irregularity over time -- rigid stays at (numerically) zero throughout;
## unconstrained grows as the six different divergence rates pull it away from regular.
p3 = plot(title="Shape Irregularity (std. of radii)", xlabel="time (s)", ylabel="std(agent distance to centroid) (m)")
plot!(p3, collect(time_grid), rigid_irregularity, color=colors.rigid, label="rigid", linewidth=2)
plot!(p3, collect(time_grid), unconstrained_irregularity, color=colors.unconstrained, label="unconstrained", linewidth=2)

## Panel 4: final-frame side-by-side shapes, for a crisp direct comparison.
p4 = plot(title="Final Shape Comparison (t = $(T_FINAL) s)", aspect_ratio=1, xlabel="x (m)", ylabel="y (m)")
formation_snapshot!(p4, rigid_agents[end], colors.rigid, 1.0, "rigid")
formation_snapshot!(p4, unconstrained_agents[end], colors.unconstrained, 1.0, "unconstrained")

plot(p1, p2, p3, p4, layout=(2, 2), size=(1100, 900),
    plot_title="Emergent Rescaling vs. Shear Under Divergent Targets")

# The rigid formation traces a perfect (growing) regular hexagon at every snapshot -- its shape
# irregularity panel should read as (numerically) zero throughout, confirming it never shears, only
# rescales. The unconstrained swarm's irregularity grows steadily, and its final shape is visibly
# not a regular hexagon -- the six different divergence rates show up directly as distortion.

println("Emergent rescaling formation example complete.")
