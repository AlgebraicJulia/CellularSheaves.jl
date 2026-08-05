# # Diffusion versus Direct: Two Routes to the Harmonic Extension
#
# A pinned cellular sheaf defines one harmonic formation ``q^\star``. There are two ways to
# fly a fleet onto it, and this example benchmarks them against each other.
#
# **Diffusion** is the decentralized law of
# [*Heterogeneous Multi-Agent Multi-Target Tracking using Cellular Sheaves*](https://arxiv.org/abs/2512.24886)
# (eq. 4 and 5). Each agent feeds back its own sheaf disagreement on every tick:
#
# ```math
# u = -k\, g^{+}\eta, \qquad \eta = \mathcal{H}q - Bp .
# ```
#
# **Direct** solves the harmonic extension first, then tracks the delivered reference
# locally:
#
# ```math
# u = -k\, g^{+}\!\left(q - q^\star\right), \qquad q^\star = \mathcal{H}^{-1}Bp .
# ```
#
# ## The two laws solve the same linear system
#
# It is worth being precise about what Diffusion does, because it is easy to read it as a
# heuristic alongside a solver. Substituting the residual into the control law and the
# control law into the plant gives, for one tick of step ``\Delta t``,
#
# ```math
# q \;\leftarrow\; q - k\,\Delta t\,\bigl(\mathcal{H}q - Bp\bigr)
#   \;=\; \bigl(I - k\,\Delta t\,\mathcal{H}\bigr)q + k\,\Delta t\,Bp .
# ```
#
# That is Richardson iteration on ``\mathcal{H}q^\star = Bp`` with step size
# ``k\,\Delta t``. Diffusion is an unpreconditioned first-order iterative solver for exactly
# the system that Direct factors. It converges to ``q^\star`` rather than computing it, and
# three consequences follow that the rest of this page measures:
#
# 1. The error contracts by the spectral radius of ``I - k\,\Delta t\,\mathcal{H}``, so the
#    tick count scales with ``\kappa(\mathcal{H})``.
# 2. The iteration is matrix-free. ``\mathcal{H}`` is applied by neighbour sums and is never
#    assembled, so there is nothing to store.
# 3. The iterate is the fleet's physical position. Diffusion does not hold a copy of
#    ``q^\star`` anywhere, because the formation itself is the iterate.
#
# Direct takes the other classical route: a sparse Cholesky factorization of
# ``\mathcal{H}``, reused across solves. The comparison below is therefore the familiar
# one between an iterative matrix-free method and a direct factorization, carried out on a
# control problem where the iterative method's iterate is made of aircraft.

using CellularSheaves
using CellularSheaves.ControlSheaves.CoordinationBenchmarks
using CellularSheaves.ControlSheaves.CoordinationProfiling
using LinearAlgebra
using Markdown
using Printf
using Plots

default(framestyle = :box, grid = true, gridalpha = 0.18, gridstyle = :dot,
    titlefontsize = 10, guidefontsize = 9, legendfontsize = 8, tickfontsize = 8,
    markerstrokewidth = 0, size = (720, 380))

comma(v) = replace(string(v), r"(?<=\d)(?=(\d{3})+$)" => ",")

fmt_bytes(b) = b >= 1_048_576 ? @sprintf("%.1f MB", b / 1_048_576) :
               b >= 1024 ? @sprintf("%.1f kB", b / 1024) : "$(comma(b)) B"

# Winning cells are shaded green in every table below.

function bench_table(headers, rows)
    io = IOBuffer()
    print(io, "<table class=\"bench\"><thead><tr>")
    for h in headers
        print(io, "<th>", h, "</th>")
    end
    print(io, "</tr></thead><tbody>")
    for (cells, wins) in rows
        print(io, "<tr>")
        for (c, w) in zip(cells, wins)
            print(io, w ? "<td class=\"win\">" : "<td>", c, "</td>")
        end
        print(io, "</tr>")
    end
    print(io, "</tbody></table>")
    return HTML(String(take!(io)))
end

nowin(n) = falses(n)

# ## Both laws share the harmonic fixed point
#
# Diffusion drives ``\eta \to 0``, and by Lemma 1 of the paper ``\eta`` vanishes exactly at
# the harmonic extension. The two laws therefore select the same formation, which has to
# hold before any comparison of how they get there means anything.

scenario = coordination_scenario(:chain; size_parameter = 32)

plan = direct_plan(scenario)
targets = [2.2 -2.2; 0.0 0.0]
qstar = harmonic_reference(plan, targets)

eta = zeros(length(qstar))
diffusion_residual!(eta, scenario, qstar, targets)
@printf("residual at Direct's solution: ‖η(q*)‖∞ = %.2e\n", norm(eta, Inf))

# The cached tree solve is also checked against the library's canonical
# `harmonic_extension`, which rebuilds the blocks and refactors on every call.

oracle = CoordinationBenchmarks.harmonic_reference_oracle(scenario, targets)
@printf("cached tree solve versus harmonic_extension: %.2e\n", norm(qstar - oracle, Inf))

# ## Why the transients differ
#
# Substituting each law into the tracking error gives
#
# ```math
# \dot e_{\mathrm{Diffusion}} = -k\mathcal{H}e \qquad\text{versus}\qquad
# \dot e_{\mathrm{Direct}} = -k e .
# ```
#
# Diffusion inherits the sheaf spectrum and Direct does not. Two quantities follow, and
# between them they explain every number on this page.

spec = spectral_summary(scenario)
@printf("λ_min = %.4f, λ_max = %.3f, κ(ℋ) = %.1f\n",
        spec.minimum, spec.maximum, spec.condition)

# ``\lambda_{\min}`` is set by how much of the fleet senses a target, which is the Dirichlet
# boundary, and ``\lambda_{\max}`` by graph degree. Diffusion is squeezed at both ends: its
# slowest mode decays at ``k\lambda_{\min}``, while a zero-order hold caps its gain at
# ``2/(\Delta t\,\lambda_{\max})``. Its settling time therefore scales with the ratio.

DT = 0.02
@printf("deployable gain ceiling: Diffusion %.1f, Direct %.1f\n",
        stable_gain_ceiling(:diffusion, scenario, DT),
        stable_gain_ceiling(:direct, scenario, DT))

# ## What is counted, and to whom
#
# Fairness is the whole point of this page, so the accounting is stated before the results.
#
# Both laws run from the same initial state, against the same targets, with the same
# integrator, step size, ``g^{+}`` path and plant. The only line that differs is the one
# producing ``u``. Each law runs at `safety` times **its own** discrete stability ceiling,
# because a shared gain would not be neutral: Diffusion's per-mode gain is ``k\lambda_i``.
#
# Direct's reported time includes **the entire harmonic extension**: assembling the blocks,
# the symbolic analysis, the numeric factorization, the workspace allocation, and every
# triangular solve. None of it is treated as a free precomputation. The stage tables below
# show exactly where that time goes, and it is the majority of Direct's budget rather than a
# rounding error.
#
# Diffusion has no corresponding setup, because there is no matrix to factor. Its cost is
# the ticks it takes to converge.
#
# Communication is charged separately, in half-duplex slots, since wall-clock on one machine
# says nothing about radio traffic. Diffusion pays ``\Delta`` slots on every tick; Direct
# pays one gather and one scatter over the elimination tree, once per solve.

# ## Three scenarios
#
# The two laws differ in one number, closed-loop bandwidth, but that shows up differently
# depending on the task. Three tasks on three formations, each run under both laws.

track = coordination_scenario(:chain; size_parameter = 16)
@printf("Diffusion bandwidth kλ_min = %.3f rad/s\n", tracking_bandwidth(:diffusion, track; dt = DT))
@printf("Direct bandwidth    k      = %.3f rad/s\n", tracking_bandwidth(:direct, track; dt = DT))

# Against a moving reference the two error systems are
#
# ```math
# \dot e_{\mathrm{diff}} = -k\mathcal{H}e - \dot q^\star ,
# \qquad
# \dot e_{\mathrm{dir}} = -k e - \dot q^\star ,
# ```
#
# two stable filters driven by the same ``\dot q^\star``. Neither reaches zero error while
# the targets move. Each settles to a lag of order ``\lVert\dot q^\star\rVert`` over its own
# bandwidth. Diffusion tracks perfectly well, from further behind.
#
# ### Scenario 1: formation assembly
#
# A 16-agent chain, scattered well away from its harmonic formation, with the targets held
# still. This isolates one question: how long does it take to arrive at all?

assembly = coordination_scenario(:chain; size_parameter = 16)
asm_diff = rollout_coordination(:diffusion, assembly; dt = DT, mode = :static,
                                horizon = 40.0, spread = 1.6)
asm_dir = rollout_coordination(:direct, assembly; dt = DT, mode = :static,
                               horizon = 40.0, spread = 1.6)

settle_tick(r; frac = 0.01) =
    something(findfirst(<=(frac * r.errors[1]), r.errors), length(r.errors))

@printf("ticks to within 1%% of the formation: Diffusion %d, Direct %d\n",
        settle_tick(asm_diff), settle_tick(asm_dir))

# Both panels below are at the instant Direct has finished.

plot(plot(asm_diff, settle_tick(asm_dir) + 2), plot(asm_dir, settle_tick(asm_dir) + 2);
     layout = (1, 2), size = (1280, 480))

# Filled markers are agents coloured by the target they escort, open circles are the
# harmonic references delivered to them, and stars are the targets. Direct is on station.
# Diffusion has barely left its scattered start, because information spreads one hop per
# tick along the chain.

animate_coordination(asm_diff, asm_dir;
    filename = "scenario_assembly.gif", frames = 1:700, frame_step = 8, fps = 15)
nothing # hide

# ![Formation assembly](scenario_assembly.gif)

plot(asm_diff.times, asm_diff.errors; yscale = :log10, color = :steelblue, linewidth = 2,
     label = "Diffusion", xlabel = "time (s)", ylabel = "‖q − q*‖ / √N",
     title = "Scenario 1: formation assembly (16-agent chain, static targets)",
     legend = :topright)
plot!(asm_dir.times, asm_dir.errors; color = :darkorange, linewidth = 2, label = "Direct")

# ### Scenario 2: steady orbit tracking
#
# A 4x4 grid escorting two targets that circle at a constant rate. The reference never
# stops moving, so neither law converges. Each holds a steady lag set by its bandwidth.

orbit = coordination_scenario(:grid; size_parameter = 4)
orb_diff = rollout_coordination(:diffusion, orbit; dt = DT, mode = :orbit)
orb_dir = rollout_coordination(:direct, orbit; dt = DT, mode = :orbit)

@printf("steady lag: Diffusion %.3e, Direct %.3e, ratio %.1f\n",
        orb_diff.errors[end], orb_dir.errors[end],
        orb_diff.errors[end] / orb_dir.errors[end])

plot(plot(orb_diff, length(orb_diff.times)), plot(orb_dir, length(orb_dir.times));
     layout = (1, 2), size = (1280, 480))

# The gap between each agent and its own open circle is the tracking error. Under Diffusion
# the grid trails its reference around the orbit; under Direct the two coincide.

animate_coordination(orb_diff, orb_dir;
    filename = "scenario_orbit.gif", frame_step = 6, fps = 15)
nothing # hide

# ![Steady orbit tracking](scenario_orbit.gif)

# How far behind depends on how fast the targets move. Sweeping the orbit period shows the
# quasi-static regime appearing for both laws: both lags fall in proportion to target speed
# and the ratio between them stays put.

speed_rows = map((12.0, 30.0, 90.0, 300.0, 900.0)) do period
    re = rollout_coordination(:diffusion, track; dt = DT, period)
    ra = rollout_coordination(:direct, track; dt = DT, period)
    tail(r) = maximum(r.errors[(end - 50):end])
    (period, 2pi / period, tail(re), tail(ra), tail(re) / tail(ra))
end

bench_table(["orbit period", "ω (rad/s)", "Diffusion lag", "Direct lag", "ratio"],
    [((@sprintf("%.0f s", r[1]), @sprintf("%.4f", r[2]), @sprintf("%.2e", r[3]),
       @sprintf("%.2e", r[4]), @sprintf("%.1f×", r[5])),
      (false, false, false, true, false)) for r in speed_rows])

# Direct's lag is small but not zero, and it shrinks with target speed exactly as
# Diffusion's does. Slow the targets enough and both laws hold the formation tightly. The
# difference is only ever how slow is slow enough.
#
# ### Scenario 3: target reassignment
#
# A 16-agent ring orbiting normally until, half-way through the run, the targets jump to
# opposite stations. That is a step in the reference. Both laws take the same instantaneous
# hit, so this measures recovery rate rather than steady-state offset.

swap = coordination_scenario(:ring; size_parameter = 16)
mvr_diff = rollout_coordination(:diffusion, swap; dt = DT, mode = :maneuver)
mvr_dir = rollout_coordination(:direct, swap; dt = DT, mode = :maneuver)

function recovery(r; threshold = 0.15)
    peak = argmax(r.errors)
    back = findfirst(<=(threshold), r.errors[peak:end])
    return (r.errors[peak], back === nothing ? -1 : back)
end

pk_d, back_d = recovery(mvr_diff)
pk_r, back_r = recovery(mvr_dir)
@printf("peak error: Diffusion %.2f, Direct %.2f (the step hits both equally)\n", pk_d, pk_r)
@printf("ticks back under 0.15: Diffusion %d, Direct %d\n", back_d, back_r)

plot(plot(mvr_diff, argmax(mvr_diff.errors) + 8), plot(mvr_dir, argmax(mvr_dir.errors) + 8);
     layout = (1, 2), size = (1280, 480))

# Direct has re-formed on the new stations while Diffusion is still strung out between the
# old and new ones.

animate_coordination(mvr_diff, mvr_dir;
    filename = "scenario_maneuver.gif", frame_step = 6, fps = 15)
nothing # hide

# ![Target reassignment](scenario_maneuver.gif)

plot(mvr_diff.times, mvr_diff.errors; yscale = :log10, color = :steelblue, linewidth = 2,
     label = "Diffusion", xlabel = "time (s)", ylabel = "‖q − q*‖ / √N",
     title = "Scenario 3: target reassignment (16-agent ring)", legend = :right)
plot!(mvr_dir.times, mvr_dir.errors; color = :darkorange, linewidth = 2, label = "Direct")
vline!([mvr_diff.times[argmax(mvr_diff.errors)]]; color = :gray45, linestyle = :dash,
       label = "reassignment")

# Both spike to the same height. The difference is entirely in the decay afterwards, which
# is the ``k\lambda_{\min}`` against ``k`` story again, now visible as a transient.
#
# ### Is the comparison filtered?
#
# A filter on one law and not the other would manufacture exactly this kind of lag gap.
# Neither law is filtered above. Both are memoryless proportional feedback on the same
# current target data, and neither carries filter state.
#
# To settle it, `rollout_coordination` can apply a first-order command lag
# ``\epsilon\dot v = -v + u``, ``\dot q = v``, identically to both laws.

filter_rows = map((0.0, 0.05, 0.15, 0.30)) do eps
    re = rollout_coordination(:diffusion, track; dt = DT, epsilon = eps)
    ra = rollout_coordination(:direct, track; dt = DT, epsilon = eps)
    tail(r) = maximum(r.errors[(end - 50):end])
    (eps, tail(re), tail(ra), tail(re) / tail(ra))
end

bench_table(["command filter", "Diffusion lag", "Direct lag", "ratio"],
    [((r[1] == 0 ? "none (both unfiltered)" : @sprintf("ε = %.2f (both)", r[1]),
       @sprintf("%.3e", r[2]), @sprintf("%.3e", r[3]), @sprintf("%.1f×", r[4])),
      (false, false, true, false)) for r in filter_rows])

# The ratio barely moves. Adding the same actuator lag to both costs each of them a little
# and the gap between them survives, so it is closed-loop bandwidth rather than an artefact
# of filtering one side.

# ## Cost to reach the formation
#
# With the targets held still, the question becomes how much each law spends to arrive.
# Per-tick cost alone misleads: Diffusion's tick is cheap precisely because it does not
# solve anything, so it needs ``O(\kappa)`` of them, while Direct is at ``q^\star`` after
# one solve. `settle_to_formation` charges each law everything it spends until the fleet is
# within 0.1% of the harmonic formation.

diffusion_run = settle_to_formation(:diffusion, scenario)
direct_run = settle_to_formation(:direct, scenario)

bench_table(["law", "ticks", "compute", "slots"],
    [(("Diffusion", comma(diffusion_run.ticks), @sprintf("%.1f µs", 1e6diffusion_run.seconds),
       comma(diffusion_run.slots)), (false, false, false, false)),
     (("Direct", comma(direct_run.ticks), @sprintf("%.1f µs", 1e6direct_run.seconds),
       comma(direct_run.slots)), (false, true, true, true))])

# The settling transients plot directly.

plot(diffusion_run; title = "Settling on a 32-agent chain (κ ≈ 441)")
plot!(direct_run)

# ## Every stage of both algorithms
#
# This section takes both algorithms apart from the first setup step to the last plant
# integration, and measures each stage along every axis that decides whether a law can be
# deployed: time and its spread, heap allocation, arithmetic, matrix entries touched, and
# the memory an agent must hold between calls.
#
# Stage counts are driven by each law's actual settling tick count on this scenario, so the
# per-stage costs sum to the end-to-end cost rather than to an arbitrary horizon. With the
# targets fixed, Direct solves once.

profile_scenario = coordination_scenario(:grid; size_parameter = 6)
diff_settle = settle_to_formation(:diffusion, profile_scenario)
dir_settle = settle_to_formation(:direct, profile_scenario)

diff_stages = stage_profile(:diffusion, profile_scenario; ticks = diff_settle.ticks, solves = 0)
dir_stages = stage_profile(:direct, profile_scenario; ticks = dir_settle.ticks, solves = 1)

@printf("6x6 grid: Diffusion converges in %d ticks, Direct in %d\n",
        diff_settle.ticks, dir_settle.ticks)

# Timing is reported as a minimum and a 90th percentile. The spread between them is the
# honest measure of how much garbage collection and scheduler noise a stage suffers, which a
# bare minimum hides. Allocation is measured inside a function, because at global scope the
# boxing of non-constant globals is attributed to the call and reports allocations that the
# stage never performs.

stage_rows(stages) = [((string(st.stage), st.description, st.cadence, comma(st.calls),
        @sprintf("%.3f µs", 1e6st.t_min), @sprintf("%.3f µs", 1e6st.t_p90),
        @sprintf("%.1f µs", 1e6 * st.calls * st.t_min),
        fmt_bytes(st.bytes), comma(st.allocs), comma(st.flops), comma(st.nonzeros),
        fmt_bytes(st.resident)), nowin(12)) for st in stages]

stage_headers = ["stage", "what it does", "cadence", "calls", "min", "p90", "total",
    "alloc/call", "allocs", "flop/call", "nonzeros", "per-agent memory"]

# **Diffusion**, top to bottom. There is no setup: the loop is a neighbour sum and two
# vector updates.

bench_table(stage_headers, stage_rows(diff_stages))

# **Direct**, top to bottom. The first four stages are the harmonic extension, and they are
# charged in full.

bench_table(stage_headers, stage_rows(dir_stages))

diff_tot = stage_totals(diff_stages)
dir_tot = stage_totals(dir_stages)
setup_share = 100 * sum(st.t_min for st in dir_stages if st.cadence == "once") /
              sum(st.calls * st.t_min for st in dir_stages)

bench_table(["law", "total compute", "heap allocated", "allocations", "arithmetic",
             "per-agent memory", "fleet memory"],
    [(("Diffusion", @sprintf("%.1f µs", 1e6diff_tot.seconds), fmt_bytes(diff_tot.bytes),
       comma(diff_tot.allocs), @sprintf("%.2f Mflop", diff_tot.flops / 1e6),
       fmt_bytes(diff_tot.resident), fmt_bytes(diff_tot.resident_fleet)),
      (false, false, true, true, false, true, true)),
     (("Direct", @sprintf("%.1f µs", 1e6dir_tot.seconds), fmt_bytes(dir_tot.bytes),
       comma(dir_tot.allocs), @sprintf("%.2f Mflop", dir_tot.flops / 1e6),
       fmt_bytes(dir_tot.resident), fmt_bytes(dir_tot.resident_fleet)),
      (false, true, false, false, true, false, false))])

@printf("The harmonic extension is %.0f%% of Direct's total compute on this scenario.\n",
        setup_share)

# Four things fall out of these tables, and they do not all point the same way.
#
# **The harmonic extension dominates Direct's cost.** It is the large majority of Direct's
# budget here, not a negligible precomputation. Direct wins the end-to-end comparison in
# spite of paying for the factorization, not by having it excused. A per-tick view would
# report the opposite, because it amortises a one-time cost over a horizon Direct never
# runs for: with the targets fixed, Direct converges in a handful of ticks.
#
# **Diffusion performs far more arithmetic and is still slower in wall-clock.** The flop
# counts and the timings tell opposite stories. Diffusion's operations are cheap neighbour
# sums, of which it needs ``O(\kappa)`` rounds; Direct's are dense panel operations on
# frontal matrices, of which it needs one pass.
#
# **Diffusion allocates nothing at all.** Zero bytes across every tick, because the
# iteration is matrix-free and its iterate is the fleet's own position. It can run in a
# hard-real-time loop with the collector disabled, as written. Direct's allocation happens
# in the setup stages, which is exactly where a deployment would want it: off the control
# path, before the mission starts.
#
# **Per-agent memory is the treewidth currency again.** A diffusion agent holds its own
# state and its neighbours', which is ``O(\text{degree})`` and flat in fleet size. A Direct
# worker holds its own frontal panel, which is ``O(\text{treewidth}^2)``. Neither holds a
# fleet-wide object, and reporting the whole factor against one agent's working set would
# compare quantities that live on different machines.

memory_rows = map(((:chain, 32), (:chain, 128), (:grid, 6), (:grid, 10), (:rgg, 128))) do (nm, sp)
    sc = coordination_scenario(nm; size_parameter = sp)
    ds = stage_profile(:diffusion, sc; ticks = 1, solves = 0, reps = 5)
    rs = stage_profile(:direct, sc; ticks = 1, solves = 1, reps = 5)
    st = tree_statistics(direct_plan(sc))
    (sc, st, stage_totals(ds), stage_totals(rs))
end

bench_table(["topology", "agents", "treewidth", "Diffusion per agent", "Direct per worker",
             "Direct fleet total"],
    [((sc.label, comma(sc.nagents), string(st.treewidth), fmt_bytes(dt.resident),
       fmt_bytes(rt.resident), fmt_bytes(rt.resident_fleet)),
      (false, false, false, true, false, false))
     for (sc, st, dt, rt) in memory_rows])

# On a chain both are flat in fleet size, and the difference between them is a few dozen
# bytes. On a random geometric graph the separators are wide and a Direct worker's panel
# grows with their square. That is the same structural axis that governs its fill and its
# arithmetic, so it is a real cost of the method rather than an artefact of accounting.
# The absolute magnitudes are small enough that this axis is unlikely to decide anything at
# these fleet sizes; it would begin to matter on microcontroller-class nodes or at fleets
# well beyond the range benchmarked here.

# ### What it would cost on real hardware
#
# Every timing above is single-core, which is symmetric but describes neither deployment.
# On a real fleet both laws parallelize, and they parallelize very differently.
#
# Diffusion parallelizes almost perfectly across agents: every agent forms its own
# disagreement from data it already holds, so a tick costs the busiest agent rather than the
# sum over the fleet. Direct parallelizes across the elimination tree, but only down to its
# critical path, since independent subtrees factor and substitute concurrently while a chain
# of supernodes must run in sequence.
#
# Both numbers below are derived from the sheaf and the factor structure, so neither depends
# on this machine.

cp_diff = critical_path_profile(:diffusion, profile_scenario; ticks = diff_settle.ticks)
cp_dir = critical_path_profile(:direct, profile_scenario; ticks = dir_settle.ticks, solves = 1)

bench_table(["law", "serial arithmetic", "critical-path arithmetic", "available speedup",
             "synchronous depth"],
    [(("Diffusion", comma(cp_diff.serial_flops), comma(cp_diff.parallel_flops),
       @sprintf("%.1f×", cp_diff.speedup), comma(cp_diff.depth)),
      (false, false, false, true, false)),
     (("Direct", comma(cp_dir.serial_flops), comma(cp_dir.parallel_flops),
       @sprintf("%.1f×", cp_dir.speedup), comma(cp_dir.depth)),
      (false, true, true, false, true))])

# This is the most important correction on the page, and it cuts against Direct. Diffusion
# has far more parallelism available, because its work is spread across every agent. Direct
# has much less, because the elimination tree serializes. Single-core timing therefore
# overstates Direct's advantage.
#
# After parallelizing both, Direct still leads on arithmetic along the critical path, and it
# leads by a much wider margin on synchronous depth: a few tree levels traversed twice,
# against one round per tick for the whole settling time. Latency, not arithmetic, is what
# separates them on real hardware.

# ## The formation families
#
# The benchmark sweeps twenty topology families, chosen to span the two structural axes that
# matter: trees and lattices at one end, dense blocks and expanders at the other. Every
# scenario is a real `EuclideanSheaf` with identity restriction maps and stalk dimension 2.
# Blue circles are free agents, red diamonds are pinned targets, and dashed links are the
# agent-to-target sensing edges that supply the Dirichlet boundary.

atlas_group(specs) = plot(
    [plot(coordination_scenario(n; size_parameter = sp)) for (n, sp) in specs]...;
    layout = (cld(length(specs), 2), 2), size = (900, 450 * cld(length(specs), 2)),
    titlefontsize = 11, top_margin = 4Plots.mm, bottom_margin = 4Plots.mm)

# Trees have treewidth 1 and no off-diagonal fill, but very different degree profiles.

atlas_group([(:chain, 14), (:star, 14), (:tree, 4), (:caterpillar, 7)])

# Sparse lattices: treewidth climbs slowly with the width of the strip or sheet.

atlas_group([(:ring, 14), (:ladder, 7), (:prism, 7), (:grid, 4)])

# Geometric and networked graphs: proximity models and the classic random networks.

atlas_group([(:torus, 4), (:grid3d, 3), (:rgg, 16), (:smallworld, 16),
             (:scalefree, 16), (:expander, 16)])

# Hubs and dense blocks, where separators get wide and the factor fills in.

atlas_group([(:wheel, 14), (:bipartite, 7), (:barbell, 6), (:lollipop, 6),
             (:complete, 10), (:twoclique, 8)])

# Each panel is titled with its agent count and treewidth. Chain, star, binary tree and
# caterpillar are all trees despite very different degree profiles. Complete bipartite,
# complete and the two-clique construction carry the widest separators. That axis turns out
# to be independent of ``\kappa``, which is the point developed below.

# ## Sweeping the topology families
#
# `benchmark_coordination` runs both laws across all twenty families, three to four sizes
# each, and both target-coverage regimes. Sparse coverage pins two agents at the extremes;
# full coverage gives every agent a target edge. Coverage sets ``\lambda_{\min}``, so
# sweeping it as well as the graph turns the comparison from a point into a map. Every
# cached solve is verified against `harmonic_extension` before timing.

result = benchmark_coordination()
@printf("%d scenarios, worst oracle residual %.2e\n",
        length(result), maximum(r.oracle_residual for r in result))

# ### Ticks to reach the formation

plot(result, :settling)

# ### Communication
#
# Diffusion pays ``\Delta`` slots on every tick until information has crossed the formation,
# so its bill is ``\text{ticks}\times\Delta`` and inherits the ``O(\kappa)``. Direct pays one
# gather up the elimination tree and one scatter back down, ``2\,\mathrm{cp}(T)`` slots,
# once.
#
# One caveat is enforced in the model rather than left to the reader.
# ``2\,\mathrm{cp}(T)`` counts only messages crossing between supernodes, which presumes each
# supernode sits on its own worker. On dense graphs the elimination tree collapses to one or
# two supernodes, meaning a single machine holds the whole factor. That is a centralized
# solve, not a distributed one, and would otherwise be charged a near-zero slot count. Those
# cases are charged the centralized ``2n`` pigeonhole bound instead.

plot(result, :communication)

# ### The advantage, in both currencies

plot(result, :speedup)

# ### The map
#
# Plotting every scenario against both structural axes at once, with treewidth on the
# horizontal and ``\kappa`` on the vertical, shows where each law owns the space.

plot(result, :landscape; size = (860, 560))

# The separation is horizontal rather than diagonal. Diffusion's wins form a band at low
# ``\kappa`` running across the entire treewidth range. Which law is faster is decided by
# conditioning; treewidth only sets how much Direct's win costs it when it does win.

# ### The full sweep

summary_rows = map(unique(r.name for r in result)) do fam
    fr = filter(r -> r.name == fam && r.coverage == :sparse, result.rows)
    last(sort(fr; by = r -> r.agents))
end

bench_table(["topology", "agents", "κ(ℋ)", "treewidth", "Diffusion ticks", "Direct ticks",
             "Diffusion compute", "Direct compute", "Direct faster", "Direct fewer slots"],
    [((r.label, comma(r.agents), @sprintf("%.1f", r.condition), string(r.treewidth),
       comma(r.diffusion.ticks), comma(r.direct.ticks),
       @sprintf("%.1f µs", 1e6r.diffusion.seconds), @sprintf("%.1f µs", 1e6r.direct.seconds),
       @sprintf("%.1f×", r.speedup), @sprintf("%.0f×", r.slot_ratio)),
      (false, false, false, false, false, r.speedup >= 1, r.speedup < 1, r.speedup >= 1,
       false, false)) for r in summary_rows])

# One row per family at its largest sparsely-covered size, with the faster law's tick count
# shaded. Here is how the two coverage regimes split.

sparse_rows = filter(r -> r.coverage == :sparse, result.rows)
full_rows = filter(r -> r.coverage == :full, result.rows)
@printf("sparse coverage: Direct faster in %d of %d, median %.1f×\n",
        count(r -> r.speedup >= 1, sparse_rows), length(sparse_rows),
        sort([r.speedup for r in sparse_rows])[end ÷ 2])
@printf("full coverage:   Direct faster in %d of %d, median %.2f×\n",
        count(r -> r.speedup >= 1, full_rows), length(full_rows),
        sort([r.speedup for r in full_rows])[end ÷ 2])

# ## Two different currencies: conditioning and treewidth
#
# The two laws are billed by different structural properties of the same graph, and
# measuring both is what turns the comparison into a map.
#
# Diffusion's cost is spectral, governed by
# ``\kappa(\mathcal H) = \lambda_{\max}/\lambda_{\min}``: degree at one end, target coverage
# at the other. Direct's cost is combinatorial, governed by the treewidth of the formation
# graph, which sets the width of the widest frontal matrix and therefore the fill, the
# arithmetic per solve, and the memory each worker holds.

structure_rows = map(((:chain, 64), (:star, 64), (:grid, 8), (:rgg, 64),
                      (:twoclique, 32), (:complete, 36), (:expander, 64))) do (name, sp)
    s = coordination_scenario(name; size_parameter = sp)
    t = tree_statistics(direct_plan(s))
    (s.label, s.nagents, spectral_summary(s).condition, t.treewidth, t.fill,
     t.offdiagonal_fill, t.depth)
end

bench_table(["topology", "agents", "κ(ℋ)", "treewidth", "factor fill", "off-diagonal",
             "tree depth"],
    [((r[1], comma(r[2]), @sprintf("%.1f", r[3]), string(r[4]), comma(r[5]),
       comma(r[6]), string(r[7])), nowin(7)) for r in structure_rows])

# Chain and star are both trees, with treewidth 1, yet their ``\kappa`` differs by a factor
# of five. The random geometric graph carries a wide separator at a middling ``\kappa``, and
# the two-clique graph has the worst treewidth in the table with a ``\kappa`` comparable to
# the chain's.
#
# The off-diagonal column is worth reading against the fill column. On the complete graph
# the elimination tree is a single supernode, so every stored entry sits in the diagonal
# block and the off-diagonal count is zero, which would look like a free factorization if
# the diagonal blocks were not counted. They are.
#
# That independence is what makes the comparison interesting rather than circular:
#
# | | low ``\kappa`` | high ``\kappa`` |
# |---|---|---|
# | **low treewidth** | both cheap | Direct's best case: trivial factor, starved Diffusion |
# | **high treewidth** | Direct's worst case: the expander below | expensive for Direct, worse for Diffusion |
#
# Tree depth matters separately: it sets ``\mathrm{cp}(T)`` and therefore Direct's slot
# count. The star has depth 2 and clears in a handful of slots, but its factor lives on a
# single worker, so a shallow tree means little parallelism rather than a fast distributed
# solve.

# ## Where Diffusion wins
#
# Every table so far reports sparsely covered formations, and Diffusion loses all of them.
# That is not the whole picture. Across the full sweep Diffusion takes a substantial
# minority of scenarios, and they are not scattered at random.

diffusion_wins = sort(filter(r -> r.speedup < 1, result.rows); by = r -> r.speedup)
direct_wins = filter(r -> r.speedup >= 1, result.rows)

@printf("Diffusion faster in %d of %d scenarios\n", length(diffusion_wins), length(result))
@printf("  Diffusion's wins: median κ %.1f, median treewidth %d, %d%% fully covered\n",
        sort([r.condition for r in diffusion_wins])[end ÷ 2],
        sort([r.treewidth for r in diffusion_wins])[end ÷ 2],
        round(Int, 100count(r -> r.coverage === :full, diffusion_wins) / length(diffusion_wins)))
@printf("  Direct's wins:    median κ %.1f, median treewidth %d, %d%% fully covered\n",
        sort([r.condition for r in direct_wins])[end ÷ 2],
        sort([r.treewidth for r in direct_wins])[end ÷ 2],
        round(Int, 100count(r -> r.coverage === :full, direct_wins) / length(direct_wins)))

# The eight widest margins:

bench_table(["topology", "agents", "coverage", "κ(ℋ)", "treewidth", "Diffusion", "Direct",
             "Diffusion faster by"],
    [((r.label, comma(r.agents), string(r.coverage), @sprintf("%.1f", r.condition),
       string(r.treewidth), @sprintf("%.1f µs", 1e6r.diffusion.seconds),
       @sprintf("%.1f µs", 1e6r.direct.seconds), @sprintf("%.2f×", 1 / r.speedup)),
      (false, false, false, false, false, true, false, false))
     for r in diffusion_wins[1:min(8, end)]])

# Every one is fully covered and every one has ``\kappa`` in single digits. That is the
# rule: give every agent its own target edge, ``\lambda_{\min}`` snaps to 1, ``\kappa``
# collapses to ``\lambda_{\max}(L) + 1``, and Diffusion converges in a handful of rounds, at
# which point Direct's factorization is no longer amortised over enough solves to pay for
# itself. Note what that regime is, though: every agent already senses its own target, so
# there is barely a coordination problem left for the sheaf to solve.
#
# ### Constructing the hard case deliberately
#
# From the two-currency table, the way to beat Direct on structure rather than on coverage
# is high treewidth at bounded degree: keep ``\lambda_{\max}`` small so Diffusion's gain
# ceiling stays high, while denying the elimination tree any small separators. A
# bounded-degree expander does exactly that. Its spectral gap does not decay, so ``\kappa``
# is flat in ``n`` while treewidth and fill grow with the fleet.

adversarial_rows = map(((:expander, 32), (:expander, 64), (:expander, 128))) do (name, sp)
    sc = coordination_scenario(name; size_parameter = sp, coverage = :full)
    st = tree_statistics(direct_plan(sc))
    dif = settle_to_formation(:diffusion, sc)
    dir = settle_to_formation(:direct, sc)
    (sc, st, dif, dir)
end

bench_table(["n", "κ(ℋ)", "treewidth", "fill", "Diffusion ticks", "Direct ticks",
             "Diffusion", "Direct", "Diffusion slots", "Direct slots", "winner"],
    [((comma(sc.nagents), @sprintf("%.2f", spectral_summary(sc).condition),
       string(st.treewidth), comma(st.fill), comma(dif.ticks), comma(dir.ticks),
       @sprintf("%.1f µs", 1e6dif.seconds), @sprintf("%.1f µs", 1e6dir.seconds),
       comma(dif.slots), comma(dir.slots),
       dif.seconds < dir.seconds ? @sprintf("Diffusion, %.2f×", dir.seconds / dif.seconds) :
                                   @sprintf("Direct, %.1f×", dif.seconds / dir.seconds)),
      (false, false, false, false, false, false,
       dif.seconds < dir.seconds, dif.seconds >= dir.seconds, false, true,
       dif.seconds < dir.seconds))
     for (sc, st, dif, dir) in adversarial_rows])

# The construction works, and it is the only family in the sweep where Diffusion wins at
# every size rather than only under full coverage. The margins are low single digits, and at
# this scale wall-clock is noisy enough that these are minimum-of-runs, so the direction is
# the result and the exact factor is indicative. The communication ledger still runs the
# other way, because a diffusion pays ``\Delta`` slots every round no matter how few rounds
# it needs.
#
# For contrast, two overlapping cliques is the intuitive candidate, dense enough that
# information crosses in one hop, and it does not work despite carrying the widest
# separators in the study. Density drives ``\lambda_{\max}`` up, which throttles Diffusion's
# gain ceiling, so a clique ends up as badly conditioned as a chain.

clique = coordination_scenario(:twoclique; size_parameter = 32)
clique_spec = spectral_summary(clique)
@printf("two overlapping cliques, n = %d: λ_min = %.3f, λ_max = %.1f, κ = %.1f\n",
        clique.nagents, clique_spec.minimum, clique_spec.maximum, clique_spec.condition)

# The distinction is between mixing and convergence. Information does reach everyone in one
# hop, but the iteration's rate is set by ``\kappa``, and ``\lambda_{\min}`` in a pinned
# Dirichlet problem is governed by how much boundary there is rather than by how connected
# the interior is.

# ## How much to trust these numbers
#
# Tick counts are exact. The trajectories are deterministic and the counts agree with the
# Richardson rate predicted from the spectrum: the slowest mode contracts by
# ``2\,\mathrm{safety}/\kappa`` per tick for Diffusion and by ``2\,\mathrm{safety}`` for
# Direct. Measured counts land below the worst-case prediction, the shortfall being initial
# conditions that do not fully excite the slowest mode, so they are a property of this
# initial condition as well as of the sheaf.
#
# Timings are indicative. Direct's cost is dominated by the factorization, which allocates
# heavily, and a garbage-collection pause landing inside it has been measured to inflate a
# single run by over 100 times. Every timing is a minimum of seven runs for that reason, and
# the settling loop runs twice, once instrumented to establish the tick count and once clean
# to take the clock, so the convergence check is never charged to the law that executes it
# more often. Read a reported speedup as one significant figure.
#
# Slot counts, flop counts, critical-path lengths and per-agent memory are exact and
# machine-independent, being properties of the graph and the elimination tree.
#
# ## Scope
#
# Direct's advantage is contingent on the formation graph having bounded treewidth. Physical
# fleets fly proximity graphs, which have good separators; expanders deliberately do not and
# require long-range links that range-limited radios cannot form. Timings are single-machine
# and single-threaded, so the slot counts and the critical-path model rather than the
# wall-clock are the deployment currencies. Plants are single integrators, and each law runs
# at its own stability ceiling, so Direct commands a more aggressive transient that a
# saturating actuator would clip.
