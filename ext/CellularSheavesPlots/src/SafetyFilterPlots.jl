module SafetyFilterPlots

using Plots

using LinearAlgebra
using Printf
using CellularSheaves.ControlSheaves.SafetyFilters: FilterRollout, barrier_window,
                                                   saturation_window,
                                                   relaxation_window

# Shared with LayeredEscortPlots so the safety figures read as part of the same series.
const AGENT_COLOR = :steelblue
const TARGET_COLOR = :black
const GUIDE_COLOR = :gray80
const REFERENCE_COLOR = :purple
const BOUND_COLOR = :black
const BREACH_COLOR = :firebrick
const FILTERED_COLOR = :seagreen
const UNFILTERED_COLOR = :darkorange

export comparison_frame, trace_panel, circle, GUIDE_COLOR, FILTERED_COLOR, UNFILTERED_COLOR

# ==========================
# Metric primitives
# ==========================

# By default, plot the whole rollout if `k` is not provided.
@recipe function f(res::FilterRollout, metric::Symbol; k = length(res.times))
    ts = res.times

    if metric == :separation
        xlabel --> "time (s)"
        ylabel --> "separation (m)"
        title --> "Closest Pair Separation"
        xlims --> (0, ts[end])
        ylims --> (0, maximum(filter(isfinite, res.separation)) * 1.15 + 0.01)
        @series begin
            color := GUIDE_COLOR
            linestyle := :dash
            linewidth := 1
            label := "keep-out"
            [0.0, ts[end]], [res.separation_floor, res.separation_floor]
        end
        linewidth --> 1.5
        color --> AGENT_COLOR
        label --> false
        return ts[1:k], res.separation[1:k]

    elseif metric == :clearance
        res.clearance === nothing && error("this rollout has no obstacle to clear")
        xlabel --> "time (s)"
        ylabel --> "clearance (m)"
        title --> "Clearance to Obstacle"
        xlims --> (0, ts[end])
        ylims --> (0, maximum(res.clearance) * 1.15 + 0.01)
        @series begin
            color := GUIDE_COLOR
            linestyle := :dash
            linewidth := 1
            label := "keep-out"
            [0.0, ts[end]], [res.keep_out, res.keep_out]
        end
        linewidth --> 1.5
        color --> AGENT_COLOR
        label --> false
        return ts[1:k], res.clearance[1:k]

    elseif metric == :command
        xlabel --> "time (s)"
        ylabel --> "peak ‖u‖"
        title --> "Commanded Magnitude"
        xlims --> (0, ts[end])
        top = max(maximum(res.command), isfinite(res.command_cap) ? res.command_cap : 0.0)
        ylims --> (0, top * 1.15 + 0.01)
        if isfinite(res.command_cap)
            @series begin
                color := GUIDE_COLOR
                linestyle := :dash
                linewidth := 1
                label := "u_max"
                [0.0, ts[end]], [res.command_cap, res.command_cap]
            end
        end
        linewidth --> 1.5
        color --> :orange
        label --> false
        return ts[1:k], res.command[1:k]

    elseif metric == :lyapunov
        xlabel --> "time (s)"
        ylabel --> "V(t)"
        title --> "Tracking Lyapunov Function"
        xlims --> (0, ts[end])
        # Shade where the relaxation was open. The certified inequality is
        # Vdot <= -c3 V + delta, so that is exactly the interval over which V is permitted to
        # rise — and attributing a rise to any other constraint would be guesswork.
        relaxed = relaxation_window(res)
        if relaxed !== nothing
            @series begin
                seriestype := :vspan
                color := :goldenrod
                alpha := 0.15
                label := "relaxation open"
                [relaxed[1], relaxed[2]]
            end
        end
        linewidth --> 1.5
        color --> REFERENCE_COLOR
        label --> false
        return ts[1:k], res.lyapunov[1:k]

    elseif metric == :relaxation
        xlabel --> "time (s)"
        ylabel --> "δ(t)"
        title --> "CLF Relaxation"
        xlims --> (0, ts[end])
        ylims --> (0, maximum(res.slack) * 1.2 + 1e-6)
        linewidth --> 1.5
        color --> :purple
        label --> false
        return ts[1:k], res.slack[1:k]

    elseif metric == :top_down
        agents = size(res.positions, 2)
        xs = vec(res.positions[:, :, 1])
        ys = vec(res.positions[:, :, 2])
        cx, cy = (minimum(xs) + maximum(xs)) / 2, (minimum(ys) + maximum(ys)) / 2
        span = max(maximum(xs) - minimum(xs), maximum(ys) - minimum(ys)) / 2
        span += max(span * 0.15, 0.3)

        aspect_ratio --> 1
        xlims --> (cx - span, cx + span)
        ylims --> (cy - span, cy + span)
        xlabel --> "x position (m)"
        ylabel --> "y position (m)"
        title --> "World Top-Down View (x-y Plane)"

        for i in 1:agents
            @series begin
                seriestype := :path
                marker := :circle
                markersize := 3
                alpha := 0.6
                linewidth := 1.4
                color := AGENT_COLOR
                label := false
                res.positions[1:k, i, 1], res.positions[1:k, i, 2]
            end
        end
        if agents > 2
            ring_x = [res.positions[k, i, 1] for i in 1:agents]
            ring_y = [res.positions[k, i, 2] for i in 1:agents]
            push!(ring_x, res.positions[k, 1, 1])
            push!(ring_y, res.positions[k, 1, 2])
            @series begin
                color := GUIDE_COLOR
                linestyle := :dot
                linewidth := 1
                label := false
                ring_x, ring_y
            end
        end
    else
        error("unknown metric $(metric); expected one of :top_down, :separation, " *
              ":clearance, :command, :lyapunov, :relaxation")
    end
end

# ==========================
# Comparison animation
# ==========================

function circle(centre, radius; n = 90)
    θ = range(0, 2π; length = n)
    return centre[1] .+ radius .* cos.(θ), centre[2] .+ radius .* sin.(θ)
end

"""
    frame_status(res, k, mode) -> (breached, caption)

Whether the quantity the given constraint family governs is currently violated, and the
caption reporting it. Returns `(false, "")` when the family has nothing to report.
"""
function frame_status(res::FilterRollout, k::Int, mode::Symbol)
    if mode === :safety
        margin = res.clearance === nothing ? res.separation[k] : res.clearance[k]
        bound = res.clearance === nothing ? res.separation_floor : res.keep_out
        (isfinite(margin) && isfinite(bound) && bound > 0) || return (false, "")
        return (margin < bound, @sprintf("%.2f m  (need %.2f)", margin, bound))
    elseif mode === :actuator
        isfinite(res.command_cap) ||
            return (false, @sprintf("‖u‖ = %.2f  (uncapped)", res.command[k]))
        return (res.command[k] > res.command_cap * (1 + 1e-6),
                @sprintf("‖u‖ = %.2f  (cap %.2f)", res.command[k], res.command_cap))
    elseif mode === :stability
        growing = k > 1 && res.lyapunov[k] > res.lyapunov[k - 1]
        return (growing, @sprintf("V = %.2e", res.lyapunov[k]))
    elseif mode === :contended
        # Both hard bounds at once, since this is the frame where they compete.
        margin = res.separation[k]
        breached = isfinite(margin) && margin < res.separation_floor
        capped = isfinite(res.command_cap) ?
                 @sprintf("  ‖u‖ %.2f/%.2f", res.command[k], res.command_cap) :
                 @sprintf("  ‖u‖ %.2f", res.command[k])
        return (breached, @sprintf("sep %.2f (need %.2f)%s", margin,
                                   res.separation_floor, capped))
    end
    return (false, "")
end

"""
    comparison_frame(res, k, limits, obstacle; mode = :safety)

One frame of a side-by-side comparison. What is drawn is decided by `mode`, the constraint
family the demonstration is about, so that a figure never asserts something its scenario
does not enforce: an actuator demonstration draws the commanded vector against its cap and
no keep-out region, a stability demonstration draws the airframe attitude, and only a safety
demonstration draws keep-out discs and flags a breach of them. `:contended` is the exception
that draws two families at once, because it exists to show them binding together.
"""
function comparison_frame(res::FilterRollout, k::Int, limits, obstacle;
                          mode::Symbol = :safety, command_reference::Float64 = NaN)
    breached, caption = frame_status(res, k, mode)
    edge = breached ? BREACH_COLOR : BOUND_COLOR
    p = plot(aspect_ratio = 1, xlims = limits[1], ylims = limits[2],
             xlabel = "x position (m)", ylabel = "y position (m)",
             title = res.label, legend = false)
    span = limits[1][2] - limits[1][1]

    if mode === :safety && obstacle !== nothing
        o = obstacle(res.times[k])
        plot!(p, circle(o, res.keep_out)...; color = edge, linewidth = 2,
              fill = (0, 0.10, edge))
        scatter!(p, [o[1]], [o[2]]; marker = :xcross, markersize = 7, color = TARGET_COLOR)
    end

    agents = size(res.positions, 2)
    tail = max(1, k - 60)
    for i in 1:agents
        plot!(p, res.positions[tail:k, i, 1], res.positions[tail:k, i, 2];
              color = AGENT_COLOR, linewidth = 1.4, alpha = 0.6)
        centre = res.positions[k, i, :]

        if mode === :safety && obstacle === nothing && isfinite(res.separation_floor) &&
           res.separation_floor > 0
            # Each agent carries half the pairwise keep-out, so the discs touch at contact.
            plot!(p, circle(centre, res.separation_floor / 2)...; color = edge,
                  linewidth = 2, fill = (0, 0.12, edge))
            scatter!(p, [centre[1]], [centre[2]]; marker = :circle, markersize = 4,
                     color = AGENT_COLOR)

        elseif mode === :stability && res.attitude !== nothing
            # The airframe as a bar tilted by its roll angle, so an attitude excursion shows
            # as the frame turning over rather than as a dot that barely moves.
            φ = res.attitude[k, i]
            arm = 0.10 * span
            d = (arm * cos(φ), arm * sin(φ))
            tips = ([centre[1] - d[1], centre[1] + d[1]],
                    [centre[2] - d[2], centre[2] + d[2]])
            plot!(p, tips...; color = edge, linewidth = 4)
            scatter!(p, tips...; marker = :circle, markersize = 5, color = edge)

        elseif mode === :actuator
            # The commanded vector itself, drawn against a dashed circle of radius u_max, so
            # that exceeding the cap is visible as the arrow leaving the circle.
            reference = isfinite(command_reference) ? command_reference :
                        (isfinite(res.command_cap) ? res.command_cap : maximum(res.command))
            scale = 0.20 * span / max(reference, eps())
            heading = k > 1 ? res.positions[k, i, :] - res.positions[k - 1, i, :] :
                      [1.0, 0.0]
            direction = norm(heading) > 1e-12 ? heading ./ norm(heading) : [1.0, 0.0]
            tip = centre .+ (scale * res.command[k]) .* direction
            if isfinite(res.command_cap)
                plot!(p, circle(centre, scale * res.command_cap)...; color = GUIDE_COLOR,
                      linestyle = :dash, linewidth = 1)
            end
            plot!(p, [centre[1], tip[1]], [centre[2], tip[2]]; color = edge,
                  linewidth = 3, arrow = true)
            scatter!(p, [centre[1]], [centre[2]]; marker = :circle, markersize = 4,
                     color = AGENT_COLOR)

        elseif mode === :contended
            # Both families drawn together: the keep-out disc each agent carries, and the
            # command it is allowed against the cap. This is the only frame where a reader
            # can see the two hard bounds binding at the same instant.
            if isfinite(res.separation_floor) && res.separation_floor > 0
                plot!(p, circle(centre, res.separation_floor / 2)...; color = edge,
                      linewidth = 2, fill = (0, 0.12, edge))
            end
            if isfinite(res.command_cap)
                reference = isfinite(command_reference) ? command_reference : res.command_cap
                scale = 0.16 * span / max(reference, eps())
                heading = k > 1 ? res.positions[k, i, :] - res.positions[k - 1, i, :] :
                          [1.0, 0.0]
                direction = norm(heading) > 1e-12 ? heading ./ norm(heading) : [1.0, 0.0]
                tip = centre .+ (scale * res.command[k]) .* direction
                plot!(p, circle(centre, scale * res.command_cap)...; color = :orange,
                      linestyle = :dash, linewidth = 1)
                plot!(p, [centre[1], tip[1]], [centre[2], tip[2]]; color = :orange,
                      linewidth = 3, arrow = true)
            end
            scatter!(p, [centre[1]], [centre[2]]; marker = :circle, markersize = 4,
                     color = AGENT_COLOR)

        else
            scatter!(p, [centre[1]], [centre[2]]; marker = :circle, markersize = 4,
                     color = AGENT_COLOR)
        end
    end

    if agents > 2
        ring_x = [res.positions[k, i, 1] for i in 1:agents]
        ring_y = [res.positions[k, i, 2] for i in 1:agents]
        push!(ring_x, res.positions[k, 1, 1])
        push!(ring_y, res.positions[k, 1, 2])
        plot!(p, ring_x, ring_y; color = GUIDE_COLOR, linestyle = :dot, linewidth = 1)
    end

    if !isempty(caption)
        annotate!(p, limits[1][1] + 0.5 * (limits[1][2] - limits[1][1]),
                  limits[2][2] - 0.08 * (limits[2][2] - limits[2][1]),
                  text(caption, 8, breached ? BREACH_COLOR : TARGET_COLOR))
    end
    return p
end

"""
    trace_panel(unfiltered, filtered, k, mode)

The scalar quantity the given family governs, drawn for both runs up to step `k` against its
bound. A stability trace is logarithmic, because a diverging run and a converging one separate
by many orders of magnitude and on a linear axis the converging one collapses onto zero.
"""
function trace_panel(unfiltered::FilterRollout, filtered::FilterRollout, k::Int,
                     mode::Symbol)
    ts = filtered.times
    if mode === :contended
        # A raster rather than a trace: with three families the question is not how one
        # scalar evolves but which constraints are acting, and when they coincide.
        cbf = filtered.barrier_active[1:k]
        cone = filtered.saturated[1:k]
        relax = filtered.slack[1:k] .> 1e-2
        together = cbf .& cone .& relax
        tk = ts[1:k]
        p = plot(xlabel = "time (s)", xlims = (0, ts[end]), ylims = (0.4, 3.6),
                 yticks = ([1.0, 2.0, 3.0], ["δ open", "cone", "barrier"]),
                 legend = false, title = "Which constraint family is acting")
        any(together) && vspan!(p, [tk[findfirst(together)], tk[findlast(together)]];
                                color = :goldenrod, alpha = 0.25)
        any(relax) && scatter!(p, tk[relax], fill(1.0, count(relax));
                               color = REFERENCE_COLOR, markersize = 2)
        any(cone) && scatter!(p, tk[cone], fill(2.0, count(cone));
                              color = :orange, markersize = 2)
        any(cbf) && scatter!(p, tk[cbf], fill(3.0, count(cbf));
                             color = AGENT_COLOR, markersize = 2)
        return p
    end
    if mode === :actuator
        series = (unfiltered.command, filtered.command)
        bound = filtered.command_cap
        label, scale = "peak ‖u‖", :identity
    elseif mode === :stability
        series = (max.(unfiltered.lyapunov, 1e-8), max.(filtered.lyapunov, 1e-8))
        bound = NaN
        label, scale = "V(t)", :log10
    elseif filtered.clearance !== nothing
        series = (unfiltered.clearance, filtered.clearance)
        bound = filtered.keep_out
        label, scale = "clearance (m)", :identity
    else
        series = (unfiltered.separation, filtered.separation)
        bound = filtered.separation_floor
        label, scale = "separation (m)", :identity
    end

    finite = filter(isfinite, vcat(series[1], series[2]))
    top = scale === :log10 ? maximum(finite) * 4 : maximum(finite) * 1.15 + 0.01
    base = scale === :log10 ? max(minimum(filter(>(0), finite)) / 4, 1e-9) : 0.0

    p = plot(ts[1:k], series[1][1:k]; color = UNFILTERED_COLOR, linewidth = 1.8,
             label = unfiltered.label, xlabel = "time (s)", ylabel = label,
             xlims = (0, ts[end]), ylims = (base, top), yscale = scale,
             legend = :topright)
    plot!(p, ts[1:k], series[2][1:k]; color = FILTERED_COLOR, linewidth = 1.8,
          label = filtered.label)
    if isfinite(bound) && bound > 0
        plot!(p, [0.0, ts[end]], [bound, bound]; color = GUIDE_COLOR, linestyle = :dash,
              linewidth = 1.5, label = "bound")
    end
    return p
end

end # module SafetyFilterPlots
