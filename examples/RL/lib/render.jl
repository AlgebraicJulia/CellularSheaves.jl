#
# render.jl — shared, Flux-free rendering for the sheaf-track demonstrations.
#
# Assumes `sheaf_rl.jl` is already included (`using .SheafRL`). Reuses the library's
# env (sheaf, targets, drift, base control, ball projection) so a rollout here is identical to
# the one used in training/eval — it just additionally records positions for plotting.
#
# Produces, for a trained model: a top-down 3-panel animation (analytic | drift-oracle | sheaf+RL),
# a tracking-error plot, and static top-down + 3-D trajectory plots. `render_all` does all of them.
#
using Plots, JLD2, Random, Printf, Statistics, LinearAlgebra
gr()

# A plain feed-forward evaluation of the saved weights (no Flux dependency).
function load_actor(model)
    Ws, bs = model["Ws"], model["bs"]
    function actor(o)
        h = Float32.(o)
        for k in 1:length(Ws)-1
            h = max.(0f0, Ws[k]*h .+ bs[k])
        end
        tanh.(Ws[end]*h .+ bs[end])
    end
end

# Rebuild the exact training Config from a saved model, with optional overrides for the maps/cap
# that older model files may not carry.
function config_from_model(model; rotate_consensus=nothing, planar_agents=nothing, actuator_cap=nothing)
    Config(
        drift_amp        = model["DRIFT"],
        drift_freq       = model["DRIFTW"],
        residual_scale   = model["RESS"],
        gain             = model["KGAIN"],
        dt               = model["DT"],
        horizon          = model["K"],
        actuator_cap     = actuator_cap     === nothing ? get(model, "UCAP", Inf) : actuator_cap,
        rotate_consensus = rotate_consensus === nothing ? get(model, "rotate_consensus", false) : rotate_consensus,
        planar_agents    = planar_agents    === nothing ? Set(get(model, "planar_agents", Int[])) : planar_agents,
    )
end

# Closed-loop rollout that records agent positions. `mode` ∈ (:analytic, :oracle, :rl).
function rollout(cfg::Config, actor, seed::Int, mode::Symbol)
    env = reset!(SheafEnv(cfg), MersenneTwister(seed))
    K = cfg.horizon
    pos = Array{Float64}(undef, STALK_DIM, K, N_AGENTS)
    for t in 1:K
        base = base_control(env)
        for i in 1:N_AGENTS
            raw = mode === :analytic ? base[i] :
                  mode === :oracle   ? base[i] .- drift_field(cfg, env.x[i], env.t) :
                                       base[i] .+ cfg.residual_scale .* Float64.(vec(actor(observation(env, i))))
            u = project_ball(cfg, raw)
            pos[:, t, i] = env.x[i]
            push_history!(env, i, u)
            env.x[i] = env.x[i] .+ cfg.dt .* (u .+ drift_field(cfg, env.x[i], env.t))
        end
        env.t += 1
    end
    pos, copy(env.reference), env.targets
end

tracking_error(pos, reference) =
    [mean(norm(pos[:, t, i] .- reference[:, t, i]) for i in 1:N_AGENTS) for t in 1:size(pos, 2)]

# Marker shape per agent: square = planar (x-y-only), diamond = interior, circle = pinned.
marker(cfg, i) = i in cfg.planar_agents ? :rect : pinned_target(i) == 0 ? :diamond : :circle
agent_color(i) = pinned_target(i) == 0 ? :gray : pinned_target(i)

# One top-down (x-y) frame at time `t`.
function topdown_panel(cfg, pos, targets, t, title)
    plt = plot(; title, titlefontsize=10, legend=false, aspect_ratio=:equal,
               xlabel="x", ylabel="y", xlims=(-7, 7), ylims=(-7, 7))
    for (i, j) in AGENT_EDGES
        plot!(plt, [pos[1,t,i], pos[1,t,j]], [pos[2,t,i], pos[2,t,j]], color=:gray, alpha=0.2)
    end
    for k in 1:N_TARGETS
        lo = max(1, t-40)
        plot!(plt, targets[1, lo:t, k], targets[2, lo:t, k], color=:black, alpha=0.3)
        scatter!(plt, [targets[1,t,k]], [targets[2,t,k]], color=:black, markershape=:star5, markersize=8)
    end
    for i in 1:N_AGENTS
        scatter!(plt, [pos[1,t,i]], [pos[2,t,i]], color=agent_color(i), markersize=5, markershape=marker(cfg, i))
    end
    plt
end

# One 3-D frame at time `t` (the heterogeneity — planar agents off their target's z — is only visible here).
function frame3d(cfg, pos, targets, t, title)
    plt = plot(; title, titlefontsize=10, legend=false, xlabel="x", ylabel="y", zlabel="z",
               xlims=(-7,7), ylims=(-7,7), zlims=(-4,4), camera=(35, 18))
    for (i, j) in AGENT_EDGES
        plot!(plt, [pos[1,t,i],pos[1,t,j]], [pos[2,t,i],pos[2,t,j]], [pos[3,t,i],pos[3,t,j]], color=:gray, alpha=0.2)
    end
    for k in 1:N_TARGETS
        scatter!(plt, [targets[1,t,k]], [targets[2,t,k]], [targets[3,t,k]], color=:black, markershape=:star5, markersize=7)
    end
    for i in 1:N_AGENTS
        lo = max(1, t-25)
        plot!(plt, pos[1,lo:t,i], pos[2,lo:t,i], pos[3,lo:t,i], color=agent_color(i), alpha=0.35)
        scatter!(plt, [pos[1,t,i]], [pos[2,t,i]], [pos[3,t,i]], color=agent_color(i), markersize=4, markershape=marker(cfg, i))
    end
    plt
end

"""
Render every figure for a saved `model_path` to `outdir`, tagged `tag`. Set `three_d=true` for the
heterogeneous demo (adds a 3-D animation). `cfg_overrides` are forwarded to `config_from_model`.
"""
function render_all(model_path; tag, outdir, seed=7001, three_d=false, header="multi-agent structure", cfg_overrides...)
    mkpath(outdir)
    model = load(model_path)
    cfg   = config_from_model(model; cfg_overrides...)
    actor = load_actor(model)
    capstr = isinf(cfg.actuator_cap) ? "" : @sprintf("  ‖u‖≤%.2g", cfg.actuator_cap)

    pa, ref, targets = rollout(cfg, actor, seed, :analytic)
    po, _, _         = rollout(cfg, actor, seed, :oracle)
    pr, refr, _      = rollout(cfg, actor, seed, :rl)
    ea, eo, er = tracking_error(pa, ref), tracking_error(po, ref), tracking_error(pr, refr)
    @printf("%s (%s): analytic %.4f  drift-oracle %.4f  sheaf+RL %.4f\n", header, tag, mean(ea), mean(eo), mean(er))
    # machine-readable metrics for the docs number-autoupdate (refresh_viz.sh reads these)
    open(joinpath(outdir, "results_$tag.toml"), "w") do io
        @printf(io, "analytic = %.4f\noracle = %.4f\nrl = %.4f\n", mean(ea), mean(eo), mean(er))
    end

    stride = max(1, cfg.horizon ÷ 150)
    anim = @animate for t in 1:stride:cfg.horizon
        plot(topdown_panel(cfg, pa, targets, t, @sprintf("analytic  %.3f", mean(ea))),
             topdown_panel(cfg, po, targets, t, @sprintf("drift-oracle  %.3f", mean(eo))),
             topdown_panel(cfg, pr, targets, t, @sprintf("sheaf+RL  %.3f", mean(er))),
             layout=(1,3), size=(1500,500),
             plot_title="$header$capstr, t=$t   [★ target  ◇ interior  ▪ planar]")
    end
    gif(anim, joinpath(outdir, "anim_$tag.gif"), fps=20)

    ep = plot(title="tracking error ‖x−q*‖ — $header$capstr", xlabel="step", ylabel="mean ‖x−q*‖")
    plot!(ep, ea, label="analytic", lw=2); plot!(ep, eo, label="drift-oracle", lw=2, ls=:dot); plot!(ep, er, label="sheaf+RL", lw=2, ls=:dash)
    png(ep, joinpath(outdir, "trkerr_$tag.png"))

    p2 = plot(title="top-down trajectories — sheaf+RL ($tag)", xlabel="x", ylabel="y", legend=false, aspect_ratio=:equal)
    for k in 1:N_TARGETS; plot!(p2, targets[1,1:cfg.horizon,k], targets[2,1:cfg.horizon,k], color=:black, ls=:dash, alpha=0.6, lw=2); end
    for i in 1:N_AGENTS; plot!(p2, pr[1,:,i], pr[2,:,i], color=agent_color(i), alpha=0.7); end
    png(p2, joinpath(outdir, "traj2d_$tag.png"))

    p3 = plot(title="3-D trajectories — sheaf+RL ($tag)", xlabel="x", ylabel="y", zlabel="z", legend=false)
    for i in 1:N_AGENTS; plot!(p3, pr[1,:,i], pr[2,:,i], pr[3,:,i], color=agent_color(i), alpha=0.7); end
    for k in 1:N_TARGETS; plot!(p3, targets[1,1:cfg.horizon,k], targets[2,1:cfg.horizon,k], targets[3,1:cfg.horizon,k], color=:black, ls=:dash, alpha=0.6); end
    png(p3, joinpath(outdir, "traj3d_$tag.png"))

    if three_d
        anim3 = @animate for t in 1:stride:cfg.horizon
            plot(frame3d(cfg, pa, targets, t, @sprintf("analytic  %.3f", mean(ea))),
                 frame3d(cfg, pr, targets, t, @sprintf("sheaf+RL  %.3f", mean(er))),
                 layout=(1,2), size=(1300,650),
                 plot_title="$header in 3-D, t=$t   [▪ planar agents track x-y only → off-z]")
        end
        gif(anim3, joinpath(outdir, "anim3d_$tag.gif"), fps=18)
    end
    (; ea, eo, er, outdir)
end
