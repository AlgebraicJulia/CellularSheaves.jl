#
# Visualize the learned BC policy against the paper's analytic control law -k(Lz) on fresh
# scenarios. For each seed, both controllers run from the IDENTICAL initial conditions, topology,
# and target trajectories (drawn before the control loop), so the only difference is the controller.
#
# Produces, per seed:
#   anim_bc_vs_analytic_<seed>.gif  — side-by-side: analytic (left) vs BC (right), agents tracking
#                                     their targets, trajectory trails, consensus edges.
#   track_err_<seed>.png            — tracking error vs time, BC vs analytic (the "how close" plot).
# And a summary printed/saved: mean tracking error BC vs analytic across seeds.
#
#   julia --project=examples/RL examples/RL/viz/bc_vs_analytic.jl
# Env: BCS_MODEL (bc_sheaf.jld2), BCS_VIZ_SEEDS (csv), BCS_VIZ_N, BCS_VIZ_K, BCS_OUTDIR
#
using CellularSheaves
import CellularSheaves.NetworkSheaves.EuclideanSheaves: harmonic_extension
using LinearAlgebra, Random, JLD2, Printf, Statistics, Plots
# Flux-free: the policy is reconstructed from raw weight arrays (Ws,bs) saved by the trainer, so
# this renders locally with no Flux/CUDA and no cross-julia-version coupling. Runs headless too.
get(ENV, "DISPLAY", "") == "" && (ENV["GKSwstype"] = "100")
gr()

const MODEL  = get(ENV,"BCS_MODEL", joinpath(@__DIR__,"..","cache","bc_sheaf.jld2"))
const SEEDS  = parse.(Int, split(get(ENV,"BCS_VIZ_SEEDS","9001,9002,9003"), ","))
const VN     = parse(Int, get(ENV,"BCS_VIZ_N","4"))
const VK     = parse(Int, get(ENV,"BCS_VIZ_K","400"))
const OUTDIR = get(ENV,"BCS_OUTDIR", joinpath(@__DIR__,"..","cache","viz"))
mkpath(OUTDIR)

@info "Loading BC model" MODEL
M = load(MODEL)
const ND = M["ND"]; const OBS = M["OBS"]; const HORIZON = M["HORIZON"]
const XM = Float32.(M["XM"]); const XS = Float32.(M["XS"])
const YM = Float32.(M["YM"]); const YS = Float32.(M["YS"])
meta = M["meta"]; const DT = meta["DT"]; const KGAIN = meta["KGAIN"]
const CWEIGHT = get(meta, "CWEIGHT", 1.0); const DRIFT = get(meta, "DRIFT", 0.0); const XEDGES = get(meta, "XEDGES", 1)
@info "Model meta" ND OBS HORIZON DT KGAIN CWEIGHT DRIFT
drift(q) = DRIFT == 0.0 ? zero(q) : DRIFT .* vcat(-sin(q[2]), sin(q[1]), zeros(length(q) - 2))

const Ws = M["Ws"]; const bs = M["bs"]      # raw weights: Dense(relu)×(n-1) → Dense
function mlp(x)
    h = x
    for k in 1:length(Ws)-1
        h = max.(0f0, Ws[k] * h .+ bs[k])    # relu hidden layers
    end
    return Ws[end] * h .+ bs[end]            # linear output
end
policy(obs) = (on = (Float32.(obs) .- XM) ./ XS; Float64.(mlp(on) .* YS .+ YM))

# --- helpers (match collect_bc_data.jl) ------------------------------------------------------
eye(nd) = Matrix{Float64}(I, nd, nd)
function build_tracking_sheaf(N, M_, nd, agent_edges, pin)
    s = EuclideanSheaf{Float64}(fill(nd, N + M_))
    cw = sqrt(CWEIGHT) * eye(nd)
    for (i, j) in agent_edges; add_sheaf_edge!(s, i, j, cw, cw); end
    for i in 1:N, j in 1:M_
        w = pin[i, j]; w != 0.0 && add_sheaf_edge!(s, i, N + j, sqrt(abs(w))*eye(nd), sqrt(abs(w))*eye(nd))
    end
    return s, collect(1:N), collect(N+1:N+M_)
end
function make_topology(N, rng)
    edges = Tuple{Int,Int}[(i, i % N + 1) for i in 1:N]
    for _ in 1:XEDGES
        i = rand(rng, 1:N); j = rand(rng, 1:N); i != j && push!(edges, (min(i,j), max(i,j)))
    end
    pin = zeros(N, N); for i in 1:N; pin[i, i] = 1.0; end
    edges, pin
end
function target_traj(rng, nd, K)
    c = (rand(rng, nd).-0.5).*6.0; A = 0.5 .+ rand(rng, nd).*2.0
    omega = (0.3 .+ rand(rng, nd).*0.7).*(2π/(K*DT)); phi = rand(rng, nd).*2π
    drift = (rand(rng, nd).-0.5).*0.3
    P = Array{Float64}(undef, nd, K+1)
    for t in 0:K; tt = t*DT; @. P[:, t+1] = c + A*sin(omega*tt+phi) + drift*tt; end
    P
end
futtgt(Tg, t, i, K) = HORIZON > 0 ? reduce(vcat, [Tg[:, min(t+τ, K+1), i] for τ in 1:HORIZON]) : Float64[]

# --- one rollout under a chosen controller (identical setup for both modes via the seed) ------
function rollout(N, seed, K; mode::Symbol)
    rng = MersenneTwister(seed)
    edges, pin = make_topology(N, rng)
    s, _av, tverts = build_tracking_sheaf(N, N, ND, edges, pin); d = coboundary_map(s)
    Tg = Array{Float64}(undef, ND, K+1, N); for j in 1:N; Tg[:, :, j] = target_traj(rng, ND, K); end
    x = [Tg[:, 1, i] .+ 0.5 .* randn(rng, ND) for i in 1:N]
    traj = Array{Float64}(undef, ND, K, N)
    for t in 1:K
        z = zeros(ND*(N+N))
        for i in 1:N; z[(i-1)*ND+1:i*ND] = x[i]; end
        for j in 1:N; z[(N+j-1)*ND+1:(N+j)*ND] = Tg[:, t, j]; end
        Lz = d' * (d * z)
        qvec = Vector(harmonic_extension(s, Dict(tverts[j] => Tg[:, t, j] for j in 1:N))[1])
        for i in 1:N
            e = qvec[(i-1)*ND+1:i*ND] .- x[i]
            u = mode === :analytic ? -KGAIN .* Lz[(i-1)*ND+1:i*ND] :
                                     policy(vcat(x[i], e, futtgt(Tg, t+1, i, K)))
            traj[:, t, i] = x[i]; x[i] = x[i] .+ DT .* (u .+ drift(x[i]))
        end
    end
    return traj, Tg, edges
end

trackerr(traj, Tg) = [mean(norm(traj[:, t, i] .- Tg[:, t, i]) for i in 1:size(traj,3)) for t in 1:size(traj,2)]

# --- render -----------------------------------------------------------------------------------
function panel(traj, Tg, edges, t, N, title)
    plt = plot(title=title, legend=false, aspect_ratio=:equal, xlabel="x", ylabel="y")
    for (i, j) in edges                                     # consensus edges (light)
        plot!(plt, [traj[1,t,i], traj[1,t,j]], [traj[2,t,i], traj[2,t,j]], color=:gray, alpha=0.3)
    end
    for i in 1:N
        lo = max(1, t-40)
        plot!(plt, traj[1,lo:t,i], traj[2,lo:t,i], color=i, alpha=0.5)               # trail
        scatter!(plt, [traj[1,t,i]], [traj[2,t,i]], color=i, markersize=6)           # agent
        scatter!(plt, [Tg[1,t,i]], [Tg[2,t,i]], color=i, markershape=:star5, markersize=8) # target
    end
    plt
end

summary = String[]
ta_all = Float64[]; tb_all = Float64[]
for seed in SEEDS
    @info "Rendering seed" seed
    ta, Tg, edges = rollout(VN, seed, VK; mode=:analytic)
    tb, _, _      = rollout(VN, seed, VK; mode=:bc)
    ea = trackerr(ta, Tg); eb = trackerr(tb, Tg)
    push!(ta_all, mean(ea)); push!(tb_all, mean(eb))

    stride = max(1, VK ÷ 150)
    anim = @animate for t in 1:stride:VK
        plot(panel(ta, Tg, edges, t, VN, "analytic  −k(Lz)"),
             panel(tb, Tg, edges, t, VN, "BC policy"),
             layout=(1,2), size=(1000,500), plot_title="t=$t  (seed $seed)")
    end
    gif(anim, joinpath(OUTDIR, "anim_bc_vs_analytic_$seed.gif"), fps=20)

    errplt = plot(title="tracking error (seed $seed)", xlabel="step", ylabel="mean |x-target|")
    plot!(errplt, ea, label="analytic", lw=2); plot!(errplt, eb, label="BC", lw=2, ls=:dash)
    png(errplt, joinpath(OUTDIR, "track_err_$seed.png"))
    push!(summary, @sprintf("seed %d:  analytic %.4f   BC %.4f   (Δ %+.1f%%)",
          seed, mean(ea), mean(eb), 100*(mean(eb)-mean(ea))/mean(ea)))
end

println("\n=== BC vs analytic — mean tracking error ===")
for line in summary; println("  ", line); end
@printf("OVERALL: analytic %.4f   BC %.4f   (Δ %+.1f%%)\n",
        mean(ta_all), mean(tb_all), 100*(mean(tb_all)-mean(ta_all))/mean(ta_all))
println("GIFs + PNGs in: $OUTDIR")
