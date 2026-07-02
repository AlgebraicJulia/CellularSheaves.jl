#
# Stage-2 RL: residual TD3 on the sheaf-tracking env with a TIME-VARYING nonlinear drift f(q,t).
# Because f changes with t, the optimal feedforward is NOT a function of the current state alone — the
# policy must INFER the drift from recent experience. We give it a HISTORY WINDOW: the last H (x,u)
# pairs, from which the drift is exactly recoverable (f(x^s,s) = (x^{s+1}−x^s)/DT − u^s) and can be
# extrapolated to now. This is genuine online adaptation to unknown, non-autonomous dynamics — the
# regime the ECC paper (f_i(q,t)) is about; removes the Stage-1 "observable drift" cheat.
#
#   u_i = −k(Lz)_i + RESS·Δ_θ([x_i; e_i; last H (x,u)]),   reward −‖x_i − q*_i‖
#   oracle (ceiling) = −k(Lz)_i − f(x_i, t)   (knows the true time-varying drift)
#
#   julia --project=examples/RL examples/RL/train/single_integrator.jl
# Env: RLS_KGAIN, RLS_DRIFT, RLS_DRIFTW (time-var freq), RLS_HIST (H), RLS_N, RLS_K, RLS_DT,
#      RLS_RESSCALE, RLS_TOTAL, RLS_NENVS, RLS_OUT, RLS_SEED
#
using CellularSheaves
import CellularSheaves.NetworkSheaves.EuclideanSheaves: harmonic_extension
using Flux, CUDA, cuDNN, Random, Statistics, LinearAlgebra, JLD2, Printf
const dev = CUDA.functional() ? gpu : cpu
CUDA.functional() && @info "GPU: $(CUDA.name(CUDA.device()))"

const ND     = 2
const KGAIN  = parse(Float64, get(ENV,"RLS_KGAIN","2.0"))
const DRIFT  = parse(Float64, get(ENV,"RLS_DRIFT","2.0"))
const DRIFTW = parse(Float64, get(ENV,"RLS_DRIFTW","3.0"))     # angular freq of the TIME variation
const HIST   = parse(Int,     get(ENV,"RLS_HIST","6"))         # history-window length (past (x,u) pairs)
const HORIZON= parse(Int,     get(ENV,"RLS_HORIZON","0"))      # future-reference window (anticipation; q*(t+1..t+H))
const N      = parse(Int,     get(ENV,"RLS_N","4"))
const K      = parse(Int,     get(ENV,"RLS_K","400"))
const DT     = parse(Float64, get(ENV,"RLS_DT","0.02"))
const RESS   = parse(Float64, get(ENV,"RLS_RESSCALE","4.0"))
const UBAR   = parse(Float64, get(ENV,"RLS_UBAR","Inf"))       # actuator ball ‖u‖≤ū (Inf = unconstrained)
const TOTAL  = parse(Int,     get(ENV,"RLS_TOTAL","300000"))
const NENVS  = parse(Int,     get(ENV,"RLS_NENVS","32"))
const XEDGES = 1
const OBS    = 2ND + 2ND*HIST + ND*HORIZON   # [x_i ; e_i ; histx ; histu ; future-ref(H·ND)]
const SEED   = parse(Int,     get(ENV,"RLS_SEED","1"))
const OUT    = get(ENV,"RLS_OUT", joinpath(@__DIR__,"..","cache","rl_sheaf_tv.jld2"))
Random.seed!(SEED)

# TIME-VARYING nonlinear drift: a swirl whose phase rotates with time → f(q,t) ≠ f(q).
drift(q, t) = DRIFT == 0.0 ? zero(q) :
    DRIFT .* vcat(-sin(q[2] + DRIFTW*t*DT), sin(q[1] + DRIFTW*t*DT), zeros(length(q)-2))

# actuator ball: hard-clamp the physical control ‖u‖ ≤ ū (applied to every controller, fairly)
function clampu(u)
    isfinite(UBAR) || return u
    n = norm(u); n > UBAR ? u .* (UBAR / n) : u
end

# --- env helpers (match collect_bc_data.jl) --------------------------------------------------
eye(nd) = Matrix{Float64}(I, nd, nd)
function build_sheaf(n, m, nd, edges, pin)
    s = EuclideanSheaf{Float64}(fill(nd, n + m))
    for (i, j) in edges; add_sheaf_edge!(s, i, j, eye(nd), eye(nd)); end
    for i in 1:n, j in 1:m; pin[i, j] != 0 && add_sheaf_edge!(s, i, n + j, eye(nd), eye(nd)); end
    s, collect(n+1:n+m)
end
function topo(n, rng)
    e = Tuple{Int,Int}[(i, i % n + 1) for i in 1:n]
    for _ in 1:XEDGES; a = rand(rng,1:n); b = rand(rng,1:n); a != b && push!(e,(min(a,b),max(a,b))); end
    pin = zeros(n, n); for i in 1:n; pin[i,i] = 1.0; end
    e, pin
end
function targ(rng, nd)
    c = (rand(rng,nd).-0.5).*6; A = 0.5 .+ rand(rng,nd).*2
    om = (0.3 .+ rand(rng,nd).*0.7).*(2π/(K*DT)); ph = rand(rng,nd).*2π; dr = (rand(rng,nd).-0.5).*0.3
    P = Array{Float64}(undef, nd, K+1); for t in 0:K; tt=t*DT; @. P[:,t+1]=c+A*sin(om*tt+ph)+dr*tt; end; P
end

# --- per-scenario env (N agents) with per-agent (x,u) history --------------------------------
mutable struct SheafEnv
    d
    Tg :: Array{Float64,3}
    Q  :: Array{Float64,3}
    x  :: Vector{Vector{Float64}}
    hx :: Array{Float64,3}        # ND × HIST × N   recent states (oldest..newest)
    hu :: Array{Float64,3}        # ND × HIST × N   recent controls
    t  :: Int
end
function reset_env!(e::SheafEnv, rng)
    edges, pin = topo(N, rng)
    s, tv = build_sheaf(N, N, ND, edges, pin)
    e.d = coboundary_map(s)
    e.Tg = Array{Float64}(undef, ND, K+1, N); for j in 1:N; e.Tg[:,:,j] = targ(rng, ND); end
    e.Q = Array{Float64}(undef, ND, K, N)
    for t in 1:K
        qv = Vector(harmonic_extension(s, Dict(tv[j] => e.Tg[:,t,j] for j in 1:N))[1])
        for i in 1:N; e.Q[:,t,i] = qv[(i-1)*ND+1:i*ND]; end
    end
    e.x = [e.Tg[:,1,i] .+ 0.5.*randn(rng, ND) for i in 1:N]
    e.hx = zeros(ND, HIST, N); e.hu = zeros(ND, HIST, N)
    e.t = 1
    e
end
function anchor(e::SheafEnv)
    z = zeros(ND*2N)
    for i in 1:N; z[(i-1)*ND+1:i*ND] = e.x[i]; end
    for j in 1:N; z[(N+j-1)*ND+1:(N+j)*ND] = e.Tg[:,e.t,j]; end
    Lz = e.d' * (e.d * z)
    [-KGAIN .* Lz[(i-1)*ND+1:i*ND] for i in 1:N]
end
function obs_i(e, i)
    base = vcat(e.x[i], e.Q[:,e.t,i] .- e.x[i], vec(e.hx[:,:,i]), vec(e.hu[:,:,i]))
    fut = HORIZON > 0 ? reduce(vcat, [e.Q[:, min(e.t+τ, K), i] .- e.x[i] for τ in 1:HORIZON]) : Float64[]
    Float32.(vcat(base, fut))
end
function step!(e::SheafEnv, Δ::Matrix{Float64})
    ua = anchor(e); r = zeros(Float32, N)
    for i in 1:N
        u = clampu(ua[i] .+ RESS .* Δ[:,i])
        r[i] = Float32(-norm(e.x[i] .- e.Q[:,e.t,i]))
        if HIST > 0   # record (x,u) into history BEFORE advancing (drift recoverable from consecutive states)
            e.hx[:,1:HIST-1,i] = e.hx[:,2:HIST,i]; e.hx[:,HIST,i] = e.x[i]
            e.hu[:,1:HIST-1,i] = e.hu[:,2:HIST,i]; e.hu[:,HIST,i] = u
        end
        e.x[i] = e.x[i] .+ DT .* (u .+ drift(e.x[i], e.t))
    end
    e.t += 1
    r, e.t > K
end
newenv() = SheafEnv(nothing, zeros(0,0,0), zeros(0,0,0), Vector{Float64}[], zeros(0,0,0), zeros(0,0,0), 1)

# --- nets ------------------------------------------------------------------------------------
mk_actor()  = Chain(Dense(OBS,256,relu), Dense(256,256,relu), Dense(256,ND,tanh))
mk_critic() = Chain(Dense(OBS+ND,256,relu), Dense(256,256,relu), Dense(256,1))

# --- TD3 (mirrors qp_td3_core.jl) ------------------------------------------------------------
mutable struct Replay; s; a; r; sp; d; cap; idx; full; end
Replay(cap) = Replay(zeros(Float32,OBS,cap), zeros(Float32,ND,cap), zeros(Float32,cap),
                     zeros(Float32,OBS,cap), zeros(Float32,cap), cap, 0, false)
function store!(rb, s,a,r,sp,d)
    for i in 1:length(r)
        p=(rb.idx % rb.cap)+1; rb.s[:,p]=s[:,i]; rb.a[:,p]=a[:,i]; rb.r[p]=r[i]; rb.sp[:,p]=sp[:,i]; rb.d[p]=d[i]
        rb.idx+=1; rb.idx>=rb.cap && (rb.full=true)
    end
end
nstored(rb) = rb.full ? rb.cap : rb.idx
function sample(rb,B); n=nstored(rb); idx=rand(1:n,B); (rb.s[:,idx],rb.a[:,idx],rb.r[idx],rb.sp[:,idx],rb.d[idx]); end

function train(; rng, total, n_envs, γ=0.99f0, τ=0.005f0, alr=3f-4, clr=1f-3, batch=256,
                 expl=0.15f0, tgt_noise=0.2f0, tgt_clip=0.5f0, pdelay=2, warmup=2000, res_anchor=0.05f0)
    actor=dev(mk_actor()); critic1=dev(mk_critic()); critic2=dev(mk_critic())
    at=deepcopy(actor); c1t=deepcopy(critic1); c2t=deepcopy(critic2)
    oa=Flux.setup(Adam(alr),actor); oc1=Flux.setup(Adam(clr),critic1); oc2=Flux.setup(Adam(clr),critic2)
    rb=Replay(300_000)
    polyak!(t,s)=for (pt,ps) in zip(Flux.params(t),Flux.params(s)); pt.=(1-τ).*pt.+τ.*ps; end
    envs=[reset_env!(newenv(), rng) for _ in 1:n_envs]
    cols() = reduce(hcat,[reduce(hcat,[obs_i(e,i) for i in 1:N]) for e in envs])
    obs=cols(); losses=Float32[]; steps=0; t0=time()
    while steps < total
        a = cpu(clamp.(actor(dev(obs)) .+ dev(expl.*randn(Float32,ND,size(obs,2))), -1f0,1f0))
        rvec=Float32[]; dvec=Float32[]; ai=1
        for e in envs
            Δ = Float64.(a[:, ai:ai+N-1]); ai+=N
            r,done = step!(e, Δ); append!(rvec,r); append!(dvec, fill(done ? 1f0 : 0f0, N))
            done && reset_env!(e, rng)
        end
        sp=cols(); store!(rb, obs, a, rvec, sp, dvec); obs=sp; steps += size(a,2)
        if nstored(rb) >= max(batch,warmup)
            s,ac,r,s2,d=sample(rb,batch)
            sd,acd,rd,s2d,dd = dev(s),dev(ac),dev(r),dev(s2),dev(d)
            ε=dev(clamp.(tgt_noise.*randn(Float32,ND,batch), -tgt_clip,tgt_clip))
            a2=clamp.(at(s2d).+ε, -1f0,1f0)
            y=rd .+ γ.*(1f0 .- dd).*min.(vec(c1t(vcat(s2d,a2))), vec(c2t(vcat(s2d,a2))))
            l1,gc1=Flux.withgradient(m->mean((vec(m(vcat(sd,acd))).-y).^2), critic1); Flux.update!(oc1,critic1,gc1[1])
            l2,gc2=Flux.withgradient(m->mean((vec(m(vcat(sd,acd))).-y).^2), critic2); Flux.update!(oc2,critic2,gc2[1])
            push!(losses, Float32(l1))
            if (steps ÷ size(a,2)) % pdelay == 0
                _,ga=Flux.withgradient(actor) do m
                    aa=m(sd); -mean(critic1(vcat(sd,aa))) + res_anchor*mean(aa.^2)
                end
                Flux.update!(oa,actor,ga[1]); polyak!(at,actor); polyak!(c1t,critic1); polyak!(c2t,critic2)
            end
        end
        if steps % (n_envs*N*50) < n_envs*N
            @printf("  steps %7d  buff %7d  critic_l %.4f  elapsed %.0fs\n",
                    steps, nstored(rb), isempty(losses) ? NaN : losses[end], time()-t0); flush(stdout)
        end
    end
    actor
end

# --- eval: analytic vs oracle(time-varying) vs RL(history) on fresh scenarios -----------------
function eval_policy(actor; rng, nscn=40)
    function run(mode)
        tot=0.0; cnt=0
        for _ in 1:nscn
            e=reset_env!(newenv(), rng)
            for _ in 1:K
                ua=anchor(e)
                for i in 1:N
                    u = clampu(mode===:analytic ? ua[i] :
                        mode===:oracle   ? ua[i] .- drift(e.x[i], e.t) :
                                           ua[i] .+ RESS .* Float64.(vec(actor(obs_i(e,i)))))
                    tot += norm(e.x[i] .- e.Q[:,e.t,i]); cnt += 1
                    if HIST > 0
                        e.hx[:,1:HIST-1,i]=e.hx[:,2:HIST,i]; e.hx[:,HIST,i]=e.x[i]
                        e.hu[:,1:HIST-1,i]=e.hu[:,2:HIST,i]; e.hu[:,HIST,i]=u
                    end
                    e.x[i] = e.x[i] .+ DT .* (u .+ drift(e.x[i], e.t))
                end
                e.t += 1
            end
        end
        tot/cnt
    end
    (analytic=run(:analytic), oracle=run(:oracle), rl=run(:rl))
end

@info "Stage-2 RL config (time-varying drift + history)" KGAIN DRIFT DRIFTW HIST N K RESS OBS TOTAL
rng=MersenneTwister(SEED)
actor=train(; rng=rng, total=TOTAL, n_envs=NENVS)
res=eval_policy(cpu(actor); rng=MersenneTwister(SEED+777))
@printf("\n=== EVAL (mean ‖x−q*‖, 40 fresh scenarios, TIME-VARYING drift) ===\n")
@printf("  analytic = %.4f\n  oracle   = %.4f  (knows f(q,t))\n  RL       = %.4f\n", res.analytic, res.oracle, res.rl)
recovered = 100*(res.analytic-res.rl)/(res.analytic-res.oracle)
@printf("  RL recovers %.1f%% of the analytic→oracle headroom\n", recovered)
cm=cpu(actor); Ws=[Array(l.weight) for l in cm.layers]; bs=[Array(l.bias) for l in cm.layers]
jldsave(OUT; Ws=Ws, bs=bs, KGAIN=KGAIN, DRIFT=DRIFT, DRIFTW=DRIFTW, HIST=HIST, ND=ND, DT=DT, RESS=RESS,
        UBAR=UBAR, HORIZON=HORIZON, K=K, XEDGES=XEDGES, OBS=OBS,
        eval=Dict("analytic"=>res.analytic,"oracle"=>res.oracle,"rl"=>res.rl,"recovered_pct"=>recovered))
@info "Saved Stage-2 RL policy" OUT
