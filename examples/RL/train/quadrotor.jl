#
# Underactuated extension: each agent is a planar QUADROTOR (state ℝ⁶ = [y,z,φ,vy,vz,ω], control ℝ²
# = two rotor-thrust deviations). The quad is UNDERACTUATED (2 controls, 6 states; can't translate
# without rolling), so g is NOT full row rank → the ECC paper's analytic law u=−k·g⁺·η is INAPPLICABLE
# (its Assumption 2 fails). There is no analytic sheaf controller here — a learned/optimal one is
# REQUIRED. That is the point of this experiment.
#
#   Outer (sheaf, unchanged): harmonic extension q*_i ∈ ℝ² (desired y,z) from the targets.
#   Baseline: LQR about hover tracking q*_i (textbook, trustworthy — NOT the paper's law).
#   RL: residual on LQR,  u = clamp(u_LQR + RESS·Δ_θ([x_i; q*_i−p_i])),  reward −‖p_i−q*_i‖.
#   Δ≈0 at init ⇒ policy IS the (stabilizing) LQR — safe start on an open-loop-unstable plant.
#
#   julia --project=examples/RL examples/RL/train/quadrotor.jl
# Env: RLQ_N, RLQ_K, RLQ_DT, RLQ_RESSCALE, RLQ_TOTAL, RLQ_NENVS, RLQ_DRAG, RLQ_WIND, RLQ_OUT, RLQ_SEED
#
using CellularSheaves
import CellularSheaves.NetworkSheaves.EuclideanSheaves: harmonic_extension
using Flux, Random, Statistics, LinearAlgebra, JLD2, Printf
include(joinpath(@__DIR__, "..", "arch", "nonlinear_dynamics.jl"))   # QuadParams, nl_rk4_step, nl_linearize_c, c2d_affine

const NX = 6; const NU = 2; const ND = 2     # quad state / control / position-sheaf dim
const N      = parse(Int,     get(ENV,"RLQ_N","4"))
const K      = parse(Int,     get(ENV,"RLQ_K","300"))
const DT     = parse(Float64, get(ENV,"RLQ_DT","0.05"))
const RESS   = parse(Float64, get(ENV,"RLQ_RESSCALE","2.0"))
const TOTAL  = parse(Int,     get(ENV,"RLQ_TOTAL","400000"))
const NENVS  = parse(Int,     get(ENV,"RLQ_NENVS","32"))
const DRAG   = parse(Float64, get(ENV,"RLQ_DRAG","0.0"))     # unknown aero drag (cd)
const WIND   = parse(Float64, get(ENV,"RLQ_WIND","0.0"))     # unknown lateral wind accel
const CLASS  = parse(Bool,    get(ENV,"RLQ_CLASS","false"))  # sample a CLASS of plants per episode (E9 regime)
const XEDGES = 1
const OBS    = NX + ND                        # [x_i(6) ; q*_i − p_i(2)]
const SEED   = parse(Int,     get(ENV,"RLQ_SEED","1"))
const OUT    = get(ENV,"RLQ_OUT", joinpath(@__DIR__,"..","cache","rl_quad.jld2"))
Random.seed!(SEED)

const PNOM = QuadParams()                                  # nominal plant (LQR designed on this)
const PTRUE = QuadParams(; cd=DRAG, wy=WIND)               # true plant (unknown drag/wind = the drift)
const UMAX = PNOM.m*PNOM.g/2                               # thrust deviation bound: T_i = mg/2 + u_i ≥ 0

# --- LQR about hover (the trustworthy baseline; replaces the paper's inapplicable analytic law) ----
function dlqr(A, B, Q, R; iters=2000)
    P = copy(Q)
    for _ in 1:iters
        K = (R .+ B'*P*B) \ (B'*P*A)
        Pn = Q .+ A'*P*A .- A'*P*B*K
        norm(Pn .- P) < 1e-10 && (P = Pn; break); P = Pn
    end
    (R .+ B'*P*B) \ (B'*P*A)
end
const Ac, Bc, _c = nl_linearize_c(zeros(NX), zeros(NU), PNOM)
const Ad, Bd, _cd = c2d_affine(Ac, Bc, zeros(NX), DT)
const QLQR = diagm([12.0, 12.0, 1.0, 1.0, 1.0, 0.5])
const RLQR = diagm([0.5, 0.5])
const KLQR = dlqr(Ad, Bd, QLQR, RLQR)
@info "LQR designed at hover" Kspectral=round.(KLQR; digits=2) underactuated="g rank $(rank(Bc)) < $NX"

# --- env helpers: position sheaf (nd=2) + moving targets -------------------------------------
eye(nd) = Matrix{Float64}(I, nd, nd)
function build_sheaf(n, m, nd, edges, pin)
    s = EuclideanSheaf{Float64}(fill(nd, n + m))
    for (i, j) in edges; add_sheaf_edge!(s, i, j, eye(nd), eye(nd)); end
    for i in 1:n, j in 1:m; pin[i, j] != 0 && add_sheaf_edge!(s, i, n + j, eye(nd), eye(nd)); end
    s, collect(n+1:n+m)
end
function topo(n, rng)
    e = Tuple{Int,Int}[(i, i % n + 1) for i in 1:n]
    for _ in 1:XEDGES; a=rand(rng,1:n); b=rand(rng,1:n); a!=b && push!(e,(min(a,b),max(a,b))); end
    pin = zeros(n,n); for i in 1:n; pin[i,i]=1.0; end; e, pin
end
function targ(rng)   # smooth (y,z) target trajectory, gentle so the quad can follow
    c = (rand(rng,ND).-0.5).*4; A = 0.3 .+ rand(rng,ND).*0.8
    om = (0.2 .+ rand(rng,ND).*0.4).*(2π/(K*DT)); ph = rand(rng,ND).*2π
    P = Array{Float64}(undef, ND, K+1); for t in 0:K; tt=t*DT; @. P[:,t+1]=c+A*sin(om*tt+ph); end; P
end

mutable struct QEnv
    d; Tg::Array{Float64,3}; Q::Array{Float64,3}; x::Vector{Vector{Float64}}; t::Int; p::QuadParams
end
newenv() = QEnv(nothing, zeros(0,0,0), zeros(0,0,0), Vector{Float64}[], 1, PNOM)
# sample the TRUE plant for an episode: fixed PTRUE, or a per-episode CLASS (unknown to the policy)
sample_plant(rng) = CLASS ?
    QuadParams(; m=PNOM.m*(0.8+0.4rand(rng)), cd=0.8rand(rng), wy=(rand(rng)-0.5)*6, wz=(rand(rng)-0.5)*2) : PTRUE
function reset_env!(e::QEnv, rng)
    e.p = sample_plant(rng)
    edges, pin = topo(N, rng); s, tv = build_sheaf(N, N, ND, edges, pin); e.d = coboundary_map(s)
    e.Tg = Array{Float64}(undef, ND, K+1, N); for j in 1:N; e.Tg[:,:,j] = targ(rng); end
    e.Q = Array{Float64}(undef, ND, K, N)
    for t in 1:K
        qv = Vector(harmonic_extension(s, Dict(tv[j] => e.Tg[:,t,j] for j in 1:N))[1])
        for i in 1:N; e.Q[:,t,i] = qv[(i-1)*ND+1:i*ND]; end
    end
    e.x = [vcat(e.Q[:,1,i] .+ 0.2.*randn(rng,ND), zeros(NX-ND)) for i in 1:N]   # start at hover near q*
    e.t = 1; e
end
xref(e, i) = vcat(e.Q[:,e.t,i], zeros(NX-ND))             # desired state: pos=q*, hover attitude, zero vel
ulqr(e, i) = -KLQR * (e.x[i] .- xref(e, i))
clampu(u)  = clamp.(u, -UMAX, UMAX)
obs_i(e, i) = Float32.(vcat(e.x[i], e.Q[:,e.t,i] .- e.x[i][1:ND]))
function step!(e::QEnv, Δ::Matrix{Float64})
    r = zeros(Float32, N)
    for i in 1:N
        u = clampu(ulqr(e, i) .+ RESS .* Δ[:,i])
        perr = norm(e.x[i][1:ND] .- e.Q[:,e.t,i])
        r[i] = Float32(-perr - 0.05*e.x[i][3]^2)         # track position; small attitude penalty
        e.x[i] = nl_rk4_step(e.x[i], u, e.p, DT)
    end
    e.t += 1; r, e.t > K
end

# --- nets + TD3 (same machinery as train_rl_sheaf.jl) ----------------------------------------
mk_actor()  = Chain(Dense(OBS,256,relu), Dense(256,256,relu), Dense(256,NU,tanh))
mk_critic() = Chain(Dense(OBS+NU,256,relu), Dense(256,256,relu), Dense(256,1))
mutable struct Replay; s; a; r; sp; d; cap; idx; full; end
Replay(cap) = Replay(zeros(Float32,OBS,cap), zeros(Float32,NU,cap), zeros(Float32,cap),
                     zeros(Float32,OBS,cap), zeros(Float32,cap), cap, 0, false)
function store!(rb,s,a,r,sp,d)
    for i in 1:length(r)
        p=(rb.idx%rb.cap)+1; rb.s[:,p]=s[:,i]; rb.a[:,p]=a[:,i]; rb.r[p]=r[i]; rb.sp[:,p]=sp[:,i]; rb.d[p]=d[i]
        rb.idx+=1; rb.idx>=rb.cap && (rb.full=true)
    end
end
nstored(rb)=rb.full ? rb.cap : rb.idx
sample(rb,B)=(n=nstored(rb); idx=rand(1:n,B); (rb.s[:,idx],rb.a[:,idx],rb.r[idx],rb.sp[:,idx],rb.d[idx]))
function train(; rng, total, n_envs, γ=0.99f0, τ=0.005f0, alr=3f-4, clr=1f-3, batch=256,
                 expl=0.1f0, tgt_noise=0.2f0, tgt_clip=0.5f0, pdelay=2, warmup=2000, res_anchor=0.05f0)
    actor=mk_actor(); critic1=mk_critic(); critic2=mk_critic()
    at=deepcopy(actor); c1t=deepcopy(critic1); c2t=deepcopy(critic2)
    oa=Flux.setup(Adam(alr),actor); oc1=Flux.setup(Adam(clr),critic1); oc2=Flux.setup(Adam(clr),critic2)
    rb=Replay(300_000); polyak!(t,s)=for (pt,ps) in zip(Flux.params(t),Flux.params(s)); pt.=(1-τ).*pt.+τ.*ps; end
    envs=[reset_env!(newenv(),rng) for _ in 1:n_envs]
    cols()=reduce(hcat,[reduce(hcat,[obs_i(e,i) for i in 1:N]) for e in envs])
    obs=cols(); losses=Float32[]; steps=0; t0=time()
    while steps < total
        a=clamp.(actor(obs).+expl.*randn(Float32,NU,size(obs,2)),-1f0,1f0)
        rvec=Float32[]; dvec=Float32[]; ai=1
        for e in envs
            Δ=Float64.(a[:,ai:ai+N-1]); ai+=N; r,done=step!(e,Δ)
            append!(rvec,r); append!(dvec, fill(done ? 1f0 : 0f0, N)); done && reset_env!(e,rng)
        end
        sp=cols(); store!(rb,obs,a,rvec,sp,dvec); obs=sp; steps+=size(a,2)
        if nstored(rb)>=max(batch,warmup)
            s,ac,r,s2,d=sample(rb,batch)
            ε=clamp.(tgt_noise.*randn(Float32,NU,batch),-tgt_clip,tgt_clip); a2=clamp.(at(s2).+ε,-1f0,1f0)
            y=r .+ γ.*(1f0 .- d).*min.(vec(c1t(vcat(s2,a2))),vec(c2t(vcat(s2,a2))))
            _,gc1=Flux.withgradient(m->mean((vec(m(vcat(s,ac))).-y).^2),critic1); Flux.update!(oc1,critic1,gc1[1])
            l2,gc2=Flux.withgradient(m->mean((vec(m(vcat(s,ac))).-y).^2),critic2); Flux.update!(oc2,critic2,gc2[1])
            push!(losses,Float32(l2))
            if (steps÷size(a,2))%pdelay==0
                _,ga=Flux.withgradient(actor) do m; aa=m(s); -mean(critic1(vcat(s,aa)))+res_anchor*mean(aa.^2); end
                Flux.update!(oa,actor,ga[1]); polyak!(at,actor); polyak!(c1t,critic1); polyak!(c2t,critic2)
            end
        end
        if steps % (n_envs*N*50) < n_envs*N
            @printf("  steps %7d  buff %7d  critic_l %.4f  elapsed %.0fs\n",
                    steps,nstored(rb), isempty(losses) ? NaN : losses[end], time()-t0); flush(stdout)
        end
    end
    actor
end

# --- eval: LQR vs RL (and the analytic-law failure) on fresh scenarios ------------------------
function eval_policy(actor; rng, nscn=40)
    function run(mode)
        tot=0.0; cnt=0; div=0
        for _ in 1:nscn
            e=reset_env!(newenv(),rng); bad=false
            for _ in 1:K
                for i in 1:N
                    u = mode===:lqr ? clampu(ulqr(e,i)) : clampu(ulqr(e,i).+RESS.*Float64.(vec(actor(obs_i(e,i)))))
                    tot += norm(e.x[i][1:ND].-e.Q[:,e.t,i]); cnt += 1
                    e.x[i]=nl_rk4_step(e.x[i],u,e.p,DT)
                    any(!isfinite, e.x[i]) && (bad=true)
                end
                e.t += 1; bad && break
            end
            bad && (div += 1)
        end
        (err=tot/cnt, diverged=div)
    end
    (lqr=run(:lqr), rl=run(:rl))
end

@info "Underactuated quad RL config" N K DT RESS DRAG WIND UMAX OBS TOTAL
rng=MersenneTwister(SEED)
actor=train(; rng=rng, total=TOTAL, n_envs=NENVS)
res=eval_policy(actor; rng=MersenneTwister(SEED+777))
@printf("\n=== EVAL (mean ‖p−q*‖, 40 fresh scenarios; underactuated quad) ===\n")
@printf("  LQR baseline = %.4f   (%d diverged)\n", res.lqr.err, res.lqr.diverged)
@printf("  RL (resid)   = %.4f   (%d diverged)\n", res.rl.err, res.rl.diverged)
@printf("  RL vs LQR    = %+.1f%%\n", 100*(res.rl.err-res.lqr.err)/res.lqr.err)
cm=cpu(actor); Ws=[Array(l.weight) for l in cm.layers]; bs=[Array(l.bias) for l in cm.layers]
jldsave(OUT; Ws=Ws, bs=bs, N=N, K=K, DT=DT, RESS=RESS, DRAG=DRAG, WIND=WIND, UMAX=UMAX, OBS=OBS, NX=NX, NU=NU,
        KLQR=KLQR, eval=Dict("lqr"=>res.lqr.err,"rl"=>res.rl.err))
@info "Saved underactuated-quad RL policy" OUT
