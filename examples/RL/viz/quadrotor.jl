#
# Render the underactuated planar-quad result: sheaf coordination (q* per agent) + learned residual
# control vs the LQR baseline, on the NONLINEAR underactuated quad (the regime where the paper's
# analytic law is inapplicable). Side-by-side GIF + tracking-error plot. Flux-free local render.
#
#   julia --project=examples/RL examples/RL/viz/quadrotor.jl
# Env: AQ_MODEL, AQ_SEEDS, AQ_DRAG, AQ_WIND (viz disturbance), AQ_OUTDIR
#
get(ENV, "DISPLAY", "") == "" && (ENV["GKSwstype"] = "100")
using CellularSheaves
import CellularSheaves.NetworkSheaves.EuclideanSheaves: harmonic_extension
using LinearAlgebra, Random, JLD2, Printf, Statistics, Plots
include(joinpath(@__DIR__, "..", "arch", "nonlinear_dynamics.jl"))
gr()

const MODEL  = get(ENV,"AQ_MODEL", joinpath(@__DIR__,"..","cache","rl_quad_class.jld2"))
const SEEDS  = parse.(Int, split(get(ENV,"AQ_SEEDS","7001,7002,7003"), ","))
const OUTDIR = get(ENV,"AQ_OUTDIR", joinpath(@__DIR__,"..","cache","viz_quad"))
mkpath(OUTDIR)
M = load(MODEL)
const Ws=M["Ws"]; const bs=M["bs"]; const KLQR=M["KLQR"]
const N=M["N"]; const K=M["K"]; const DT=M["DT"]; const RESS=M["RESS"]; const UMAX=M["UMAX"]
const NX=M["NX"]; const NU=M["NU"]; const ND=2; const XEDGES=1
# viz plant: a representative unknown disturbance (the policy was trained over a CLASS of these)
const PVIZ = QuadParams(; cd=parse(Float64,get(ENV,"AQ_DRAG","0.4")), wy=parse(Float64,get(ENV,"AQ_WIND","3.0")))
@info "Quad viz" N K DT RESS eval=get(M,"eval",nothing) plant=(PVIZ.cd, PVIZ.wy)

actor(o) = (h=Float32.(o); for k in 1:length(Ws)-1; h=max.(0f0, Ws[k]*h .+ bs[k]); end; tanh.(Ws[end]*h .+ bs[end]))
eye(nd)=Matrix{Float64}(I,nd,nd)
function build_sheaf(n,m,nd,edges,pin)
    s=EuclideanSheaf{Float64}(fill(nd,n+m))
    for (i,j) in edges; add_sheaf_edge!(s,i,j,eye(nd),eye(nd)); end
    for i in 1:n,j in 1:m; pin[i,j]!=0 && add_sheaf_edge!(s,i,n+j,eye(nd),eye(nd)); end
    s, collect(n+1:n+m)
end
function topo(n,rng)
    e=Tuple{Int,Int}[(i,i%n+1) for i in 1:n]
    for _ in 1:XEDGES; a=rand(rng,1:n);b=rand(rng,1:n);a!=b&&push!(e,(min(a,b),max(a,b))); end
    pin=zeros(n,n); for i in 1:n; pin[i,i]=1.0; end; e,pin
end
function targ(rng)
    c=(rand(rng,ND).-0.5).*4; A=0.3 .+ rand(rng,ND).*0.8; om=(0.2 .+ rand(rng,ND).*0.4).*(2π/(K*DT)); ph=rand(rng,ND).*2π
    P=Array{Float64}(undef,ND,K+1); for t in 0:K; tt=t*DT; @. P[:,t+1]=c+A*sin(om*tt+ph); end; P
end
xref(Q,t,i)=vcat(Q[:,t,i],zeros(NX-ND))
ulqr(x,Q,t,i)=-KLQR*(x.-xref(Q,t,i))
clampu(u)=clamp.(u,-UMAX,UMAX)
obs_i(x,Q,t,i)=Float32.(vcat(x, Q[:,t,i].-x[1:ND]))

function rollout(N,seed,mode)
    rng=MersenneTwister(seed); edges,pin=topo(N,rng); s,tv=build_sheaf(N,N,ND,edges,pin); d=coboundary_map(s)
    Tg=Array{Float64}(undef,ND,K+1,N); for j in 1:N; Tg[:,:,j]=targ(rng); end
    Q=Array{Float64}(undef,ND,K,N)
    for t in 1:K; qv=Vector(harmonic_extension(s,Dict(tv[j]=>Tg[:,t,j] for j in 1:N))[1]); for i in 1:N; Q[:,t,i]=qv[(i-1)*ND+1:i*ND]; end; end
    x=[vcat(Q[:,1,i].+0.2.*randn(rng,ND), zeros(NX-ND)) for i in 1:N]
    pos=Array{Float64}(undef,ND,K,N); phi=Array{Float64}(undef,K,N)
    for t in 1:K
        for i in 1:N
            u = mode===:lqr ? clampu(ulqr(x[i],Q,t,i)) : clampu(ulqr(x[i],Q,t,i).+RESS.*Float64.(vec(actor(obs_i(x[i],Q,t,i)))))
            pos[:,t,i]=x[i][1:ND]; phi[t,i]=x[i][3]; x[i]=nl_rk4_step(x[i],u,PVIZ,DT)
        end
    end
    pos,Q,Tg,phi,edges
end
trkerr(pos,Q)=[mean(norm(pos[:,t,i].-Q[:,t,i]) for i in 1:N) for t in 1:K]

function panel(pos,Q,Tg,edges,t,title)
    plt=plot(title=title,legend=false,aspect_ratio=:equal,xlabel="y",ylabel="z")
    for (i,j) in edges; plot!(plt,[pos[1,t,i],pos[1,t,j]],[pos[2,t,i],pos[2,t,j]],color=:gray,alpha=0.25); end
    for i in 1:N
        lo=max(1,t-30); plot!(plt,pos[1,lo:t,i],pos[2,lo:t,i],color=i,alpha=0.4)
        scatter!(plt,[Q[1,t,i]],[Q[2,t,i]],color=i,markershape=:cross,markersize=5,alpha=0.6)  # q* reference
        scatter!(plt,[pos[1,t,i]],[pos[2,t,i]],color=i,markersize=6)                            # quad
    end
    plt
end

summary=String[]
for seed in SEEDS
    @info "Rendering seed" seed
    pl,Ql,Tg,_,edges = rollout(N,seed,:lqr)
    pr,Qr,_,_,_      = rollout(N,seed,:rl)
    el=trkerr(pl,Ql); er=trkerr(pr,Qr)
    open(joinpath(OUTDIR,"results_$seed.toml"),"w") do io   # metrics for the docs number-autoupdate (analytic slot = LQR baseline)
        @printf(io,"analytic = %.4f\nrl = %.4f\n", mean(el), mean(er))
    end
    stride=max(1,K÷150)
    anim=@animate for t in 1:stride:K
        plot(panel(pl,Ql,Tg,edges,t,@sprintf("LQR  %.3f",mean(el))),
             panel(pr,Qr,Tg,edges,t,@sprintf("sheaf+RL  %.3f",mean(er))),
             layout=(1,2),size=(1000,500),plot_title="underactuated quad, t=$t (seed $seed) [+ = q*]")
    end
    gif(anim, joinpath(OUTDIR,"anim_quad_$seed.gif"), fps=20)
    ep=plot(title="tracking error ‖p−q*‖ (seed $seed)",xlabel="step",ylabel="mean ‖p−q*‖")
    plot!(ep,el,label="LQR",lw=2); plot!(ep,er,label="sheaf+RL",lw=2,ls=:dash); png(ep,joinpath(OUTDIR,"trkerr_$seed.png"))
    push!(summary,@sprintf("seed %d:  LQR %.4f   sheaf+RL %.4f",seed,mean(el),mean(er)))
end
println("\n=== underactuated quad: LQR vs sheaf+RL (mean ‖p−q*‖) ===")
for l in summary; println("  ",l); end
println("GIFs + PNGs in: $OUTDIR")
