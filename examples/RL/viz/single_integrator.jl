#
# Render the time-varying-drift result: learned u(x,b*) with a HISTORY window, tracking under a
# genuinely unknown (non-autonomous) drift f(q,t), vs the analytic law and the drift-oracle.
# Flux-free local render (raw weight matrices), CUDA-optional. Handles history obs + time-varying
# drift (reads everything from meta).
#
#   julia --project=examples/RL examples/RL/viz/single_integrator.jl
# Env: ATV_MODEL, ATV_SEEDS, ATV_OUTDIR
#
get(ENV, "DISPLAY", "") == "" && (ENV["GKSwstype"] = "100")
using CellularSheaves
import CellularSheaves.NetworkSheaves.EuclideanSheaves: harmonic_extension
using LinearAlgebra, Random, JLD2, Printf, Statistics, Plots, CUDA
gr()
const dev = CUDA.functional() ? CuArray : identity
CUDA.functional() && @info "GPU: $(CUDA.name(CUDA.device()))"

const MODEL  = get(ENV,"ATV_MODEL", joinpath(@__DIR__,"..","cache","rl_sheaf_tv.jld2"))
const SEEDS  = parse.(Int, split(get(ENV,"ATV_SEEDS","9101,9102,9103"), ","))
const OUTDIR = get(ENV,"ATV_OUTDIR", joinpath(@__DIR__,"..","cache","viz_tv"))
mkpath(OUTDIR)
M = load(MODEL)
const Ws=dev.(M["Ws"]); const bs=dev.(M["bs"])
const ND=M["ND"]; const DT=M["DT"]; const KGAIN=M["KGAIN"]; const DRIFT=M["DRIFT"]
const DRIFTW=M["DRIFTW"]; const HIST=M["HIST"]; const RESS=M["RESS"]; const XEDGES=get(M,"XEDGES",1)
const K = M["K"]
@info "TV viz" ND DT KGAIN DRIFT DRIFTW HIST RESS eval=get(M,"eval",nothing)

actor(o)=(h=dev(Float32.(o)); for k in 1:length(Ws)-1; h=max.(0f0,Ws[k]*h.+bs[k]); end; Array(tanh.(Ws[end]*h.+bs[end])))
drift(q,t)= DRIFT==0.0 ? zero(q) : DRIFT.*vcat(-sin(q[2]+DRIFTW*t*DT), sin(q[1]+DRIFTW*t*DT), zeros(length(q)-2))
eye(nd)=Matrix{Float64}(I,nd,nd)
function build_sheaf(n,m,nd,edges,pin)
    s=EuclideanSheaf{Float64}(fill(nd,n+m))
    for (i,j) in edges; add_sheaf_edge!(s,i,j,eye(nd),eye(nd)); end
    for i in 1:n,j in 1:m; pin[i,j]!=0 && add_sheaf_edge!(s,i,n+j,eye(nd),eye(nd)); end
    s,collect(n+1:n+m)
end
function topo(n,rng)
    e=Tuple{Int,Int}[(i,i%n+1) for i in 1:n]
    for _ in 1:XEDGES; a=rand(rng,1:n);b=rand(rng,1:n);a!=b&&push!(e,(min(a,b),max(a,b))); end
    pin=zeros(n,n); for i in 1:n; pin[i,i]=1.0; end; e,pin
end
function targ(rng,nd)
    c=(rand(rng,nd).-0.5).*6; A=0.5 .+ rand(rng,nd).*2; om=(0.3 .+ rand(rng,nd).*0.7).*(2π/(K*DT)); ph=rand(rng,nd).*2π; dr=(rand(rng,nd).-0.5).*0.3
    P=Array{Float64}(undef,nd,K+1); for t in 0:K; tt=t*DT; @. P[:,t+1]=c+A*sin(om*tt+ph)+dr*tt; end; P
end

function rollout(N,seed,mode)
    rng=MersenneTwister(seed); edges,pin=topo(N,rng); s,tv=build_sheaf(N,N,ND,edges,pin); d=coboundary_map(s)
    Tg=Array{Float64}(undef,ND,K+1,N); for j in 1:N; Tg[:,:,j]=targ(rng,ND); end
    Q=Array{Float64}(undef,ND,K,N)
    for t in 1:K; qv=Vector(harmonic_extension(s,Dict(tv[j]=>Tg[:,t,j] for j in 1:N))[1]); for i in 1:N; Q[:,t,i]=qv[(i-1)*ND+1:i*ND]; end; end
    x=[Tg[:,1,i].+0.5.*randn(rng,ND) for i in 1:N]
    hx=zeros(ND,HIST,N); hu=zeros(ND,HIST,N)
    pos=Array{Float64}(undef,ND,K,N)
    for t in 1:K
        z=zeros(ND*2N); for i in 1:N; z[(i-1)*ND+1:i*ND]=x[i]; end; for j in 1:N; z[(N+j-1)*ND+1:(N+j)*ND]=Tg[:,t,j]; end
        Lz=d'*(d*z)
        for i in 1:N
            ua=-KGAIN.*Lz[(i-1)*ND+1:i*ND]; e=Q[:,t,i].-x[i]
            obs=vcat(x[i],e,vec(hx[:,:,i]),vec(hu[:,:,i]))
            u = mode===:analytic ? ua : mode===:oracle ? ua.-drift(x[i],t) : ua.+RESS.*Float64.(vec(actor(obs)))
            pos[:,t,i]=x[i]
            if HIST>0; hx[:,1:HIST-1,i]=hx[:,2:HIST,i]; hx[:,HIST,i]=x[i]; hu[:,1:HIST-1,i]=hu[:,2:HIST,i]; hu[:,HIST,i]=u; end
            x[i]=x[i].+DT.*(u.+drift(x[i],t))
        end
    end
    pos,Q,Tg,edges
end
trkerr(pos,Q,N)=[mean(norm(pos[:,t,i].-Q[:,t,i]) for i in 1:N) for t in 1:K]

function panel(pos,Q,edges,t,N,title)
    plt=plot(title=title,legend=false,aspect_ratio=:equal,xlabel="x",ylabel="y")
    for (i,j) in edges; plot!(plt,[pos[1,t,i],pos[1,t,j]],[pos[2,t,i],pos[2,t,j]],color=:gray,alpha=0.25); end
    for i in 1:N
        lo=max(1,t-40); plot!(plt,pos[1,lo:t,i],pos[2,lo:t,i],color=i,alpha=0.4)
        scatter!(plt,[Q[1,t,i]],[Q[2,t,i]],color=i,markershape=:cross,markersize=5,alpha=0.6)
        scatter!(plt,[pos[1,t,i]],[pos[2,t,i]],color=i,markersize=6)
    end
    plt
end

const VN=4
summary=String[]
for seed in SEEDS
    @info "Rendering seed" seed
    pa,Qa,Tg,edges=rollout(VN,seed,:analytic)
    po,Qo,_,_=rollout(VN,seed,:oracle)
    pr,Qr,_,_=rollout(VN,seed,:rl)
    ea=trkerr(pa,Qa,VN); eo=trkerr(po,Qo,VN); er=trkerr(pr,Qr,VN)
    open(joinpath(OUTDIR,"results_$seed.toml"),"w") do io   # metrics for the docs number-autoupdate
        @printf(io,"analytic = %.4f\noracle = %.4f\nrl = %.4f\n", mean(ea), mean(eo), mean(er))
    end
    stride=max(1,K÷150)
    anim=@animate for t in 1:stride:K
        plot(panel(pa,Qa,edges,t,VN,@sprintf("analytic  %.3f",mean(ea))),
             panel(pr,Qr,edges,t,VN,@sprintf("sheaf+RL (history)  %.3f",mean(er))),
             layout=(1,2),size=(1000,500),plot_title="time-varying drift, t=$t (seed $seed) [+ = q*]")
    end
    gif(anim, joinpath(OUTDIR,"anim_tv_$seed.gif"), fps=20)
    ep=plot(title="tracking error ‖x−q*‖ (seed $seed)",xlabel="step",ylabel="mean ‖x−q*‖")
    plot!(ep,ea,label="analytic",lw=2); plot!(ep,eo,label="oracle",lw=2,ls=:dot); plot!(ep,er,label="sheaf+RL",lw=2,ls=:dash)
    png(ep,joinpath(OUTDIR,"trkerr_$seed.png"))
    push!(summary,@sprintf("seed %d:  analytic %.4f  oracle %.4f  sheaf+RL %.4f",seed,mean(ea),mean(eo),mean(er)))
end
println("\n=== time-varying drift: analytic vs oracle vs sheaf+RL (mean ‖x−q*‖) ===")
for l in summary; println("  ",l); end
println("GIFs + PNGs in: $OUTDIR")
