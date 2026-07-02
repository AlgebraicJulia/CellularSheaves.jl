#
# 2x4 zero-shot drift-generalization animation: the policy trained ONLY on the ω=3 sinusoidal
# swirl, evaluated (no retraining) on four drift fields. Rows = controller (analytic / sheaf+RL),
# columns = drift field. Same scenario seed across every panel. Flux-free local render.
#
#   julia --project=examples/RL examples/RL/viz/generalization.jl
# Env: EG_MODEL, EG_SEED, EG_OUTDIR
#
get(ENV,"DISPLAY","")=="" && (ENV["GKSwstype"]="100")
using CellularSheaves
import CellularSheaves.NetworkSheaves.EuclideanSheaves: harmonic_extension
using LinearAlgebra, Random, JLD2, Printf, Statistics, Plots
gr()

const MODEL = get(ENV,"EG_MODEL", joinpath(@__DIR__,"..","cache","rl_multiagent_f2.jld2"))
const SEED  = parse(Int, get(ENV,"EG_SEED","7001"))
const OUTDIR= get(ENV,"EG_OUTDIR", joinpath(@__DIR__,"..","cache","viz_multiagent")); mkpath(OUTDIR)
M=load(MODEL)
const Ws=M["Ws"]; const bs=M["bs"]; const EDGES=M["EDGES"]; const PIN=M["PIN"]
const ND=M["ND"]; const N=M["N"]; const MT=M["M"]; const KGAIN=M["KGAIN"]; const DT=M["DT"]
const K=M["K"]; const RESS=M["RESS"]; const A=M["DRIFT"]; const W=M["DRIFTW"]; const HIST=M["HIST"]

actor(o)=(h=Float32.(o); for k in 1:length(Ws)-1; h=max.(0f0,Ws[k]*h.+bs[k]); end; tanh.(Ws[end]*h.+bs[end]))
eye(nd)=Matrix{Float64}(I,nd,nd)
function build_sheaf()
    s=EuclideanSheaf{Float64}(fill(ND,N+MT))
    for (i,j) in EDGES; add_sheaf_edge!(s,i,j,eye(ND),eye(ND)); end
    for i in 1:N,j in 1:MT; PIN[i,j]!=0 && add_sheaf_edge!(s,i,N+j,eye(ND),eye(ND)); end
    s, collect(N+1:N+MT)
end
const SHEAF,TVERTS=build_sheaf(); const DMAT=coboundary_map(SHEAF)
function targets(rng)
    centers=[[-3.0,-3,0.0],[3.0,-3,0.0],[3.0,3,0.0],[-3.0,3,0.0]]; Tg=Array{Float64}(undef,ND,K+1,MT)
    for j in 1:MT
        c=centers[j].+0.3.*randn(rng,ND); om=(0.6+0.4*rand(rng))*(2π/(K*DT)); ph=rand(rng)*2π
        Ax=1.8+0.5*rand(rng); Ay=1.2+0.4*rand(rng); Az=1.3+0.5*rand(rng)
        for t in 0:K; tt=t*DT
            Tg[1,t+1,j]=c[1]+Ax*sin(om*tt+ph); Tg[2,t+1,j]=c[2]+Ay*sin(2*om*tt+ph); Tg[3,t+1,j]=c[3]+Az*sin(om*tt+ph+π/2)
        end
    end; Tg
end
pinof(i)=(j=findfirst(!=(0.0),PIN[i,:]); j===nothing ? 0 : j)

const DRIFTS = [
  ("sin swirl ω=3 (train)", (q,t)-> A.*[-sin(q[2]+W*t*DT),      sin(q[1]+W*t*DT),      0.0]),
  ("constant wind",         (q,t)-> A.*[0.7, -0.5, 0.0]),
  ("tanh 'square' swirl",   (q,t)-> A.*[-tanh(2*(q[2]+W*t*DT)), tanh(2*(q[1]+W*t*DT)), 0.0]),
  ("linear shear",          (q,t)-> [0.4*q[2], -0.4*q[1], 0.0]),
]

function rollout(drift,mode)
    rng=MersenneTwister(SEED); Tg=targets(rng); Q=Array{Float64}(undef,ND,K,N)
    for t in 1:K; qv=Vector(harmonic_extension(SHEAF,Dict(TVERTS[j]=>Tg[:,t,j] for j in 1:MT))[1]); for i in 1:N; Q[:,t,i]=qv[(i-1)*ND+1:i*ND]; end; end
    x=[Q[:,1,i].+0.4.*randn(rng,ND) for i in 1:N]; hx=zeros(ND,HIST,N); hu=zeros(ND,HIST,N); pos=Array{Float64}(undef,ND,K,N); tot=0.0
    for t in 1:K
        z=zeros(ND*(N+MT)); for i in 1:N; z[(i-1)*ND+1:i*ND]=x[i]; end; for j in 1:MT; z[(N+j-1)*ND+1:(N+j)*ND]=Tg[:,t,j]; end
        Lz=DMAT'*(DMAT*z)
        for i in 1:N
            ua=-KGAIN.*Lz[(i-1)*ND+1:i*ND]; obs=vcat(x[i],Q[:,t,i].-x[i],vec(hx[:,:,i]),vec(hu[:,:,i]))
            u = mode===:analytic ? ua : ua.+RESS.*Float64.(vec(actor(obs)))
            pos[:,t,i]=x[i]; tot+=norm(x[i].-Q[:,t,i])
            if HIST>0; hx[:,1:HIST-1,i]=hx[:,2:HIST,i];hx[:,HIST,i]=x[i];hu[:,1:HIST-1,i]=hu[:,2:HIST,i];hu[:,HIST,i]=u; end
            x[i]=x[i].+DT.*(u.+drift(x[i],t))
        end
    end
    pos,Tg,tot/(K*N)
end

# roll all 8 cells (4 drifts × {analytic, rl}) once; keep trajectories + mean errors
cells=[]   # (row, col) → (pos, Tg, err, title)
for (di,(name,d)) in enumerate(DRIFTS)
    pa,Tg,ea = rollout(d,:analytic); pr,_,er = rollout(d,:rl)
    push!(cells, (1,di,pa,Tg,ea,@sprintf("analytic · %s  %.2f",name,ea)))
    push!(cells, (2,di,pr,Tg,er,@sprintf("sheaf+RL · %s  %.2f",name,er)))
end

function panel(pos,Tg,t,title)
    plt=plot(title=title,titlefontsize=9,legend=false,aspect_ratio=:equal,xlims=(-7,7),ylims=(-7,7),xticks=false,yticks=false)
    for (i,j) in EDGES; plot!(plt,[pos[1,t,i],pos[1,t,j]],[pos[2,t,i],pos[2,t,j]],color=:gray,alpha=0.2); end
    for j in 1:MT; lo=max(1,t-40); plot!(plt,Tg[1,lo:t,j],Tg[2,lo:t,j],color=:black,alpha=0.3); scatter!(plt,[Tg[1,t,j]],[Tg[2,t,j]],color=:black,markershape=:star5,markersize=7); end
    for i in 1:N; c=pinof(i)==0 ? :gray : pinof(i); scatter!(plt,[pos[1,t,i]],[pos[2,t,i]],color=c,markersize=4,markershape=(pinof(i)==0 ? :diamond : :circle)); end
    plt
end

# order cells row-major for layout=(2,4): row1 = analytic×4, row2 = rl×4
ordered = vcat([c for c in cells if c[1]==1], [c for c in cells if c[1]==2])
stride=max(1,K÷120)
anim=@animate for t in 1:stride:K
    plot([panel(c[3],c[4],t,c[6]) for c in ordered]..., layout=(2,4), size=(2000,1000),
         plot_title="zero-shot drift generalization (trained on ω=3 swirl only), t=$t")
end
gif(anim, joinpath(OUTDIR,"anim_generalization.gif"), fps=18)
@printf("\nzero-shot generalization (seed %d):\n",SEED)
for (di,(name,_)) in enumerate(DRIFTS)
    ea=[c[5] for c in cells if c[1]==1 && c[2]==di][1]; er=[c[5] for c in cells if c[1]==2 && c[2]==di][1]
    @printf("  %-22s  analytic %.3f   sheaf+RL %.3f\n",name,ea,er)
end
println("→ ", joinpath(OUTDIR,"anim_generalization.gif"))
