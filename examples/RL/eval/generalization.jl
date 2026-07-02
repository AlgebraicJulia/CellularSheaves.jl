#
# Zero-shot drift generalization: take the policy trained on the sinusoidal swirl and evaluate it
# — WITHOUT retraining — on drift fields of different form (constant, other frequencies, tanh
# "square" swirl, multi-frequency). Tests whether the history window learned a field-agnostic
# drift estimator (f ≈ Δx/dt − u, fed forward) or just memorized the training sinusoid.
#
#   julia --project=examples/RL examples/RL/eval/generalization.jl
#
get(ENV,"DISPLAY","")=="" && (ENV["GKSwstype"]="100")
using CellularSheaves
import CellularSheaves.NetworkSheaves.EuclideanSheaves: harmonic_extension
using LinearAlgebra, Random, JLD2, Printf, Statistics

const MODEL = get(ENV,"EG_MODEL", joinpath(@__DIR__,"..","cache","rl_multiagent_f2.jld2"))
const NSCN  = parse(Int, get(ENV,"EG_NSCN","25"))
M=load(MODEL)
const Ws=M["Ws"]; const bs=M["bs"]; const EDGES=M["EDGES"]; const PIN=M["PIN"]
const ND=M["ND"]; const N=M["N"]; const MT=M["M"]; const KGAIN=M["KGAIN"]; const DT=M["DT"]
const K=M["K"]; const RESS=M["RESS"]; const A=M["DRIFT"]; const W=M["DRIFTW"]; const HIST=M["HIST"]
@info "Drift-generalization eval" MODEL trained_drift=(A,W) HIST NSCN

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

# drift fields to test (all bounded ~A); the policy only ever saw `sin_swirl_w3`
const DRIFTS = [
  ("sin swirl ω=3  (TRAIN)", (q,t)-> A.*[-sin(q[2]+W*t*DT),       sin(q[1]+W*t*DT),       0.0]),
  ("sin swirl ω=1",          (q,t)-> A.*[-sin(q[2]+1.0*t*DT),     sin(q[1]+1.0*t*DT),     0.0]),
  ("sin swirl ω=6",          (q,t)-> A.*[-sin(q[2]+6.0*t*DT),     sin(q[1]+6.0*t*DT),     0.0]),
  ("constant wind",          (q,t)-> A.*[0.7, -0.5, 0.0]),
  ("tanh 'square' swirl ω=3",(q,t)-> A.*[-tanh(2*(q[2]+W*t*DT)),  tanh(2*(q[1]+W*t*DT)),  0.0]),
  ("multi-freq swirl",       (q,t)-> A.*[-0.6sin(q[2]+W*t*DT)-0.4sin(2.3*t*DT), 0.6sin(q[1]+W*t*DT)+0.4cos(1.7*t*DT), 0.0]),
  ("linear shear",           (q,t)-> [0.4*q[2], -0.4*q[1], 0.0]),
]

function rollout(seed, drift, mode)
    rng=MersenneTwister(seed); Tg=targets(rng); Q=Array{Float64}(undef,ND,K,N)
    for t in 1:K; qv=Vector(harmonic_extension(SHEAF,Dict(TVERTS[j]=>Tg[:,t,j] for j in 1:MT))[1]); for i in 1:N; Q[:,t,i]=qv[(i-1)*ND+1:i*ND]; end; end
    x=[Q[:,1,i].+0.4.*randn(rng,ND) for i in 1:N]; hx=zeros(ND,HIST,N); hu=zeros(ND,HIST,N); tot=0.0
    for t in 1:K
        z=zeros(ND*(N+MT)); for i in 1:N; z[(i-1)*ND+1:i*ND]=x[i]; end; for j in 1:MT; z[(N+j-1)*ND+1:(N+j)*ND]=Tg[:,t,j]; end
        Lz=DMAT'*(DMAT*z)
        for i in 1:N
            ua=-KGAIN.*Lz[(i-1)*ND+1:i*ND]
            obs=vcat(x[i],Q[:,t,i].-x[i],vec(hx[:,:,i]),vec(hu[:,:,i]))
            u = mode===:analytic ? ua : mode===:oracle ? ua.-drift(x[i],t) : ua.+RESS.*Float64.(vec(actor(obs)))
            tot+=norm(x[i].-Q[:,t,i])
            if HIST>0; hx[:,1:HIST-1,i]=hx[:,2:HIST,i];hx[:,HIST,i]=x[i];hu[:,1:HIST-1,i]=hu[:,2:HIST,i];hu[:,HIST,i]=u; end
            x[i]=x[i].+DT.*(u.+drift(x[i],t))
        end
    end
    tot/(K*N)
end
meanerr(drift,mode)=mean(rollout(s,drift,mode) for s in 1:NSCN)

@printf("\n%-26s  %8s  %8s  %8s   %s\n","drift field","analytic","d-oracle","sheaf+RL","RL recovers")
println("-"^78)
for (name,d) in DRIFTS
    a=meanerr(d,:analytic); o=meanerr(d,:oracle); r=meanerr(d,:rl)
    rec = a>o ? 100*(a-r)/(a-o) : NaN     # % of analytic→oracle gap closed
    @printf("%-26s  %8.4f  %8.4f  %8.4f   %5.0f%%\n",name,a,o,r,rec)
end
