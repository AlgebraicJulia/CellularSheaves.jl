#
# Is the RL margin "smarter control" or just "more authority"? Measure mean control magnitude
# ‖u‖ alongside mean tracking error for analytic / drift-oracle / sheaf+RL on the trained drift.
# If RL's ‖u‖ is much larger than analytic's, the sub-oracle gap is partly extra gain.
#
get(ENV,"DISPLAY","")=="" && (ENV["GKSwstype"]="100")
using CellularSheaves
import CellularSheaves.NetworkSheaves.EuclideanSheaves: harmonic_extension
using LinearAlgebra, Random, JLD2, Printf, Statistics
const MODEL=get(ENV,"DM_MODEL",joinpath(@__DIR__,"..","cache","rl_multiagent_f2.jld2")); const NSCN=parse(Int,get(ENV,"DM_NSCN","25"))
M=load(MODEL)
const Ws=M["Ws"];const bs=M["bs"];const EDGES=M["EDGES"];const PIN=M["PIN"]
const ND=M["ND"];const N=M["N"];const MT=M["M"];const KGAIN=M["KGAIN"];const DT=M["DT"]
const K=M["K"];const RESS=M["RESS"];const A=M["DRIFT"];const W=M["DRIFTW"];const HIST=M["HIST"]
actor(o)=(h=Float32.(o);for k in 1:length(Ws)-1;h=max.(0f0,Ws[k]*h.+bs[k]);end;tanh.(Ws[end]*h.+bs[end]))
drift(q,t)=A.*[-sin(q[2]+W*t*DT),sin(q[1]+W*t*DT),0.0]
eye(n)=Matrix{Float64}(I,n,n)
function build();s=EuclideanSheaf{Float64}(fill(ND,N+MT));for (i,j) in EDGES;add_sheaf_edge!(s,i,j,eye(ND),eye(ND));end;for i in 1:N,j in 1:MT;PIN[i,j]!=0&&add_sheaf_edge!(s,i,N+j,eye(ND),eye(ND));end;s,collect(N+1:N+MT);end
const S,TV=build();const D=coboundary_map(S)
function targets(rng);C=[[-3.0,-3,0.0],[3.0,-3,0.0],[3.0,3,0.0],[-3.0,3,0.0]];Tg=Array{Float64}(undef,ND,K+1,MT)
 for j in 1:MT;c=C[j].+0.3.*randn(rng,ND);om=(0.6+0.4*rand(rng))*(2π/(K*DT));ph=rand(rng)*2π;Ax=1.8+0.5*rand(rng);Ay=1.2+0.4*rand(rng);Az=1.3+0.5*rand(rng)
  for t in 0:K;tt=t*DT;Tg[1,t+1,j]=c[1]+Ax*sin(om*tt+ph);Tg[2,t+1,j]=c[2]+Ay*sin(2*om*tt+ph);Tg[3,t+1,j]=c[3]+Az*sin(om*tt+ph+π/2);end;end;Tg;end
function roll(seed,mode);rng=MersenneTwister(seed);Tg=targets(rng);Q=Array{Float64}(undef,ND,K,N)
 for t in 1:K;qv=Vector(harmonic_extension(S,Dict(TV[j]=>Tg[:,t,j] for j in 1:MT))[1]);for i in 1:N;Q[:,t,i]=qv[(i-1)*ND+1:i*ND];end;end
 x=[Q[:,1,i].+0.4.*randn(rng,ND) for i in 1:N];hx=zeros(ND,HIST,N);hu=zeros(ND,HIST,N);se=0.0;su=0.0
 for t in 1:K;z=zeros(ND*(N+MT));for i in 1:N;z[(i-1)*ND+1:i*ND]=x[i];end;for j in 1:MT;z[(N+j-1)*ND+1:(N+j)*ND]=Tg[:,t,j];end;Lz=D'*(D*z)
  for i in 1:N;ua=-KGAIN.*Lz[(i-1)*ND+1:i*ND];obs=vcat(x[i],Q[:,t,i].-x[i],vec(hx[:,:,i]),vec(hu[:,:,i]))
   u = mode===:analytic ? ua : mode===:oracle ? ua.-drift(x[i],t) : ua.+RESS.*Float64.(vec(actor(obs)))
   se+=norm(x[i].-Q[:,t,i]);su+=norm(u)
   if HIST>0;hx[:,1:HIST-1,i]=hx[:,2:HIST,i];hx[:,HIST,i]=x[i];hu[:,1:HIST-1,i]=hu[:,2:HIST,i];hu[:,HIST,i]=u;end
   x[i]=x[i].+DT.*(u.+drift(x[i],t));end;end;(se/(K*N),su/(K*N));end
mn(mode)=(e=0.0;u=0.0;for s in 1:NSCN;(ee,uu)=roll(s,mode);e+=ee;u+=uu;end;(e/NSCN,u/NSCN))
@printf("\n%-12s  %10s  %10s\n","controller","mean‖x−q*‖","mean‖u‖"); println("-"^36)
for m in (:analytic,:oracle,:rl); (e,u)=mn(m); @printf("%-12s  %10.4f  %10.4f\n",string(m),e,u); end
