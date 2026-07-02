#
# Headroom gate (no trained model needed). For each (DRIFT, KGAIN), compare:
#   analytic  u = −k·(Lz)_i              — the paper's law, IGNORES the drift
#   oracle    u = −k·(Lz)_i − f(x_i)     — analytic + EXACT drift cancellation (the RL upper bound)
# The gap = mean tracking error(analytic) − error(oracle) = how much a policy that LEARNS f could win.
# If the gap is ~0 the regime has no headroom (RL would just tie); pick a regime where it's clearly > 0.
#
#   julia -t auto --project=examples/RL examples/RL/eval/drift_headroom.jl
# Env: HR_NSEED, HR_N, HR_DRIFTS (csv), HR_GAINS (csv)
#
using CellularSheaves                       # EuclideanSheaf, add_sheaf_edge!, coboundary_map
using LinearAlgebra, Random, Statistics, Printf
import Base.Threads: @threads, nthreads

const ND = 2; const DT = 0.02; const K = 400; const XEDGES = 1
const NSEED  = parse(Int, get(ENV,"HR_NSEED","40"))
const N      = parse(Int, get(ENV,"HR_N","4"))
const DRIFTS = parse.(Float64, split(get(ENV,"HR_DRIFTS","0.0,0.5,1.0,2.0,4.0"), ","))
const GAINS  = parse.(Float64, split(get(ENV,"HR_GAINS","8.0,4.0,2.0"), ","))

eye(nd) = Matrix{Float64}(I, nd, nd)
function build_sheaf(N, M, nd, edges, pin)
    s = EuclideanSheaf{Float64}(fill(nd, N + M))
    for (i, j) in edges; add_sheaf_edge!(s, i, j, eye(nd), eye(nd)); end
    for i in 1:N, j in 1:M; pin[i, j] != 0 && add_sheaf_edge!(s, i, N + j, eye(nd), eye(nd)); end
    s
end
function topo(N, rng)
    e = Tuple{Int,Int}[(i, i % N + 1) for i in 1:N]
    for _ in 1:XEDGES; a = rand(rng, 1:N); b = rand(rng, 1:N); a != b && push!(e, (min(a,b), max(a,b))); end
    pin = zeros(N, N); for i in 1:N; pin[i, i] = 1.0; end
    e, pin
end
function targ(rng, nd)
    c = (rand(rng,nd).-0.5).*6; A = 0.5 .+ rand(rng,nd).*2
    om = (0.3 .+ rand(rng,nd).*0.7).*(2π/(K*DT)); ph = rand(rng,nd).*2π; dr = (rand(rng,nd).-0.5).*0.3
    P = Array{Float64}(undef, nd, K+1); for t in 0:K; tt = t*DT; @. P[:,t+1] = c + A*sin(om*tt+ph) + dr*tt; end; P
end
drift(q, amp) = amp == 0.0 ? zero(q) : amp .* vcat(-sin(q[2]), sin(q[1]), zeros(length(q)-2))

# mean tracking error |x − target| under a controller; mode = :analytic | :oracle
function rollout(N, seed, kgain, damp, mode)
    rng = MersenneTwister(seed); edges, pin = topo(N, rng); s = build_sheaf(N, N, ND, edges, pin); d = coboundary_map(s)
    Tg = Array{Float64}(undef, ND, K+1, N); for j in 1:N; Tg[:,:,j] = targ(rng, ND); end
    x = [Tg[:,1,i] .+ 0.5.*randn(rng, ND) for i in 1:N]; err = 0.0; cnt = 0
    for t in 1:K
        z = zeros(ND*2N); for i in 1:N; z[(i-1)*ND+1:i*ND] = x[i]; end
        for j in 1:N; z[(N+j-1)*ND+1:(N+j)*ND] = Tg[:,t,j]; end
        Lz = d' * (d * z)
        for i in 1:N
            f = drift(x[i], damp); ua = -kgain .* Lz[(i-1)*ND+1:i*ND]
            u = mode === :oracle ? ua .- f : ua
            err += norm(x[i] .- Tg[:,t,i]); cnt += 1
            x[i] = x[i] .+ DT .* (u .+ f)
        end
    end
    err / cnt
end

println("Headroom sweep — N=$N, $NSEED seeds, K=$K, DT=$DT, threads=$(nthreads())"); flush(stdout)
@printf("%7s %7s %10s %10s %10s %8s\n", "DRIFT", "KGAIN", "analytic", "oracle", "gap", "gap%")
for kg in GAINS, dm in DRIFTS
    seeds = collect(1:NSEED)
    a = Vector{Float64}(undef, NSEED); o = Vector{Float64}(undef, NSEED)
    @threads for q in 1:NSEED
        a[q] = rollout(N, 1000+seeds[q], kg, dm, :analytic)
        o[q] = rollout(N, 1000+seeds[q], kg, dm, :oracle)
    end
    ma = mean(a); mo = mean(o)
    @printf("%7.2f %7.1f %10.4f %10.4f %10.4f %7.1f%%\n", dm, kg, ma, mo, ma-mo, 100*(ma-mo)/mo)
    flush(stdout)
end
println("Done. Pick the (DRIFT,KGAIN) with a clearly positive gap% for the RL stage.")
