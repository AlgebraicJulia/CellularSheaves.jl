#
# Phase B (sheaf-native) — behavior-clone the decentralized policy π(x_i, e_i) from the
# no-drift dataset produced by collect_bc_data.jl.
#
#   obs_i = [ x_i ; e_i = q*_i - x_i ; (optional) future_target(1..HORIZON) ]
#   label = u_i   (the paper's analytic control)
#
# HORIZON=0 is the deliberately myopic policy (reactive feedback → chases the after-image, matches
# the analytic law). HORIZON>0 hands the policy its own future target window — anticipation, the
# model-free way to LEAD rather than chase (no [A B] involved).
#
# Mirrors train_bc_conic.jl conventions: standardize obs+labels, 3×256 MLP, Adam+MSE, save
# model_state + normalization stats so eval can invert. GPU via Flux/CUDA. Run on the CLUSTER.
#
#   julia --project=examples/RL examples/RL/train/behaviour_cloning.jl
# Env: BCS_DATA, BCS_EPOCHS, BCS_BATCH, BCS_LR, BCS_HORIZON, BCS_SEED, BCS_OUT
#
using LinearAlgebra, Random, Statistics, Printf, Flux, CUDA, cuDNN, JLD2
const dev = CUDA.functional() ? gpu : cpu
CUDA.functional() && @info "GPU: $(CUDA.name(CUDA.device()))"

const DATA    = get(ENV,"BCS_DATA", joinpath(@__DIR__,"..","cache","bc_sheaf_dataset.jld2"))
const EPOCHS  = parse(Int,     get(ENV,"BCS_EPOCHS","60"))
const BATCH   = parse(Int,     get(ENV,"BCS_BATCH","2048"))
const LR      = parse(Float64, get(ENV,"BCS_LR","3e-4"))
const HORIZON = parse(Int,     get(ENV,"BCS_HORIZON","0"))    # 0 = reactive/myopic; >0 = anticipation
const SEED    = parse(Int,     get(ENV,"BCS_SEED","1"))
const OUT     = get(ENV,"BCS_OUT", "/blue/fairbanksj/$(get(ENV,"USER","itaykadosh"))/sheaf_rl/bc_sheaf.jld2")
Random.seed!(SEED); isdir(dirname(OUT)) || mkpath(dirname(OUT))

@info "Loading dataset" DATA
D  = load(DATA)
Ns = D["Ns"]; Xs = D["X"]; Qs = D["Q"]; Us = D["U"]; Tgs = D["Tg"]
ND = D["meta"]["ND"]; K = D["meta"]["K"]
const OBS = 2 * ND + ND * HORIZON                          # [x_i ; e_i ; future_tgt(1..HORIZON)]
futtgt(Tg, t, i) = reduce(vcat, [Tg[:, min(t + τ, K + 1), i] for τ in 1:HORIZON])

@info "Assembling (obs,label) pairs..." OBS nu=ND HORIZON
xbuf = Vector{Float32}[]; ybuf = Vector{Float32}[]
for q in 1:length(Ns)
    N = Ns[q]; X = Xs[q]; Q = Qs[q]; U = Us[q]; Tg = Tgs[q]; ks = size(X, 2)
    for t in 1:ks, i in 1:N
        e   = Q[:, t, i] .- X[:, t, i]
        obs = HORIZON > 0 ? vcat(X[:, t, i], e, futtgt(Tg, t + 1, i)) : vcat(X[:, t, i], e)
        push!(xbuf, Float32.(obs))
        push!(ybuf, Float32.(U[:, t, i]))
    end
end
Xtr = reduce(hcat, xbuf); Ytr = reduce(hcat, ybuf)
@info "Dataset assembled" nsamples=size(Xtr, 2) OBS nu=ND
# data-health checks visible in the log
nbad = count(!isfinite, Xtr) + count(!isfinite, Ytr)
@printf("  obs  range [%.3f, %.3f]   label range [%.3f, %.3f]   non-finite=%d\n",
        minimum(Xtr), maximum(Xtr), minimum(Ytr), maximum(Ytr), nbad)
@printf("  sample obs[:,1] = %s\n", string(round.(Xtr[:, 1]; digits=3)))
@printf("  sample lbl[:,1] = %s\n", string(round.(Ytr[:, 1]; digits=3)))
flush(stdout)
nbad == 0 || error("non-finite values in training data — fix datagen before training")

# standardize obs + labels (store stats to invert at eval)
const XM = vec(mean(Xtr; dims=2)); const XS = vec(std(Xtr; dims=2)) .+ 1f-4
const YM = vec(mean(Ytr; dims=2)); const YS = vec(std(Ytr; dims=2)) .+ 1f-4
Xn = dev((Xtr .- XM) ./ XS); Yn = dev((Ytr .- YM) ./ YS)

model = dev(Chain(Dense(OBS,256,relu), Dense(256,256,relu), Dense(256,256,relu), Dense(256,ND)))
opt   = Flux.setup(Adam(LR), model)
n = size(Xn, 2)
@info "Training BC" EPOCHS BATCH n
t0 = time()
for ep in 1:EPOCHS
    perm = dev(randperm(n)); tot = 0f0; nb = 0
    for sidx in 1:BATCH:n
        idx = perm[sidx:min(sidx + BATCH - 1, n)]
        l, gr = Flux.withgradient(m -> Flux.mse(m(Xn[:, idx]), Yn[:, idx]), model)
        isfinite(l) || error("non-finite loss at epoch $ep — diverged (lower BCS_LR?)")
        Flux.update!(opt, model, gr[1]); tot += l; nb += 1
    end
    if ep == 1 || ep % 5 == 0
        @printf("  epoch %3d  mse(norm)=%.5f  elapsed=%.1fs\n", ep, tot / nb, time() - t0)
        flush(stdout)
    end
end

# Save BOTH the Flux state (for cluster reload) AND raw weight arrays (for Flux-free / cross-version
# local viz — the MLP is Dense(relu)×3 → Dense, reconstructable with plain matmuls).
cmodel = cpu(model)
Ws = [Array(l.weight) for l in cmodel.layers]
bs = [Array(l.bias)   for l in cmodel.layers]
jldsave(OUT; model_state=Flux.state(cmodel), Ws=Ws, bs=bs,
    XM=Array(XM), XS=Array(XS), YM=Array(YM), YS=Array(YS),
    OBS, ND, HORIZON, data=DATA, meta=D["meta"])
@info "Saved BC policy" OUT layers=length(Ws)
