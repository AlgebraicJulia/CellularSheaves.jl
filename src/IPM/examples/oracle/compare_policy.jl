# A/B driver: run the live solver (controller ON) on an array of problems under the current α policy
# (policy=0) and Tier 1 / Policy 1 (policy=1), and report the headline metric — cumulative F^-1
# applications over the whole solve — plus iters, final status, and how many iterations paid refinement.
#
# Usage:  julia --project=examples examples/oracle/compare_policy.jl [key1 key2 ...]
#   default keys: e01 e02 e04 07 ; both solvers each.

using CellularSheaves
using CellularSheaves.IPM: HSDSettings, IPMSettings, init
import CellularSheaves.IPM as IPM
using LinearAlgebra, SparseArrays

const EX = dirname(@__DIR__)   # examples/
gp(x) = x isa Tuple ? x[1] : x

for k in ("e01", "e02", "e04", "07", "e08", "e15", "e11", "e13")
    include(joinpath(EX, k == "07" ? "e07.jl" : "$k.jl"))
end

function buildprob(key)
    key == "e01" && return gp(build_merton_chain(merton_instance(; nx=21, olap=5, nsteps=10)))
    key == "e02" && return gp(build_qm_problem(generate_qm_instance(; N=10)))
    key == "e04" && return gp(build_sos_spline(sos_instance(; n=13)))
    key == "07"  && return gp(build_poisson_tv(poisson_instance(; N=128, Tsz=16, m=16, k=-1, K=6, R=12, q=3, seed=3)))
    # e11 (dense-Q SOC FIR equalizer) and e13 (dense-Q SOC fuel landing): endgame window-closers on IPM
    key == "e11" && (i = equalizer_instance(); return gp(build_equalizer(make_system(i, 8), i; form=:exact, white=true)))
    key == "e13" && return gp(build_landing(landing_instance(), 60.0, 20))
    if key == "e08"
        dt = 4.0 / 64; inst = exec_instance(; n=1, T=64)
        return gp(build_exec(inst; Σ=fill(dt, 1, 1), X0=[1.5], η=[1.0 / sqrt(dt)], γ=1.0))
    end
    if key == "e15"
        train = [[0.0, 60.0 * (v - 1), 0.0, 0.0, 0.0, 0.0] for v in 1:V_BASE]
        return gp(build_cw(train, [zeros(6) for _ in 1:V_BASE], N_BASE; gam=1e-1, kap=1e-1))
    end
    error("no recipe for $key")
end

# total applications of F^-1 over the solve: per role, base + refinement CRAIG, summed over all iterations.
function totalapps(hist)
    a = sum(hist.pbase) + sum(hist.prefn) + sum(hist.cbase) + sum(hist.crefn)
    hasproperty(hist, :wbase) && (a += sum(hist.wbase) + sum(hist.wrefn))
    return a
end

# iterations that paid at least one refinement pass in some role (the md's "overshoot" indicator).
function refinefired(hist)
    n = length(hist.ppass)
    fired = 0
    for i in 1:n
        p = hist.ppass[i] + hist.cpass[i]
        hasproperty(hist, :wpass) && (p += hist.wpass[i])
        p > 0 && (fired += 1)
    end
    return fired
end

function runone(prob, solv, policy; tol=1e-10)
    settings = solv == "ipm" ?
        IPMSettings{Float64}(; feas_tol=tol, gap_tol=tol, itmax=300, policy=policy) :
        HSDSettings{Float64}(; feas_tol=tol, gap_tol=tol, itmax=300, policy=policy)
    s = init(prob, settings)
    res = IPM.solve!(s)
    return (apps=totalapps(s.hist), iters=length(s.hist), fired=refinefired(s.hist),
            status=string(res.status))
end

keys = length(ARGS) ≥ 1 ? ARGS : ["e01", "e02", "e04", "07", "e11", "e13"]

println(rpad("problem", 10), rpad("solver", 7),
        rpad("apps(cur)", 11), rpad("apps(T1)", 11), rpad("Δapps", 8),
        rpad("it(cur)", 8), rpad("it(T1)", 8), rpad("fired c/1", 11),
        rpad("status(cur)", 16), "status(T1)")
println("-"^110)

for key in keys
    prob = buildprob(key)
    for solv in ("ipm", "hsd")
        c = runone(prob, solv, 0)
        t = runone(prob, solv, 1)
        Δ = t.apps - c.apps
        println(rpad(key, 10), rpad(solv, 7),
                rpad(c.apps, 11), rpad(t.apps, 11), rpad(Δ, 8),
                rpad(c.iters, 8), rpad(t.iters, 8), rpad("$(c.fired)/$(t.fired)", 11),
                rpad(c.status, 16), t.status)
    end
end
