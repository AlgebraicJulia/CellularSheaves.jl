# GK kmax diagnostic: run solve_logged at kmax ∈ {10,20,40} on one organic problem, dump sigma2min per
# (iter,α) so the within-sweep slope of log σ²min vs log α can be measured per kmax. Registered test:
# slope flattens with kmax at low/moderate α·σ² if the σ²min perch bias is GK TRUNCATION; stays broken at
# high α if it's the 1/(1-μ) RECOVERY amplification. Sets σ²min's documentation permanently.
using CellularSheaves.IPM: IPMSettings, init, solve_logged, write_oracle_csv
import CellularSheaves.IPM as IPM
const EX = dirname(@__DIR__)
include("$EX/e02.jl")
gp(x) = x isa Tuple ? x[1] : x
prob = gp(build_qm_problem(generate_qm_instance(; N = 10)))
settings = IPMSettings{Float64}(feas_tol = 1e-10, gap_tol = 1e-10, itmax = 300)
OUT = get(ENV, "ORACLE_OUT", @__DIR__)
for k in (10, 20, 40)
    s0 = init(prob, settings)
    _, recs = solve_logged(s0, collect(IPM.DEFAULT_ALPHA_GRID); kmax = k)
    write_oracle_csv(joinpath(OUT, "kmaxdiag_e02_k$(k).csv"), recs)
    println("kmax=$k -> $(length(recs)) records")
end
println("KMAX DIAG DONE")
