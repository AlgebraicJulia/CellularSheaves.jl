# Independent audit of every quantity published by the Diffusion versus Direct
# example. Each check re-derives the number from the sheaf rather than from the
# benchmark's own bookkeeping: H and B come from `sheaf_laplacian_matrix`, the
# condition number from a fresh dense eigendecomposition, and the gain ceilings
# from first principles. The analytic tick prediction is the important one --
# internal consistency alone will happily certify a wrong headline number.
#
#     julia --project=docs docs/scripts/verify_coordination_benchmarks.jl

using CellularSheaves, Printf, LinearAlgebra, Graphs, SparseArrays
using CellularSheaves.ControlSheaves.CoordinationBenchmarks
const CB = CoordinationBenchmarks
fails = String[]; chk(ok,msg) = ok ? nothing : push!(fails,msg)
predict(rate,tol) = ceil(Int, log(tol)/log1p(-rate))
r = benchmark_coordination()
@printf("sweep: %d scenarios\n", length(r))
for row in r.rows
    s = coordination_scenario(row.name; size_parameter=row.size_parameter, coverage=row.coverage)
    tag = "$(row.label) n=$(row.agents) $(row.coverage)"
    L = Matrix(sheaf_laplacian_matrix(s.sheaf)); nf = s.dim*s.nagents
    chk(norm(Matrix(s.H)-L[1:nf,1:nf],Inf)<1e-10, "$tag: H")
    chk(norm(Matrix(s.Bmat)+L[1:nf,nf+1:end],Inf)<1e-10, "$tag: B")
    ev = eigvals(Symmetric(Matrix(s.H)))
    chk(isapprox(row.condition, last(ev)/first(ev); rtol=1e-8), "$tag: kappa")
    chk(first(ev)>0, "$tag: not posdef")
    chk(isapprox(row.direct.gain, r.safety*2/r.dt; rtol=1e-8), "$tag: direct gain")
    chk(isapprox(row.diffusion.gain, r.safety*2/(r.dt*last(ev)); rtol=1e-8), "$tag: diff gain")
    chk(row.diffusion.slots == row.diffusion.ticks*Graphs.Δ(s.agent_graph), "$tag: diff slots")
    chk(iseven(row.direct.slots), "$tag: direct slots odd")
    st = tree_statistics(direct_plan(s))
    st.supernodes<=2 && chk(row.direct.slots >= 2s.nagents, "$tag: degenerate not 2n")
    chk(row.treewidth == max(st.bag_width÷s.dim-1,0), "$tag: treewidth")
    chk(row.fill == st.fill && st.fill >= st.offdiagonal_fill, "$tag: fill")
    chk(isapprox(row.speedup, row.diffusion.seconds/row.direct.seconds; rtol=1e-9), "$tag: speedup")
    chk(row.diffusion.converged && row.direct.converged, "$tag: not converged")
    chk(row.oracle_residual < 1e-9, "$tag: oracle $(row.oracle_residual)")
    chk(row.direct.ticks <= row.diffusion.ticks, "$tag: direct more ticks")
    chk(row.direct.peak_command > 0 && row.diffusion.peak_command > 0, "$tag: no effort")
    pd = predict(2r.safety/row.condition, r.tolerance)
    chk(0.2*pd <= row.diffusion.ticks <= 1.3*pd, "$tag: ticks $(row.diffusion.ticks) vs analytic $pd")
end
for (nm,sp) in ((:chain,32),(:grid,4),(:twoclique,16),(:expander,32),(:star,32),(:complete,24))
    s = coordination_scenario(nm; size_parameter=sp); pl = direct_plan(s)
    p = CB._target_state(s); q = harmonic_reference(pl,p)
    eta = zeros(length(q)); diffusion_residual!(eta,s,q,p)
    chk(norm(eta,Inf)<1e-9, "$(s.label): eta(q*)")
    chk(norm(s.H*q - s.Bmat*vec(p),Inf)<1e-9, "$(s.label): Hq*-Bp")
    chk(norm(q - CB.harmonic_reference_oracle(s,p),Inf)<1e-9, "$(s.label): vs oracle")
    v = randn(size(s.H,1))
    chk(CB._permute!(similar(v),v,pl.perm) == pl.factor.P' \ v, "$(s.label): permute")
    chk(CB._invpermute!(similar(v),v,pl.perm) == pl.factor.P \ v, "$(s.label): invpermute")
end
println(isempty(fails) ? "ALL CHECKS PASS" : "FAILURES ($(length(fails))):")
for f in unique(fails); println("  ",f); end
