# Data for the policy-iteration demo page: the planar double integrator with a
# disc thrust bound, solved on an n⁴ grid by GridPolicyIteration, recording the
# value function and the policy after every iteration (the callback).
#
#   sysimage/julia.sh docs/scripts/hjb_policy_iteration_demo.jl [n] [output]
#
# Writes a little-endian binary file (default docs/figures/hjb/policy_iteration.bin,
# next to the page docs/figures/hjb/policy_iteration.html, which reads it; the
# file is not committed):
#   Int32  magic 0x48_4A_42_31 ("HJB1"), n, D = 4, d = 2, K snapshots
#   Float32 lower[4], spacing[4], discount, control bound, control weight,
#           position weight, velocity weight
#   K × { Int32 iteration, Float32 change, Float16 values[N], Float16 u₁[N], Float16 u₂[N] }
# with N = n⁴ in column-major grid order (q₁ fastest, then q₂, v₁, v₂). Half
# precision keeps the page's download small; the plot needs 3 digits at most.
using CellularSheaves
using CellularSheaves.ControlSheaves.DoubleIntegratorHJB
using Printf

n = length(ARGS) >= 1 ? parse(Int, ARGS[1]) : 17
out = length(ARGS) >= 2 ? ARGS[2] : joinpath(@__DIR__, "..", "figures", "hjb", "policy_iteration.bin")
prob = HJBProblem(StateGrid(fill(-2.0, 4), fill(2.0, 4), fill(n, 4)); control_bound = 1.0, constraint = :disc)

snapshots = []
t = @elapsed sol = solve(prob, GridPolicyIteration(callback = s -> push!(snapshots, s)))
@printf("n=%d: %d policy iterations in %.2f s, converged %s\n", n, sol.iterations, t, sol.converged)
for s in snapshots
    @printf("  iteration %2d  value change %.3e  V range [%.3f, %.3f]\n", s.iteration, s.change, extrema(s.values)...)
end

mkpath(dirname(out))
open(out, "w") do io
    write(io, htol.(Int32[0x484A4231, n, 4, prob.axes, length(snapshots)]))
    g = prob.grid
    write(io, htol.(Float32[g.lower; g.spacing; prob.discount; prob.control_bound; prob.control_weight;
                            prob.position_weight; prob.velocity_weight]))
    for s in snapshots
        write(io, htol(Int32(s.iteration)), htol(Float32(s.change)))
        write(io, htol.(Float16.(s.values)))
        for k in 1:prob.axes
            write(io, htol.(Float16.(s.controls[k, :])))
        end
    end
end
println("wrote ", out, " (", round(filesize(out) / 2^20; digits = 1), " MiB)")
