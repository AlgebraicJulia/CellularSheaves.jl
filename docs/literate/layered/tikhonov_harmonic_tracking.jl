# # Tikhonov-Filtered Harmonic Extension
#
# This executable example separates the harmonic solve from its normalized
# dynamic realization:
#
#     H * qstar(t) = -L_IB * p(t)
#     epsilon * xdot = -x + qstar(t).
#
# The feature guide develops the full 12-agent pipeline. This compact example
# is intended as a runnable API check.

get!(ENV, "GKSwstype", "100")

using CellularSheaves
using CellularSheaves.ControlSheaves.Tikhonov
using LinearAlgebra, Plots, Statistics
import CellularSheaves.NetworkSheaves.EuclideanSheaves:
    _harmonic_extension_restricted_laplacian

# ## A fixed sheaf with moving boundary data

const NA, D, TARGET = 4, 2, 5
const I2 = Matrix{Float64}(I, D, D)

sheaf = EuclideanSheaf{Float64}(fill(D, NA + 1))
for (i, j) in ((1, 2), (2, 3), (3, 4))
    add_sheaf_edge!(sheaf, i, j, I2, I2)
end
for i in (1, 4)
    add_sheaf_edge!(sheaf, i, TARGET, I2, I2)
end

target(t) = [1.2cos(0.7t), 0.9sin(0.7t)]
target_rate(t) = [-0.84sin(0.7t), 0.63cos(0.7t)]
boundary0 = Dict(TARGET => target(0.0))
_, _, Hraw, LIBraw = _harmonic_extension_restricted_laplacian(sheaf, boundary0)
H, LIB = Matrix(Hraw), Matrix(LIBraw)
qstar(t) = tikhonov_equilibrium(H, -LIB * target(t))
qstar_rate(t) = tikhonov_reference_rate(H, -LIB * target_rate(t))

# ## Normalized filtering and feedforward

function rollout(epsilon; feedforward=false)
    times = collect(0.0:0.02:12.0)
    filter = TikhonovFilter(qstar(0.0); epsilon)
    error = zeros(length(times))
    for k in eachindex(times)
        error[k] = norm(filter.x - qstar(times[k]))
        k == length(times) && continue
        q0, q1 = qstar(times[k]), qstar(times[k + 1])
        if feedforward
            u0 = tikhonov_feedforward_reference(q0, qstar_rate(times[k]), epsilon)
            u1 = tikhonov_feedforward_reference(q1, qstar_rate(times[k + 1]), epsilon)
            tikhonov_step!(filter, u0, u1, times[k + 1] - times[k])
        else
            tikhonov_step!(filter, q0, q1, times[k + 1] - times[k])
        end
    end
    (; times, error)
end

baseline = rollout(0.2)
corrected = rollout(0.2; feedforward=true)
tail = baseline.times .>= 6.0
@assert sqrt(mean(abs2, corrected.error[tail])) <
        1e-3 * sqrt(mean(abs2, baseline.error[tail]))

plot(baseline.times, baseline.error;
    xlabel="time (s)", ylabel="planner error", label="uncompensated",
    linewidth=2, color=:crimson, title="Normalized harmonic-reference filter")
plot!(corrected.times, corrected.error;
    label="analytic feedforward", linewidth=2, color=:seagreen)
