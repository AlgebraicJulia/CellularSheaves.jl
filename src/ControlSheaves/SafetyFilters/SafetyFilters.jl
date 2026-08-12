"""
Local safety-and-stability execution filters for the layered control architecture.

The high-level sheaf solve supplies each agent with a reference `r_i` (a slot of the
harmonic extension) and a nominal command `u_nom`. This module implements the local,
one-way execution map of the ISS cascade: it minimally perturbs `u_nom` so that the agent
makes provable progress toward `r_i` while never violating a collision barrier or its
actuator limits.

For a control-affine agent

```math
\\dot{x}_i = f_i(x_i) + g_i(x_i) u_i + d_i(t, x_i), \\qquad \\|d_i\\| \\le \\Delta_i,
```

with tracking error ``e_i = x_i - r_i`` and ISS-CLF ``V_i = \\tfrac{1}{2} e_i^\\top M e_i``,
the filter solves the second-order cone program

```math
\\begin{aligned}
\\min_{u_i, \\delta}\\quad & \\tfrac{1}{2}\\|u_i - u_{\\text{nom}}\\|^2 + p\\,\\delta^2\\\\
\\text{s.t.}\\quad
& (g_i^\\top M e_i)^\\top u_i \\le -c_3 V_i - e_i^\\top M f_i + e_i^\\top M \\dot{r}_i - \\|M e_i\\| \\Delta_i + \\delta\\\\
& \\text{one hard row per neighbour and obstacle enforcing } \\dot{h} \\ge -\\gamma h\\\\
& \\|W u_i\\|_2 \\le u_{\\max}.
\\end{aligned}
```

# Structure

The three constraint families are independent and can be built, tested, and composed
separately:

| family | type | contributes |
|:--|:--|:--|
| stability | [`StabilityTerm`](@ref) | one relaxable row, reports ``V`` |
| safety | [`SafetyTerm`](@ref) over an [`AbstractBarrier`](@ref) | one hard row per target, reports ``h`` |
| actuator | [`ActuatorTerm`](@ref) | the cone ``\\|W u\\|_2 \\le u_{\\max}`` |

A term reads a [`FilterContext`](@ref) and emits [`LinearRow`](@ref)s through
[`constraints`](@ref); it never inspects another term. [`safety_filter`](@ref) concatenates
them and [`filter_program`](@ref) encodes whatever it is given, so **correctness composes**:
if each term's rows are right and the encoder assembles any row list correctly — which the
test suite establishes against an exact active-set enumeration — the composed program is
right by construction. **Feasibility does not compose**: hard barrier rows and a hard cone
can have empty intersection, which is why the stability row is the relaxable one and why a
filter can legitimately return no command.

[`CBFCLFParams`](@ref) bundles the three families for the common case.

# Notation

| symbol | meaning | code |
|:--|:--|:--|
| ``u_i`` | control input | `u_nom` on input, `result.command` on output |
| ``x_i``, ``r_i``, ``\\dot{r}_i`` | state, reference, reference velocity | `ctx.x`, `ctx.ref`, `ctx.ref_velocity` |
| ``f_i``, ``g_i`` | drift and input map | `model.drift`, `model.input` |
| ``P``, ``P_v`` | configuration and velocity selectors | `model.position`, `model.velocity` |
| ``M``, ``c_3``, ``p``, ``\\delta`` | CLF metric, rate, penalty, relaxation | `StabilityTerm` fields, `result.slack` |
| ``\\gamma``, ``d_{\\text{safe}}`` | class-K gain, separation | `SafetyTerm.gamma`, `barrier.d_safe` |
| ``\\Delta``, ``\\bar{\\Delta}`` | CLF and barrier disturbance bounds | the `disturbance_bound` of each term |
| ``W``, ``u_{\\max}`` | actuator weight and bound | `ActuatorTerm.metric`, `ActuatorTerm.bound` |
| ``h_{ij}`` | barrier value | `result.cbf_margin` |

Rows are reconstructed from these formulas and checked against the filter's own output in
the test suite, so the correspondence is executable rather than asserted.

The program is solved with `CellularSheaves.IPM` as a conic quadratic program in the primal
form ``\\min \\tfrac{1}{2} p^\\top Q p + c^\\top p`` subject to ``Bp = g``, ``p \\in K``.
"""
module SafetyFilters

using LinearAlgebra
using ArgCheck: @argcheck

using ...IPM: IPMProblem, IPMSettings, UzawaSettings, solve, allocblockdiag,
              AbstractCone, CofreeCone, PositiveCone, SecondOrderCone,
              OPTIMAL, NEAR_OPTIMAL
using ...BlockSparseArrays: blocksparse, block, colrange
using ..AgentControllers: AbstractAgentDynamics, AbstractAgentController, AgentState,
                          LQRController, continuous_matrices, position_indices,
                          velocity_indices
using ..Tikhonov: tikhonov_step!, JointTikhonovFilter

import ..AgentControllers: step_agent!

export ControlAffineModel, double_integrator_model, RelativeDegreeError,
       laplacian_block_metric, lyapunov_metric, lyapunov_rate,
       LinearRow, ConeRow, filter_program,
       AbstractBarrier, DistanceBarrier, BrakingBarrier, barrier_row, barrier_value,
       AbstractFilterTerm, FilterContext, constraints,
       StabilityTerm, DiscreteStabilityTerm, SafetyTerm, ActuatorTerm,
       SafetyFilterResult, safety_filter, safety_filter_all, safety_filter_settings,
       CBFCLFParams, filter_terms, CBFCLFController, barrier_values,
       FilterRollout, barrier_window, saturation_window, relaxation_window, animate_filter_rollout

"""
    safety_filter_settings(; itmax=200, raug=1e3)

Default `IPMSettings` for the local filters. Active barrier sets need more central-path
iterations than an unconstrained tracking solve, hence the explicit budget.

The Uzawa augmentation `raug` is deliberately moderate: a large one makes the KKT system
fail numerically whenever the CLF relaxation is driven positive, which is exactly the
saturated regime the relaxation exists to handle.
"""
safety_filter_settings(; itmax::Int = 200, raug::Float64 = 1e3) =
    IPMSettings{Float64}(kkt = UzawaSettings{Float64}(raug = raug), itmax = itmax)

include("model.jl")
include("program.jl")
include("barriers.jl")
include("terms.jl")
include("filter.jl")
include("compat.jl")
include("controller.jl")
include("rollout.jl")

end # module SafetyFilters
