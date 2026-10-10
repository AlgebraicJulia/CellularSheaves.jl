# Rediscretized against Galerkin (T A E) coarse operators of the HJB policy
# evaluations, under the optimal policy, per kind of coefficient: the diagonal,
# and the couplings across faces of position and of momentum (or velocity)
# dimensions. Reports the fraction of coarse points where they differ (by more
# than 1e-10 relative to the diagonal) and the median and largest relative
# difference, at two resolutions, to see the first-order agreement.
#
#   julia --project=<env with PATHPACK> docs/scripts/hjb_coarse_operators.jl
using CellularSheaves
using CellularSheaves.ControlSheaves.DoubleIntegratorHJB
using CellularSheaves.ControlSheaves.MechanicalHJB
using CellularSheaves.ControlSheaves.DoubleIntegratorHJB: solve
using KernelAbstractions
using PATHPACK
using Printf
using Statistics

const DIH = CellularSheaves.ControlSheaves.DoubleIntegratorHJB
const GS = CellularSheaves.NetworkSheaves.GridSchwarz

function coarse_operators(prob, U)
    n = Tuple(DIH._state_grid(prob).points)
    D = length(n)
    layout = GS.BoxLayout(n, ntuple(_ -> 1, D), 0; periodic = DIH._periodic(prob))
    op = GS.box_operator(layout, DIH._UpwindStencil(U, DIH._KernelProblem(prob, ntuple(_ -> 0, D))))
    nc = coarse_points(n)
    galerkin = zeros(nc..., 2D + 1)
    galerkin_coefficients!(galerkin, op)
    Uc = zeros(nc..., size(U, D + 1))
    DIH._restrict_controls_kernel!(KernelAbstractions.CPU())(Uc, U, Val(D), Val(size(U, D + 1)); ndrange = nc)
    opc = GS.box_operator(coarse_layout(layout), DIH._UpwindStencil(Uc, DIH._KernelProblem(DIH._coarsen(prob), ntuple(_ -> 0, D))))
    redisc = zeros(nc..., 2D + 1)
    stencil_coefficients!(redisc, opc)
    return redisc, galerkin
end

function compare(label, prob)
    sol = solve(prob, GridPolicyIteration(multigrid = Multigrid(levels = 2)))
    n = DIH._state_grid(prob).points
    D = length(n)
    U = reshape(permutedims(sol.controls), n..., size(sol.controls, 1))
    redisc, galerkin = coarse_operators(prob, U)
    R, G = reshape(redisc, :, 2D + 1), reshape(galerkin, :, 2D + 1)
    scale = abs.(R[:, 1])
    d = D ÷ 2
    kinds = [("diagonal", [1]), ("position faces", vcat(1 .+ (1:d), 1 + D .+ (1:d))),
             ("momentum faces", vcat(1 + d .+ (1:d), 1 + D + d .+ (1:d)))]
    for (kind, columns) in kinds
        rel = vec(maximum(abs.(R[:, columns] .- G[:, columns]) ./ scale; dims = 2))
        @printf("%-28s %-15s differ at %5.1f%% of coarse points; relative difference median %.2e, max %.2e\n",
            label, kind, 100 * count(>(1e-10), rel) / length(rel), median(rel), maximum(rel))
    end
    flush(stdout)
end

for n in (12, 24)
    compare("double integrator $(n)^4", HJBProblem(StateGrid(fill(-2.0, 4), fill(2.0, 4), fill(n, 4)); control_bound = 1.0,
        constraint = :disc))
end
for (na, np) in ((16, 12), (32, 24))
    compare("arm $(na)^2x$(np)^2", MechanicalHJBProblem(TwoLinkArm(); angle_points = na, momentum_points = np,
        torque_bound = (6.0, 3.0)))
end
