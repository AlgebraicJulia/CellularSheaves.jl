# Needs PATHPACK (AlgebraicJulia/PathPack.jl, unregistered: Pkg.develop a checkout) in the environment;
# "cuda" as an argument also solves on the GPU. Feasibility: PathPack's semiring LU as the policy-evaluation solver. Policy
# evaluation A_u V = b, A_u = D - N (an M-matrix), is the fixpoint V = P V + c with
# P = D⁻¹N substochastic and c = D⁻¹b, so V = P* c, the PlusProd closure.
using CellularSheaves, CellularSheaves.ControlSheaves.DoubleIntegratorHJB, CellularSheaves.ControlSheaves.MechanicalHJB
using PATHPACK, PATHPACK.CPU
using PATHPACK.CPU: PlusProd, mlu
using SparseArrays, LinearAlgebra, Printf
const DIH = CellularSheaves.ControlSheaves.DoubleIntegratorHJB
const GS = CellularSheaves.NetworkSheaves.GridSchwarz
BLAS.set_num_threads(1)
gpu = "cuda" in ARGS
gpu && @eval using CUDA
gpu && @eval (println("device: ", CUDA.name(CUDA.device())); flush(stdout))     # the watchdog waits for this

# The sparse matrix of the implicit stencil (periodic wrap), its right-hand side,
# and the grid-graph pattern: every axis neighbour stored, zeros included.
function assembled(prob, U)
    g = DIH._state_grid(prob); n = Tuple(g.points); D = length(n); periodic = DIH._periodic(prob)
    layout = GS.BoxLayout(n, ntuple(_ -> 1, D), 0; periodic)
    kp = DIH._KernelProblem(prob, ntuple(_ -> 0, D))
    op = GS.box_operator(layout, DIH._UpwindStencil(U, kp))
    b = GS.grid_zeros(op)
    DIH._launch!(DIH._grid_rhs_kernel!, op, b, U, kp, op.origin, ntuple(_ -> 0, D))
    L = LinearIndices(n)
    I, J, V = Int[], Int[], Float64[]
    diag = zeros(length(L))
    for C in CartesianIndices(n)
        c0, cm, cp = GS._coefficients(op.stencil, C)
        diag[L[C]] = c0
        for j in 1:D, (side, c) in ((-1, cm[j]), (1, cp[j]))
            k = C[j] + side
            if !(1 <= k <= n[j])
                periodic[j] || continue
                k = mod1(k, n[j])
            end
            push!(I, L[C]); push!(J, L[Base.setindex(Tuple(C), k, j)...]); push!(V, c)   # c may be 0: pattern
        end
    end
    N = sparse(I, J, V, length(L), length(L))
    return N, diag, vec(GS.interior(b, op))
end

function check(label, prob, U)
    N, d, b = assembled(prob, U)
    A = spdiagm(d) - N
    t_ref = @elapsed Vref = A \ b                                   # UMFPACK, the reference
    P = spdiagm(1 ./ d) * N                                         # keeps the stored zeros
    c = b ./ d
    t_f = @elapsed F = mlu(PlusProd(), P)
    t_s = @elapsed V = lmul!(F, copy(c))
    @printf("%-26s n=%8d nnz=%9d  UMFPACK %6.2f s | PATHPACK CPU factor %6.2f s solve %.3f s  rel.err %.1e\n",
        label, size(A, 1), nnz(P), t_ref, t_f, t_s, norm(V - Vref, Inf) / norm(Vref, Inf))
    if gpu
        Ext = PATHPACK.GPU.extension()
        try
            Pt = sparse(transpose(P))                                   # rows ⇄ right multiplication: Vᵀ = cᵀ (Pᵀ)*
            Ft = mlu(PlusProd(), Pt)
            G = Ext.GPUSLU(Ft)
            B = CuArray(reshape(c, 1, :))
            Ext.rmul_gpu!(B, G); CUDA.synchronize()
            B = CuArray(reshape(c, 1, :))
            t_g = @elapsed (Ext.rmul_gpu!(B, G); CUDA.synchronize())
            Vg = vec(Array(B))
            @printf("%-26s   GPU solve %.4f s  rel.err %.1e\n", "", t_g, norm(Vg - Vref, Inf) / norm(Vref, Inf))
        catch e
            println("   GPU: ", sprint(showerror, e)[1:min(end, 300)])
        end
    end
    flush(stdout)
end

for n in (9, 13, 17)
    prob = HJBProblem(StateGrid(fill(-2.0, 4), fill(2.0, 4), fill(n, 4)); control_bound = 1.0, constraint = :disc)
    U = reshape(permutedims(DIH._lqr_controls(prob)), n, n, n, n, 2)
    check("double integrator $(n)^4", prob, U)
end
for (na, np) in ((12, 9), (16, 13), (20, 17))
    prob = MechanicalHJBProblem(TwoLinkArm(); angle_points = na, momentum_points = np, torque_bound = (6.0, 3.0))
    sol = solve(prob, GridPolicyIteration(forcing = 0.1))
    U = reshape(permutedims(sol.controls), prob.grid.points..., 2)
    check("arm $(na)^2x$(np)^2 (optimal u)", prob, U)
end
