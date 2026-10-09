# Implicit (matrix-free) policy iteration on boxes of the grid: right-hand side, policy
# improvement and the initial policy as KernelAbstractions kernels over a
# GridOperator (GridSchwarz), so one code path runs on CPU threads, on GPUs
# and, box by box, on the ranks of a distributed solve. Included in
# DoubleIntegratorHJB.

# Everything a kernel needs to know about the problem and the box, as plain
# bits. Point I of the box (1-based) is the grid point offset .+ I.
struct _KernelProblem{D,M}
    lower::NTuple{D,Float64}
    spacing::NTuple{D,Float64}
    points::NTuple{D,Int}
    offset::NTuple{D,Int}
    position_weight::Float64
    velocity_weight::Float64
    control_weight::Float64
    discount::Float64
    control_bound::Float64
    disc::Bool
    riccati::NTuple{M,Float64}
end

function _KernelProblem(prob::HJBProblem, offset::NTuple{D,Int}) where {D}
    g = prob.grid
    @argcheck ndims(g) == D
    return _KernelProblem{D,D * D}(Tuple(g.lower), Tuple(g.spacing), Tuple(g.points), offset,
        prob.position_weight, prob.velocity_weight, prob.control_weight, prob.discount, prob.control_bound,
        prob.constraint === :disc, Tuple(vec(prob.riccati)))
end

@inline _coordinates(kp::_KernelProblem{D}, I) where {D} =
    ntuple(k -> kp.lower[k] + (kp.offset[k] + I[k] - 1) * kp.spacing[k], Val(D))

@inline function _riccati(kp::_KernelProblem{D}, x) where {D}
    acc = 0.0
    for i in 1:D, j in 1:D
        acc += x[i] * kp.riccati[i + D * (j - 1)] * x[j]
    end
    return acc / 2
end

@inline function _cost(kp::_KernelProblem{D}, x, u) where {D}
    d = D ÷ 2
    q = 0.0
    v = 0.0
    for k in 1:d
        q += x[k]^2
        v += x[d + k]^2
    end
    uu = 0.0
    for k in 1:d
        uu += u[k]^2
    end
    return (kp.position_weight * q + kp.velocity_weight * v + kp.control_weight * uu) / 2
end

@inline _moved(x::NTuple{D}, j, step) where {D} = ntuple(k -> k == j ? x[k] + step : x[k], Val(D))

@inline _controls(U, I, ::Val{d}) where {d} = ntuple(k -> U[I, k], Val(d))

@inline _drift(x, u, j, d) = j <= d ? x[d + j] : u[j - d]

# The upwind stencil of a policy, evaluated on the fly from the controls and
# the grid coordinates: the same discretization as `_assemble`, never stored.
# A coefficient towards a neighbour outside the grid is zero; its Riccati
# boundary value is moved into the right-hand side (`_grid_rhs_kernel!`).
struct _UpwindStencil{D,A,K} <: AbstractStencil{D}
    U::A
    kp::K
end

_UpwindStencil(U::AbstractArray, kp::_KernelProblem{D}) where {D} = _UpwindStencil{D,typeof(U),typeof(kp)}(U, kp)

function Adapt.adapt_structure(to, s::_UpwindStencil{D}) where {D}
    U = Adapt.adapt(to, s.U)
    return _UpwindStencil{D,typeof(U),typeof(s.kp)}(U, s.kp)
end
KernelAbstractions.get_backend(s::_UpwindStencil) = get_backend(s.U)

@inline function GridSchwarz._coefficients(s::_UpwindStencil{D}, I) where {D}
    kp = s.kp
    d = D ÷ 2
    x = _coordinates(kp, I)
    u = _controls(s.U, I, Val(D ÷ 2))
    c0 = kp.discount
    for j in 1:D
        c0 += abs(_drift(x, u, j, d)) / kp.spacing[j]
    end
    cm = ntuple(Val(D)) do j
        f = _drift(x, u, j, d)
        f < 0 && kp.offset[j] + I[j] > 1 ? -f / kp.spacing[j] : 0.0
    end
    cp = ntuple(Val(D)) do j
        f = _drift(x, u, j, d)
        f > 0 && kp.offset[j] + I[j] < kp.points[j] ? f / kp.spacing[j] : 0.0
    end
    return c0, cm, cp
end

# The right-hand side of the policy evaluation: running cost plus the Riccati
# values of upwind neighbours outside the grid.
@kernel function _grid_rhs_kernel!(b, @Const(U), kp::_KernelProblem{D}, g::Int) where {D}
    I = @index(Global, Cartesian)
    d = D ÷ 2
    x = _coordinates(kp, I)
    u = _controls(U, I, Val(D ÷ 2))
    rhs = _cost(kp, x, u)
    for j in 1:D
        fj = _drift(x, u, j, d)
        G = kp.offset[j] + I[j]
        if fj > 0 && G == kp.points[j]
            rhs += fj / kp.spacing[j] * _riccati(kp, _moved(x, j, kp.spacing[j]))
        elseif fj < 0 && G == 1
            rhs -= fj / kp.spacing[j] * _riccati(kp, _moved(x, j, -kp.spacing[j]))
        end
    end
    b[I + _shift(g, Val(D))] = rhs
end

# Howard's improvement: minimize the discrete Hamiltonian pointwise (see
# `_minimize_hamiltonian!`), reading V on the ghost layer for neighbours in
# other boxes and the Riccati value outside the grid.
@kernel function _grid_improve_kernel!(U, @Const(V), kp::_KernelProblem{D}, g::Int) where {D}
    I = @index(Global, Cartesian)
    d = D ÷ 2
    x = _coordinates(kp, I)
    J = I + _shift(g, Val(D))
    v0 = V[J]
    Dp = ntuple(Val(D ÷ 2)) do k
        j = d + k
        vp = kp.offset[j] + I[j] < kp.points[j] ? V[J + _unit(j, Val(D))] : _riccati(kp, _moved(x, j, kp.spacing[j]))
        (vp - v0) / kp.spacing[j]
    end
    Dm = ntuple(Val(D ÷ 2)) do k
        j = d + k
        vm = kp.offset[j] + I[j] > 1 ? V[J - _unit(j, Val(D))] : _riccati(kp, _moved(x, j, -kp.spacing[j]))
        (v0 - vm) / kp.spacing[j]
    end
    u = _argmin_hamiltonian(kp, Dp, Dm)
    for k in 1:d
        U[I, k] = u[k]
    end
end

@inline function _argmin_hamiltonian(kp::_KernelProblem, Dp::NTuple{d}, Dm::NTuple{d}) where {d}
    r, ū = kp.control_weight, kp.control_bound
    if !kp.disc || d == 1
        return ntuple(Val(d)) do k
            up = clamp(-Dp[k] / r, 0.0, ū)
            um = clamp(-Dm[k] / r, -ū, 0.0)
            r * up^2 / 2 + up * Dp[k] <= r * um^2 / 2 + um * Dm[k] ? up : um
        end
    end
    best = Inf
    u1 = u2 = 0.0
    for s1 in (1, -1), s2 in (1, -1)
        g1 = s1 > 0 ? Dp[1] : Dm[1]
        g2 = s2 > 0 ? Dp[2] : Dm[2]
        z1 = s1 > 0 ? max(-g1 / r, 0.0) : min(-g1 / r, 0.0)
        z2 = s2 > 0 ? max(-g2 / r, 0.0) : min(-g2 / r, 0.0)
        nz = sqrt(z1^2 + z2^2)
        if nz > ū
            z1 *= ū / nz
            z2 *= ū / nz
        end
        H = r * (z1^2 + z2^2) / 2 + z1 * g1 + z2 * g2
        if H < best
            best = H
            u1, u2 = z1, z2
        end
    end
    return ntuple(k -> k == 1 ? u1 : u2, Val(d))
end

# The initial policy (clipped LQR feedback) and value (Riccati) on the box.
@kernel function _grid_initialize_kernel!(U, V, kp::_KernelProblem{D}, g::Int) where {D}
    I = @index(Global, Cartesian)
    d = D ÷ 2
    x = _coordinates(kp, I)
    V[I + _shift(g, Val(D))] = _riccati(kp, x)
    u = ntuple(Val(D ÷ 2)) do k
        acc = 0.0
        for j in 1:D
            acc -= kp.riccati[(d + k) + D * (j - 1)] * x[j]
        end
        acc / kp.control_weight
    end
    ū = kp.control_bound
    if kp.disc && d == 2
        nu = sqrt(u[1]^2 + u[2]^2)
        scale = nu > ū ? ū / nu : 1.0
        u = ntuple(k -> u[k] * scale, Val(D ÷ 2))
    else
        u = ntuple(k -> clamp(u[k], -ū, ū), Val(D ÷ 2))
    end
    for k in 1:d
        U[I, k] = u[k]
    end
end

function _launch!(kernel, op::GridOperator, args...)
    backend = get_backend(op)
    kernel(backend)(args...; ndrange=size(op))
    KernelAbstractions.synchronize(backend)
end

"""
    GridPolicyIteration(; backend = KernelAbstractions.CPU(), preconditioner = :red_black_sgs,
                        tol = 1e-8, maxiter = 50, linear_tol = 1e-10, linear_maxiter = 2000)

Howard's policy iteration (as [`PolicyIteration`](@ref)) without assembling
sparse matrices or storing stencil coefficients: the policy evaluation operator is
an implicit `GridOperator` whose upwind stencil is computed inside each kernel
from the controls and the grid coordinates, and right-hand side, policy
evaluation and improvement are
KernelAbstractions kernels on `backend` (`CPU()` runs them on Julia threads;
a GPU backend such as `CUDABackend()` from CUDA.jl runs them on the device).
Each evaluation is solved with BiCGStab (`grid_bicgstab!`) right-preconditioned
by one red–black symmetric Gauss–Seidel sweep (`red_black_sgs!`,
`preconditioner = :red_black_sgs`) or unpreconditioned (`:none`), warm-started
from the previous value function. Gives the same discrete solution as
[`PolicyIteration`](@ref), whose `KrylovPolicyEvaluation(method = :bicgstab)`
uses the identical preconditioner on the assembled matrices.
"""
Base.@kwdef struct GridPolicyIteration{B}
    backend::B = KernelAbstractions.CPU()
    preconditioner::Symbol = :red_black_sgs
    tol::Float64 = 1e-8
    maxiter::Int = 50
    linear_tol::Float64 = 1e-10
    linear_maxiter::Int = 2000
end

function CommonSolve.solve(prob::HJBProblem, alg::GridPolicyIteration)
    @argcheck alg.maxiter >= 1 && alg.tol > 0
    @argcheck alg.preconditioner in (:red_black_sgs, :none) "preconditioner must be :red_black_sgs or :none"
    g = prob.grid
    D, d = ndims(g), prob.axes
    n = Tuple(g.points)
    kp = _KernelProblem(prob, ntuple(_ -> 0, D))
    backend = alg.backend
    U = KernelAbstractions.zeros(backend, Float64, n..., d)
    op = GridOperator(_UpwindStencil(U, kp), n; ghost=1)
    V, Vnew, b = grid_zeros(op), grid_zeros(op), grid_zeros(op)
    ws = GridWorkspace(op)
    precondition! = alg.preconditioner === :none ? copyto! : (y, v) -> red_black_sgs!(y, op, v)
    _launch!(_grid_initialize_kernel!, op, U, V, kp, op.ghost)
    changes, linear_iterations = Float64[], Int[]
    converged = false
    t_assembly = t_linear = t_improvement = 0.0
    for _ in 1:alg.maxiter
        t_assembly += @elapsed _launch!(_grid_rhs_kernel!, op, b, U, kp, op.ghost)
        t_linear += @elapsed begin
            copyto!(Vnew, V)
            its, ok, _ = grid_bicgstab!(Vnew, op, precondition!, b, ws; tol=alg.linear_tol, maxiter=alg.linear_maxiter)
        end
        ok || @warn "policy evaluation did not converge"
        push!(changes, grid_reduce((a, c) -> abs(a - c), max, Vnew, V, op, ws.partial))
        push!(linear_iterations, its)
        V, Vnew = Vnew, V
        t_improvement += @elapsed _launch!(_grid_improve_kernel!, op, U, V, kp, op.ghost)
        if length(changes) > 1 && changes[end] <= alg.tol * max(1.0, grid_reduce((a, c) -> abs(a), max, V, V, op, ws.partial))
            converged = true
            break
        end
    end
    values = vec(Array(interior(V, op)))
    controls = permutedims(reshape(Array(U), length(g), d))
    seconds = (assembly = t_assembly, setup = 0.0, linear = t_linear, improvement = t_improvement)
    return HJBSolution(prob, values, controls, length(changes), changes, linear_iterations, converged, seconds)
end
