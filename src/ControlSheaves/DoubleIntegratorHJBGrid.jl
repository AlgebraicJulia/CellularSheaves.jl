# Implicit (matrix-free) policy iteration on boxes of the grid: right-hand side, policy
# improvement and the initial policy as KernelAbstractions kernels over a
# GridOperator (GridSchwarz), so one code path runs on CPU threads, on GPUs
# and, box by box, on the ranks of a distributed solve. Included in
# DoubleIntegratorHJB.
#
# The kernels are written for any controlled mechanical system whose state is
# ordered (positions, then momenta or velocities), d = D ÷ 2 of each, with one
# control per momentum. A problem type plugs in by providing
#
#   _state_grid(prob) -> StateGrid        _periodic(prob) -> NTuple{D,Bool}
#   _control_count(prob) -> d             _discount(prob) -> ρ
#   _kernel_model(prob) -> m              (a plain-bits struct, so it runs on GPUs)
#
# and methods of the model hooks for typeof(m):
#
#   _drifts(m, x, u) -> NTuple{D}         the vector field f(x, u)
#   _cost(m, x, u)                        the running cost ℓ(x, u)
#   _dirichlet(m) -> Bool                 at the edge of a non-periodic dimension,
#   _boundary_value(m, x)                   true: V = _boundary_value outside the grid;
#                                           false: the outward flux is dropped (a
#                                           reflecting boundary, a state constraint)
#   _argmin_hamiltonian(m, x, Dp, Dm)     the control minimizing the upwind discrete
#                                           Hamiltonian, given the forward and backward
#                                           differences of V along the d momenta
#   _initial_value(m, x), _initial_control(m, x)   the starting guess

# Everything a kernel needs to know about the problem and the box, as plain
# bits. Point I of the box (1-based) is the grid point offset .+ I.
struct _KernelProblem{D,S}
    lower::NTuple{D,Float64}
    spacing::NTuple{D,Float64}
    inverse_spacing::NTuple{D,Float64}
    points::NTuple{D,Int}
    offset::NTuple{D,Int}
    periodic::NTuple{D,Bool}
    discount::Float64
    model::S
end

function _KernelProblem(prob, offset::NTuple{D,Int}) where {D}
    g = _state_grid(prob)
    @argcheck ndims(g) == D
    m = _kernel_model(prob)
    return _KernelProblem{D,typeof(m)}(Tuple(g.lower), Tuple(g.spacing), Tuple(1 ./ g.spacing), Tuple(g.points),
        offset, _periodic(prob), _discount(prob), m)
end

@inline _coordinates(kp::_KernelProblem{D}, I) where {D} =
    ntuple(k -> kp.lower[k] + (kp.offset[k] + I[k] - 1) * kp.spacing[k], Val(D))

# Whether grid point I has a neighbour below / above it in dimension j.
@inline _has_lower(kp::_KernelProblem, I, j) = kp.periodic[j] || kp.offset[j] + I[j] > 1
@inline _has_upper(kp::_KernelProblem, I, j) = kp.periodic[j] || kp.offset[j] + I[j] < kp.points[j]

@inline _moved(x::NTuple{D}, j, step) where {D} = ntuple(k -> k == j ? x[k] + step : x[k], Val(D))

@inline _controls(U, I, ::Val{d}) where {d} = ntuple(k -> U[I, k], Val(d))

# The minimizer over u ∈ [-ū, ū] of the upwind Hamiltonian of one momentum,
#   φ(u) = r u² / 2 + (g + u)⁺ Dp + (g + u)⁻ Dm,
# whose drift g + u is the force g without control plus the control u: on each
# side of the kink u = -g, φ is a quadratic with its own one-sided difference.
@inline function _box_argmin(r, ū, g, Dp, Dm)
    φ(u) = r * u^2 / 2 + max(g + u, 0.0) * Dp + min(g + u, 0.0) * Dm
    up = clamp(-Dp / r, max(-g, -ū), ū)               # drift ≥ 0: the forward difference
    um = clamp(-Dm / r, -ū, min(-g, ū))               # drift ≤ 0: the backward difference
    -g > ū && return um                               # the drift cannot be made nonnegative
    -g < -ū && return up                              # nor nonpositive
    return φ(up) <= φ(um) ? up : um
end

# ---------------------------------------------------------------------------
# the double integrator: q̇ = v, v̇ = u, with the Riccati value outside the grid
# ---------------------------------------------------------------------------

struct _DoubleIntegratorModel{M}
    position_weight::Float64
    velocity_weight::Float64
    control_weight::Float64
    control_bound::Float64
    disc::Bool
    riccati::NTuple{M,Float64}
end

_state_grid(prob::HJBProblem) = prob.grid
_periodic(prob::HJBProblem) = ntuple(_ -> false, ndims(prob.grid))
_control_count(prob::HJBProblem) = prob.axes
_discount(prob::HJBProblem) = prob.discount
_kernel_model(prob::HJBProblem) =
    _DoubleIntegratorModel{length(prob.riccati)}(prob.position_weight, prob.velocity_weight, prob.control_weight,
        prob.control_bound, prob.constraint === :disc, Tuple(vec(prob.riccati)))

@inline function _riccati(m::_DoubleIntegratorModel, x::NTuple{D}) where {D}
    acc = 0.0
    for i in 1:D, j in 1:D
        acc += x[i] * m.riccati[i + D * (j - 1)] * x[j]
    end
    return acc / 2
end

@inline _drifts(::_DoubleIntegratorModel, x::NTuple{D}, u) where {D} =
    ntuple(j -> j <= D ÷ 2 ? x[D ÷ 2 + j] : u[j - D ÷ 2], Val(D))

@inline function _cost(m::_DoubleIntegratorModel, x::NTuple{D}, u) where {D}
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
    return (m.position_weight * q + m.velocity_weight * v + m.control_weight * uu) / 2
end

_dirichlet(::_DoubleIntegratorModel) = true
@inline _boundary_value(m::_DoubleIntegratorModel, x) = _riccati(m, x)
@inline _initial_value(m::_DoubleIntegratorModel, x) = _riccati(m, x)

# The clipped LQR feedback.
@inline function _initial_control(m::_DoubleIntegratorModel, x::NTuple{D}) where {D}
    d = D ÷ 2
    # (No variable a closure captures is reassigned: GPUs can't run the boxed
    # captures that would create.)
    ulqr = ntuple(Val(D ÷ 2)) do k
        acc = 0.0
        for j in 1:D
            acc -= m.riccati[(d + k) + D * (j - 1)] * x[j]
        end
        acc / m.control_weight
    end
    ū = m.control_bound
    if m.disc && d == 2
        nu = sqrt(ulqr[1]^2 + ulqr[2]^2)
        scale = nu > ū ? ū / nu : 1.0
        return ntuple(k -> ulqr[k] * scale, Val(D ÷ 2))
    end
    return ntuple(k -> clamp(ulqr[k], -ū, ū), Val(D ÷ 2))
end

@inline function _argmin_hamiltonian(m::_DoubleIntegratorModel, x, Dp::NTuple{d}, Dm::NTuple{d}) where {d}
    r, ū = m.control_weight, m.control_bound
    if !m.disc || d == 1
        return ntuple(k -> _box_argmin(r, ū, 0.0, Dp[k], Dm[k]), Val(d))
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
    return let a = u1, b = u2                       # u1, u2 are reassigned above: bind them
        ntuple(k -> k == 1 ? a : b, Val(d))
    end
end

# ---------------------------------------------------------------------------
# kernels
# ---------------------------------------------------------------------------

# The upwind stencil of a policy, evaluated on the fly from the controls and
# the grid coordinates: the same discretization as `_assemble`, never stored.
# Towards a neighbour outside the grid the coefficient is zero: with Dirichlet
# data its value moves into the right-hand side (`_grid_rhs_kernel!`) and the
# flux stays on the diagonal; at a reflecting boundary the flux is dropped.
# Across the seam of a periodic dimension the neighbour is in the ghost layer.
# The controls of box point I are U[I + ushift, :] (the controls are stored on
# the overlap-extended box, which the owned box sits inside).
struct _UpwindStencil{D,A,K} <: AbstractStencil{D}
    U::A
    kp::K
    ushift::NTuple{D,Int}
end

_UpwindStencil(U::AbstractArray, kp::_KernelProblem{D}, ushift::NTuple{D,Int}=ntuple(_ -> 0, D)) where {D} =
    _UpwindStencil{D,typeof(U),typeof(kp)}(U, kp, ushift)

function Adapt.adapt_structure(to, s::_UpwindStencil{D}) where {D}
    U = Adapt.adapt(to, s.U)
    return _UpwindStencil{D,typeof(U),typeof(s.kp)}(U, s.kp, s.ushift)
end
KernelAbstractions.get_backend(s::_UpwindStencil) = get_backend(s.U)

Base.@propagate_inbounds function GridSchwarz._coefficients(s::_UpwindStencil{D}, I) where {D}
    kp = s.kp
    x = _coordinates(kp, I)
    u = _controls(s.U, I + CartesianIndex(s.ushift), Val(D ÷ 2))
    f = _drifts(kp.model, x, u)
    a = ntuple(j -> f[j] * kp.inverse_spacing[j], Val(D))                 # signed upwind rates
    dirichlet = _dirichlet(kp.model)
    c0 = kp.discount
    for j in 1:D
        if dirichlet || (a[j] > 0 ? _has_upper(kp, I, j) : _has_lower(kp, I, j))
            c0 += abs(a[j])
        end
    end
    cm = ntuple(j -> a[j] < 0 && _has_lower(kp, I, j) ? -a[j] : 0.0, Val(D))
    cp = ntuple(j -> a[j] > 0 && _has_upper(kp, I, j) ? a[j] : 0.0, Val(D))
    return c0, cm, cp
end

# The right-hand side of the policy evaluation: running cost plus, with
# Dirichlet data, the boundary values of upwind neighbours outside the grid.
@kernel function _grid_rhs_kernel!(b, @Const(U), kp::_KernelProblem{D}, g, ushift) where {D}
    I = @index(Global, Cartesian)
    m = kp.model
    x = _coordinates(kp, I)
    u = _controls(U, I + CartesianIndex(ushift), Val(D ÷ 2))
    rhs = _cost(m, x, u)
    if _dirichlet(m)
        f = _drifts(m, x, u)
        for j in 1:D
            if f[j] > 0 && !_has_upper(kp, I, j)
                rhs += f[j] / kp.spacing[j] * _boundary_value(m, _moved(x, j, kp.spacing[j]))
            elseif f[j] < 0 && !_has_lower(kp, I, j)
                rhs -= f[j] / kp.spacing[j] * _boundary_value(m, _moved(x, j, -kp.spacing[j]))
            end
        end
    end
    b[I + _shift(g, Val(D))] = rhs
end

# The value one step from x along dimension j (step = ±h), seen from grid point
# I: from V (or its ghost layer), the boundary data, or, at a reflecting
# boundary, v0 itself (a zero difference: no flux through the boundary).
@inline function _neighbour_value(kp::_KernelProblem{D}, V, I, J, x, v0, j, upper::Bool) where {D}
    if upper ? _has_upper(kp, I, j) : _has_lower(kp, I, j)
        return upper ? V[J + _unit(j, Val(D))] : V[J - _unit(j, Val(D))]
    end
    _dirichlet(kp.model) || return v0
    return _boundary_value(kp.model, _moved(x, j, upper ? kp.spacing[j] : -kp.spacing[j]))
end

# Howard's improvement: minimize the discrete Hamiltonian pointwise (see
# `_minimize_hamiltonian!`), reading V on the ghost layer for neighbours in
# other boxes or across a periodic seam.
@kernel function _grid_improve_kernel!(U, @Const(V), kp::_KernelProblem{D}, g) where {D}
    I = @index(Global, Cartesian)
    d = D ÷ 2
    x = _coordinates(kp, I)
    J = I + _shift(g, Val(D))
    v0 = V[J]
    Dp = ntuple(k -> (_neighbour_value(kp, V, I, J, x, v0, d + k, true) - v0) / kp.spacing[d + k], Val(D ÷ 2))
    Dm = ntuple(k -> (v0 - _neighbour_value(kp, V, I, J, x, v0, d + k, false)) / kp.spacing[d + k], Val(D ÷ 2))
    u = _argmin_hamiltonian(kp.model, x, Dp, Dm)
    for k in 1:d
        U[I, k] = u[k]
    end
end

# The initial policy and value on the box.
@kernel function _grid_initialize_kernel!(U, V, kp::_KernelProblem{D}, g) where {D}
    I = @index(Global, Cartesian)
    x = _coordinates(kp, I)
    V[I + _shift(g, Val(D))] = _initial_value(kp.model, x)
    u = _initial_control(kp.model, x)
    for k in 1:(D ÷ 2)
        U[I, k] = u[k]
    end
end

function _launch!(kernel, op::GridOperator, args...)
    backend = get_backend(op)
    kernel(backend)(args...; ndrange=size(op))
    KernelAbstractions.synchronize(backend)
end

"""
    GridPolicyIteration(; backend = KernelAbstractions.CPU(), communicator = SerialBoxes(),
                        ranks = nothing, overlap = 1, preconditioner = :red_black_sgs,
                        tol = 1e-8, maxiter = 50, linear_tol = 1e-10, linear_maxiter = 2000,
                        forcing = 0, coarse_blocks = nothing, coarse_correction = :multiplicative,
                        callback = nothing)

Howard's policy iteration (as [`PolicyIteration`](@ref)) without assembling
sparse matrices or storing stencil coefficients: the policy evaluation operator
is an implicit `GridOperator` whose upwind stencil is computed inside each
kernel from the controls and the grid coordinates, and the right-hand side,
policy evaluation and improvement are KernelAbstractions kernels on `backend`
(`CPU()` runs them on Julia threads; a GPU backend such as `CUDABackend()` from
CUDA.jl runs them on the device).

The grid is split into one box per rank of `communicator` (`SerialBoxes()`: a
single box; `mpi_boxes(comm)` with MPI.jl: one box per MPI rank), `ranks`
boxes per dimension (default `balanced_ranks`), each storing its owned points
plus ghost cells. Each evaluation is solved with BiCGStab (`grid_bicgstab!`,
inner products reduced over all ranks), right-preconditioned by

- `:red_black_sgs`: one global red–black symmetric Gauss–Seidel sweep, the
  same operator for any number of ranks (two halo exchanges per sweep);
- `:ras`: restricted additive Schwarz on the boxes extended by `overlap`
  points, with one local red–black sweep per box (one halo exchange per
  application; the boxes are the vertices of the overlap sheaf and the
  exchanges run along its edges). Equal to `:red_black_sgs` on one box;
- `:none`.

With `coarse_blocks` (aggregates per dimension, e.g. `[4, 4, 2, 2]`) the
preconditioner gets a second level: a piecewise-constant coarse space on that
grid of aggregates ([`AggregateCoarseSpace`](@ref CellularSheaves.NetworkSheaves.GridSchwarz.AggregateCoarseSpace)),
independent of the boxes, whose Galerkin operator ``A_0 = R_0 A R_0^\\mathsf{T}``
is computed from the stencil, summed over all ranks and factorized by dense LU
on every rank once per policy evaluation. `coarse_correction = :multiplicative`
applies the coarse correction ``z = R_0^\\mathsf{T} A_0^{-1} R_0 v`` first and the
one-level preconditioner ``M`` to the remaining residual,
``z + M(v - Az)``; `:additive` returns ``M v + z``. The coarse LU uses BLAS: with
Julia threads spinning (`JULIA_THREAD_SLEEP_THRESHOLD=infinite`) run BLAS on
one thread (`BLAS.set_num_threads(1)`), or the competing threads slow the
small factorization down by orders of magnitude.

Warm started from the previous value function. With `forcing = η > 0` the
evaluations are inexact, as in an inexact Newton method (policy iteration is
Newton's method on the discrete HJB equation): evaluation ``k`` is solved to
the relative residual ``\\min(0.1, \\max(\\mathrm{linear\\_tol}, η\\, δ_{k-1} / \\max(1, \\|V\\|_∞)))``,
``δ_{k-1}`` the previous change of the values, so the first evaluations, of
policies still far from optimal, take few Krylov iterations; convergence is
only declared after an evaluation solved to `linear_tol`. With `forcing = 0`
(the default) every evaluation is solved to `linear_tol`.

The policy is improved on each extended box after a halo exchange of the values. Every rank returns the whole
solution. Also solves [`MechanicalHJBProblem`](@ref CellularSheaves.ControlSheaves.MechanicalHJB.MechanicalHJBProblem)s,
whose joint angles are periodic dimensions of the grid. On one box this gives the same discrete solution as
[`PolicyIteration`](@ref), whose `KrylovPolicyEvaluation(method = :bicgstab)`
uses the identical preconditioner on the assembled matrices.

`callback(snapshot)`, if given, is called on every rank with the state of the
iteration as whole-grid host arrays: once with the initial guess
(`iteration = 0`: the initial guess, for `HJBProblem` the Riccati values and the
clipped LQR policy), then after
each policy iteration ``k`` with the values ``V^k`` of the policy just evaluated
and the policy improved from them. `snapshot` is a named tuple
`(iteration, values, controls, change)` with `values` a vector and `controls`
a `d × N` matrix in the layout of [`HJBSolution`](@ref), and `change` the
largest change of the values in that iteration (`NaN` at the start).
"""
Base.@kwdef struct GridPolicyIteration{B,C<:BoxCommunicator,F}
    backend::B = KernelAbstractions.CPU()
    communicator::C = SerialBoxes()
    ranks::Union{Nothing,Vector{Int}} = nothing
    overlap::Int = 1
    preconditioner::Symbol = :red_black_sgs
    tol::Float64 = 1e-8
    maxiter::Int = 50
    linear_tol::Float64 = 1e-10
    linear_maxiter::Int = 2000
    forcing::Float64 = 0.0
    coarse_blocks::Union{Nothing,Vector{Int}} = nothing
    coarse_correction::Symbol = :multiplicative
    callback::F = nothing
end

CommonSolve.solve(prob::HJBProblem, alg::GridPolicyIteration) = _grid_solve(prob, alg)

# The solve for any problem type that implements the interface at the top of this file.
function _grid_solve(prob, alg::GridPolicyIteration)
    @argcheck alg.maxiter >= 1 && alg.tol > 0
    @argcheck alg.preconditioner in (:red_black_sgs, :ras, :none) "preconditioner must be :red_black_sgs, :ras or :none"
    g = _state_grid(prob)
    D, d = ndims(g), _control_count(prob)
    comm = alg.communicator
    points = Tuple(g.points)
    periodic = _periodic(prob)
    # Red–black ordering around a periodic dimension needs an even number of points.
    @argcheck alg.preconditioner === :none || all(k -> !periodic[k] || iseven(points[k]), 1:D) "red–black sweeps need an even number of points in periodic dimensions"
    ranks = alg.ranks === nothing ? balanced_ranks(box_count(comm), points) : Tuple(alg.ranks)
    @argcheck length(ranks) == D && prod(ranks) == box_count(comm) "ranks must have one entry per dimension and multiply to the number of boxes"
    layout = BoxLayout(points, ranks, box_rank(comm); overlap=alg.overlap, periodic)
    exchange(x) = exchange!(comm, layout, x)
    allreduce(x, op) = box_allreduce(comm, x, op)
    backend = alg.backend
    # Controls live on the extended box; the owned box reads them shifted.
    ushift = first.(layout.owned) .- first.(layout.extended)
    U = KernelAbstractions.zeros(backend, Float64, length.(layout.extended)..., d)
    kp_owned = _KernelProblem(prob, first.(layout.owned) .- 1)
    kp_extended = _KernelProblem(prob, first.(layout.extended) .- 1)
    op = box_operator(layout, _UpwindStencil(U, kp_owned, ushift), :owned)
    op_extended = box_operator(layout, _UpwindStencil(U, kp_extended), :extended)
    V, Vnew, b = grid_zeros(op), grid_zeros(op), grid_zeros(op)
    ws = GridWorkspace(op)
    smoother! = if alg.preconditioner === :none
        copyto!
    elseif alg.preconditioner === :ras
        (y, v) -> (exchange(v); red_black_sgs!(y, op_extended, v))
    else
        (y, v) -> red_black_sgs!(y, op, v; exchange)
    end
    # The optional second level: an aggregation coarse space, its dense LU
    # refreshed for each policy (factor[]).
    @argcheck alg.coarse_correction in (:multiplicative, :additive) "coarse_correction must be :multiplicative or :additive"
    coarse = alg.coarse_blocks === nothing ? nothing : AggregateCoarseSpace(layout, alg.coarse_blocks; backend)
    factor = Ref(lu(ones(1, 1)))
    z, Az = grid_zeros(op), grid_zeros(op)
    coarse_solve(v) = factor[] \ coarse_restrict(coarse, v, op, comm)
    precondition! = if coarse === nothing
        smoother!
    elseif alg.coarse_correction === :additive
        (y, v) -> (smoother!(y, v); coarse_prolong_add!(y, coarse, coarse_solve(v), op))
    else
        function (y, v)
            fill!(z, 0)
            coarse_prolong_add!(z, coarse, coarse_solve(v), op)      # z = R₀ᵀ A₀⁻¹ R₀ v
            exchange(z)
            apply!(Az, op, z)
            _lincomb!(Az, op, 1, v, -1, Az, 0, Az)                      # v - A z
            smoother!(y, Az)
            _lincomb!(y, op, 1, y, 1, z, 0, z)                          # z + M (v - A z)
        end
    end
    _launch!(_grid_initialize_kernel!, op_extended, U, V, kp_extended, op_extended.origin)
    # The whole grid's values and controls, on every rank (for the callback and the result).
    function gathered(V)                 # V passed in: the loop swaps it, and a captured V would be boxed
        values = vec(gather_boxes(comm, layout, V, op))
        owned_controls = Array(copy(view(U, ntuple(k -> ushift[k] .+ (1:length(layout.owned[k])), D)..., :)))
        controls = reduce(vcat, (permutedims(vec(gather_boxes(comm, layout, selectdim(owned_controls, D + 1, k))))
                                 for k in 1:d))
        return values, controls
    end
    function report(iteration, change, V)
        alg.callback === nothing && return
        values, controls = gathered(V)
        alg.callback((; iteration, values, controls, change))
    end
    report(0, NaN, V)
    changes, linear_iterations = Float64[], Int[]
    converged = false
    t_assembly = t_setup = t_linear = t_improvement = 0.0
    scale = 1.0
    tight = alg.forcing == 0                 # solve every evaluation to linear_tol
    for _ in 1:alg.maxiter
        δ = isempty(changes) ? 1.0 : changes[end] / max(1.0, scale)
        tolk = tight ? alg.linear_tol : clamp(alg.forcing * δ, alg.linear_tol, 0.1)
        t_assembly += @elapsed _launch!(_grid_rhs_kernel!, op, b, U, kp_owned, op.origin, ushift)
        coarse === nothing || (t_setup += @elapsed factor[] = lu(coarse_matrix(coarse, op, comm)))
        t_linear += @elapsed begin
            copyto!(Vnew, V)
            its, ok, _ = grid_bicgstab!(Vnew, op, precondition!, b, ws; tol=tolk,
                maxiter=alg.linear_maxiter, reduce=x -> allreduce(x, +), exchange)
        end
        ok || @warn "policy evaluation did not converge"
        push!(changes, allreduce(grid_reduce((a, c) -> abs(a - c), max, Vnew, V, op, ws.partial), max))
        push!(linear_iterations, its)
        V, Vnew = Vnew, V
        t_improvement += @elapsed begin
            exchange(V)
            _launch!(_grid_improve_kernel!, op_extended, U, V, kp_extended, op_extended.origin)
        end
        report(length(changes), changes[end], V)
        scale = allreduce(grid_reduce((a, c) -> abs(a), max, V, V, op, ws.partial), max)
        if length(changes) > 1 && changes[end] <= alg.tol * max(1.0, scale)
            if tolk <= alg.linear_tol
                converged = true
                break
            end
            tight = true                     # settled on inexact evaluations: confirm with an exact one
        end
    end
    values, controls = gathered(V)
    seconds = (assembly = t_assembly, setup = t_setup, linear = t_linear, improvement = t_improvement)
    return HJBSolution(prob, values, controls, length(changes), changes, linear_iterations, converged, seconds)
end
