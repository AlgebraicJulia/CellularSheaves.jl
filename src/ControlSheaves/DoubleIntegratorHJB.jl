"""
    DoubleIntegratorHJB

The Hamilton–Jacobi–Bellman equation of a double integrator steering to a
target, solved on a grid by policy iteration, with each policy evaluation
solved by overlapping Schwarz domain decomposition on the overlap sheaf of the
grid (see [`SchwarzMethods`](@ref CellularSheaves.NetworkSheaves.SchwarzMethods)).

**The control problem.** One agent with position ``q \\in \\mathbb R^d`` and
velocity ``v \\in \\mathbb R^d`` (``d = 1`` for one axis, ``d = 2`` for the
planar double integrator) follows ``\\dot q = v``, ``\\dot v = u``, with the
thrust limited to a box ``|u_k| \\le \\bar u`` or a disc ``\\|u\\| \\le \\bar u``.
It minimizes the discounted cost

```math
V(x) = \\min_{u(\\cdot)} \\int_0^\\infty e^{-\\rho t}\\Big(\\tfrac w2\\|q\\|^2
    + \\tfrac c2\\|v\\|^2 + \\tfrac r2\\|u\\|^2\\Big)\\,dt ,
```

with ``q`` measured from the target (the problem is translation invariant, so
one target at ``p`` is the problem with ``q - p`` in place of ``q``). The value
function solves the HJB equation

```math
\\rho V(x) = \\min_{u \\in U}\\Big[\\ell(x, u) + \\nabla V(x)\\cdot f(x, u)\\Big],
\\qquad f(x, u) = (v, u),
```

whose minimizer is the feedback ``u^\\star(x) = -\\operatorname{proj}_U(\\nabla_v V / r)``.
Without a bound this is the discounted linear–quadratic regulator,
``V(x) = \\tfrac12 x^\\mathsf{T} P x`` with ``P`` from a Riccati equation
([`riccati_value_matrix`](@ref)); that value is the exact reference, and the
Dirichlet data outside the grid. With a box bound the planar problem splits
into two independent one-axis problems; with a disc bound it is genuinely
four-dimensional.

**Discretization.** A monotone upwind scheme on a tensor grid
([`StateGrid`](@ref)): for a fixed policy ``u`` the equation becomes the
linear system ``A_u V = b_u`` with

```math
\\Big(\\rho + \\sum_k \\tfrac{|f_k|}{h_k}\\Big) V_x - \\sum_k \\tfrac{|f_k|}{h_k} V_{x + \\operatorname{sign}(f_k) h_k e_k}
    = \\ell(x, u(x)),
```

a nonsingular M-matrix (strictly diagonally dominant by ``\\rho``). This is the
Markov-chain approximation of Kushner and Dupuis (*Numerical Methods for
Stochastic Control Problems in Continuous Time*, 2001).

**Policy iteration** (Howard 1960; Bokanowski, Maroso and Zidani, *Some
convergence results for Howard's algorithm*, SIAM J. Numer. Anal. 2009)
alternates evaluation, ``A_u V = b_u``, with improvement, which minimizes the
discrete Hamiltonian at every grid point in closed form. Values decrease
monotonically and converge to the solution of the discrete HJB equation,
typically in a handful of iterations.

**The sheaf.** Each policy evaluation is a sparse nonsymmetric M-matrix system,
solved by [`SchwarzPolicyEvaluation`](@ref): the grid is cut into overlapping
boxes, which are the vertices of the overlap sheaf, with their overlaps as edges
and ghost layers carried by the restriction maps. A Schwarz sweep is then a
round of local HJB solves with boundary data from the neighbouring subdomains
(compare the splitting method of Falcone, Lanucara and Seghini, *A splitting
algorithm for Hamilton–Jacobi–Bellman equations*, Appl. Numer. Math. 1994).
"""
module DoubleIntegratorHJB

using ArgCheck
using CommonSolve
using CommonSolve: solve
using LinearAlgebra
using SparseArrays
using CellularSheaves.NetworkSheaves.SchwarzMethods: SchwarzDecomposition, SchwarzProblem,
    SchwarzIteration, MulticolorSweep, refactor!, _stable_lu, _threaded_mul!,
    LocalSolver, ExactLocalSolve, SymmetricGaussSeidel, _bicgstab!
using Krylov: gmres

export StateGrid, HJBProblem, riccati_value_matrix, riccati_value,
    PolicyIteration, DirectPolicyEvaluation, KrylovPolicyEvaluation, SchwarzPolicyEvaluation, HJBSolution,
    grid_partition, grid_subdomains, value_at, control_at

# ===========================================================================
# grid
# ===========================================================================

"""
    StateGrid(lower, upper, points)

A tensor grid on the box ``\\prod_k [\\text{lower}_k, \\text{upper}_k]`` with
`points[k]` equally spaced points in dimension `k`, numbered in column-major
order (the first dimension varies fastest). For the double integrator the
dimensions are ordered ``(q_1, \\dots, q_d, v_1, \\dots, v_d)``.
"""
struct StateGrid
    lower::Vector{Float64}
    upper::Vector{Float64}
    points::Vector{Int}
    spacing::Vector{Float64}
    strides::Vector{Int}
end

function StateGrid(lower::AbstractVector{<:Real}, upper::AbstractVector{<:Real}, points::AbstractVector{<:Integer})
    @argcheck length(lower) == length(upper) == length(points)
    @argcheck all(upper .> lower) "need lower < upper in every dimension"
    @argcheck all(points .>= 3) "need at least 3 points per dimension"
    spacing = (Float64.(upper) .- Float64.(lower)) ./ (points .- 1)
    strides = [prod(points[1:(k - 1)]; init = 1) for k in eachindex(points)]
    return StateGrid(Float64.(lower), Float64.(upper), Vector{Int}(points), spacing, strides)
end

Base.ndims(g::StateGrid) = length(g.points)
Base.length(g::StateGrid) = prod(g.points)
_cartesian(g::StateGrid) = CartesianIndices(Tuple(g.points))

function _coordinates!(x, g::StateGrid, I::CartesianIndex)
    @inbounds for k in eachindex(x)
        x[k] = g.lower[k] + (I[k] - 1) * g.spacing[k]
    end
    return x
end

"""
    grid_partition(grid, blocks) -> Vector{Int}

Label each grid point with the box of a `blocks[1] × blocks[2] × …` partition
of the grid into nearly equal index ranges.
"""
function grid_partition(g::StateGrid, blocks::AbstractVector{<:Integer})
    @argcheck length(blocks) == ndims(g)
    @argcheck all(1 .<= blocks .<= g.points)
    shape = Tuple(blocks)
    labels = LinearIndices(shape)
    return vec([labels[ntuple(k -> 1 + ((I[k] - 1) * blocks[k]) ÷ g.points[k], ndims(g))...] for I in _cartesian(g)])
end

"""
    grid_subdomains(grid, blocks; overlap = 1) -> Vector{Vector{Int}}

The boxes of [`grid_partition`](@ref), each extended by `overlap` grid points
in every direction: the overlapping subdomains of the overlap sheaf.
"""
function grid_subdomains(g::StateGrid, blocks::AbstractVector{<:Integer}; overlap::Integer = 1)
    @argcheck length(blocks) == ndims(g)
    @argcheck overlap >= 0
    ranges(k, b) = begin
        first_ = 1 + cld((b - 1) * g.points[k], blocks[k])
        last_ = cld(b * g.points[k], blocks[k])
        max(1, first_ - overlap):min(g.points[k], last_ + overlap)
    end
    L = LinearIndices(Tuple(g.points))
    return vec([vec(collect(L[CartesianIndices(ntuple(k -> ranges(k, B[k]), ndims(g)))]))
            for B in CartesianIndices(Tuple(blocks))])
end

# ===========================================================================
# problem
# ===========================================================================

"""
    riccati_value_matrix(axes, position_weight, velocity_weight, control_weight, discount) -> Matrix

The matrix ``P`` of the discounted LQR value ``V(x) = \\tfrac12 x^\\mathsf{T} P x``
of the unbounded double integrator in `axes` dimensions: the stabilizing
solution of

```math
Q + (A - \\tfrac\\rho2 I)^\\mathsf{T} P + P (A - \\tfrac\\rho2 I) - P B R^{-1} B^\\mathsf{T} P = 0,
```

computed from the stable invariant subspace of the Hamiltonian matrix.
"""
function riccati_value_matrix(axes::Integer, w::Real, c::Real, r::Real, ρ::Real)
    d = Int(axes)
    n = 2d
    A = [zeros(d, d) Matrix(1.0I, d, d); zeros(d, n)] - (ρ / 2) * I
    B = [zeros(d, d); Matrix(1.0I, d, d)]
    Q = Matrix(Diagonal([fill(Float64(w), d); fill(Float64(c), d)]))
    H = [A -(B * B') / r; -Q -A']
    F = schur(H)
    ordschur!(F, real.(F.values) .< 0)
    X1, X2 = F.Z[1:n, 1:n], F.Z[(n + 1):end, 1:n]
    P = X2 / X1
    return (P + P') / 2
end

"""
    HJBProblem(grid; position_weight = 1, velocity_weight = 0.1, control_weight = 0.1,
               discount = 0.5, control_bound = Inf, constraint = :box)

The discounted HJB problem of the module docstring on `grid`, which must have
2 dimensions ``(q, v)`` or 4 dimensions ``(q_1, q_2, v_1, v_2)``.
`constraint` is `:box` (``|u_k| \\le \\bar u``) or `:disc`
(``\\|u\\| \\le \\bar u``, planar only). Outside the grid the value is taken
from the unbounded LQR value [`riccati_value`](@ref).
"""
struct HJBProblem
    grid::StateGrid
    axes::Int
    position_weight::Float64
    velocity_weight::Float64
    control_weight::Float64
    discount::Float64
    control_bound::Float64
    constraint::Symbol
    riccati::Matrix{Float64}
end

function HJBProblem(grid::StateGrid; position_weight::Real = 1.0, velocity_weight::Real = 0.1,
        control_weight::Real = 0.1, discount::Real = 0.5, control_bound::Real = Inf, constraint::Symbol = :box)
    @argcheck ndims(grid) in (2, 4) "the grid must have 2 (one axis) or 4 (planar) dimensions"
    axes = ndims(grid) ÷ 2
    @argcheck position_weight > 0 && velocity_weight >= 0 && control_weight > 0
    @argcheck discount > 0 "the discount rate must be positive"
    @argcheck control_bound > 0
    @argcheck constraint in (:box, :disc)
    @argcheck constraint === :box || axes == 2 "the disc constraint needs the planar (4-dimensional) problem"
    P = riccati_value_matrix(axes, position_weight, velocity_weight, control_weight, discount)
    return HJBProblem(grid, axes, position_weight, velocity_weight, control_weight, discount,
        control_bound, constraint, P)
end

"""
    riccati_value(problem, x) -> Float64

The value ``\\tfrac12 x^\\mathsf{T} P x`` of the unbounded problem at `x`.
"""
riccati_value(prob::HJBProblem, x::AbstractVector) = dot(x, prob.riccati, x) / 2

_running_cost(prob::HJBProblem, x, u) =
    (prob.position_weight * sum(abs2, view(x, 1:prob.axes)) +
     prob.velocity_weight * sum(abs2, view(x, (prob.axes + 1):(2prob.axes))) +
     prob.control_weight * sum(abs2, u)) / 2

function _project!(u, prob::HJBProblem)
    ū = prob.control_bound
    isfinite(ū) || return u
    if prob.constraint === :box
        u .= clamp.(u, -ū, ū)
    else
        nu = norm(u)
        nu > ū && (u .*= ū / nu)
    end
    return u
end

# Value at the neighbour of grid point n (with coordinates x) one step of size
# ±h in dimension j: from V inside the grid, from the Riccati value outside.
@inline function _neighbour_value(prob::HJBProblem, V, x, I::CartesianIndex, n::Int, j::Int, step::Int)
    g = prob.grid
    if 1 <= I[j] + step <= g.points[j]
        return V[n + step * g.strides[j]]
    end
    x[j] += step * g.spacing[j]
    value = riccati_value(prob, x)
    x[j] -= step * g.spacing[j]
    return value
end

# ===========================================================================
# policy evaluation: assembly
# ===========================================================================

# The full first-order stencil of the grid (each point and its ±1 neighbours
# in every dimension) as the rows of a sparse matrix: the union of the
# patterns of every upwind matrix A_u. Rows are stored as the columns of Aᵀ
# (`colptr`, `rowval`), with the slot of each entry so that rows can be filled
# in parallel; a slot is 0 where the neighbour is outside the grid.
struct _Stencil
    colptr::Vector{Int}
    rowval::Vector{Int}
    diagonal::Vector{Int}
    plus::Matrix{Int}
    minus::Matrix{Int}
end

function _stencil(g::StateGrid)
    D, N = ndims(g), length(g)
    C = _cartesian(g)
    colptr = ones(Int, N + 1)
    for n in 1:N
        I = C[n]
        colptr[n + 1] = colptr[n] + 1 + count(j -> I[j] > 1, 1:D) + count(j -> I[j] < g.points[j], 1:D)
    end
    rowval = Vector{Int}(undef, colptr[end] - 1)
    diagonal = zeros(Int, N)
    plus, minus = zeros(Int, D, N), zeros(Int, D, N)
    for n in 1:N
        I = C[n]
        p = colptr[n]
        for j in D:-1:1                   # n − s_D < … < n − s_1 < n < n + s_1 < … < n + s_D
            I[j] > 1 || continue
            rowval[p] = n - g.strides[j]; minus[j, n] = p; p += 1
        end
        rowval[p] = n; diagonal[n] = p; p += 1
        for j in 1:D
            I[j] < g.points[j] || continue
            rowval[p] = n + g.strides[j]; plus[j, n] = p; p += 1
        end
    end
    return _Stencil(colptr, rowval, diagonal, plus, minus)
end

# The full stencil as a symmetric pattern matrix (for the Schwarz cover).
_stencil_pattern(st::_Stencil) =
    SparseMatrixCSC(length(st.diagonal), length(st.diagonal), copy(st.colptr), copy(st.rowval), ones(length(st.rowval)))

# The upwind system A_u V = b_u for the controls U (d × N), assembled row by row
# in parallel into the fixed stencil. Returns A, Aᵀ (A by rows) and b; entries
# of the stencil the policy does not use are stored zeros.
function _assemble(prob::HJBProblem, U::AbstractMatrix, st::_Stencil)
    g = prob.grid
    d, D, N = prob.axes, ndims(g), length(g)
    C = _cartesian(g)
    vals = zeros(length(st.rowval))
    b = zeros(N)
    chunks = collect(Iterators.partition(1:N, cld(N, 8 * Threads.nthreads())))
    Threads.@threads for chunk in chunks
        x = zeros(D)
        @inbounds for n in chunk
            I = C[n]
            _coordinates!(x, g, I)
            u = view(U, :, n)
            diagonal = prob.discount
            rhs = _running_cost(prob, x, u)
            for j in 1:D
                fj = j <= d ? x[d + j] : u[j - d]
                iszero(fj) && continue
                a = abs(fj) / g.spacing[j]
                diagonal += a
                slot = fj > 0 ? st.plus[j, n] : st.minus[j, n]
                if slot != 0
                    vals[slot] = -a
                else
                    rhs += a * _neighbour_value(prob, nothing, x, I, n, j, fj > 0 ? 1 : -1)
                end
            end
            vals[st.diagonal[n]] = diagonal
            b[n] = rhs
        end
    end
    At = SparseMatrixCSC(N, N, st.colptr, st.rowval, vals)
    return copy(transpose(At)), At, b
end

function _assemble(prob::HJBProblem, U::AbstractMatrix)
    A, _, b = _assemble(prob, U, _stencil(prob.grid))
    return A, b
end

# ===========================================================================
# policy improvement
# ===========================================================================

# Minimize the discrete Hamiltonian ½r‖u‖² + Σ_k (u_k⁺ D⁺_k − u_k⁻ D⁻_k) over the
# constraint set, writing the minimizer into u. The upwind difference D± used
# for component k is the one the assembly uses for that sign of u_k, so this is
# exactly Howard's improvement step for the discrete problem.
function _minimize_hamiltonian!(u, prob::HJBProblem, Dp, Dm)
    r, ū = prob.control_weight, prob.control_bound
    if prob.constraint === :box || prob.axes == 1
        @inbounds for k in eachindex(u)
            up = clamp(-Dp[k] / r, 0.0, ū)
            um = clamp(-Dm[k] / r, -ū, 0.0)
            Hp = r * up^2 / 2 + up * Dp[k]
            Hm = r * um^2 / 2 + um * Dm[k]
            u[k] = Hp <= Hm ? up : um
        end
        return u
    end
    # Disc: in each quadrant the problem is the projection of −g/r onto the
    # quadrant ∩ disc, i.e. clip the signs, then scale into the disc.
    best = Inf
    @inbounds for s1 in (1, -1), s2 in (1, -1)
        g1 = s1 > 0 ? Dp[1] : Dm[1]
        g2 = s2 > 0 ? Dp[2] : Dm[2]
        z1 = s1 > 0 ? max(-g1 / r, 0.0) : min(-g1 / r, 0.0)
        z2 = s2 > 0 ? max(-g2 / r, 0.0) : min(-g2 / r, 0.0)
        nz = hypot(z1, z2)
        if nz > ū
            z1 *= ū / nz
            z2 *= ū / nz
        end
        H = r * (z1^2 + z2^2) / 2 + z1 * g1 + z2 * g2
        if H < best
            best = H
            u[1], u[2] = z1, z2
        end
    end
    return u
end

# Pointwise, so the grid is split into chunks improved concurrently.
function _improve!(U::AbstractMatrix, prob::HJBProblem, V::AbstractVector)
    g = prob.grid
    d, D = prob.axes, ndims(g)
    C = _cartesian(g)
    chunks = collect(Iterators.partition(1:length(g), cld(length(g), 8 * Threads.nthreads())))
    Threads.@threads for chunk in chunks
        x = zeros(D)
        Dp, Dm, u = zeros(d), zeros(d), zeros(d)
        @inbounds for n in chunk
            I = C[n]
            _coordinates!(x, g, I)
            for k in 1:d
                j = d + k
                h = g.spacing[j]
                Dp[k] = (_neighbour_value(prob, V, x, I, n, j, 1) - V[n]) / h
                Dm[k] = (V[n] - _neighbour_value(prob, V, x, I, n, j, -1)) / h
            end
            _minimize_hamiltonian!(u, prob, Dp, Dm)
            U[:, n] .= u
        end
    end
    return U
end

# The unbounded LQR feedback u = −R⁻¹BᵀPx, projected onto the constraint set.
function _lqr_controls(prob::HJBProblem)
    g = prob.grid
    d = prob.axes
    K = prob.riccati[(d + 1):(2d), :] ./ prob.control_weight
    U = zeros(d, length(g))
    x = zeros(ndims(g))
    for (n, I) in enumerate(_cartesian(g))
        _coordinates!(x, g, I)
        u = -K * x
        U[:, n] .= _project!(u, prob)
    end
    return U
end

# ===========================================================================
# policy evaluation: linear solvers
# ===========================================================================

"""
    DirectPolicyEvaluation()

Solve each policy evaluation ``A_u V = b_u`` with a sparse LU factorization of
the whole grid (partial pivoting; see `SchwarzMethods`). The reference; its
fill-in grows quickly with the dimension
(on a ``n^4`` grid the separators have ``n^3`` points).
"""
struct DirectPolicyEvaluation end

"""
    SchwarzPolicyEvaluation(blocks; overlap = 1, algorithm = SchwarzIteration(...),
                            local_solver = ExactLocalSolve())

Solve each policy evaluation by Schwarz domain decomposition: the grid is cut
into `blocks[1] × blocks[2] × …` boxes ([`grid_partition`](@ref)), extended by
`overlap` points ([`grid_subdomains`](@ref)); each box owns its points. Every
policy iteration refactors one `SchwarzDecomposition` (built at the first
step, then updated in place with `refactor!`) for the new
``A_u`` on the same subdomains and solves it with `algorithm`
(a `SchwarzIteration` or `SchwarzGMRES`), warm-started from the previous
value function. With `local_solver = SymmetricGaussSeidelLocalSolve()` each
subdomain does the same symmetric Gauss–Seidel pass as the
[`KrylovPolicyEvaluation`](@ref) smoothers (one shared kernel) instead of an
exact sparse LU, so setup only gathers the new matrix values.
"""
struct SchwarzPolicyEvaluation{A}
    blocks::Vector{Int}
    overlap::Int
    algorithm::A
    local_solver::LocalSolver
end

SchwarzPolicyEvaluation(blocks::AbstractVector{<:Integer}; overlap::Integer = 1,
        algorithm = SchwarzIteration(sweep = MulticolorSweep(), tol = 1e-10, maxiter = 10_000),
        local_solver::LocalSolver = ExactLocalSolve()) =
    SchwarzPolicyEvaluation(Vector{Int}(blocks), Int(overlap), algorithm, local_solver)

"""
    KrylovPolicyEvaluation(; method = :gmres, preconditioner = :symmetric_gauss_seidel,
                           blocks = nothing, tol = 1e-10, maxiter = 5000, memory = 50)

Solve each policy evaluation on the whole grid with restarted GMRES
(`method = :gmres`, Krylov.jl) or with a BiCGStab whose vector operations are
all threaded (`method = :bicgstab`, the solver of `SchwarzBiCGStab`), right
preconditioned by a point smoother of ``A_u = D + L + U``: `:gauss_seidel`
(``(D + L)^{-1}``), `:symmetric_gauss_seidel`
(``(D + U)^{-1} D (D + L)^{-1}``), `:jacobi` (``D^{-1}``), `:none`, or
`:block_symmetric_gauss_seidel`: block Jacobi over the boxes of
[`grid_partition`](@ref)`(grid, blocks)` with a symmetric Gauss–Seidel sweep
inside each box, the boxes swept concurrently. Warm started from the previous
value function; matrix–vector products are threaded.

Policy iteration is a semismooth Newton method for the discrete HJB equation
(Bokanowski, Maroso and Zidani 2009). With a point smoother this is the serial
Newton–Krylov baseline; the block smoother is its parallel counterpart, equal
to it with one block, and the nonoverlapping, inexact relative of
[`SchwarzPolicyEvaluation`](@ref).
"""
Base.@kwdef struct KrylovPolicyEvaluation
    method::Symbol = :gmres
    preconditioner::Symbol = :symmetric_gauss_seidel
    blocks::Union{Nothing, Vector{Int}} = nothing
    tol::Float64 = 1e-10
    maxiter::Int = 5000
    memory::Int = 50
end

struct _Jacobi
    inverse_diagonal::Vector{Float64}
end
LinearAlgebra.ldiv!(y::AbstractVector, P::_Jacobi, x::AbstractVector) = (y .= P.inverse_diagonal .* x)

struct _GaussSeidel{L}
    lower::L
end
LinearAlgebra.ldiv!(y::AbstractVector, P::_GaussSeidel, x::AbstractVector) = ldiv!(P.lower, copyto!(y, x))

# The symmetric Gauss–Seidel pass is SchwarzMethods.SymmetricGaussSeidel, the
# same kernel SchwarzPolicyEvaluation uses inside each subdomain.
_symmetric_gauss_seidel(A::SparseMatrixCSC) = SymmetricGaussSeidel(A)

# Block Jacobi with a symmetric Gauss–Seidel sweep inside each block: the
# blocks are independent and are swept concurrently. With one block it is the
# global symmetric Gauss–Seidel preconditioner.
struct _BlockSymmetricGaussSeidel{P}
    blocks::Vector{Vector{Int}}
    smoothers::Vector{P}
end

function _block_symmetric_gauss_seidel(A::SparseMatrixCSC, blocks::Vector{Vector{Int}})
    first_smoother = _symmetric_gauss_seidel(A[blocks[1], blocks[1]])
    smoothers = Vector{typeof(first_smoother)}(undef, length(blocks))
    smoothers[1] = first_smoother
    Threads.@threads for k in 2:length(blocks)
        smoothers[k] = _symmetric_gauss_seidel(A[blocks[k], blocks[k]])
    end
    return _BlockSymmetricGaussSeidel(blocks, smoothers)
end

function LinearAlgebra.ldiv!(y::AbstractVector, P::_BlockSymmetricGaussSeidel, x::AbstractVector)
    Threads.@threads for k in eachindex(P.blocks)
        idx = P.blocks[k]
        yk = x[idx]
        ldiv!(yk, P.smoothers[k], copy(yk))
        y[idx] = yk
    end
    return y
end

function _preconditioner(A::SparseMatrixCSC, kind::Symbol, blocks)
    kind === :none && return I
    kind === :jacobi && return _Jacobi(1 ./ Vector(diag(A)))
    kind === :gauss_seidel && return _GaussSeidel(LowerTriangular(tril(A)))
    kind === :symmetric_gauss_seidel && return _symmetric_gauss_seidel(A)
    if kind === :block_symmetric_gauss_seidel
        @argcheck blocks !== nothing "the block preconditioner needs `blocks`"
        return _block_symmetric_gauss_seidel(A, blocks)
    end
    throw(ArgumentError("unknown preconditioner $kind"))
end

# A by rows (Aᵀ stored column-wise): a threaded matrix–vector product for GMRES.
struct _RowMatrix
    At::SparseMatrixCSC{Float64, Int}
end
Base.size(M::_RowMatrix) = (size(M.At, 2), size(M.At, 1))
Base.size(M::_RowMatrix, d::Integer) = size(M)[d]
Base.eltype(::_RowMatrix) = Float64
LinearAlgebra.mul!(y::AbstractVector, M::_RowMatrix, x::AbstractVector) = _threaded_mul!(y, M.At, x)

function _block_indices(labels::Vector{Int})
    blocks = [Int[] for _ in 1:maximum(labels)]
    for (k, b) in enumerate(labels)
        push!(blocks[b], k)
    end
    return blocks
end

_evaluation_data(::HJBProblem, ::DirectPolicyEvaluation, st) = nothing
_evaluation_data(prob::HJBProblem, e::KrylovPolicyEvaluation, st) =
    e.blocks === nothing ? nothing : _block_indices(grid_partition(prob.grid, e.blocks))
_evaluation_data(prob::HJBProblem, e::SchwarzPolicyEvaluation, st) =
    (grid_subdomains(prob.grid, e.blocks; overlap = e.overlap), grid_partition(prob.grid, e.blocks),
     _stencil_pattern(st), Ref{Any}(nothing))

# Each evaluation returns (V, iterations, converged, setup seconds).
function _evaluate(::DirectPolicyEvaluation, data, A, At, b, V0)
    setup = @elapsed F = _stable_lu(A)
    V = F \ b
    residual = norm(A * V - b) / max(norm(b), eps())
    return V, 0, residual <= 1e-8, setup
end

function _evaluate(e::KrylovPolicyEvaluation, blocks, A, At, b, V0)
    @argcheck e.method in (:gmres, :bicgstab) "unknown Krylov method $(e.method)"
    setup = @elapsed P = _preconditioner(A, e.preconditioner, blocks)
    op = _RowMatrix(At)
    r0 = b - mul!(similar(b), op, V0)
    iszero(norm(r0)) && return copy(V0), 0, true, setup
    rtol = e.tol * max(norm(b), eps()) / norm(r0)
    if e.method === :bicgstab
        dV = zeros(length(b))
        precondition! = P isa UniformScaling ? copyto! : (y, x) -> ldiv!(y, P, x)
        its, ok, _ = _bicgstab!(dV, (y, x) -> _threaded_mul!(y, At, x), precondition!, r0;
            tol = min(rtol, 0.5), maxiter = e.maxiter)
        return V0 + dV, its, ok, setup
    end
    dV, stats = gmres(op, r0; N = P, ldiv = !(P isa UniformScaling), memory = e.memory, restart = true,
        atol = 0.0, rtol = min(rtol, 0.5), itmax = e.maxiter)
    return V0 + dV, stats.niter, stats.solved, setup
end

# The cover, ownership and coloring are built once, from the full stencil (its
# explicit zeros kept), at the first Newton step; later steps gather the new
# values into the local problems and refactor their LUs numerically in place.
function _evaluate(e::SchwarzPolicyEvaluation, (subdomains, parts, structure, cached), A, At, b, V0)
    setup = @elapsed dd = cached[] === nothing ?
        SchwarzDecomposition(A, subdomains; owner = parts, structure, dropzeros = false,
local_solver = e.local_solver) :
        refactor!(cached[], A; At)
    cached[] = dd
    result = solve(SchwarzProblem(dd, b; u0 = V0), e.algorithm)
    return result.u, result.iterations, result.converged, setup
end

# ===========================================================================
# policy iteration
# ===========================================================================

"""
    PolicyIteration(; evaluation = DirectPolicyEvaluation(), tol = 1e-8, maxiter = 50)

Howard's policy iteration: evaluate the current policy with `evaluation`
([`DirectPolicyEvaluation`](@ref), [`KrylovPolicyEvaluation`](@ref) or
[`SchwarzPolicyEvaluation`](@ref)), then
improve it pointwise. Starts from the clipped LQR feedback and stops when the
value function changes by at most `tol` (relative to its largest value) in
one iteration.
"""
Base.@kwdef struct PolicyIteration{E}
    evaluation::E = DirectPolicyEvaluation()
    tol::Float64 = 1e-8
    maxiter::Int = 50
end

"""
    HJBSolution

The result of `solve(problem, PolicyIteration(...))`.

# Fields
- `problem`: the [`HJBProblem`](@ref).
- `values`: the value function at the grid points.
- `controls`: `d × N`, the optimal feedback at the grid points.
- `iterations`: policy iterations used.
- `value_changes`: per iteration, the largest change of the value function.
- `linear_iterations`: per iteration, Schwarz or GMRES iterations of the
  policy evaluation (zero for the direct solver).
- `converged`: whether the value change fell below the tolerance.
- `seconds`: wall time summed over the iterations, split into `assembly` of
  ``A_u, b_u``, `setup` of the linear solver (LU, Schwarz decomposition with
  its local factorizations, or preconditioner), the `linear` solve itself, and
  policy `improvement`.
"""
struct HJBSolution
    problem::HJBProblem
    values::Vector{Float64}
    controls::Matrix{Float64}
    iterations::Int
    value_changes::Vector{Float64}
    linear_iterations::Vector{Int}
    converged::Bool
    seconds::NamedTuple{(:assembly, :setup, :linear, :improvement), NTuple{4, Float64}}
end

function CommonSolve.solve(prob::HJBProblem, alg::PolicyIteration)
    @argcheck alg.maxiter >= 1 && alg.tol > 0
    g = prob.grid
    U = _lqr_controls(prob)
    x = zeros(ndims(g))
    V = vec([riccati_value(prob, _coordinates!(x, g, I)) for I in _cartesian(g)])
    st = _stencil(g)
    data = _evaluation_data(prob, alg.evaluation, st)
    changes, linear_iterations = Float64[], Int[]
    converged = false
    t_assembly = t_setup = t_linear = t_improvement = 0.0
    for _ in 1:alg.maxiter
        t_assembly += @elapsed A, At, b = _assemble(prob, U, st)
        t_eval = @elapsed Vnew, its, ok, setup = _evaluate(alg.evaluation, data, A, At, b, V)
        t_setup += setup
        t_linear += t_eval - setup
        ok || @warn "policy evaluation did not converge"
        push!(changes, maximum(abs, Vnew - V))
        push!(linear_iterations, its)
        V = Vnew
        t_improvement += @elapsed _improve!(U, prob, V)
        if length(changes) > 1 && changes[end] <= alg.tol * max(1.0, maximum(abs, V))
            converged = true
            break
        end
    end
    seconds = (assembly = t_assembly, setup = t_setup, linear = t_linear, improvement = t_improvement)
    return HJBSolution(prob, V, U, length(changes), changes, linear_iterations, converged, seconds)
end

# ===========================================================================
# interpolation
# ===========================================================================

# Multilinear interpolation of the columns of F (one column per grid point, or
# a vector) at x; nothing outside the grid.
function _interpolate(g::StateGrid, F::AbstractVecOrMat, x::AbstractVector)
    D = ndims(g)
    lo = zeros(Int, D)
    t = zeros(D)
    for k in 1:D
        s = (x[k] - g.lower[k]) / g.spacing[k]
        (s < 0 || s > g.points[k] - 1) && return nothing
        lo[k] = min(floor(Int, s), g.points[k] - 2)
        t[k] = s - lo[k]
    end
    acc = F isa AbstractVector ? 0.0 : zeros(size(F, 1))
    for corner in 0:(2^D - 1)
        weight = 1.0
        n = 1
        for k in 1:D
            bit = (corner >> (k - 1)) & 1
            weight *= bit == 1 ? t[k] : 1 - t[k]
            n += (lo[k] + bit) * g.strides[k]
        end
        iszero(weight) && continue
        if F isa AbstractVector
            acc += weight * F[n]
        else
            acc .+= weight .* view(F, :, n)
        end
    end
    return acc
end

"""
    value_at(solution, x) -> Float64

The value function at the state `x`, by multilinear interpolation on the grid;
the unbounded LQR value outside it.
"""
function value_at(sol::HJBSolution, x::AbstractVector)
    v = _interpolate(sol.problem.grid, sol.values, x)
    return v === nothing ? riccati_value(sol.problem, x) : v
end

"""
    control_at(solution, x) -> Vector

The optimal feedback at the state `x`: the grid policy interpolated
multilinearly and projected onto the constraint set; the clipped LQR feedback
outside the grid.
"""
function control_at(sol::HJBSolution, x::AbstractVector)
    prob = sol.problem
    u = _interpolate(prob.grid, sol.controls, x)
    if u === nothing
        d = prob.axes
        u = -(prob.riccati[(d + 1):(2d), :] * x) ./ prob.control_weight
    end
    return _project!(u, prob)
end

end # module DoubleIntegratorHJB
