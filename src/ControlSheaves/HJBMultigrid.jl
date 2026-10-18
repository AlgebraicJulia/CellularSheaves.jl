# Geometric multigrid for the policy evaluations of GridPolicyIteration: a
# V-cycle over coarse grids (GridMultigrid), preconditioning BiCGStab on the
# fine grid, with red–black SGS smoothing on every level but the coarsest,
# which is solved directly (dense LU, or PathPack's semiring LU with the
# PATHPACK package extension). Included in DoubleIntegratorHJB.

export CoarseCorrection, AggregateCoarseCorrection, Multigrid, CoarseOperator, Rediscretization, GalerkinOperator,
    CoarsestSolver, PathPackCoarsest, DenseCoarsest

"""
    CoarseCorrection

A second level for the preconditioner of [`GridPolicyIteration`](@ref) (its
`coarse` option): [`AggregateCoarseCorrection`](@ref) or [`Multigrid`](@ref).
"""
abstract type CoarseCorrection end

"""
    AggregateCoarseCorrection(blocks; composition = MultiplicativeComposition())

A piecewise-constant coarse space on a grid of `blocks` aggregates per dimension
(e.g. `[4, 4, 2, 2]`;
[`AggregateCoarseSpace`](@ref CellularSheaves.NetworkSheaves.GridSchwarz.AggregateCoarseSpace)),
independent of the boxes, whose Galerkin operator ``A_0 = R_0 A R_0^\\mathsf{T}``
is computed from the stencil, summed over all ranks and factorized by dense LU
on every rank once per policy evaluation. With `MultiplicativeComposition()`
the coarse correction ``z = R_0^\\mathsf{T} A_0^{-1} R_0 v`` comes first and the
one-level smoother ``M`` acts on the remaining residual, ``z + M(v - Az)``;
with `AdditiveComposition()` the preconditioner is ``M v + z``.

The coarse LU uses BLAS: with Julia threads spinning
(`JULIA_THREAD_SLEEP_THRESHOLD=infinite`) run BLAS on one thread
(`BLAS.set_num_threads(1)`), or the competing threads slow the small
factorization down by orders of magnitude.
"""
struct AggregateCoarseCorrection{C<:LevelComposition} <: CoarseCorrection
    blocks::Vector{Int}
    composition::C
end
AggregateCoarseCorrection(blocks::AbstractVector{<:Integer}; composition::LevelComposition = MultiplicativeComposition()) =
    AggregateCoarseCorrection(Vector{Int}(blocks), composition)

"""
    CoarseOperator

The operator of each coarse level of a [`Multigrid`](@ref):
[`Rediscretization`](@ref) or [`GalerkinOperator`](@ref).
"""
abstract type CoarseOperator end

"""
    Rediscretization()

Each coarse level is the same HJB problem discretized on the coarse grid, under
the policy pushed forward to it (the block average of the controls, which stays
in the convex control set).
"""
struct Rediscretization <: CoarseOperator end

"""
    GalerkinOperator()

Each coarse level is the composite ``T A E`` of the next finer level's operator.
"""
struct GalerkinOperator <: CoarseOperator end

"""
    CoarsestSolver

How a [`Multigrid`](@ref) solves its coarsest grid: [`PathPackCoarsest`](@ref) or
[`DenseCoarsest`](@ref). A new solver `S` implements `_coarsest_solver(::S, P)`,
returning a solver of ``V = P V + c`` for matrices with the pattern of `P`,
refactored by `_refactor!(solver, P)` and applied by `_solve!(solver, x, c)`.
"""
abstract type CoarsestSolver end

"""
    PathPackCoarsest()

The closure ``V = P^* (D^{-1} b)`` of the coarsest grid's discounted transition
matrix ``P = D^{-1} N`` over the (+, ×) semiring, with PathPack's chordal LU
(`using PATHPACK`: the symbolic analysis is done once and every policy
refactorizes numerically).
"""
struct PathPackCoarsest <: CoarsestSolver end

"""
    DenseCoarsest()

LAPACK's LU of the dense matrix ``I - P`` (small coarsest grids).
"""
struct DenseCoarsest <: CoarsestSolver end

"""
    Multigrid(; levels = nothing, coarsest_size = 4096, coarse_operator = Rediscretization(),
              coarsest = PathPackCoarsest(), presmooth = 1, postsmooth = 1)

A V-cycle preconditioner for [`GridPolicyIteration`](@ref): pass it as
`GridPolicyIteration(coarse = Multigrid(...))`. Each coarse grid halves
the points of every dimension that can be coarsened (an odd count leaves one
block of one point; a periodic dimension is halved while its count is a
multiple of 4; see [`GridMultigrid`](@ref CellularSheaves.NetworkSheaves.GridMultigrid)),
and the levels are linked by the pullback ``E`` (copy each coarse value onto
its block) and the pushforward transfer ``T`` (block average) of the constant
sheaf along the aggregation homomorphism.

- `levels = nothing` coarsens until the coarsest grid has at most
  `coarsest_size` points, or cannot be coarsened further: ``O(\\log n)``
  levels on a grid of ``n`` points per dimension. An integer asks for exactly
  that many levels, fine grid included.
- `coarse_operator`: [`Rediscretization`](@ref) or [`GalerkinOperator`](@ref).
  The two coincide where each drift component is constant along its own axis
  within a block, and otherwise differ at first order in the grid spacing
  (Galerkin takes the drift on the outflow face, rediscretization at the block
  centre).
- `coarsest`: the [`CoarsestSolver`](@ref) of the coarsest grid's evaluation.
- `presmooth`, `postsmooth`: red–black SGS sweeps before and after each
  coarse-grid correction.

Each application of the preconditioner is one V-cycle from a zero initial
guess. All work arrays are allocated once per solve. Single box only
(`SerialBoxes()`; CPU threads or one GPU, whose coarsest solve runs on the host).
"""
Base.@kwdef struct Multigrid{O<:CoarseOperator,S<:CoarsestSolver} <: CoarseCorrection
    levels::Union{Nothing,Int} = nothing
    coarsest_size::Int = 4096
    coarse_operator::O = Rediscretization()
    coarsest::S = PathPackCoarsest()
    presmooth::Int = 1
    postsmooth::Int = 1
end

# ---------------------------------------------------------------------------
# problems on coarse grids
# ---------------------------------------------------------------------------

_with_grid(prob::HJBProblem, grid::StateGrid) =
    HJBProblem(grid, prob.axes, prob.position_weight, prob.velocity_weight, prob.control_weight, prob.discount,
        prob.control_bound, prob.constraint, prob.riccati)

# The problem on the coarse grid (GridMultigrid's coarsening): a halved
# dimension of even count keeps the block centres (lower + h/2, spacing 2h),
# one of odd count every other point (same ends, spacing 2h); either way the
# coarse grid is uniform, so it is again a StateGrid.
function _coarsen(prob)
    g = _state_grid(prob)
    r = coarsening_factors(g.points, _periodic(prob))
    halved = map((rk, n) -> rk == 2 && iseven(n), r, g.points)
    lower = map((c, lo, h) -> c ? lo + h / 2 : lo, halved, g.lower, g.spacing)
    upper = map((c, hi, h) -> c ? hi - h / 2 : hi, halved, g.upper, g.spacing)
    return _with_grid(prob, StateGrid(lower, upper, coarse_points(g.points, r)))
end

# The policy pushed forward to the coarse grid: the average over each block.
@kernel function _restrict_controls_kernel!(Uc, @Const(U), r, nf, ::Val{D}, ::Val{d}) where {D,d}
    Ic = @index(Global, Cartesian)
    children = _children(Ic, r, nf, Val(D))
    for c in 1:d
        acc = 0.0
        for I in children
            acc += U[I, c]
        end
        Uc[Ic, c] = acc / length(children)
    end
end

# ---------------------------------------------------------------------------
# the coarsest solve
# ---------------------------------------------------------------------------

# Policy evaluation on the coarsest grid as the fixpoint V = P V + c, P = D⁻¹N
# (D the diagonal of the stencil, N its couplings), c = D⁻¹b. The pattern of P
# is the grid graph with every neighbour stored, zero or not, so it is the same
# for every policy (and strongly connected, as PathPack's GPU solver needs).
struct _CoarsestSystem{D,S,R}
    points::NTuple{D,Int}
    P::SparseMatrixCSC{Float64,Int}
    slot::Matrix{Int}                  # N × 2D: index into P.nzval of each coupling, 0 if none
    diagonal::Vector{Float64}
    coef::Array{Float64,R}             # the stencil coefficients on the host, (points..., 2D + 1)
    rhs::Vector{Float64}
    solution::Vector{Float64}
    solver::S
end

function _coarsest_pattern(points::NTuple{D,Int}, periodic::NTuple{D,Bool}) where {D}
    L = LinearIndices(points)
    I, J = Int[], Int[]
    neighbours = zeros(Int, length(L), 2D)
    for C in CartesianIndices(points), j in 1:D, (s, side) in ((1, -1), (2, 1))
        k = C[j] + side
        if !(1 <= k <= points[j])
            periodic[j] || continue
            k = mod1(k, points[j])
        end
        neighbours[L[C], (s - 1) * D + j] = L[Base.setindex(Tuple(C), k, j)...]
        push!(I, L[C]); push!(J, neighbours[L[C], (s - 1) * D + j])
    end
    P = sparse(I, J, zeros(length(I)), length(L), length(L), (a, b) -> a)
    slot = zeros(Int, length(L), 2D)
    rows = rowvals(P)
    for col in 1:length(L), p in nzrange(P, col)
        i = rows[p]
        for c in 1:(2D)
            neighbours[i, c] == col && (slot[i, c] = p)
        end
    end
    return P, slot
end

"""
    _coarsest_solver(s::CoarsestSolver, P::SparseMatrixCSC) -> solver

A solver of ``V = P V + c`` for matrices with the pattern of `P`, refactored by
`_refactor!(solver, P)` and applied by `_solve!(solver, x, c)`. `DenseCoarsest`
here; `PathPackCoarsest` in the PATHPACK package extension.
"""
function _coarsest_solver end

mutable struct _DenseCoarsest
    factor::LU{Float64,Matrix{Float64},Vector{Int}}
    matrix::Matrix{Float64}
end

_coarsest_solver(::DenseCoarsest, P::SparseMatrixCSC) =
    _DenseCoarsest(lu(Matrix{Float64}(I, size(P)...)), zeros(size(P)))

function _refactor!(s::_DenseCoarsest, P::SparseMatrixCSC)
    s.matrix .= .-P
    for i in axes(s.matrix, 1)
        s.matrix[i, i] += 1
    end
    s.factor = lu!(s.matrix)                     # I - P, an M-matrix: no pivoting needed, LAPACK's is harmless
    return s
end

_solve!(s::_DenseCoarsest, x::Vector{Float64}, c::Vector{Float64}) = ldiv!(x, s.factor, c)

_coarsest_solver(::PathPackCoarsest, ::AbstractMatrix) =
    throw(ArgumentError("PathPackCoarsest needs PathPack: run `using PATHPACK` first"))
_coarsest_solver(s::CoarsestSolver, ::AbstractMatrix) =
    throw(ArgumentError("$(typeof(s)) does not implement _coarsest_solver; see the CoarsestSolver docstring"))

function _CoarsestSystem(kind::CoarsestSolver, points::NTuple{D,Int}, periodic::NTuple{D,Bool}) where {D}
    P, slot = _coarsest_pattern(points, periodic)
    N = prod(points)
    solver = _coarsest_solver(kind, P)
    return _CoarsestSystem{D,typeof(solver),D + 1}(points, P, slot, zeros(N), zeros(points..., 2D + 1), zeros(N), zeros(N),
        solver)
end

# Load the coefficients (already in sys.coef) into P and refactorize.
function _refactor!(sys::_CoarsestSystem{D}) where {D}
    coef = reshape(sys.coef, :, 2D + 1)
    nz = nonzeros(sys.P)
    for i in axes(coef, 1)
        d = coef[i, 1]
        sys.diagonal[i] = d
        for c in 1:(2D)
            p = sys.slot[i, c]
            p == 0 || (nz[p] = coef[i, 1 + c] / d)
        end
    end
    _refactor!(sys.solver, sys.P)
    return sys
end

# ---------------------------------------------------------------------------
# the hierarchy
# ---------------------------------------------------------------------------

struct _Level{O,L,X,A,U,C}
    op::O
    layout::L
    x::X                                # the level's iterate and right-hand side (the fine level uses BiCGStab's)
    b::X
    r::A
    e::A
    controls::U                         # rediscretization: this level's policy
    coef::C                             # Galerkin: this level's stored coefficients
end

# The fine level (which borrows BiCGStab's vectors) and the coarse levels, which
# share one concrete type, so that a V-cycle is type-stable at every depth.
struct _Hierarchy{F<:_Level,L<:_Level,S,B,C}
    fine::F
    coarse::Vector{L}                   # levels 2, 3, …; never empty
    coarsest::S
    coarsest_buffer::B                  # contiguous device copy of the coarsest interior
    coarsest_coef::C                    # the coarsest stencil coefficients on the device
    galerkin::Bool
    presmooth::Int
    postsmooth::Int
end

# The number of levels of the hierarchy of a grid of `points` (fine grid
# included): `mg.levels`, checked, or with `nothing` as many as it takes to bring
# the coarsest grid down to mg.coarsest_size points or to a grid that cannot
# be coarsened.
function _plan_levels(mg::Multigrid, points::NTuple{D,Int}, periodic::NTuple{D,Bool}) where {D}
    @argcheck mg.presmooth >= 1 && mg.postsmooth >= 0
    @argcheck mg.levels === nothing || mg.levels >= 2 "levels must be nothing (automatic) or an integer ≥ 2"
    @argcheck mg.coarsest_size >= 1
    levels, n = 1, points
    while true
        r = coarsening_factors(n, periodic)
        done = mg.levels === nothing ? (prod(n) <= mg.coarsest_size || all(==(1), r)) : levels == mg.levels
        done && break
        @argcheck any(==(2), r) "level $levels ($(n) points) cannot be coarsened: ask for fewer levels or use levels = nothing"
        n = coarse_points(n, r)
        levels += 1
    end
    @argcheck levels >= 2 "the grid ($(points) points) is already at most coarsest_size = $(mg.coarsest_size) points"
    return levels
end

function _Hierarchy(mg::Multigrid, prob, op, layout, U, backend)
    D = length(gridsize(op))
    d = size(U, D + 1)
    nlevels = _plan_levels(mg, gridsize(op), layout.periodic)
    galerkin = mg.coarse_operator isa GalerkinOperator
    first_level = _Level(op, layout, nothing, nothing, grid_zeros(op), grid_zeros(op), U, nothing)
    levels = _Level[]
    problem = prob
    for _ in 2:nlevels
        problem = _coarsen(problem)
        layout = coarse_layout(layout)
        n = layout.points
        controls = galerkin ? nothing : KernelAbstractions.zeros(backend, Float64, n..., d)
        coef = galerkin ? KernelAbstractions.zeros(backend, Float64, n..., 2D + 1) : nothing
        stencil = galerkin ? CoefficientStencil(coef) :
                  _UpwindStencil(controls, _KernelProblem(problem, ntuple(_ -> 0, D)))
        opc = box_operator(layout, stencil)
        push!(levels, _Level(opc, layout, grid_zeros(opc), grid_zeros(opc), grid_zeros(opc), grid_zeros(opc),
            controls, coef))
    end
    coarse = [level for level in levels]                # narrowed to the coarse levels' common type
    coarsest_layout = coarse[end].layout
    system = _CoarsestSystem(mg.coarsest, coarsest_layout.points, coarsest_layout.periodic)
    buffer = KernelAbstractions.zeros(backend, Float64, prod(coarsest_layout.points))
    coarsest_coef = KernelAbstractions.zeros(backend, Float64, coarsest_layout.points..., 2D + 1)
    return _Hierarchy(first_level, coarse, system, buffer, coarsest_coef, galerkin, mg.presmooth, mg.postsmooth)
end

# For a new policy (the fine controls changed): the coarse operators, then the
# coarsest factorization (symbolic analysis kept).
function _update!(H::_Hierarchy)
    _update_level!(H, H.fine, H.coarse[1])
    for k in 2:length(H.coarse)
        _update_level!(H, H.coarse[k - 1], H.coarse[k])
    end
    last = H.coarse[end]
    if H.galerkin
        copyto!(H.coarsest.coef, last.coef)
    else
        copyto!(H.coarsest.coef, stencil_coefficients!(H.coarsest_coef, last.op))
    end
    _refactor!(H.coarsest)
    return H
end

# The operator of `coarse` from that of the next finer level.
function _update_level!(H::_Hierarchy, fine::_Level, coarse::_Level)
    D = length(gridsize(fine.op))
    if H.galerkin
        galerkin_coefficients!(coarse.coef, fine.op)
    else
        backend = get_backend(coarse.op)
        d = size(coarse.controls, D + 1)
        _restrict_controls_kernel!(backend)(coarse.controls, fine.controls,
            ntuple(k -> gridsize(coarse.op, k) == gridsize(fine.op, k) ? 1 : 2, D), gridsize(fine.op), Val(D), Val(d);
            ndrange=gridsize(coarse.op))
        KernelAbstractions.synchronize(backend)
    end
    return coarse
end

_exchange!(level::_Level, x) = exchange!(SerialBoxes(), level.layout, x)

function _residual!(r, level::_Level, b, x)
    _exchange!(level, x)
    apply!(r, level.op, x)
    _lincomb!(r, level.op, 1, b, -1, r, 0, r)
    return r
end

# Smooth A x = b by red–black SGS sweeps in correction form, x updated.
function _smooth!(level::_Level, x, b, sweeps)
    for _ in 1:sweeps
        _residual!(level.r, level, b, x)
        red_black_sgs!(level.e, level.op, level.r; exchange=y -> _exchange!(level, y))
        _lincomb!(x, level.op, 1, x, 1, level.e, 0, x)
    end
    return x
end

function _coarsest_solve!(H::_Hierarchy, level::_Level)
    sys = H.coarsest
    reshape(H.coarsest_buffer, gridsize(level.op)) .= interior(level.b, level.op)
    copyto!(sys.rhs, H.coarsest_buffer)
    sys.rhs ./= sys.diagonal
    _solve!(sys.solver, sys.solution, sys.rhs)
    copyto!(H.coarsest_buffer, sys.solution)
    fill!(level.x, 0)
    interior(level.x, level.op) .= reshape(H.coarsest_buffer, gridsize(level.op))
    return level.x
end

# One V-cycle for A_ℓ x = b from x = 0, x and b padded arrays of level ℓ.
_vcycle!(H::_Hierarchy, ℓ::Int, x, b) =
    ℓ == 1 ? _vcycle!(H, H.fine, 1, x, b) : _vcycle!(H, H.coarse[ℓ - 1], ℓ, x, b)

function _vcycle!(H::_Hierarchy, level::_Level, ℓ::Int, x, b)
    ℓ == length(H.coarse) + 1 && return _coarsest_solve!(H, level)
    red_black_sgs!(x, level.op, b; exchange=y -> _exchange!(level, y))            # x = S b
    _smooth!(level, x, b, H.presmooth - 1)
    _residual!(level.r, level, b, x)
    coarse = H.coarse[ℓ]                                                          # level ℓ + 1
    restrict_average!(coarse.b, coarse.op, level.r, level.op)                     # pushforward T
    _vcycle!(H, coarse, ℓ + 1, coarse.x, coarse.b)
    prolong_add!(x, level.op, coarse.x, coarse.op)                                # pullback E
    _smooth!(level, x, b, H.postsmooth)
    return x
end
