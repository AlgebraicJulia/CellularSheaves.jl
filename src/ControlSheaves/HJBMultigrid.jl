# Geometric multigrid for the policy evaluations of GridPolicyIteration: a
# V-cycle over cell-centred coarse grids (GridMultigrid), preconditioning
# BiCGStab on the fine grid, with red–black SGS smoothing on every level but
# the coarsest, which is solved directly (dense LU, or PathPack's semiring LU
# with the PATHPACK package extension). Included in DoubleIntegratorHJB.

"""
    Multigrid(; levels = 3, coarse_operator = :rediscretize, coarsest = :pathpack,
              presmooth = 1, postsmooth = 1)

A V-cycle preconditioner for [`GridPolicyIteration`](@ref): pass it as
`GridPolicyIteration(multigrid = Multigrid(...))`. The fine grid and
`levels - 1` cell-centred coarse grids, each with half the points per
dimension (every dimension of every level but the coarsest needs an even
number of points), are linked by the pullback ``E`` (copy each coarse value
onto its block) and the pushforward transfer ``T`` (block average) of the
constant sheaf along the aggregation homomorphism (see
[`GridMultigrid`](@ref CellularSheaves.NetworkSheaves.GridMultigrid)).

- `coarse_operator = :rediscretize`: each coarse level is the same HJB problem
  discretized on the coarse grid, under the policy pushed forward to it (the
  block average of the controls, which stays in the convex control set);
  `:galerkin`: the composite ``T A E`` of the next finer level's operator.
  The two agree except on blocks straddling a switching surface of the drift.
- `coarsest`: `:pathpack` factorizes the coarsest grid's evaluation, as the
  closure ``V = P^* (D^{-1} b)`` of its discounted transition matrix
  ``P = D^{-1} N`` over the (+, ×) semiring, with PathPack's chordal LU
  (`using PATHPACK`: the symbolic analysis is done once and every policy
  refactorizes numerically); `:dense` uses LAPACK on a dense matrix (small
  coarsest grids).
- `presmooth`, `postsmooth`: red–black SGS sweeps before and after each
  coarse-grid correction.

Each application of the preconditioner is one V-cycle from a zero initial
guess. All work arrays are allocated once per solve. Single box only
(`SerialBoxes()`; CPU threads or one GPU, whose coarsest solve runs on the host).
"""
Base.@kwdef struct Multigrid
    levels::Int = 3
    coarse_operator::Symbol = :rediscretize
    coarsest::Symbol = :pathpack
    presmooth::Int = 1
    postsmooth::Int = 1
end

# ---------------------------------------------------------------------------
# problems on coarse grids
# ---------------------------------------------------------------------------

_with_grid(prob::HJBProblem, grid::StateGrid) =
    HJBProblem(grid, prob.axes, prob.position_weight, prob.velocity_weight, prob.control_weight, prob.discount,
        prob.control_bound, prob.constraint, prob.riccati)

# The problem on the cell-centred coarse grid: points at the block centres,
# spacing 2h.
function _coarsen(prob)
    g = _state_grid(prob)
    h = g.spacing
    return _with_grid(prob, StateGrid(g.lower .+ h ./ 2, g.upper .- h ./ 2, g.points .÷ 2))
end

@kernel function _restrict_controls_kernel!(Uc, @Const(U), ::Val{D}, ::Val{d}) where {D,d}
    Ic = @index(Global, Cartesian)
    for c in 1:d
        acc = 0.0
        for k in 0:(2^D - 1)
            acc += U[_child(Ic, k, Val(D)), c]
        end
        Uc[Ic, c] = acc / 2^D
    end
end

# ---------------------------------------------------------------------------
# the coarsest solve
# ---------------------------------------------------------------------------

# Policy evaluation on the coarsest grid as the fixpoint V = P V + c, P = D⁻¹N
# (D the diagonal of the stencil, N its couplings), c = D⁻¹b. The pattern of P
# is the grid graph with every neighbour stored, zero or not, so it is the same
# for every policy (and strongly connected, as PathPack's GPU solver needs).
struct _CoarsestSystem{D,S}
    points::NTuple{D,Int}
    P::SparseMatrixCSC{Float64,Int}
    slot::Matrix{Int}                  # N × 2D: index into P.nzval of each coupling, 0 if none
    diagonal::Vector{Float64}
    coef::Array{Float64}               # the stencil coefficients on the host, (points..., 2D + 1)
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
    _coarsest_solver(Val(kind), P::SparseMatrixCSC) -> solver

A solver of ``V = P V + c`` for matrices with the pattern of `P`, refactored by
`_refactor!(solver, P)` and applied by `_solve!(solver, x, c)`. `:dense`
here; `:pathpack` in the PATHPACK package extension.
"""
function _coarsest_solver end

mutable struct _DenseCoarsest
    factor::LU{Float64,Matrix{Float64},Vector{Int}}
    matrix::Matrix{Float64}
end

_coarsest_solver(::Val{:dense}, P::SparseMatrixCSC) =
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

function _coarsest_solver(::Val{kind}, ::SparseMatrixCSC) where {kind}
    kind === :pathpack && error("Multigrid(coarsest = :pathpack) needs PathPack: run `using PATHPACK` first")
    throw(ArgumentError("unknown coarsest solver $kind; use :pathpack or :dense"))
end

function _CoarsestSystem(kind::Symbol, points::NTuple{D,Int}, periodic::NTuple{D,Bool}) where {D}
    P, slot = _coarsest_pattern(points, periodic)
    N = prod(points)
    solver = _coarsest_solver(Val(kind), P)
    return _CoarsestSystem{D,typeof(solver)}(points, P, slot, zeros(N), zeros(points..., 2D + 1), zeros(N), zeros(N),
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

struct _Hierarchy{V,S,B,C}
    levels::V
    coarsest::S
    coarsest_buffer::B                  # contiguous device copy of the coarsest interior
    coarsest_coef::C                    # the coarsest stencil coefficients on the device
    galerkin::Bool
    presmooth::Int
    postsmooth::Int
end

function _check_multigrid(mg::Multigrid, points::NTuple{D,Int}) where {D}
    @argcheck mg.levels >= 2 "a multigrid hierarchy needs at least 2 levels"
    @argcheck mg.coarse_operator in (:rediscretize, :galerkin) "coarse_operator must be :rediscretize or :galerkin"
    @argcheck mg.coarsest in (:pathpack, :dense) "coarsest must be :pathpack or :dense"
    @argcheck mg.presmooth >= 1 && mg.postsmooth >= 0
    n = points
    for ℓ in 1:(mg.levels - 1)
        @argcheck all(iseven, n) "level $ℓ has $(n) points: every level but the coarsest needs an even number of points per dimension"
        n = n .÷ 2
    end
    @argcheck all(n .>= 2) "the coarsest grid ($(n) points) is too small"
end

function _Hierarchy(mg::Multigrid, prob, op, layout, U, backend)
    D = ndims(op)
    d = size(U, D + 1)
    _check_multigrid(mg, size(op))
    galerkin = mg.coarse_operator === :galerkin
    first_level = _Level(op, layout, nothing, nothing, grid_zeros(op), grid_zeros(op), U, nothing)
    levels = Any[first_level]
    problem = prob
    for _ in 2:mg.levels
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
    coarsest_layout = levels[end].layout
    system = _CoarsestSystem(mg.coarsest, coarsest_layout.points, coarsest_layout.periodic)
    buffer = KernelAbstractions.zeros(backend, Float64, prod(coarsest_layout.points))
    coarsest_coef = KernelAbstractions.zeros(backend, Float64, coarsest_layout.points..., 2D + 1)
    return _Hierarchy(Tuple(levels), system, buffer, coarsest_coef, galerkin, mg.presmooth, mg.postsmooth)
end

# For a new policy (the fine controls changed): the coarse operators, then the
# coarsest factorization (symbolic analysis kept).
function _update!(H::_Hierarchy)
    levels = H.levels
    D = ndims(levels[1].op)
    for ℓ in 2:length(levels)
        fine, coarse = levels[ℓ - 1], levels[ℓ]
        if H.galerkin
            galerkin_coefficients!(coarse.coef, fine.op)
        else
            backend = get_backend(coarse.op)
            d = size(coarse.controls, D + 1)
            _restrict_controls_kernel!(backend)(coarse.controls, fine.controls, Val(D), Val(d);
                ndrange=size(coarse.op))
            KernelAbstractions.synchronize(backend)
        end
    end
    last = levels[end]
    if H.galerkin
        copyto!(H.coarsest.coef, last.coef)
    else

        copyto!(H.coarsest.coef, stencil_coefficients!(H.coarsest_coef, last.op))
    end
    _refactor!(H.coarsest)
    return H
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
    reshape(H.coarsest_buffer, size(level.op)) .= interior(level.b, level.op)
    copyto!(sys.rhs, H.coarsest_buffer)
    sys.rhs ./= sys.diagonal
    _solve!(sys.solver, sys.solution, sys.rhs)
    copyto!(H.coarsest_buffer, sys.solution)
    fill!(level.x, 0)
    interior(level.x, level.op) .= reshape(H.coarsest_buffer, size(level.op))
    return level.x
end

# One V-cycle for A_ℓ x = b from x = 0, x and b padded arrays of level ℓ.
function _vcycle!(H::_Hierarchy, ℓ::Int, x, b)
    levels = H.levels
    level = levels[ℓ]
    ℓ == length(levels) && return _coarsest_solve!(H, level)
    red_black_sgs!(x, level.op, b; exchange=y -> _exchange!(level, y))            # x = S b
    _smooth!(level, x, b, H.presmooth - 1)
    _residual!(level.r, level, b, x)
    coarse = levels[ℓ + 1]
    restrict_average!(coarse.b, coarse.op, level.r, level.op)                     # pushforward T
    _vcycle!(H, ℓ + 1, coarse.x, coarse.b)
    prolong_add!(x, level.op, coarse.x, coarse.op)                                # pullback E
    _smooth!(level, x, b, H.postsmooth)
    return x
end
