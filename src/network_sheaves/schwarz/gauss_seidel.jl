# Inexact local solves for SchwarzMethods: symmetric Gauss–Seidel sweeps of the
# local matrices, in natural or multicolor order. Included in SchwarzMethods.

export SymmetricGaussSeidelLocalSolve, SymmetricGaussSeidel, SweepOrdering, NaturalOrdering, MulticolorOrdering

"""
    SweepOrdering

The order in which a Gauss–Seidel sweep visits the unknowns:
[`NaturalOrdering`](@ref) or [`MulticolorOrdering`](@ref).
"""
abstract type SweepOrdering end

"""
    NaturalOrdering()

Sweep the unknowns in index order, column by column of the matrix: inherently
sequential.
"""
struct NaturalOrdering <: SweepOrdering end

"""
    MulticolorOrdering(colors = nothing)

Sweep the unknowns color by color for a coloring of the graph of
``A + A^\\mathsf{T}``: `colors` is a vector of index vectors, or `nothing` for a
greedy coloring in index order. Unknowns of one color are not coupled, so each
color can be updated concurrently. On a grid whose stencil couples only axis
neighbours the greedy coloring is the red–black (checkerboard) one.
"""
struct MulticolorOrdering <: SweepOrdering
    colors::Union{Nothing,Vector{Vector{Int}}}
end
MulticolorOrdering() = MulticolorOrdering(nothing)
MulticolorOrdering(colors::AbstractVector{<:AbstractVector{<:Integer}}) =
    MulticolorOrdering([Vector{Int}(c) for c in colors])

"""
    SymmetricGaussSeidelLocalSolve(; sweeps = 1, ordering = MulticolorOrdering())

Solve each local problem inexactly with `sweeps` symmetric Gauss–Seidel sweeps
of ``\\tilde A_i = D + L + U`` ([`SymmetricGaussSeidel`](@ref), in the given
`ordering`), i.e. the stationary iteration with
``M_i^{-1} = (D + U)^{-1} D (D + L)^{-1}``, applied in correction form
``x_i \\leftarrow x_i + M_i^{-1}(b_i - \\tilde A_i x_i)`` from the current
iterate. Setup needs no factorization. For a nonsingular M-matrix the
Gauss–Seidel splitting is regular, so the Schwarz methods with inexact local
solves still converge (Frommer and Szyld, *Weighted max norm estimates for
additive Schwarz methods*, 1999). With `sweeps = 1` each subdomain does
exactly the symmetric Gauss–Seidel pass a global point smoother would.
"""
struct SymmetricGaussSeidelLocalSolve{O<:SweepOrdering} <: LocalSolver
    sweeps::Int
    ordering::O
    function SymmetricGaussSeidelLocalSolve(sweeps::Integer, ordering::O=MulticolorOrdering()) where {O<:SweepOrdering}
        @argcheck sweeps >= 1 "need at least one sweep"
        @argcheck !(ordering isa MulticolorOrdering) || ordering.colors === nothing "a local solver colors each local matrix itself"
        return new{O}(Int(sweeps), ordering)
    end
end
SymmetricGaussSeidelLocalSolve(; sweeps::Integer=1, ordering::SweepOrdering=MulticolorOrdering()) =
    SymmetricGaussSeidelLocalSolve(sweeps, ordering)

_local_factor(s::SymmetricGaussSeidelLocalSolve, Ai::SparseMatrixCSC, symmetric::Bool, i::Int) =
    SymmetricGaussSeidel(Ai; sweeps=s.sweeps, ordering=s.ordering, threaded=false)

"""
    SymmetricGaussSeidel(A::SparseMatrixCSC; sweeps = 1, ordering = MulticolorOrdering(),
                         threaded = size(A, 1) >= 50_000)

The symmetric Gauss–Seidel preconditioner of ``A = D + L + U`` (``L``, ``U``
strictly lower and upper in the chosen ordering of the unknowns):
`ldiv!(y, M, r)` runs `sweeps` sweeps of the stationary iteration
``y \\leftarrow y + (D + U)^{-1} D (D + L)^{-1}(r - A y)`` from ``y = 0``; one
sweep applies ``(D + U)^{-1} D (D + L)^{-1}``, a forward Gauss–Seidel pass
followed by a backward one. The diagonal must be stored and nonzero.

- [`NaturalOrdering`](@ref) sweeps the unknowns in index order, column by column
  of `A`: inherently sequential.
- [`MulticolorOrdering`](@ref) sweeps them color by color, each color updated
  concurrently (when `threaded`), on CPU threads or a GPU alike.

New values of `A` on the same pattern take effect after [`refactor!`](@ref)
gathers them (or `nonzeros(M.A)` is overwritten followed by `_refresh!(M)`),
without any other setup.
"""
struct SymmetricGaussSeidel{T}
    A::SparseMatrixCSC{T,Int}
    diagonal::Vector{Int}           # position of A[j, j] in nonzeros(A)
    sweeps::Int
    rows::SparseMatrixCSC{T,Int}    # Aᵀ, i.e. A by rows (multicolor only)
    transpose_map::Vector{Int}      # nonzeros(rows) == nonzeros(A)[transpose_map]
    row_diagonal::Vector{Int}       # position of A[i, i] in nonzeros(rows)
    colors::Vector{Vector{Int}}     # empty for the natural ordering
    threaded::Bool
end

function SymmetricGaussSeidel(A::SparseMatrixCSC; sweeps::Integer=1, ordering::SweepOrdering=MulticolorOrdering(),
                              threaded::Bool=size(A, 1) >= 50_000)
    @argcheck size(A, 1) == size(A, 2) "A must be square"
    @argcheck sweeps >= 1 "need at least one sweep"
    S = _float_sparse(A)
    diagonal = _diagonal_positions(S)
    return _symmetric_gauss_seidel(S, diagonal, Int(sweeps), ordering, threaded)
end

function _symmetric_gauss_seidel(S, diagonal, sweeps, ::NaturalOrdering, threaded)
    empty = SparseMatrixCSC(0, 0, [1], Int[], eltype(S)[])
    return SymmetricGaussSeidel(S, diagonal, sweeps, empty, Int[], Int[], Vector{Int}[], threaded)
end

function _symmetric_gauss_seidel(S, diagonal, sweeps, ordering::MulticolorOrdering, threaded)
    n = size(S, 1)
    groups = ordering.colors === nothing ? _matrix_coloring(S) : ordering.colors
    _check_coloring(S, groups)
    index = SparseMatrixCSC(n, n, S.colptr, rowvals(S), collect(1:nnz(S)))
    transposed = copy(transpose(index))
    transpose_map = nonzeros(transposed)
    rows = SparseMatrixCSC(n, n, transposed.colptr, rowvals(transposed), nonzeros(S)[transpose_map])
    return SymmetricGaussSeidel(S, diagonal, sweeps, rows, transpose_map, _diagonal_positions(rows),
        groups, threaded)
end

function _diagonal_positions(S::SparseMatrixCSC)
    return map(1:size(S, 2)) do j
        p = findfirst(q -> rowvals(S)[q] == j, nzrange(S, j))
        @argcheck p !== nothing && !iszero(nonzeros(S)[nzrange(S, j)[p]]) "the diagonal entry $j is zero"
        nzrange(S, j)[p]
    end
end

# Greedy coloring in index order of the graph of S + Sᵀ: each unknown takes the
# smallest color none of its earlier neighbours has.
function _matrix_coloring(S::SparseMatrixCSC)
    P = _structure(S)
    n = size(P, 1)
    color = zeros(Int, n)
    stamp = Int[]
    for i in 1:n
        for p in nzrange(P, i)
            j = rowvals(P)[p]
            (j != i && color[j] > 0) && (stamp[color[j]] = i)
        end
        c = findfirst(!=(i), stamp)
        c === nothing && (push!(stamp, 0); c = length(stamp))
        color[i] = c
    end
    groups = [Int[] for _ in 1:maximum(color; init=0)]
    for i in 1:n
        push!(groups[color[i]], i)
    end
    return groups
end

function _check_coloring(S::SparseMatrixCSC, groups::Vector{Vector{Int}})
    n = size(S, 1)
    color = zeros(Int, n)
    for (c, g) in enumerate(groups), i in g
        @argcheck 1 <= i <= n && color[i] == 0 "the colors must partition 1:$n"
        color[i] = c
    end
    @argcheck all(>(0), color) "the colors must partition 1:$n"
    for j in 1:n, p in nzrange(S, j)
        i = rowvals(S)[p]
        @argcheck i == j || color[i] != color[j] "unknowns $i and $j are coupled but have the same color"
    end
end

# Bring the row copy up to date after new values were written into nonzeros(M.A).
function _refresh!(M::SymmetricGaussSeidel)
    isempty(M.colors) || (nonzeros(M.rows) .= view(nonzeros(M.A), M.transpose_map))
    return M
end

Base.size(M::SymmetricGaussSeidel) = size(M.A)
Base.size(M::SymmetricGaussSeidel, d::Integer) = size(M.A, d)
Base.eltype(::SymmetricGaussSeidel{T}) where {T} = T

_with_sweeps(M::SymmetricGaussSeidel, s::Int) = SymmetricGaussSeidel(M.A, M.diagonal, s, M.rows,
    M.transpose_map, M.row_diagonal, M.colors, M.threaded)

# One sweep from zero: y = (D + U)⁻¹ D (D + L)⁻¹ r.
function _sgs_sweep!(y::AbstractVector, M::SymmetricGaussSeidel, r::AbstractVector)
    isempty(M.colors) || return _multicolor_sweep!(y, M, r)
    A, d = M.A, M.diagonal
    rows, vals = rowvals(A), nonzeros(A)
    copyto!(y, r)
    @inbounds for j in 1:size(A, 2)
        yj = y[j] /= vals[d[j]]
        for p in nzrange(A, j)
            i = rows[p]
            i > j && (y[i] -= vals[p] * yj)
        end
    end
    @inbounds for j in 1:size(A, 2)
        y[j] *= vals[d[j]]
    end
    @inbounds for j in size(A, 2):-1:1
        yj = y[j] /= vals[d[j]]
        for p in nzrange(A, j)
            i = rows[p]
            i < j && (y[i] -= vals[p] * yj)
        end
    end
    return y
end

# Forward pass over the colors from y = 0, then a backward pass in reverse color
# order, each unknown updated from its whole row with the latest values:
# exactly (D + U)⁻¹ D (D + L)⁻¹ r in the color-by-color ordering.
function _multicolor_sweep!(y::AbstractVector, M::SymmetricGaussSeidel, r::AbstractVector)
    fill!(y, zero(eltype(y)))
    for c in M.colors
        _color_pass!(y, M, r, c)
    end
    for c in Iterators.drop(Iterators.reverse(M.colors), 1)   # the last color would only repeat itself
        _color_pass!(y, M, r, c)
    end
    return y
end

function _color_pass!(y::AbstractVector, M::SymmetricGaussSeidel, r::AbstractVector, color::Vector{Int})
    R, d = M.rows, M.row_diagonal
    cols, vals = rowvals(R), nonzeros(R)
    @inline function update!(i)
        acc = r[i]
        @inbounds for p in nzrange(R, i)
            j = cols[p]
            j != i && (acc -= vals[p] * y[j])
        end
        @inbounds y[i] = acc / vals[d[i]]
    end
    if M.threaded
        Threads.@threads for t in eachindex(color)
            update!(@inbounds color[t])
        end
    else
        for i in color
            update!(i)
        end
    end
    return y
end

# `sweeps` sweeps of the stationary iteration, correcting y in place.
function _sgs_correct!(y::AbstractVector, M::SymmetricGaussSeidel, b::AbstractVector)
    r, z = similar(y), similar(y)
    for _ in 1:M.sweeps
        isempty(M.colors) || !M.threaded ? mul!(r, M.A, y) : _threaded_mul!(r, M.rows, y)
        r .= b .- r
        y .+= _sgs_sweep!(z, M, r)
    end
    return y
end

function LinearAlgebra.ldiv!(y::AbstractVector, M::SymmetricGaussSeidel, r::AbstractVector)
    @argcheck length(y) == length(r) == size(M, 1) "y and r must have one entry per unknown"
    @argcheck !Base.mightalias(y, r) "y and r must not alias; use ldiv!(M, r) in place"
    _sgs_sweep!(y, M, r)
    M.sweeps == 1 && return y
    return _sgs_correct!(y, _with_sweeps(M, M.sweeps - 1), r)
end

LinearAlgebra.ldiv!(M::SymmetricGaussSeidel, r::AbstractVector) = copyto!(r, ldiv!(similar(r), M, r))
Base.:\(M::SymmetricGaussSeidel, r::AbstractVector) = ldiv!(similar(r, promote_type(eltype(M), eltype(r))), M, r)

_factor_solve(M::SymmetricGaussSeidel, b::AbstractVector) = ldiv!(similar(b, eltype(M)), M, b)
_factor_solve!(z::AbstractVector, M::SymmetricGaussSeidel, b::AbstractVector) = ldiv!(z, M, b)
_local_update(M::SymmetricGaussSeidel, b::AbstractVector, x::AbstractVector) = _sgs_correct!(copy(x), M, b)

# refactor! for Gauss–Seidel local factors: a fresh smoother, or new values on its pattern.
_refactored(old::SymmetricGaussSeidel, solver, Ai, i) = SymmetricGaussSeidel(Ai; sweeps=old.sweeps,
    ordering=isempty(old.colors) ? NaturalOrdering() : MulticolorOrdering(old.colors), threaded=old.threaded)
_refactor_values!(dd, i, M::SymmetricGaussSeidel, values) = (nonzeros(M.A) .= values; _refresh!(M))
