"""
Module for trajectory-valued sheaves encoding linear dynamical system state
evolution, together with a colocation boundary-value solver via harmonic
extension.

A `TrajectorySheaf` wraps an `EuclideanSheaf` on a path graph of length `k+1`
whose restriction maps encode the autonomous discrete-time dynamics
``x_{t+1} = A x_t``.  Global sections of the inner sheaf correspond exactly
to valid `k`-step trajectories.

This module provides:

- [`TrajectorySheaf`](@ref) — struct holding `k`, `A`, and the inner sheaf.
- [`build_trajectory_sheaf`](@ref) — construct the path-graph sheaf for given
  `F`, `A`, and `k`.
- [`colocation_trajectory`](@ref) — solve the two-point boundary-value problem
  ``x_0 = x_{\\text{given}}``, ``x_k = x_{\\text{given}}`` via harmonic
  extension.
"""
module TrajectorySheaves

export TrajectorySheaf, build_trajectory_sheaf, colocation_trajectory

using ArgCheck: @argcheck
using LinearAlgebra
using BlockArrays

using ..EuclideanSheaves: EuclideanSheaf, add_sheaf_edge!, harmonic_extension
using ..SheafInterface: vertex_stalks

# ---------------------------------------------------------------------------
# Struct
# ---------------------------------------------------------------------------

"""
    TrajectorySheaf{T}

A trajectory-valued sheaf encoding `k` steps of the discrete-time linear
dynamics ``x_{t+1} = A x_t``.

Fields:
- `k  :: Int`                  — number of time steps.
- `A  :: AbstractMatrix{T}`    — state-transition matrix (``d \\times d``).
- `sheaf :: EuclideanSheaf{T}` — inner path-graph sheaf; its global sections
  are exactly the valid trajectories ``[x_0, x_1, \\ldots, x_k]``.
"""
struct TrajectorySheaf{T}
    k::Int
    A::AbstractMatrix{T}
    sheaf::EuclideanSheaf{T}
end

# ---------------------------------------------------------------------------
# Construction helpers
# ---------------------------------------------------------------------------

"""
    build_trajectory_sheaf(F::EuclideanSheaf{T}, A::AbstractMatrix{T}, k::Int)
        -> EuclideanSheaf{T}

Construct the inner path-graph sheaf for a `TrajectorySheaf`.

The state dimension `d` equals the total dimension of the 0-cochain space of
`F` (i.e. `sum(vertex_stalks(F))`).  The returned sheaf has `k+1` vertices,
each with stalk ``\\mathbb{R}^d``, connected by `k` edges with an induced
orientation in the coboundary map. For each edge between vertices `t` and
`t+1`, using the induced orientation ``t \\to t+1``:

- The restriction map from vertex `t`   is ``A``    (dynamics matrix).
- The restriction map from vertex `t+1` is ``I_d``  (identity).

A 0-cochain ``[x_0, \\ldots, x_k]`` is a global section iff
``A x_t - x_{t+1} = 0``, i.e. ``x_{t+1} = A x_t`` (forward dynamics),

# Arguments
- `F` — base `EuclideanSheaf` whose total 0-cochain dimension gives `d`.
- `A` — ``d \\times d`` state-transition matrix.
- `k` — number of time steps (must be ≥ 1).

# Throws
- `ArgumentError` if `k < 1`.
- `ArgumentError` if `A` is not square.
- `ArgumentError` if `size(A, 1) ≠ d`.
"""
function build_trajectory_sheaf(F::EuclideanSheaf{T}, A::AbstractMatrix{T}, k::Int) where T
    d = sum(vertex_stalks(F))

    @argcheck k >= 1 "k must be at least 1, got $k"
    @argcheck size(A, 1) == size(A, 2) "A must be a square matrix, got size $(size(A))"
    @argcheck size(A, 1) == d "A must be $d×$d (total 0-cochain dimension of F), got $(size(A, 1))×$(size(A, 2))"

    sheaf = EuclideanSheaf{T}(fill(d, k + 1))
    Id = Matrix{T}(I, d, d)
    Am = Matrix{T}(A)
    for t in 1:k
        # Edge (t, t+1):
        #   restriction from xₜ   → edge stalk: A   (so d_e x = A*x_t - I*x_{t+1})
        #   restriction from xₜ₊₁ → edge stalk: I
        # A 0-cochain is a global section iff A * x_t - I * x_{t+1} = 0,
        # i.e. x_{t+1} = A * x_t (forward dynamics).
        add_sheaf_edge!(sheaf, t, t + 1, Am, Id)
    end
    return sheaf
end

"""
    TrajectorySheaf(F::EuclideanSheaf{T}, A::AbstractMatrix{T}, k::Int)
        -> TrajectorySheaf{T}

Construct a `TrajectorySheaf` from a base sheaf `F`, dynamics matrix `A`, and
number of steps `k`.

The inner `sheaf` is built eagerly by [`build_trajectory_sheaf`](@ref).

# Throws
- `ArgumentError` if `k < 1`.
- `ArgumentError` if `A` is not square or has the wrong size relative to `F`.
"""
function TrajectorySheaf(F::EuclideanSheaf{T}, A::AbstractMatrix{T}, k::Int) where T
    sheaf = build_trajectory_sheaf(F, A, k)
    return TrajectorySheaf{T}(k, A, sheaf)
end

# ---------------------------------------------------------------------------
# Colocation
# ---------------------------------------------------------------------------

"""
    colocation_trajectory(tsheaf::TrajectorySheaf{T}, x0, xk)
        -> BlockVector

Solve the two-point boundary-value problem for the trajectory sheaf by
harmonic extension.

Given initial state `x0` (assigned to time step 1) and final state `xk`
(assigned to time step `k+1`), returns the minimum-energy 0-cochain of
`tsheaf.sheaf` with those boundary values. If the boundary data are
consistent with the dynamics over `k` steps (that is, `xk ≈ A^k * x0`),
this coincides with the exact trajectory determined by those endpoints.
If the boundary data are inconsistent, there is no exact trajectory with
those endpoints; in that case the returned value is the harmonic
minimum-energy / least-squares 0-cochain and generally will not satisfy
`x[t+1] = A * x[t]` exactly at every step.

The returned `BlockVector` has `k+1` blocks where `trajectory[Block(t)]` gives
the state at time step `t-1` (so `Block(1)` is `x_0`, `Block(2)` is `x_1`,
…, `Block(k+1)` is `x_k`).

# Arguments
- `tsheaf` — a `TrajectorySheaf{T}`.
- `x0`     — initial state vector of length `d = size(tsheaf.A, 1)`.
- `xk`     — final state vector of length `d`.

# Throws
- `ArgumentError` if `length(x0) ≠ d` or `length(xk) ≠ d`.
"""
function colocation_trajectory(tsheaf::TrajectorySheaf{T}, x0, xk) where T
    d = size(tsheaf.A, 1)
    @argcheck length(x0) == d "x0 must have length $d, got $(length(x0))"
    @argcheck length(xk) == d "xk must have length $d, got $(length(xk))"

    boundary = Dict{Int,Vector{T}}(
        1          => convert(Vector{T}, x0),
        tsheaf.k + 1 => convert(Vector{T}, xk),
    )
    x_p, _null_basis = harmonic_extension(tsheaf.sheaf, boundary)
    return x_p
end

end # TrajectorySheaves
