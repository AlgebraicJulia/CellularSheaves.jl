"""
Module for trajectory-valued sheaves encoding linear dynamical system state
evolution, together with colocation boundary-value solvers via harmonic
extension.

This module provides two trajectory-sheaf constructions:

**Autonomous (`TrajectorySheaf`)**

A `TrajectorySheaf` wraps an `EuclideanSheaf` on a path graph of length `k+1`
whose restriction maps encode the autonomous discrete-time dynamics
``x_{t+1} = A x_t``.  Global sections of the inner sheaf correspond exactly
to valid `k`-step trajectories.

- [`TrajectorySheaf`](@ref) — struct holding `k`, `A`, and the inner sheaf.
- [`build_trajectory_sheaf`](@ref) — construct the path-graph sheaf for given
  `F`, `A`, and `k`.
- [`colocation_trajectory`](@ref) — solve the two-point boundary-value problem
  ``x_0 = x_{\\text{given}}``, ``x_k = x_{\\text{given}}`` via harmonic
  extension.

**Controlled (`ControlledTrajectorySheaf`)**

A `ControlledTrajectorySheaf` wraps an `EuclideanSheaf` on a path graph of
length `k+2` whose restriction maps encode the zero-order-hold discretization
of the continuous-time LTI system ``\\dot{x}(t) = A_c x(t) + B_c u(t)``.
Global sections of the inner sheaf correspond exactly to feasible sampled
state-control trajectories.

- [`continuous_to_discrete_zoh`](@ref) — discretize a continuous-time LTI
  system using the zero-order-hold block-exponential formula.
- [`ControlledTrajectorySheaf`](@ref) — struct and constructor; builds the
  inner path-graph sheaf eagerly.
- [`feasible_control_trajectory_basis`](@ref) — return a particular feasible
  trajectory and a null-space basis for endpoint-preserving perturbations.
"""
module TrajectorySheaves

export TrajectorySheaf, build_trajectory_sheaf, colocation_trajectory,
    continuous_to_discrete_zoh, ControlledTrajectorySheaf,
    feasible_control_trajectory_basis

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

# ---------------------------------------------------------------------------
# Zero-order-hold discretization
# ---------------------------------------------------------------------------

"""
    continuous_to_discrete_zoh(Ac::AbstractMatrix{T},
                               Bc::AbstractMatrix{T},
                               h::Real) where T
        -> (Ad::Matrix{T}, Bd::Matrix{T})

Discretize the continuous-time LTI system `ẋ(t) = Ac*x(t) + Bc*u(t)` using
zero-order hold on a uniform step size `h`.

The returned matrices satisfy `x[k+1] = Ad*x[k] + Bd*u[k]`.

The computation uses the block-matrix exponential identity

```math
\\exp\\!\\left(h \\begin{bmatrix} A_c & B_c \\\\ 0 & 0 \\end{bmatrix}\\right)
=
\\begin{bmatrix} A_d & B_d \\\\ 0 & I_m \\end{bmatrix}.
```

# Arguments
- `Ac` — ``n \\times n`` continuous-time state matrix.
- `Bc` — ``n \\times m`` continuous-time input matrix.
- `h`  — sample period (must be positive).

# Throws
- `ArgumentError` if `h <= 0`.
- `ArgumentError` if `Ac` is not square.
- `ArgumentError` if `size(Bc, 1) ≠ size(Ac, 1)`.
"""
function continuous_to_discrete_zoh(Ac::AbstractMatrix{T},
                                     Bc::AbstractMatrix{T},
                                     h::Real) where T
    n = size(Ac, 1)
    m = size(Bc, 2)
    @argcheck h > 0 "h must be positive, got $h"
    @argcheck size(Ac, 1) == size(Ac, 2) "Ac must be square, got size $(size(Ac))"
    @argcheck size(Bc, 1) == n "Bc must have $n rows (same as Ac), got $(size(Bc, 1))"

    # Build the (n+m)×(n+m) augmented matrix M = h*[Ac Bc; 0 0]
    M = zeros(T, n + m, n + m)
    M[1:n, 1:n]     .= h .* Ac
    M[1:n, n+1:n+m] .= h .* Bc
    # Bottom block remains zero

    E  = exp(M)
    Ad = E[1:n, 1:n]
    Bd = E[1:n, n+1:n+m]
    return Ad, Bd
end

# ---------------------------------------------------------------------------
# ControlledTrajectorySheaf struct
# ---------------------------------------------------------------------------

"""
    ControlledTrajectorySheaf{T}

Trajectory sheaf for the zero-order-hold discretization of the continuous-time
controlled system `ẋ(t) = Ac*x(t) + Bc*u(t)`.

The inner sheaf is defined on a path graph with `k+2` vertices:
- vertex 1 is a dummy initial vertex storing only `x₁`,
- vertices `2, …, k+1` store `(xₜ, uₜ)` for `t = 1, …, k`,
- vertex `k+2` stores only `x_{k+1}`.

Global sections of `sheaf` correspond exactly to feasible sampled state-control
trajectories of the discretized system.

Fields:
- `k          :: Int`                  — number of time steps.
- `h          :: T`                    — sample period.
- `Ac         :: AbstractMatrix{T}`    — continuous-time state matrix.
- `Bc         :: AbstractMatrix{T}`    — continuous-time input matrix.
- `Ad         :: Matrix{T}`            — discrete-time state matrix.
- `Bd         :: Matrix{T}`            — discrete-time input matrix.
- `sheaf      :: EuclideanSheaf{T}`    — inner path-graph sheaf.
- `state_dim  :: Int`                  — state dimension `n`.
- `control_dim :: Int`                 — control dimension `m`.
"""
struct ControlledTrajectorySheaf{T}
    k::Int
    h::T
    Ac::AbstractMatrix{T}
    Bc::AbstractMatrix{T}
    Ad::Matrix{T}
    Bd::Matrix{T}
    sheaf::EuclideanSheaf{T}
    state_dim::Int
    control_dim::Int
end

# ---------------------------------------------------------------------------
# Construction
# ---------------------------------------------------------------------------

"""
    ControlledTrajectorySheaf(F::EuclideanSheaf{T},
                              Ac::AbstractMatrix{T},
                              Bc::AbstractMatrix{T},
                              h::Real,
                              k::Int) where T
        -> ControlledTrajectorySheaf{T}

Construct a `ControlledTrajectorySheaf` from a base sheaf `F`, continuous-time
system matrices `(Ac, Bc)`, sample period `h`, and number of steps `k`.

The state dimension `n` is the total 0-cochain dimension of `F`
(i.e. `sum(vertex_stalks(F))`).  The inner sheaf has `k+2` vertices:

- vertex 1: stalk ``\\mathbb{R}^n`` (dummy initial state),
- vertices ``2, \\ldots, k+1``: stalk ``\\mathbb{R}^{n+m}`` (state and control),
- vertex ``k+2``: stalk ``\\mathbb{R}^n`` (terminal state).

The ``k+1`` edges enforce the zero-order-hold discrete dynamics:

- dummy initial edge ``(1,2)``: ``x_1^{\\mathrm{dummy}} = x_1^{\\mathrm{traj}}``,
- dynamics edges ``(t+1, t+2)`` for ``t = 1, \\ldots, k-1``:
  ``A_d x_t + B_d u_t = x_{t+1}``,
- final dynamics edge ``(k+1, k+2)``: ``A_d x_k + B_d u_k = x_{k+1}``.

# Throws
- `ArgumentError` if `k < 1`.
- `ArgumentError` if `h <= 0`.
- `ArgumentError` if `Ac` is not square.
- `ArgumentError` if `size(Bc, 1) ≠ sum(vertex_stalks(F))`.
"""
function ControlledTrajectorySheaf(F::EuclideanSheaf{T},
                                    Ac::AbstractMatrix{T},
                                    Bc::AbstractMatrix{T},
                                    h::Real,
                                    k::Int) where T
    n = sum(vertex_stalks(F))
    m = size(Bc, 2)

    @argcheck k >= 1 "k must be at least 1, got $k"
    @argcheck h > 0 "h must be positive, got $h"
    @argcheck size(Ac, 1) == size(Ac, 2) "Ac must be square, got size $(size(Ac))"
    @argcheck size(Bc, 1) == n "Bc must have $n rows (same state dimension as F), got $(size(Bc, 1))"

    Ad, Bd = continuous_to_discrete_zoh(Ac, Bc, h)

    # Build inner sheaf
    stalks = [n; fill(n + m, k); n]
    sheaf  = EuclideanSheaf{T}(stalks)

    In  = Matrix{T}(I, n, n)
    In0 = hcat(In, zeros(T, n, m))   # [I_n | 0_{n×m}] — projects out state
    AB  = hcat(Ad, Bd)               # [A_d | B_d]

    # Dummy initial edge (1, 2):
    #   ρ_{1→e} = I_n,   ρ_{2→e} = [I_n | 0]
    add_sheaf_edge!(sheaf, 1, 2, In, In0)

    # Dynamics edges (t+1, t+2) for t = 1, …, k-1:
    #   ρ_{t+1→e} = [A_d | B_d],  ρ_{t+2→e} = [I_n | 0]
    for t in 1:k-1
        add_sheaf_edge!(sheaf, t + 1, t + 2, AB, In0)
    end

    # Final dynamics edge (k+1, k+2):
    #   ρ_{k+1→e} = [A_d | B_d],  ρ_{k+2→e} = I_n
    add_sheaf_edge!(sheaf, k + 1, k + 2, AB, In)

    return ControlledTrajectorySheaf{T}(k, T(h), Ac, Bc, Ad, Bd, sheaf, n, m)
end

# ---------------------------------------------------------------------------
# Coordinate extraction
# ---------------------------------------------------------------------------

# Extract the public trajectory vector z = (x₁, …, x_{k+1}, u₁, …, u_k)
# from the internal 0-cochain of the controlled sheaf.
#
# Internal layout (0-indexed block offsets):
#   Block 1        : x₁ (n entries, dummy vertex)
#   Block t+1      : (xₜ, uₜ)  for t = 1, …, k  (n+m entries each)
#   Block k+2      : x_{k+1}   (n entries, terminal vertex)
#
# Public layout: (k+1)*n + k*m entries
#   Indices 1 : (k+1)*n                → states  x₁, x₂, …, x_{k+1}
#   Indices (k+1)*n+1 : (k+1)*n + k*m → controls u₁, u₂, …, u_k
function _internal_to_public(x_int::AbstractVector{T}, n::Int, m::Int, k::Int) where T
    z = Vector{T}(undef, (k + 1) * n + k * m)

    # x₁ from Block 1 (indices 1:n in internal)
    z[1:n] .= x_int[1:n]

    # xₜ for t = 2, …, k from Block t+1 (first n entries)
    for t in 2:k
        int_start = n + (t - 1) * (n + m) + 1
        pub_start = (t - 1) * n + 1
        z[pub_start:pub_start+n-1] .= x_int[int_start:int_start+n-1]
    end

    # x_{k+1} from Block k+2
    int_start = n + k * (n + m) + 1
    z[k*n+1:(k+1)*n] .= x_int[int_start:int_start+n-1]

    # uₜ for t = 1, …, k from Block t+1 (last m entries)
    for t in 1:k
        int_start = n + (t - 1) * (n + m) + n + 1
        pub_start = (k + 1) * n + (t - 1) * m + 1
        z[pub_start:pub_start+m-1] .= x_int[int_start:int_start+m-1]
    end

    return z
end

# Apply _internal_to_public column-wise to a matrix.
function _internal_to_public_matrix(N_int::AbstractMatrix{T}, n::Int, m::Int, k::Int) where T
    n_pub  = (k + 1) * n + k * m
    n_cols = size(N_int, 2)
    N_pub  = Matrix{T}(undef, n_pub, n_cols)
    for j in 1:n_cols
        N_pub[:, j] = _internal_to_public(view(N_int, :, j), n, m, k)
    end
    return N_pub
end

# ---------------------------------------------------------------------------
# Feasible trajectory basis
# ---------------------------------------------------------------------------

"""
    feasible_control_trajectory_basis(ts::ControlledTrajectorySheaf{T},
                                      x1::AbstractVector{T},
                                      xk1::AbstractVector{T}) where T
        -> (z_p::Vector{T}, null_basis::Matrix{T})

Return one feasible trajectory hitting the endpoint states together with a basis
matrix whose columns span all feasible perturbations that preserve those
endpoint states.

The boundary conditions fix vertex 1 to the initial state `x1` and vertex
`k+2` to the terminal state `xk1` of the inner sheaf.  The particular solution
and each null-basis column are then converted to the **public trajectory
coordinates**

```math
z = (x_1, x_2, \\dots, x_{k+1}, u_1, u_2, \\dots, u_k).
```

The complete feasible set is `{ z_p + null_basis * α : α ∈ Rᵈ }`.

# Arguments
- `ts`  — a `ControlledTrajectorySheaf{T}`.
- `x1`  — initial state vector of length `ts.state_dim`.
- `xk1` — terminal state vector of length `ts.state_dim`.

# Throws
- `ArgumentError` if `length(x1) ≠ ts.state_dim`.
- `ArgumentError` if `length(xk1) ≠ ts.state_dim`.
"""
function feasible_control_trajectory_basis(ts::ControlledTrajectorySheaf{T},
                                            x1::AbstractVector,
                                            xk1::AbstractVector) where T
    n = ts.state_dim
    m = ts.control_dim
    k = ts.k

    @argcheck length(x1) == n "x1 must have length $n, got $(length(x1))"
    @argcheck length(xk1) == n "xk1 must have length $n, got $(length(xk1))"

    boundary = Dict{Int,Vector{T}}(
        1       => convert(Vector{T}, x1),
        k + 2   => convert(Vector{T}, xk1),
    )

    x_p_internal, null_internal = harmonic_extension(ts.sheaf, boundary)

    z_p        = _internal_to_public(Array(x_p_internal), n, m, k)
    null_basis = _internal_to_public_matrix(null_internal, n, m, k)

    return z_p, null_basis
end

end # TrajectorySheaves
