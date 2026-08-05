# 002 — Controlled trajectory sheaf and feasible-trajectory affine space

## Mathematical Background

This issue extends the autonomous trajectory-sheaf construction to a continuous-time linear control system.

```math
\dot{x}(t) = A_c x(t) + B_c u(t),
```

The state dimension is `n` and the control dimension is `m`.

The intended discretization is zero-order hold on a uniform time grid, with the
control held constant on each sample interval. The discrete dynamics are

```math
x_{k+1} = A_d x_k + B_d u_k,
```

where `(A_d, B_d)` are computed from `(A_c, B_c, h)` by the standard block
matrix exponential formula

```math
\exp\!\left(
  h \begin{bmatrix}
    A_c & B_c \\
    0   & 0
  \end{bmatrix}
\right)
=
\begin{bmatrix}
  A_d & B_d \\
  0   & I_m
\end{bmatrix}.
```

The sampled trajectory coordinates should use Julia-style 1-based indexing and
list **all states first, then all controls**:

```math
z = (x_1, x_2, \dots, x_{k+1}, u_1, u_2, \dots, u_k).
```

Internally, the sheaf should use a path graph with `k+2` vertices:

- vertex 1 is a dummy initial vertex and stores only `x_1`,
- vertices `2, \dots, k+1` store `(x_t, u_t)` for `t = 1, \dots, k`,
- vertex `k+2` stores only `x_{k+1}`,
- each edge stalk has the same dimension as the state.

This adds only one extra copy of the state to the problem size: the initial
state appears once in the dummy boundary vertex and once in the first true
trajectory vertex. The advantage is that the harmonic-extension problem can
still impose whole-vertex boundary conditions while leaving the initial control
and terminal control unspecified.

The dummy initial edge should enforce equality of the duplicated initial state:

```math
x_1^{\mathrm{dummy}} - x_1^{\mathrm{traj}} = 0.
```

The remaining `k` edges should enforce the discretized dynamics:

```math
A_d x_t + B_d u_t - x_{t+1} = 0.
```

Concretely:

- for the dummy initial edge `e_0`,
  ```math
  \rho_{1 \to e_0} = I_n, \qquad
  \rho_{2 \to e_0} = [I_n \; 0],
  ```
  so the coboundary equation is
  ```math
  x_1^{\mathrm{dummy}} - x_1^{\mathrm{traj}} = 0;
  ```
- for dynamics edges `e_t` with `1 \le t \le k-1`,
  ```math
  \rho_{t+1 \to e_t} = [A_d \; B_d], \qquad
  \rho_{t+2 \to e_t} = [I_n \; 0],
  ```
  since vertices `t+1` and `t+2` store `(x_t, u_t)` and `(x_{t+1}, u_{t+1})`;
- for the final edge `e_k`,
  ```math
  \rho_{k+1 \to e_k} = [A_d \; B_d], \qquad
  \rho_{k+2 \to e_k} = I_n,
  ```
  since vertex `k+2` stores only `x_{k+1}`.

A global section of this sheaf is therefore exactly a feasible sampled
state-control trajectory for the discretized control system.

The boundary data for this issue should fix only the state endpoints: vertex 1
is fixed to `x_1` and vertex `k+2` is fixed to `x_{k+1}`. All controls remain
free. The desired output is the affine space

```math
\mathcal{T}(x_1, x_{k+1}) = z_p + N \alpha,
```

where:

- `z_p` is one feasible trajectory hitting the endpoints;
- the columns of `N` form a basis for the homogeneous feasible trajectories
  with zero endpoint-state variation.

This should reuse the same sheaf-Laplacian / harmonic-extension viewpoint
already used in `harmonic_extension` and `colocation_trajectory`: solve a
particular boundary-value problem on the internal sheaf coordinates and then
convert the result into the public state-first trajectory coordinates shown
above.

## Codebase Orientation

| File | Why it matters |
| --- | --- |
| `src/network_sheaves/TrajectorySheaf.jl` | Current autonomous discrete-time trajectory sheaf; the new controlled version should follow the same overall API shape. |
| `src/network_sheaves/EuclideanSheaves.jl` | Contains `harmonic_extension`, `nullspace_ldlt`, and the LDLt pseudoinverse/nullspace logic that should be reused or lightly generalized. |
| `src/network_sheaves/SheafInterface.jl` | If any new public accessors are introduced, declare them here first. |
| `src/network_sheaves/NetworkSheaves.jl` | Wire any new include / re-export here. |
| `test/network_sheaves/TrajectorySheaf.jl` | Existing tests for the autonomous trajectory sheaf give the expected style and invariants. |
| `docs/literate/trajectory_sheaf.jl` | Current colocation documentation; useful as the conceptual predecessor for the controlled version. |

## Requested Implementation

Extend `src/network_sheaves/TrajectorySheaf.jl` rather than creating a separate
top-level module. The controlled construction belongs with the existing
trajectory-sheaf API.

### Exact API

```julia
"""
    continuous_to_discrete_zoh(Ac::AbstractMatrix{T},
                               Bc::AbstractMatrix{T},
                               h::Real) where T
        -> (Ad::Matrix{T}, Bd::Matrix{T})

Discretize the continuous-time LTI system `x'(t) = Ac*x(t) + Bc*u(t)` using
zero-order hold on a uniform step size `h`.

The returned matrices satisfy `x[k+1] = Ad*x[k] + Bd*u[k]`.
"""
function continuous_to_discrete_zoh(Ac, Bc, h)

"""
    ControlledTrajectorySheaf{T}

Trajectory sheaf for the zero-order-hold discretization of the continuous-time
controlled system `x'(t) = Ac*x(t) + Bc*u(t)`.
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

"""
    build_controlled_trajectory_sheaf(F::EuclideanSheaf{T},
                                      Ac::AbstractMatrix{T},
                                      Bc::AbstractMatrix{T},
                                      h::Real,
                                      k::Int) where T
        -> EuclideanSheaf{T}

Construct the path-graph sheaf whose global sections are feasible sampled
state-control trajectories of the zero-order-hold discretization.
"""
function build_controlled_trajectory_sheaf(F, Ac, Bc, h, k)

"""
    ControlledTrajectorySheaf(F::EuclideanSheaf{T},
                              Ac::AbstractMatrix{T},
                              Bc::AbstractMatrix{T},
                              h::Real,
                              k::Int) where T
        -> ControlledTrajectorySheaf{T}
"""
function ControlledTrajectorySheaf(F, Ac, Bc, h, k)

"""
    feasible_control_trajectory_basis(ts::ControlledTrajectorySheaf{T},
                                      x1::AbstractVector{T},
                                      xk1::AbstractVector{T}) where T
        -> (z_p::Vector{T}, null_basis::AbstractMatrix{T})

Return one feasible trajectory hitting the endpoint states together with a basis
matrix whose columns span all feasible perturbations preserving those endpoint
states.

The returned coordinates use the public state-first ordering described in the
Mathematical Background section.
"""
function feasible_control_trajectory_basis(ts, x1, xk1)
```

### Implementation sketch

1. **Zero-order-hold discretization.**
   Implement `continuous_to_discrete_zoh` using the block-exponential identity
   above. Use only `LinearAlgebra.exp`; do not add control-toolbox dependencies.

2. **Controlled trajectory sheaf.**
   Let the state dimension be the total 0-cochain dimension of `F`, and let the
   control dimension be the column count of `Bc`.
   Build a path graph with `k+2` vertices:
   - vertex 1 has stalk dimension `n`,
   - vertices `2:(k+1)` have stalk dimension `n + m`,
   - vertex `k+2` has stalk dimension `n`.

   Use the restriction maps from the Mathematical Background section exactly.

3. **Endpoint-state fixing by whole-vertex boundary data.**
   Because both endpoint stalks are state-only, the public
   `harmonic_extension(s, boundary)` interface can be used directly:
   fix vertex 1 to the initial state and vertex `k+2` to the terminal state.

   No fixed-coordinate generalization is required in this issue.

4. **Coordinate extraction.**
   The internal harmonic-extension result lives in the sheaf’s 0-cochain
   coordinates. Add a deterministic extraction map from the internal sheaf
   coordinates to the public trajectory coordinates

   ```math
   z = (x_1, x_2, \dots, x_{k+1}, u_1, u_2, \dots, u_k).
   ```

   Apply that same extraction map to both the particular solution and each
   null-basis column before returning them from
   `feasible_control_trajectory_basis`.

5. **Validation.**
   Use `@argcheck` for:
   - `k >= 1`,
   - `h > 0`,
   - `Ac` square,
   - `size(Bc, 1) == n`,
   - `length(x1) == n`,
   - `length(xk1) == n`.

## Tests to Write

Create `test/network_sheaves/ControlledTrajectorySheaf.jl` and include it from
`test/runtests.jl`.

Use the scalar integrator as the primary exact test case:

```math
\dot{x}(t) = u(t), \qquad A_c = [0], \quad B_c = [1].
```

With `h = 1` this gives `A_d = [1]`, `B_d = [1]`.

Recommended tests:

1. `continuous_to_discrete_zoh([0.0], [1.0], 1.0)` returns `Ad ≈ [1.0]`,
   `Bd ≈ [1.0]`.
2. `ControlledTrajectorySheaf` stores the expected `k`, `h`, `Ad`, `Bd`,
   `state_dim`, and `control_dim`.
3. `vertex_stalks(ts.sheaf)` equals `[1; fill(2, k); 1]` in the scalar
   integrator case.
4. Every edge stalk has dimension `1`.
5. The dummy initial edge restriction maps are `[1]` and `[1 0]`.
6. The dynamics edge restriction maps are `[1 1]` and `[1 0]` on non-terminal
   dynamics edges, and `[1 1]` and `[1]` on the final edge.
7. `feasible_control_trajectory_basis(ts, [0.0], [1.0])` returns a particular
   trajectory in the public ordering
   ```math
   z = (x_1, x_2, \dots, x_{k+1}, u_1, \dots, u_k)
   ```
   with initial state `0` and terminal state `1`.
8. For `k = 2`, the nullspace dimension is `1` in the scalar integrator case.
9. Every returned basis column corresponds to a global section of the internal
   sheaf after applying the inverse coordinate extraction.
10. Every null-basis column preserves the endpoint states.
11. Wrong-sized `Ac`, `Bc`, `x1`, or `xk1` throw `ArgumentError`.
12. `h <= 0` and `k < 1` throw `ArgumentError`.

## Verification Checklist

- [ ] `continuous_to_discrete_zoh` uses the block-exponential zero-order-hold formula.
- [ ] `ControlledTrajectorySheaf` extends the existing trajectory-sheaf module rather than creating an unrelated API surface.
- [ ] The controlled sheaf’s global sections are exactly the feasible sampled trajectories.
- [ ] The dummy initial vertex adds only one extra copy of the state while keeping the endpoint boundary conditions whole-vertex.
- [ ] `feasible_control_trajectory_basis` returns the public trajectory coordinates with all states before all controls.
- [ ] `feasible_control_trajectory_basis` returns both a particular feasible trajectory and a basis for all feasible perturbations.
- [ ] `test/network_sheaves/ControlledTrajectorySheaf.jl` passes within the full suite.

## Out of Scope

- Nonlinear dynamics or control-affine systems with state-dependent input matrix.
- Time-varying `A_c(t)` or `B_c(t)`.
- Inequality constraints on states or controls.
- Feedback laws `u = Kx`; this issue is about open-loop feasible trajectories.
- Cost optimization over the feasible affine space.
- Stochastic control or model-predictive control.
- Non-uniform time grids.

## References

- Justin Curry, *Sheaves, Cosheaves and Applications* — for the cochain/global-section language used to encode trajectory constraints.
- Robert Ghrist and collaborators on sheaf-theoretic signal processing and sheaf Laplacians — for the sheaf-Laplacian viewpoint behind the harmonic-extension formulation.
- The current `harmonic_extension` and `TrajectorySheaf` implementations in CellularSheaves.jl — for the numerical and API conventions this issue extends.
- Standard zero-order-hold discretization for linear systems, especially the block-matrix exponential identity used to compute `(A_d, B_d)` without extra dependencies.

## Prerequisite / Follow-up

This issue is the prerequisite for **Issue 003**, which optimizes a convex
quadratic cost over the affine space returned by
`feasible_control_trajectory_basis`.
