# 003 — Convex quadratic optimal control on the feasible trajectory space

## Mathematical Background

This issue depends on **Issue 002**, which provides:

- a `ControlledTrajectorySheaf`,
- a feasible endpoint-hitting trajectory `z_p`,
- a basis matrix `N` for all feasible perturbations preserving the endpoint
  states.

The admissible trajectories therefore form the affine space

```math
\mathcal{T}(x_1, x_{k+1}) = \{ z_p + N\alpha : \alpha \in \mathbb{R}^r \}.
```

The optimization target for this issue is the class of **convex quadratic
objectives**

```math
J(z) = \tfrac12 z^\top H z + f^\top z + c,
```

with positive-semidefinite Hessian. Restricting to the feasible affine space
gives the reduced problem

```math
\min_{\alpha \in \mathbb{R}^r}
\tfrac12 \alpha^\top (N^\top H N)\alpha
+ \alpha^\top N^\top (H z_p + f).
```

Define

```math
R = N^\top H N, \qquad
r = N^\top (H z_p + f).
```

Then the optimizer satisfies

```math
R \alpha^\star = -r.
```

Since the full Hessian is positive semidefinite, the reduced Hessian is also
positive semidefinite. When the reduced Hessian is singular, the implementation
should return the **minimum-norm optimizer**
using the same LDLt pseudoinverse logic already used in `harmonic_extension`.

The first concrete objective family should be a discrete **LQR/LQT-style**
quadratic cost over the sampled state-control trajectory

```math
z = (x_1, x_2, \dots, x_{k+1}, u_1, u_2, \dots, u_k),
```

namely

```math
J(z)
=
\tfrac12 \sum_{t=1}^{k}
\Bigl((x_t - \bar{x}_t)^\top Q (x_t - \bar{x}_t)
     + (u_t - \bar{u}_t)^\top R_u (u_t - \bar{u}_t)\Bigr)
+
\tfrac12 (x_{k+1} - \bar{x}_{k+1})^\top Q_f (x_{k+1} - \bar{x}_{k+1}),
```

with `Q \succeq 0`, `Q_f \succeq 0`, and `R_u \succ 0`.

This is a good fit for the current dependency envelope:

- no inequality constraints,
- no external QP solver,
- no new optimization packages.

Purely linear objectives should **not** be the primary target here: on an
unconstrained affine space they are often unbounded below. The implementation
may accept `H = 0` only when it can prove the reduced problem is bounded; the
default documented path should remain convex quadratic objectives.

## Codebase Orientation

| File | Why it matters |
| --- | --- |
| `docs/issues/002-controlled-trajectory-affine-space.md` | This issue’s required prerequisite; defines the feasible affine space and controlled trajectory object. |
| `src/network_sheaves/TrajectorySheaf.jl` | This is where the controlled trajectory functions from Issue 002 should live, so the optimization layer should extend the same module. |
| `src/network_sheaves/EuclideanSheaves.jl` | Reuse the LDLt pseudoinverse / nullspace pattern for singular reduced Hessians. |
| `docs/literate/trajectory_sheaf.jl` | Existing documentation example showing colocation and the `M^\top M = L_{II}` relationship; useful template for the new LQR-style example. |
| `test/network_sheaves/TrajectorySheaf.jl` | Current style for trajectory-related tests and examples. |
| `src/network_sheaves/NetworkSheaves.jl` | Any new public exports must be re-exported here. |

## Requested Implementation

Extend the trajectory-sheaf module from Issue 002 with:

1. a generic convex-quadratic optimizer over the feasible trajectory affine
   space;
2. a helper that assembles the standard stagewise LQR/LQT quadratic cost into
   `(H, f, c)`;
3. a literate documentation example using an LQR-style objective.

### Exact API

```julia
"""
    lqr_objective(ts::ControlledTrajectorySheaf{T},
                  Q::AbstractMatrix{T},
                  Ru::AbstractMatrix{T};
                  Qf::Union{Nothing,AbstractMatrix{T}}=nothing,
                  x_ref::Union{Nothing,AbstractMatrix{T}}=nothing,
                  u_ref::Union{Nothing,AbstractMatrix{T}}=nothing) where T
        -> (H::SparseMatrixCSC{T,Int}, f::Vector{T}, c::T)

Assemble the quadratic objective

    0.5*z'*H*z + f'*z + c

for the controlled sampled trajectory `z`.

The reference arrays are optional; when omitted, use zero references. The state
reference has one column for each sampled state, and the control reference has
one column for each sampled control.
"""
function lqr_objective(ts, Q, Ru; Qf=nothing, x_ref=nothing, u_ref=nothing)

"""
    optimal_control_trajectory(ts::ControlledTrajectorySheaf{T},
                               x1::AbstractVector{T},
                               xk1::AbstractVector{T},
                               H::AbstractMatrix{T},
                               f::AbstractVector{T}=zeros(T, size(H, 1))) where T
        -> (z_opt::BlockVector, α_opt::Vector{T}, z_p::BlockVector, null_basis::AbstractMatrix{T})

Compute the minimum-cost feasible trajectory for the convex quadratic objective
subject to fixed initial and terminal states.
"""
function optimal_control_trajectory(ts, x1, xk1, H, f=zeros(...))
```

### Implementation sketch

1. **Build on Issue 002.**
   Call `feasible_control_trajectory_basis(ts, x1, xk1)` to obtain
   `(z_p, N)`.

2. **Reduced quadratic program.**
   Form:

   ```julia
   q_p = Array(z_p)
   Rred = Symmetric(N' * H * N)
   rred = N' * (H * q_p + f)
   ```

   Then solve

   ```math
   R_{\mathrm{red}} \alpha^\star = -r_{\mathrm{red}}
   ```

   using `ChordalLDLt` when `Rred` is sparse / symmetric PSD.

3. **Singular reduced Hessian.**
   Reuse the pseudoinverse/nullspace logic already present in
   `EuclideanSheaves.jl` to compute the minimum-norm optimizer when `Rred` is
   singular.

   If `Rred` is singular and `rred` has a component outside the range of
   `Rred`, the objective is unbounded below on the feasible affine space;
   throw `ArgumentError` rather than silently returning a bogus result.

4. **Empty nullspace case.**
   If `size(N, 2) == 0`, then the feasible trajectory is unique and the
   optimizer is just `z_p`.

5. **LQR objective assembly.**
   `lqr_objective` should construct `H`, `f`, and `c` in the public trajectory
   coordinates

   ```math
   z = (x_1, x_2, \dots, x_{k+1}, u_1, u_2, \dots, u_k),
   ```

   with all state blocks placed before all control blocks.

   The assembled `H` should be sparse and block diagonal in time.

6. **Validation.**
   Use `@argcheck` for:
   - `Q`, `Qf` square of size `n × n`,
   - `Ru` square of size `m × m`,
   - `Ru` symmetric positive definite for the documented LQR path,
   - `x_ref` shape `n × (k+1)` when provided,
   - `u_ref` shape `m × k` when provided,
   - `size(H, 1) == size(H, 2) == n*(k+1) + m*k`,
   - `length(f) == size(H, 1)`.

7. **Documentation example.**
   Add a literate example, e.g. `docs/literate/controlled_trajectory_lqr.jl`,
   using a simple continuous-time system such as the double integrator

   ```math
   \dot{x} = \begin{bmatrix}0 & 1\\0 & 0\end{bmatrix}x
           + \begin{bmatrix}0\\1\end{bmatrix}u,
   ```

   discretized internally by `continuous_to_discrete_zoh`, with endpoint
   conditions and an LQR-style running cost. The example should plot the state
   and control and explicitly describe the relationship between:
   - the feasible affine space from Issue 002,
   - the reduced quadratic program in `\alpha`,
   - the returned optimal trajectory.

## Tests to Write

Create `test/network_sheaves/ControlledOptimalControl.jl` and include it from
`test/runtests.jl`.

Use the scalar integrator as the exact regression case:

```math
\dot{x} = u, \qquad h = 1, \qquad x_1 = 0, \qquad x_{k+1} = 1.
```

With `Q = [0]`, `Ru = [1]`, `Qf = [0]`, the minimum-energy solution should have
constant controls and linearly interpolating states.

Recommended tests:

1. `lqr_objective` returns `H`, `f`, `c` of the correct sizes and sparse type.
2. With zero references, `f` is zero for the standard LQR cost.
3. `optimal_control_trajectory` returns a trajectory satisfying the endpoint
   states exactly.
4. The returned trajectory satisfies the discrete controlled dynamics implied by
   the sheaf.
5. In the scalar-integrator minimum-energy case, all controls equal `1/k`
   (within tolerance) and the states are linear in time.
6. The optimizer lies in the affine space from Issue 002:
   `Array(z_opt) ≈ Array(z_p) + null_basis * α_opt`.
7. If `size(null_basis, 2) == 0`, the function returns `z_p`.
8. Wrong-sized `Q`, `Ru`, `Qf`, `x_ref`, `u_ref`, `H`, or `f` throw
   `ArgumentError`.
9. A deliberately unbounded reduced problem throws `ArgumentError`.
10. The reduced first-order optimality condition holds:
    `norm((null_basis' * H * null_basis) * α_opt + null_basis' * (H * Array(z_p) + f)) < 1e-10`.

## Verification Checklist

- [ ] The implementation depends on Issue 002 and does not duplicate its feasible-space logic.
- [ ] `optimal_control_trajectory` optimizes over the affine space `z_p + Nα`, not over all cochains directly.
- [ ] The reduced Hessian solve uses existing sparse linear-algebra tools (`CliqueTrees.Multifrontal`) only.
- [ ] Singular reduced Hessians return minimum-norm minimizers when bounded and throw when unbounded.
- [ ] `lqr_objective` assembles the standard stagewise quadratic cost without adding dependencies.
- [ ] A literate LQR-style example is added to the docs and wired into `docs/make.jl`.
- [ ] The full test suite passes.

## Out of Scope

- Inequality-constrained optimal control.
- Box constraints, saturation, or state path constraints.
- Riccati-recursion implementations or dynamic-programming solvers.
- Non-convex costs.
- Nonlinear control systems.
- Time-varying weights `Q_t`, `R_t`, `Q_{f,t}` unless they drop out for free from the chosen implementation.
- MPC or receding-horizon controllers.
- Kalman filtering / estimation.

## References

- Justin Curry, *Sheaves, Cosheaves and Applications* — for the sheaf/cochain perspective.
- Robert Ghrist and collaborators on sheaf Laplacians and sheaf signal processing — for the harmonic-extension/global-section viewpoint behind the feasible-space construction.
- The current CellularSheaves.jl `harmonic_extension` and `TrajectorySheaf` code paths — for the LDLt-based numerical conventions this issue extends.
- Standard continuous-time LQR / linear-quadratic tracking formulations, used here only as a quadratic objective template over the feasible trajectory affine space.

## Prerequisite

This issue depends on **Issue 002** being completed first.
