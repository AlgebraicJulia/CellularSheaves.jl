# 005 — Visualize feasible control nullspace basis as a trajectory family

## Mathematical Background

In the controlled trajectory workflow, fixed endpoints produce an affine feasible
set

```math
\mathcal{T}(x_1, x_{k+1}) = \{ z_p + N\alpha : \alpha \in \mathbb{R}^r \},
```

where:

- `z_p` is a particular feasible trajectory,
- `N = [n_1 \; n_2 \; \cdots \; n_r]` is a basis for endpoint-preserving feasible
  perturbations.

Each basis direction `n_j` gives a canonical one-parameter family

```math
z^{(j)}(\lambda) = z_p + \lambda n_j.
```

For the control examples in `docs/literate/control/`, the trajectory vector uses
the public coordinate ordering

```math
z = (x_1, x_2, \dots, x_{k+1}, u_1, u_2, \dots, u_k),
```

so every `z^{(j)}(\lambda)` simultaneously specifies:

1. a full sampled state trajectory `x_1,\dots,x_{k+1}`,
2. a full sampled control sequence `u_1,\dots,u_k`.

The requested visualization is a **basis-indexed family plot**: one rendered
trajectory per basis vector, optionally with `\lambda = \pm \sigma` to show
direction and scale. This makes the nullspace geometrically interpretable:
“which endpoint-preserving control mode does this basis vector represent?”

This viewpoint aligns with sheaf-harmonic extension interpretations in the
CellularSheaves literature and provides a direct controls-facing interpretation
for users coming from LTI/LQR examples.

## Codebase Orientation

| File | Why it matters |
| --- | --- |
| `src/network_sheaves/TrajectorySheaf.jl` | Contains `feasible_control_trajectory_basis`, public trajectory coordinate layout, and controlled optimization APIs that expose the nullspace basis. |
| `src/network_sheaves/SheafInterface.jl` | Add abstract declarations for any new public helper introduced for basis-family extraction/visualization support. |
| `src/network_sheaves/NetworkSheaves.jl` | Export surface for any new public utility function added in this issue. |
| `docs/literate/control/controlled_double_integrator.jl` | Smallest control example; best place to introduce the basis-family visualization concept first. |
| `docs/literate/control/controlled_vehicle_platoon.jl` | Multi-agent setting where basis modes are especially interpretable as collective patterns. |
| `docs/literate/control/controlled_planar_quadrotor.jl` | Multi-input single-agent example to confirm visualization works beyond scalar control. |
| `docs/literate/control/controlled_mass_spring_damper_chain.jl` | Coupled mechanical benchmark where nullspace modes can be interpreted as coordinated forcing patterns. |
| `test/network_sheaves/TrajectorySheaf.jl` | Existing trajectory test conventions; add focused regression tests for basis-family extraction and reconstruction. |

## Requested Implementation

Add a small public utility that converts nullspace basis columns into
trajectory-family members in block form, then use it in the control literate
examples to plot one family member per basis vector.

### Exact API

```julia
"""
    nullspace_trajectory_family(
        ts::ControlledTrajectorySheaf{T},
        z_p::AbstractVector{T},
        null_basis::AbstractMatrix{T};
        amplitude::T = one(T),
        include_negative::Bool = false
    ) where T
        -> Vector{BlockVector{T}}

Construct a family of feasible trajectories from the affine parameterization
`z = z_p + N*α` by activating one nullspace basis direction at a time.

For each basis column `n_j` in `null_basis`, return `z_p + amplitude*n_j`.
If `include_negative=true`, also return `z_p - amplitude*n_j` for each `j`.

Each returned element is a `BlockVector` in the public trajectory block layout:
`k+1` state blocks (size `ts.state_dim`) followed by `k` control blocks
(size `ts.control_dim`).
"""
function nullspace_trajectory_family(ts, z_p, null_basis; amplitude=one(T), include_negative=false)
```

### Implementation sketch

1. **Input validation (`@argcheck`).**
   - Let `p = (ts.k + 1) * ts.state_dim + ts.k * ts.control_dim`.
   - Check `length(z_p) == p`.
   - Check `size(null_basis, 1) == p`.
   - Check `amplitude >= 0`.

2. **Block layout helper.**
   Reuse the same block partition used by `optimal_control_trajectory`:
   `vcat(fill(ts.state_dim, ts.k + 1), fill(ts.control_dim, ts.k))`.

3. **Construct basis-indexed family.**
   - For each basis column `j`, compute `z_plus = z_p + amplitude * null_basis[:, j]`.
   - Convert `z_plus` to `BlockArray(z_plus, block_sizes)` and append.
   - If `include_negative`, also construct and append the `z_minus` version.

4. **Degenerate case (`r = 0`).**
   - If `size(null_basis, 2) == 0`, return an empty vector.
   - Do not invent pseudo-directions.

5. **Documentation integration.**
   In each control example file, after obtaining
   `z_opt, α_opt, z_p, null_basis = optimal_control_trajectory(...)`:
   - call `nullspace_trajectory_family(ts, Array(z_p), null_basis; amplitude=1.0)`,
   - extract control/state traces from each returned block trajectory,
   - overlay basis-family curves with muted styling,
   - keep `z_opt` highlighted as the selected optimum.

6. **Narrative standardization across examples.**
   Add a short repeated prose block explaining:
   - each basis vector yields one endpoint-preserving control mode,
   - trajectories in these overlays are feasible but not necessarily optimal,
   - the optimal solution is a linear combination in this basis.

### Plotting contract in docs

For each `docs/literate/control/*.jl` example, add:

1. a state-family panel (one curve per basis direction),
2. a control-family panel (one curve per basis direction),
3. a legend entry keyed by basis index (`n₁`, `n₂`, …).

If the basis dimension is large (`r > 8`), limit plotting to the first 8 basis
vectors and print a note indicating truncation.

## Tests to Write

Create a dedicated test file:

`test/network_sheaves/NullspaceTrajectoryFamily.jl`

and include it from `test/runtests.jl`.

Recommended tests:

1. **Shape and block layout test**
   - build a small controlled system with nontrivial nullspace,
   - compute `(z_p, N) = feasible_control_trajectory_basis(...)`,
   - call `nullspace_trajectory_family(...)`,
   - `@test` each output is a `BlockVector` with `k+1+k` blocks and expected block sizes.

2. **Basis-direction reconstruction test**
   - for each returned `z_j`, verify
     `Array(z_j) ≈ z_p + amplitude * N[:, j]` (and negative branch when enabled).

3. **Feasibility preservation test**
   - for each generated family trajectory, verify endpoint blocks match `x1`, `xk1`,
   - verify discrete dynamics residual is below tolerance at every step.

4. **Empty-basis behavior**
   - construct a case with unique feasible trajectory (`size(N,2)==0`),
   - `@test isempty(nullspace_trajectory_family(...))`.

5. **Validation failures**
   - wrong `z_p` length throws `ArgumentError`,
   - wrong `N` row count throws `ArgumentError`,
   - negative amplitude throws `ArgumentError`.

## Verification Checklist

- [ ] `nullspace_trajectory_family` is implemented in `TrajectorySheaf.jl` with `@argcheck` validation
- [ ] New public declaration/export is wired through `SheafInterface.jl` and `NetworkSheaves.jl`
- [ ] All four control literate examples include basis-family overlays
- [ ] Example prose explains basis-family meaning and relation to optimal solution
- [ ] Plot truncation behavior for large `r` is documented and implemented
- [ ] New nullspace-family tests are added and included from `test/runtests.jl`
- [ ] `julia --project=test test/runtests.jl` passes
- [ ] `julia --project=docs docs/make.jl` passes with no new doc warnings

## Out of Scope

- Nonlinear continuation/manifold visualization of feasible trajectories
- Interactive GUI widgets or animation frameworks
- Time-varying basis re-orthogonalization or PCA compression of nullspace modes
- Introducing new plotting dependencies beyond current docs stack
- Redefining the core affine-space or optimal-control APIs from Issues 002/003
- Adding inequality-constrained feasible-set visualization in this issue

## References

- Justin Curry, *Sheaves, Cosheaves and Applications* (sheaf/harmonic-extension viewpoint).
- Robert Ghrist et al., sheaf Laplacian and harmonic extension literature (interpretation of constrained extensions and null directions).
- Standard LTI trajectory parameterization references in optimal control texts (affine feasible sets under endpoint constraints).
