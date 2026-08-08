# 004 — Controlled trajectory example progression from classical control benchmarks

## Mathematical Background

This issue requests a documentation-focused follow-on to **Issue 002** and
**Issue 003**. Once the controlled trajectory sheaf and quadratic optimal-control
APIs exist, the package should include a small progression of benchmark examples
that move from the simplest single-agent linear system to graph-coupled mechanical
systems that are close to later consensus and multi-agent work.

All examples in this issue should use the same continuous-time control model

```math
\dot{x}(t) = A_c x(t) + B_c u(t),
```

with zero-order-hold discretization handled internally by the API from
Issue 002, and quadratic optimal control posed over the feasible affine
trajectory space from Issue 003.

The requested progression is:

1. **Double integrator**
2. **Vehicle platoon**
3. **Planar quadrotor**
4. **Mass-spring-damper chain**

This ordering is intentional.

- The **double integrator** is the canonical first LQR example and should be the
  simplest possible explanation of the trajectory-space workflow.
- The **vehicle platoon** keeps essentially the same control story while moving
  to a multi-agent / path-graph setting.
- The **planar quadrotor** adds a more recognizably robotic multi-input system
  with a standard hover linearization.
- The **mass-spring-damper chain** returns to a graph-coupled mechanical system
  that is especially close to later sheaf- and consensus-style network examples.

Each example should make the same conceptual pipeline visible:

```math
(A_c, B_c)
\longrightarrow
(A_d, B_d)
\longrightarrow
\mathcal{T}(x_1, x_{k+1})
\longrightarrow
\min_{z \in \mathcal{T}(x_1, x_{k+1})} J(z).
```

The examples should therefore not be just collections of matrices. They should
show how the controlled-trajectory construction, harmonic-extension viewpoint,
and convex quadratic optimization fit together on standard benchmark systems from
the controls literature.

## Codebase Orientation

| File | Why it matters |
| --- | --- |
| `docs/issues/002-controlled-trajectory-affine-space.md` | Prerequisite API and mathematical setup for feasible endpoint-hitting controlled trajectories. |
| `docs/issues/003-controlled-trajectory-quadratic-optimal-control.md` | Prerequisite API for quadratic optimal control over the feasible affine space. |
| `docs/issues/control-examples-briefing.md` | Source document for the chosen progression, literature references, and example-selection rationale. |
| `docs/literate/trajectory_sheaf.jl` | Existing literate example showing the repo’s preferred documentation style, plotting approach, and conceptual level. |
| `docs/make.jl` | Add the new literate example pages here so they are rendered by Documenter and appear in the docs navigation. |
| `docs/src/index.md` | Add a short pointer to the new controlled-trajectory example progression if the landing page is updated. |
| `test/network_sheaves/TrajectorySheaf.jl` | Style reference for trajectory-related numerical tests. |

## Requested Implementation

Add a sequence of literate documentation examples that use the controlled
trajectory API from Issues 002 and 003. The examples should live in
`docs/literate/` and be listed in `docs/make.jl` in the same order as the
progression above.

### Exact deliverables

Add these new literate files:

```text
docs/literate/controlled_double_integrator.jl
docs/literate/controlled_vehicle_platoon.jl
docs/literate/controlled_planar_quadrotor.jl
docs/literate/controlled_mass_spring_damper_chain.jl
```

Register the generated pages in `docs/make.jl` under `"Examples"` in exactly
this order:

1. `generated/controlled_double_integrator.md`
2. `generated/controlled_vehicle_platoon.md`
3. `generated/controlled_planar_quadrotor.md`
4. `generated/controlled_mass_spring_damper_chain.md`

If `docs/src/index.md` is touched, only add a short sentence pointing readers to
the controlled-trajectory examples. Do not turn the landing page into a long
tutorial.

### Common structure for all four examples

Each literate example should:

1. state the physical system and cite one or two standard references already
   listed in `docs/issues/control-examples-briefing.md`;
2. define a concrete continuous-time `(A_c, B_c)` model or standard local
   linearization;
3. choose a small, reproducible horizon length and step size;
4. fix endpoint states and leave controls free, consistent with Issue 002;
5. build the controlled trajectory object using the public API from Issue 002;
6. assemble a quadratic objective using the public API from Issue 003;
7. solve for the optimal feasible trajectory;
8. plot or otherwise summarize the sampled state and control trajectories;
9. explain, in prose, what this example adds relative to the previous one in the
   progression.

All mathematics displayed in the prose should use fenced `math` blocks. Do not
put LaTeX source in single-backtick code spans. Follow the repo’s Literate.jl
conventions: comments above code blocks, not line comments inside function
definitions.

### Example-specific requirements

#### 1. Double integrator

Use the standard point-mass model

```math
x =
\begin{bmatrix}
p \\
\dot{p}
\end{bmatrix},
\qquad
u = a,
\qquad
A_c =
\begin{bmatrix}
0 & 1 \\
0 & 0
\end{bmatrix},
\qquad
B_c =
\begin{bmatrix}
0 \\
1
\end{bmatrix}.
```

This should be the cleanest “first look” example and the one most suitable for
users to copy-paste into their own experiments.

#### 2. Vehicle platoon

Use a small longitudinal platoon on a path graph, for example two or three
vehicles with state

```math
x_i =
\begin{bmatrix}
p_i \\
\dot{p}_i
\end{bmatrix},
\qquad
u_i = a_i.
```

The example should explicitly point out that this is the first step toward the
later consensus / formation direction. Keep the model simple: a linearized
spacing-regulation or stacked double-integrator formulation is preferable to a
heavier traffic model.

#### 3. Planar quadrotor

Use a hover-linearized planar quadrotor with state and control of the form

```math
x =
\begin{bmatrix}
y \\
z \\
\phi \\
\dot{y} \\
\dot{z} \\
\dot{\phi}
\end{bmatrix},
\qquad
u =
\begin{bmatrix}
u_1 \\
u_2
\end{bmatrix}.
```

Keep the example focused on the controlled-trajectory workflow, not on deriving
the full rigid-body model from scratch. It is enough to present the standard
hover linearization and cite the source.

#### 4. Mass-spring-damper chain

Use a short chain of coupled masses with forcing input, phrased as a graph-coupled
mechanical system. This example should make the transition back toward
network-structured mechanics explicit, since that is the bridge to future sheaf
and consensus examples.

Prefer a chain with two or three masses so the plots remain readable and the docs
build stays lightweight.

### Scope and sequencing notes

- Treat this as a documentation/examples issue, not as the place to invent new
  control APIs.
- The examples should call the public interfaces introduced by Issues 002 and 003
  rather than reaching into internal helpers.
- If Issue 003 lands with a slightly different function name than currently
  drafted, update the examples to match that final public API rather than
  preserving stale names from the issue text.

## Tests to Write

Create `test/network_sheaves/ControlledTrajectoryExamples.jl` and include it from
`test/runtests.jl`.

This test file should be lightweight and should validate the core numerical setup
 behind each documented example without depending on plotting or Documenter.

Recommended coverage:

1. **Double integrator**
   - construct the example system,
   - compute an optimal trajectory for fixed endpoint states,
   - `@test` that the returned trajectory has the expected total coordinate length,
   - `@test` that the initial and terminal state blocks match the requested endpoints.

2. **Vehicle platoon**
   - construct a small stacked system,
   - `@test` that the optimization call succeeds,
   - `@test` that the endpoint constraints are satisfied for every vehicle.

3. **Planar quadrotor**
   - construct the hover-linearized model,
   - `@test` that the returned control trajectory has two control coordinates per time step,
   - `@test` that all returned entries are finite.

4. **Mass-spring-damper chain**
   - construct a small coupled chain,
   - `@test` that the optimal trajectory satisfies the requested endpoints,
   - `@test` that the example uses a state dimension larger than the single-agent
     double-integrator case, so it actually exercises a coupled model.

The tests do not need to compare against closed-form optimal controls. Their job
is to ensure the four documented workflows remain executable and numerically sane
as the public API evolves.

## Verification Checklist

- [ ] Four new literate example files added under `docs/literate/`
- [ ] Examples appear in `docs/make.jl` in the requested progression order
- [ ] Each example cites literature already identified in the briefing document
- [ ] Each example uses the public controlled-trajectory and quadratic-objective APIs
- [ ] Each example explains its role in the progression
- [ ] A lightweight regression test file covers all four examples
- [ ] `julia --project=test test/runtests.jl` passes
- [ ] `julia --project=docs docs/make.jl` passes

## Out of Scope

- Adding nonlinear control examples or nonlinear trajectory optimization
- Adding inequality constraints, MPC, or external QP / optimal-control dependencies
- Adding consensus algorithms themselves in this issue
- Replacing the four selected examples with a larger unsorted example catalog
- Writing a literature survey beyond the short citations already needed for the examples
- Adding cart-pole, bicycle, DC motor, or satellite examples in the same issue
