# Sheaf-Coordinated Reinforcement Learning

Learned control on top of the cellular-sheaf coordinator of the **ECC paper** — *Heterogeneous
Multi-Agent Multi-Target Tracking using Cellular Sheaves* (Hanks, Nino, Bou Barcelo, Copeland,
Dixon, Fairbanks; arXiv:2512.24886). That paper coordinates a fleet by the **harmonic extension**
of a sheaf Laplacian and drives each agent with the analytic Lyapunov feedback `u_i = −k·(Lz)_i`,
which converges only to a *ball* (drift is treated as bounded, never modeled).

This example keeps the paper's sheaf coordination and **replaces the fixed feedback with a learned
residual policy** (TD3) trained per dynamics type — collapsing the drift ball and extending the
paper's single-integrator agents to underactuated and constrained dynamics.

The control law each agent applies is

```
u_i = u_base(q*)  +  α · π_θ(s_i)
```

where `q* = H⁻¹ B p` is the harmonic extension of the sheaf Laplacian system (the paper's analytic
"Layer 1"), `u_base` is the corresponding feedback, and `π_θ` is the learned residual
("Layer 2"). Layer 1 is shared across the fleet; Layer 2 is a separate policy per dynamics type
(single-integrator, underactuated quadrotor), because the state/action dimensions differ.

## Layout

```
lib/     shared library — the environment, TD3, evaluation, and rendering
  sheaf_rl.jl    module SheafRL: Config, SheafEnv, base control, drift, TD3 agent, eval
  render.jl      rollout → animation/plots (analytic | drift-oracle | sheaf+RL panels)
train/   one trainer per demonstration; writes a policy to cache/*.jld2
  single_integrator.jl   time-varying drift, single-integrator agents (the paper's dynamics)
  quadrotor.jl           underactuated quadrotor residual (extends the paper to 2nd-order dynamics)
  multiagent.jl          13-agent, 4-target formation (identity restriction maps) — the paper's setup
  heterogeneous.jl       genuinely heterogeneous sheaf (mixed restriction maps)
  actuator_cap.jl        hard actuator constraint (ball projection)
  behaviour_cloning.jl   supervised clone of the paper's analytic law
  bc_data.jl             collect behaviour-cloning demonstrations
viz/     one renderer per demonstration (reads cache/*.jld2, writes figures to cache/viz_*/)
eval/    quantitative diagnostics (generalization, control effort, drift headroom)
cache/   trained policies (*.jld2, tracked); rendered figures (cache/viz_*/, git-ignored)
```

## Trained policies (in `cache/`)

| policy                  | demonstration                          |
|-------------------------|----------------------------------------|
| `rl_multiagent_f2.jld2` | multi-agent formation under drift      |
| `rl_hetero_f2.jld2`     | heterogeneous sheaf                    |
| `rl_ballcap_u15.jld2`   | actuator-constrained (‖u‖ ≤ 1.5)       |
| `rl_quad_class.jld2`    | underactuated quadrotor                |
| `rl_sheaf_tv.jld2`      | single-integrator, time-varying drift  |
| `rl_sheaf.jld2`         | single-integrator, no drift (baseline) |
| `bc_sheaf.jld2`         | behaviour-cloned analytic law          |

## Setup

The example has its own `Project.toml`. From the repository root, develop the local
`CellularSheaves` package into it and instantiate (one time):

```
julia --project=examples/RL -e 'using Pkg; Pkg.develop(PackageSpec(path=pwd())); Pkg.instantiate()'
```

Requires Julia ≥ 1.10. `CUDA`/`cuDNN` are dependencies but **optional at runtime** — every trainer
falls back to CPU when `CUDA.functional()` is false, so the demos run on a laptop without a GPU.

## Usage

Train (writes a policy to `cache/`):

```
julia --project=examples/RL examples/RL/train/multiagent.jl
```

Render one demonstration (reads the policy, writes figures to `cache/viz_*`):

```
julia --project=examples/RL examples/RL/viz/multiagent.jl
```

The viz renderers are Flux-free (they read raw weight matrices from the `.jld2`), so they run
headless and are decoupled from the training Julia version.

Regenerate every documentation figure and commit them locally (never pushes):

```
examples/RL/refresh_viz.sh            # regenerate + copy + commit
examples/RL/refresh_viz.sh --dry-run  # print the plan only
```
