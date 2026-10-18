# Schwarz Methods and Coordination Sheaves

This page relates the overlapping Schwarz solvers
([`SchwarzMethods`](../api/schwarz_methods.md)) to the package's theory of
multi-agent coordination over cellular sheaves
([Coordination theory](distributed_sheaf_solve/theory_coordination.md)). It
covers four things: where the two match exactly, where they differ, what each
side offers the other, and where the current implementation does not yet fit
the sheaf picture.

## A dictionary

| Coordination sheaf | Overlap sheaf of a Schwarz decomposition |
|---|---|
| agent ``v`` | subdomain ``\Omega_i`` |
| vertex stalk: the agent's state | ``\mathbb R^{\Omega_i}``: the subdomain's copy of the field |
| edge stalk: the language two agents compare in | ``\mathbb R^{\Omega_i\cap\Omega_j}``: the unknowns they share |
| restriction map: what an agent reports on an edge | coordinate projection onto the overlap |
| global section: agreement | one field on the whole domain, ``H^0 \cong \mathbb R^n`` |
| Dirichlet energy ``\tfrac12\lVert\delta x\rVert^2`` | squared [`overlap_disagreement`](@ref) |
| pinned cells (targets) | the physical boundary, eliminated into ``f`` |
| a team of agents, fiber of ``\varphi : G \to H`` | an aggregate of subdomains |
| pushforward ``\varphi_* F`` and the nested tower | [`ExactPushforwardCoarseSpace`](@ref), [`TruncatedPushforwardCoarseSpace`](@ref) |
| edge coloring in the message-slot model | vertex coloring of the conflict graph ([`MulticolorSweep`](@ref)) |
| per-agent state copies in asynchronous diffusion | the 0-cochain iterate of [`ParallelSweep`](@ref) |

## Where they agree exactly

**Each local solve is a pinned harmonic extension.** The coordination
reference solves ``H q^\star = -B p``: free cells take the configuration of
least energy given the pinned cells. A Schwarz local solve with Dirichlet
transmission solves

```math
A_{\Omega_i\Omega_i}\, x_i = f|_{\Omega_i} - A_{\Omega_i\Gamma_i}\, g_i ,
```

which is the same problem with ``H = A_{\Omega_i\Omega_i}``,
``B = A_{\Omega_i\Gamma_i}`` and the neighbours' current values ``g_i`` as the
pins. A Schwarz method is therefore a team of agents, each repeatedly
computing its harmonic extension with its neighbours as moving targets. Three
facts from the coordination theory carry over unchanged:

- The *Dirichlet condition* for ``H \succ 0`` corresponds to
  ``A_{\Omega_i\Omega_i} \succ 0``.
- *Moving targets only change the right-hand side*, so each local problem is
  factored once with `ChordalLDLt` ([`LocalProblem`](@ref)) and reused every
  sweep.
- The Schwarz code is algebraic, so it applies to the coordination problem
  itself. The test set "coordination problems and the pushforward tower" solves
  ``H q = -B p`` for a random coordination sheaf on a ``6\times 6`` agent grid
  with teams of agents as subdomains. It recovers `harmonic_extension`'s answer
  with multicolor sweeps, with a coarse space, and with preconditioned CG.

**The nested tower is a Galerkin coarse problem.** Let ``B`` lift each
coarse vertex's coordinates to its fiber through the fiber-section bases (the
block-diagonal lift ``q = B_v q_H[v]`` of the nested tower, issue 007).

- Every column of ``B`` is a section on its fiber, so ``\delta_F B`` vanishes
  on edges inside a fiber.
- On a cross edge, the pushforward's restriction maps are ``F``'s maps composed
  with the fiber bases.

Together these give ``\lVert\delta_{\varphi_*F}\, q\rVert = \lVert\delta_F B q\rVert``, so

```math
L_{\varphi_* F} = B^\mathsf{T} L_F\, B ,
```

which the tests check numerically. The hierarchical solve of the tower is
therefore exactly one Galerkin coarse correction with prolongation ``B``, and
the energy gap ``E_{\text{hier}} \ge E_{\text{direct}}`` (issues 007 and 010) is
the error of a coarse solve used alone. Following that coarse correction with
local sweeps, as in [`SchwarzIteration`](@ref) with a coarse space, removes the
gap: the two-level iteration converges to the direct optimum. Teams can still
deform, at the cost of a few local solves per agent team.

[`ExactPushforwardCoarseSpace`](@ref) uses the full pushforward stalks.
[`TruncatedPushforwardCoarseSpace`](@ref) keeps a few modes of each stalk,
which is the usual way to keep a coarse problem small.

## Where they differ

**What is being minimized.** In coordination the sheaf energy *is* the
objective, and the sheaf Laplacian is the operator being solved. In Schwarz
the operator is the discretized PDE ``A``. The sheaf only constrains: the
iterate must become a global section. The Schwarz problem is

```math
\min_{x \in H^0(F)} \; \tfrac12\, u^\mathsf{T} A u - f^\mathsf{T} u,
\qquad u = \operatorname{glue}(x),
```

a quadratic objective over global sections, with an objective that is not a
sum of per-agent terms: ``A`` couples ``\Omega_i`` to ``\Gamma_i``. Coordination
problems are the special case where the objective is the sheaf energy itself.

**Algorithms.**

| Method | Per round, each agent | Analogue |
|---|---|---|
| sheaf diffusion ``x \leftarrow x - \gamma L x`` | one gradient step | Richardson iteration |
| [`ParallelSweep`](@ref) | an exact local minimization from old data | block Jacobi |
| [`MultiplicativeSweep`](@ref) / [`MulticolorSweep`](@ref) | the same, one color class at a time | block Gauss–Seidel |
| multifrontal tree solve | one elimination step along the clique tree | exact, non-overlapping |

Diffusion is cheap per round and slow to converge. The tree solve is exact but
fixes the communication pattern to the clique tree. Schwarz sits between them.
Each team solves exactly inside its subdomain, and the overlap width and a
coarse space control the iteration count.

## Corners

Two kinds of corners matter.

- **Re-entrant corners of the domain.** The [`notched_rectangle`](@ref) model
  problem has two corners of angle ``3\pi/2``, where the solution behaves like
  ``r^{2/3}``. The algebraic decomposition handles them without special code:
  [`box_partition`](@ref) drops boxes that fall inside the notch, and
  subdomains next to it simply have an irregular shape.
- **Cross points, where several subdomains meet.** Here discrete optimized
  Schwarz is fragile (Gander–Kwok). [`RobinTransmission`](@ref) adds a Robin
  term for each face, reads each face's data from the neighbour across it, and
  does not push values in the alternating sweeps. On box decompositions this
  removed the divergence for small ``p`` that we saw before.

## How the implementation fits the sheaf picture

Two earlier mismatches are now resolved:

- **Messages are restriction maps.** A Schwarz decomposition uses a
  ghost-layer cover ([`ghost_layer_cover`](@ref)). Each vertex stalk is the
  closed subdomain ``\overline\Omega_i = \Omega_i \cup \Gamma_i``. The boundary
  data of ``\Omega_i`` lies in ``\Gamma_i \cap \Omega_j \subset \overline\Omega_i \cap \overline\Omega_j``,
  the edge stalk shared with its owner ``j``. So receiving the ghost layer is
  the owner's restriction map followed by the adjoint restriction of ``i``. The
  Robin faces ``(m, k)`` lie in the same edge stalk. As a side effect,
  non-overlapping subdomains work too, as block Gauss–Seidel and block Jacobi.
- **Restriction maps need not be dense.** [`AbstractRestrictionMap`](@ref) is a
  LinearMaps-style interface (`mul!` with the map and its adjoint), with
  dense, sparse, selection and matrix-free implementations
  ([Restriction Maps](../api/restriction_maps.md)). `EuclideanSheaf{T,M}` stores
  maps of type `M`, and [`overlap_sheaf`](@ref) stores index lists. The
  coboundary is assembled sparsely or applied matrix-free with
  [`coboundary_operator`](@ref).

Remaining differences:

1. **Trivial cohomology.** With selections as restriction maps, ``H^0`` is just
   ``\mathbb R^n``. The rich structure of coordination sheaves (non-identity
   maps, nontrivial sections on fibers) enters only when Schwarz is applied to a
   coordination problem, as in the test above, and not through the overlap
   sheaf itself.
2. **No asynchrony yet.** Asynchronous diffusion keeps per-agent copies and
   updates on independent clocks. The parallel Schwarz iterate is already a
   cochain of per-subdomain copies, so asynchronous Schwarz, known to converge
   for M-matrices, is a natural next step.
3. **Conventions.** The IPM settings (Mumblebee) use `max_iter`. The Schwarz algorithms use
   `maxiter` and `tol`, and could be aligned.
