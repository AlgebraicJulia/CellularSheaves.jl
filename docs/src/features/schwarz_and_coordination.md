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

## Sheaf ADMM and Robin transmission

Hanks, Riess, Cohen, Gross, Hale and Fairbanks (arXiv:2504.02049, Algorithm 1;
code in AlgebraicOptimization.jl) solve homological programs

```math
\min_{x \in C^0} \sum_i f_i(x_i) \quad\text{subject to}\quad x \in H^0
```

by ADMM in consensus form. All of ``x``, the copies ``z`` and the scaled
multipliers ``y`` live on vertex stalks:

```math
x_i \leftarrow \operatorname{argmin}_{x_i} f_i(x_i) + \tfrac{\rho}{2}\lVert x_i - z_i + y_i\rVert^2,
\qquad z \leftarrow \Pi_{H^0}(x + y),
\qquad y_i \leftarrow y_i + x_i - z_i .
```

In the paper the projection is computed by sheaf diffusion run to convergence
(their Theorem 2). The AlgebraicOptimization.jl code computes it with CG.

[`SheafADMM`](@ref) runs this iteration on the ghost-layer overlap sheaf.
[`local_objectives`](@ref) splits the PDE energy exactly into convex quadratics
``f_i(x_i) = \tfrac12 x_i^\mathsf{T} K_i x_i - b_i^\mathsf{T} x_i`` on the
closed subdomains: each matrix edge term is shared equally by the stalks that
contain it. For this sheaf ``\Pi_{H^0}`` averages the copies of each dof, which
is a single exchange with the neighbours.

**The two local solves have the same shape.**

| | Robin (optimized Schwarz) | Sheaf ADMM |
|---|---|---|
| local operator | ``A_{\Omega_i\Omega_i} - N_i``: full operator, Dirichlet coupling removed on interface rows | ``K_i``: the subdomain's share of the energy on its closed stalk |
| penalty | ``p_{ij}`` on each interface face | ``\rho`` on the whole stalk (`penalty = :stalk`) or on shared dofs (`:shared`) |
| data | neighbour's current values on ``\Gamma_i`` and at the interface | ``z_i - y_i``: average of copies minus multiplier |
| needs overlap | yes (one layer at least) | no |
| convergence | not proven for every ``p``; fast at ``p^*`` | proven for every ``\rho > 0`` with an exact projection |

Both are "Neumann-type local operator plus penalty", and in our runs the best
``\rho`` was the optimized Robin parameter ``p^*`` (in the units of `A`). That
matches the classical picture in which the multiplier plays the role of the
interface flux, and ADMM on a consensus splitting is a Douglas–Rachford method
(Gabay 1983), like Lions' Robin method (Lions–Mercier 1979).

**Optimized Schwarz converges much faster.** Iterations to a relative residual of
``10^{-8}`` for Poisson problems with overlap 1 (unit square on a ``32 \times 32``
grid; notched rectangle on a ``31 \times 15`` grid with ``4 \times 2`` boxes):

| decomposition | Dirichlet (parallel) | Robin, ``p^*`` | ADMM, best ``\rho \approx p^*`` |
|---|---|---|---|
| 4 strips | 95 | 19 | 248 |
| ``4\times 4`` boxes | 189 | 21 | 400 |
| notched rectangle | 51 | 25 | 254 |

The difference is in what a local solve knows. A Robin subdomain has the true
operator on its whole interior and reads its neighbours' values *and fluxes*
directly. An ADMM subdomain holds only part of the stiffness on shared dofs
and sees its neighbours only through averages. The flux has to be learned
iteratively in ``y``. Penalizing only shared dofs (`penalty = :shared`) changes
little.

Replacing the exact projection by a single diffusion step
(`projection_steps = 1`), as a communication-light variant of Algorithm 1,
converged on strips but diverged on box grids, where dofs carry different
numbers of copies. The convergence proof assumes the exact projection.

**What each side can borrow.**

- ADMM could take Robin-type local problems: keep the full operator on the
  interior and penalize only the interface faces. That would be optimized
  Schwarz with ADMM's multiplier update, i.e. an augmented-Lagrangian
  formulation of optimized Schwarz.
- Schwarz could use ADMM's guarantee. ADMM converges with zero overlap and for
  every ``\rho``, which makes it a safe fallback where the Robin iteration is
  delicate.

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
