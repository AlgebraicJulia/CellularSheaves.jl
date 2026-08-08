# Agent Instructions for CellularSheaves.jl

This file captures preferences and conventions for AI agents working on this repository.
It is intended for use by Claude Code, GitHub Copilot, and similar tools.

## Project Overview

CellularSheaves.jl is a Julia package for cellular (network) sheaves on graphs.
It lives in the AlgebraicJulia ecosystem. The primary abstraction is `EuclideanSheaf`,
a sheaf whose stalks are Euclidean spaces and whose restriction maps are matrices.
The main application area is sheaf Laplacians, global sections, and functorial operations
(pushforward, morphisms, limits/colimits).

Developers of this package are generally fluent in applied category theory and sheaf theory. But users are not. So docstrings should explain some of the math, especially for externally facing functions.

Key source layout:
```
src/network_sheaves/
  EuclideanSheaves.jl    # EuclideanSheaf type, Laplacian, solvers
  Morphisms.jl           # SheafMorphism, ComplexMorphism, compose, is_morphism
  Pushforwards.jl        # pushforward_sheaf, all_fiber_bases, fiber_vertices
  GraphHomomorphisms.jl  # GraphHomomorphism type and helpers
  SheafInterface.jl      # Abstract interface — new operations declared here
  NetworkSheaves.jl      # Module aggregator — new includes/exports go here
test/network_sheaves/    # One test file per source module
docs/issues/             # Feature request drafts (markdown, numbered)
```

## Building the Docs

The following instructions are copied from the `.github/workflows/docs.yml` file to help you successfully build the docs
from the package root. You only need to Install docs dependencies the first time you try to build the docs.

```yaml
# Install Julia dependencies for the docs project.  This step is
# required so that the current
# development version of CellularSheaves is properly linked.
- name: Install docs dependencies
  run: julia --project=docs -e 'using Pkg; Pkg.develop(PackageSpec(path=pwd())); Pkg.instantiate()'

# Build documentation make sure to set the project to docs and run docs/make.jl from the package root.
- name: Build and deploy docs
  run: julia --project=docs docs/make.jl
```

These instructions are consistent with the standard julia package ecosystem as documented in 
the github action 

## Dependency Philosophy

**Prefer no new dependencies.** This is a focused linear-algebra package; heavyweight
categorical infrastructure is out of scope unless the core feature cannot be implemented
without it. Before adding a dependency, check whether the needed functionality can be
implemented in ~50 lines using `Graphs.jl`, `LinearAlgebra`, or `SparseArrays`, or `CliqueTrees` (all
already present).

Specific guidance:
- **Do not add Catlab.jl** for graph-combinatorial operations. Graph pushouts,
  homomorphisms, and related constructions should be implemented directly on
  `Graphs.SimpleGraph` using union-find or similar elementary algorithms.
- **Do not add LDLFactorizations.jl** — `CliqueTrees.Multifrontal` (`ChordalLDLt`)
  is already the package's sparse LDLt solver and should be used consistently.

## Solver Preferences

For linear systems involving a sparse symmetric positive (semi)definite matrix,
prefer the `ChordalLDLt` path from `CliqueTrees.Multifrontal` rather than base Julia's
`ldlt(Symmetric(Matrix(A)))`, since the latter converts to dense and is often significantly slower.
Small examples may still use the dense path where that behavior already exists, but new sparse-system
code should generally avoid densifying the matrix unless there is a clear reason to do so.

The factorization convention is `X = P' * L * D * L' * P`. To solve `X * v = b`:

```julia
using CliqueTrees.Multifrontal

function ldlt_solve(M, b)
    c = M.P' \ b    # c = P * b
    z = M.L  \ c    # forward solve
    w = M.D  \ z    # diagonal scaling
    y = M.L' \ w    # back solve
    return M.P \ y  # v = P⁻¹ * y
end

M   = ldlt!(ChordalLDLt(A), RowMaximum())
x   = ldlt_solve(M, b)
```

This pattern mirrors `nullspace_ldlt` in `EuclideanSheaves.jl`,
which is the authoritative reference for LDLt usage in this codebase.

## Code Conventions

- **Sparse throughout.** Keep matrices sparse from construction to solve. Avoid
  `Matrix(A)` or `Array(A)` conversions in hot paths.
- **No new top-level dependencies** without checking `Project.toml` first.
- **`@argcheck`** (from ArgCheck.jl, already a dependency) for input validation at
  public function boundaries.
- **`BlockArray`** (from BlockArrays.jl) for returning cochains; partition by
  `vertex_stalks(s)`.
- Vertex `v`'s slice in a cochain vector `x`:
  ```julia
  offsets = [0; cumsum(vertex_stalks(s))]
  x[offsets[v]+1 : offsets[v+1]]
  ```
- New exported functions go in `SheafInterface.jl` (abstract declaration) and their
  concrete implementation in the appropriate module file.
- Add `export` lines to `NetworkSheaves.jl`.
- New test files go in `test/network_sheaves/` and are included from `test/runtests.jl`.
- Do not add the package under test to the test dependencies `test/Project.toml`. This breaks the CI.
- Prefer multiple dispatch over runtime case statements and union types. 

## Documentation Structure Conformance

When updating the docs, make them conform to the package's current module structure.

- Keep the `makedocs(modules=...)` list in `docs/make.jl` aligned with the modules and
  submodules whose docstrings are referenced from the documentation.
- In `docs/src/api.md`, include a separate `@autodocs` block for every documented
  module or submodule. Do not assume a parent module's autodocs block is enough.
- Prefer this "one autodocs block per module" pattern whenever modules are added,
  renamed, or reorganized; otherwise Documenter cross-references can fail to resolve.
- When adding a new module that should appear in the docs, update `docs/make.jl` and
  `docs/src/api.md` in the same change.
- When you run the docs build, make sure that there are no missing cross references.
  Missing cross references will produce warnings in `docs/make.jl` rather than failing the build.

## Feature Request Format

Feature requests are stored as numbered markdown files in `docs/issues/`. When drafting
a new feature request, include all of the following sections — agents working from these
files should be able to implement the feature without asking follow-up questions.

### Required sections

**Mathematical Background** — precise definition of the construction, not just a name.
Include the relevant commutative diagrams or equations. Agents do not know the
project's mathematical conventions; spell them out.

**Codebase Orientation** — a table of files to read before writing any code, with a
one-line explanation of *why* each file is relevant. Include line numbers for key
functions when helpful.

**Requested Implementation** — exact function signature with a complete docstring,
plus an implementation sketch. The sketch should be concrete enough that the agent
cannot choose a wrong algorithm. For linear-algebra steps, include the matrix
expressions explicitly.

**Tests to Write** — name the test file, describe the concrete setup (specific graph,
specific stalk dimensions), and list individual `@test` statements. Tests serve as
acceptance criteria; vague test descriptions produce vague tests.

**Verification Checklist** — a bulleted checklist the agent can tick off, ending with tests passing.

**Out of Scope** — explicit list of related things the agent should *not* implement.
Without this, agents expand scope.

### Issue sizing

Each issue should cover one self-contained function or a small cluster of tightly
related helpers. If drafting an issue and it requires implementing a non-trivial
sub-problem (e.g. graph pushout as part of sheaf pushout), split into two issues
with an explicit prerequisite link. Smaller issues are easier to assign, review,
and verify.

## Mathematical Context

The owner is fluent in sheaf theory, category theory, and applied linear algebra.
Do not over-explain the mathematics, but do be precise about which variant of a
construction is intended (e.g. "pushout over a fixed graph" vs. "pushout over a span
of graph homomorphisms" are different features).

When the owner's request is mathematically ambiguous, ask for clarification before
proposing an implementation. When multiple implementations are possible (e.g. different
solver backends), briefly name the trade-off and recommend one rather than implementing
all of them.

We are working off the work of researchers such as Hansen, Ghrist, Curry, Riess, Hanks, and Fairbanks, so explicitly cite relevant papers where appropriate and use nomenclature from their papers.

## Git & Development Workflow Safety

- **Atomic Git Commits After Verification**: Commit working changes to git immediately after verifying clean test execution (`test/runtests.jl`) and clean documentation builds (`docs/make.jl`). Keeping commits atomic and frequent ensures a safe checkpoint is always available if quota limits or session interruptions occur.
- **Preserve User Content & Existing Code**: Never blindly overwrite, truncate, or delete user files, custom literate example scripts, or markdown documentation. Always view and inspect existing content first, preserve the user's comments and structure, and make targeted modifications (`replace_file_content`) rather than blanket file rewrites.

## Control Sheaves & Trajectory Tracking Conventions

- **Reference State Assembly**: When using `JointTikhonovFilter` for joint position and velocity tracking, $x_{\text{ref}}$ MUST be assembled from both `filter.x` (position target) and `filter.v` (velocity target mapped to velocity state indices `x_ref[v_idxs] = filter.v[v_idxs]`). Leaving velocity components as zero causes LQR velocity feedback to damp movement and increase tracking lag.
- **Hierarchical Reference Solves**:
  - **1-Stage (Position $q^*$)**: Solves $H q^* = -L_{IB} p(t)$. Standard LQR feedback.
  - **2-Stage (Velocity $\dot{q}^*$)**: Solves $H \dot{q}^* = -L_{IB} \dot{p}(t)$. Eliminates 1st-order formation tracking lag (~30× error reduction).
  - **3-Stage (Acceleration $\ddot{q}^*$)**: Solves $H \ddot{q}^* = -L_{IB} \ddot{p}(t)$. Pre-computes attitude tilt ($\theta_{\text{ref}} = \ddot{x}^*/g, \phi_{\text{ref}} = -\ddot{y}^*/g$) and thrust boost ($u_{1,\text{ff}} = m \ddot{z}^*$), eliminating 2nd-order centripetal lag (~sub-millimeter precision).
- **High-Level Problem Ergonomics**: Control sheaf modules (such as `Layered`) must provide high-level constructors that accept domain specifications (e.g. `LayeredEscortSpec`, `AbstractDynamicsSpec`, target trajectories) and automatically assemble the underlying `EuclideanSheaf`, `GraphHomomorphism`, pushforward sheaf, and fiber bases. Advanced users can still pass pre-built sheaves when overriding defaults.
- **Multiple Dispatch over `isa` Branching**: Avoid runtime `if spec isa Type` branching in configuration helpers. Define separate method signatures via multiple dispatch (e.g. `get_agent_dynamics_config(spec::HomogeneousDynamics, ...)`).
- **Pre-computed Factorizations in Fiber Bases**: Pre-compute matrix pseudoinverses/factorizations during initialization (e.g. `target_subbases_inv` in `LayeredFiberBases`) to replace matrix backslash divisions (`\`) with matrix-vector multiplications (`*`) in simulation loops. Avoid `Matrix(B_mat)` conversions on sparse block matrices.
- **Practical Sensor Limitations & Derivative Orders**:
  - Finite-differencing noisy sensor data amplifies noise at $O(1 / \Delta t^{2n})$. High-order finite differences ($n \ge 3$, Jerk/Snap) cause motor chatter and instability on real sensor streams.
  - For real hardware and sensor streams, cap reference derivatives at 1st/2nd order (Position, Velocity, Acceleration) filtered via low-pass Tikhonov observers ($\epsilon \approx 0.02 - 0.05\text{ s}$). Higher-order Jerk/Snap feedforward should be restricted to closed-form analytical trajectories or smooth offline splines.

## What to Avoid

- **Do not add Catlab.jl** as a dependency for combinatorial graph operations.
- **Do not add this repo** as a dependency for the test suite. This breaks the CI.
- **Do not use dense solvers** (`ldlt(Symmetric(Matrix(A)))`) where the sparse
  `ChordalLDLt` path applies.
- **Do not write a single large issue** covering multiple independent features.
  Split at natural seams and list prerequisites explicitly.
- **Do not store issue drafts in ephemeral agent scratch files.** Write them to
  `docs/issues/` in the repository so they persist across sessions.
- **Do not leave the Out of Scope section empty.** Bounding the work is as important
  as specifying it.
- **Do not put line comments inside of function definitions in literate examples.**
  This breaks the literate pipeline. Use `##` for comments inside compound blocks or put comments above the function definition.
- **Do not use Java style getters and setters.** Prefer idiomatic julia names that do
  not contain the substrings get or set.
- **Do not import modules inside of function definitions** Import and export expressions should only appear towards the top of a file.