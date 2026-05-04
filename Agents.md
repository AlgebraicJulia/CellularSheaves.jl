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

For any linear system involving a sparse symmetric positive (semi)definite matrix,
use the `ChordalLDLt` path from `CliqueTrees.Multifrontal` — **not** base Julia's
`ldlt(Symmetric(Matrix(A)))`. The latter converts to dense and is significantly slower.

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

This pattern mirrors `nullspace_ldlt` in `EuclideanSheaves.jl` (lines 314–337),
which is the authoritative reference for LDLt usage in this codebase.

## Code Conventions

- **Sparse throughout.** Keep matrices sparse from construction to solve. Avoid
  `Matrix(A)` or `Array(A)` conversions in hot paths.
- **No new top-level dependencies** without checking `Project.toml` first.
- **`@argcheck`** (from ArgCheck.jl, already a dependency) for input validation at
  public function boundaries.
- **`BlockArray`** (from BlockArrays.jl) for returning cochains; partition by
  `s.vertex_stalks`.
- Vertex `v`'s slice in a cochain vector `x`:
  ```julia
  offsets = [0; cumsum(s.vertex_stalks)]
  x[offsets[v]+1 : offsets[v+1]]
  ```
- New exported functions go in `SheafInterface.jl` (abstract declaration) and their
  concrete implementation in the appropriate module file.
- Add `export` lines to `NetworkSheaves.jl`.
- New test files go in `test/network_sheaves/` and are included from `test/runtests.jl`.
- Do not add the package under test to the test dependencies `test/Project.toml`. This breaks the CI.

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

We are working off the work of researchers such as Hansen, Ghrist, Curry, Riess, Hanks, and Fairbanks, so explicitly refer to their papers in your references section and use nomenclature from their papers.


## What to Avoid

- **Do not add Catlab.jl** as a dependency for combinatorial graph operations.
- **Do not use dense solvers** (`ldlt(Symmetric(Matrix(A)))`) where the sparse
  `ChordalLDLt` path applies.
- **Do not write a single large issue** covering multiple independent features.
  Split at natural seams and list prerequisites explicitly.
- **Do not store issue drafts in ephemeral agent scratch files.** Write them to
  `docs/issues/` in the repository so they persist across sessions.
- **Do not leave the Out of Scope section empty.** Bounding the work is as important
  as specifying it.
