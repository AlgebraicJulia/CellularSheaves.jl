# 003 — MathOptInterface wrapper for `CellularSheaves.IPM`

## Summary

Expose the interior-point conic-QP solver in `CellularSheaves.IPM` through
[MathOptInterface](https://jump.dev/MathOptInterface.jl/) so it is drivable from
JuMP as `Model(CellularSheaves.IPM.Optimizer)`. This is unblocked by issue 002:
`MOI.ObjectiveFunction{ScalarQuadraticFunction}` lowers to an arbitrary symmetric
`P`, which the solver could not faithfully represent while `Q` was constrained to
be block-diagonal. With general `Q` shipped, a faithful QP wrapper is now possible.

Two decisions (fixed):
- **Packaging:** `MathOptInterface` becomes a direct dependency of the top-level
  `CellularSheaves` package; the wrapper is a submodule of `IPM` exporting
  `IPM.Optimizer`. (Migrate to a weakdep package-extension before any registration.)
- **Scope:** all supported cones in the first pass — `Reals`, `Zeros`,
  `Nonnegatives`, `SecondOrderCone`, `PositiveSemidefiniteConeTriangle`,
  `ExponentialCone`, `PowerCone`, plus quadratic objective.

## The organizing observation: our form is already standard-form conic

The solver's native primal/dual pair (`src/IPM/src/ipm/problem.jl:6-16`) is

```
(P)  min ½ pᵀQp − cᵀp   s.t.  B p = g,   p ∈ K
(D)  max −½ pᵀQp + gᵀy   s.t.  Q p − c − Bᵀy = d ∈ K*
```

Variables `p` live **directly** in a product cone `K`, coupled by a linear
equality `Bp=g`. This is exactly the shape produced by *slacking* every conic
constraint of a JuMP model — and it is how the solver's own examples are already
built (`build_convex_clustering`: free `xᵢ` blocks + SOC slack blocks `sₑ` coupled
by `B`). So the MOI mapping is not novel; it mirrors the existing usage.

Contrast with Clarabel's geometric form (`min ½xᵀPx+qᵀx s.t. Ax+s=b, s∈K`, `x`
free). Hypatia's form (`min cᵀx s.t. Ax=b, Gx+s=h, s∈K`) is closer — it separates
equalities from cones exactly as we do — and its wrapper is the primary template.

## Result recovery is complete — no gap

`solve(prob, settings)` returns `HSDResult`/`IPMResult` (`src/IPM/src/ipm/result/`)
with `p`, `y`, `d` **de-homogenized (÷τ), un-equilibrated (`unscale!`), and
un-permuted** back to the caller's original variable order (`ldiv!(·, s.P, ·)` in
`result(...)`). `VariablePrimal` ← `p`, equality `ConstraintDual` ← `y`, conic
`ConstraintDual`/`ConstraintPrimal` ← `d`/slack, all as direct slices. The only
caveat: `d` on `CofreeCone` blocks is identically 0 by construction.

## Mathematical Background

### Slack reformulation (MOI model → standard form)

A JuMP model presents free variables `x` and constraints `fₖ(x) ∈ Sₖ`. Build the
solver's `p = [x ; s₁ ; s₂ ; …]`:

- Each MOI scalar variable → its own **1-dimensional `CofreeCone` vertex** (v1
  decision: maximally sparse; quadratic cross-terms `xᵢxⱼ` become off-diagonal
  `Q` blocks between vertices, leaning on issue 002. Grouping free variables into
  wider blocks is a deferred performance pass, not a correctness concern).
- Each `Sₖ = Zeros` constraint `fₖ(x) = 0` → equality rows of `B` directly.
- Each cone constraint `fₖ(x) ∈ Cₖ` → a slack vertex `sₖ ∈ Cₖ` and equality rows
  `fₖ(x) − sₖ = 0` in `B`.
- `Q`, `c` from the objective; `Q`'s off-diagonal blocks carry the quadratic
  cross-terms between free variables (the issue-002 capability).

Assembly convention (Hypatia `wrapper.jl:426-427`): a row of `B` stores the
**negated** VAF coefficients, and `g` stores the **added** VAF constants, so that
`Bp=g` reproduces `fₖ(x) − sₖ = 0` with the MOI function value equal to the slack.

### Cone coordinate conventions (verified against source)

| MOI set | IPM cone (`src/IPM/src/cone/`) | Mapping |
|---|---|---|
| `Reals` | `CofreeCone` (noc.jl) | direct; dual is `Zeros` (`d≡0`) |
| `Zeros` | — | equality rows of `B` |
| `Nonnegatives` | `PositiveCone` (pos.jl) | direct, self-dual, dim = n |
| `SecondOrderCone` | `SecondOrderCone` (soc.jl) | **direct** — scalar-first `[t;u]`, `t≥‖u‖` |
| `PowerCone{α}` | `PowerCone(α)` (tdc/pow.jl) | **direct**, IPM `α` = MOI exponent, order `(x,y,z)` |
| `PositiveSemidefiniteConeTriangle` | `SemidefiniteCone` (sdp.jl) | √2 off-diagonal **matches**, but **lower↔upper svec permutation** for n≥3 |
| `ExponentialCone` | `ExponentialCone` (tdc/exp.jl) | **coordinate permutation** `(u,v,w) → (w,v,u)` |

Details:

- **PSD.** IPM svec (`smat!`/`svec!` sdp.jl:80-116) is **lower-triangle,
  column-major**, off-diagonals scaled by √2:
  `(1,1),(2,1),(3,1),(2,2),(3,2),(3,3)`. MOI `PositiveSemidefiniteConeTriangle`
  is **upper-triangle, column-major** (= lower-triangle *row-major* by symmetry),
  *unscaled*: `(1,1),(2,1),(2,2),(3,1),(3,2),(3,3)`. Accepting
  `MOI.Scaled{PositiveSemidefiniteConeTriangle}` gets the √2 scaling for free
  (MOI applies it); a static per-side-dimension **index permutation** then maps
  MOI order → IPM order. For n≤2 the orders coincide. Precompute the permutation
  once per PSD block; apply it when writing that block's `B` rows and invert it
  when slicing its dual/primal out. This is the Hypatia `transform.jl`
  `needs_permute`/`permute_idxs`/`untransform_affine` pattern.
- **Exponential.** IPM `(x₁,x₂,x₃)` with `x₂·log(x₁/x₂) ≥ x₃` equals MOI
  `(u,v,w)` with `w ≥ v·exp(u/v)` under `(u,v,w) → (x₁,x₂,x₃) = (w,v,u)`. Apply
  the 3-coordinate reversal at assembly and recovery. (The internal tdc Cholesky
  pivot `1,2,3→2,3,1` is *not* an external reordering — do not apply it.)
- **`DualExponentialCone`/`DualPowerCone`** have no native IPM type; rely on MOI
  bridges or leave unsupported (declare only the primal cones).

## Codebase Orientation

| File | Why it matters |
|---|---|
| `src/IPM/src/ipm/problem.jl` | `IPMProblem(Q,B,c,g,K)` target type + ctor asserts (per-vertex dims, one cone/vertex, diagonal-Q-arc requirement) |
| `src/IPM/src/ipm/result/{hsd,ipm}.jl` | result fields `p,y,d,status,niter,pobj,dobj` — the recovery surface |
| `src/IPM/src/ipm/ipm.jl` | `@enum IPMStatus` (12 values) → `MOI.TerminationStatus` map |
| `src/IPM/src/settings/{ipm,hsd}.jl` | `IPMSettings`/`HSDSettings` fields → `RawOptimizerAttribute`s + `Silent` |
| `src/IPM/src/cone/{noc,pos,soc,sdp}.jl`, `cone/tdc/{exp,pow,tdc}.jl` | cone types + coordinate conventions |
| `src/BlockSparseArrays/src/block_sparse_matrix.jl` | `blocksparse(rid,cid,blk)`, `block`, `colrange`, `nvtxs` — how to assemble `B`/`Q` from triplets |
| `src/IPM/src/utils.jl` | `allocblockdiag(B)` — seeds `Q` with a diagonal block per vertex |
| `src/IPM/examples/convex_clustering.jl` | worked slack-reformulation example; the benchmark the wrapper should subsume |
| `.julia/packages/Clarabel/*/src/MOI_wrapper/` | copy_to / attributes / status-dict template (geometric form) |
| `.julia/packages/Hypatia/*/src/MathOptInterface/` | equality-vs-cone split + svec transform hooks (primary template) |

## Requested Implementation

New submodule `src/IPM/src/moi/MOIWrapper.jl` (`include`d from `IPM.jl`), plus
`using MathOptInterface` in `CellularSheaves`'s Project. Follow the copy_to-only
pattern (declare `MOI.supports_incremental_interface(::Optimizer) = false`; JuMP
wraps in a caching optimizer and supplies bridges).

### 1. Optimizer struct

```julia
mutable struct Optimizer{T} <: MOI.AbstractOptimizer
    settings::Union{IPMSettings{T}, HSDSettings{T}}  # HSD default (does infeasibility)
    raw::Dict{String, Any}                            # RawOptimizerAttribute overrides
    result::Union{Nothing, HSDResult{T}, IPMResult{T}}
    sense::MOI.OptimizationSense
    objconstant::T
    # copy_to bookkeeping (built fresh each solve):
    var_cols::Dict{MOI.VariableIndex, Int}            # variable → p-column
    eq_rows::Vector{UnitRange{Int}}                   # Zeros-constraint row ranges in B/g/y
    cone_blocks::Vector{UnitRange{Int}}               # slack cone column ranges in p/d
    moi_cones::Vector{MOI.AbstractVectorSet}          # original set per cone (for untransform)
    nvars::Int
end
Optimizer{T}(; kw...) where {T} = ...   # route kwargs through RawOptimizerAttribute
Optimizer(; kw...) = Optimizer{Float64}(; kw...)
```

### 2. `copy_to(dest::Optimizer, src::MOI.ModelLike)`

Algorithm:
1. Assign `var_cols` in `ListOfVariableIndices` order; `nvars = length`.
2. Read objective (`ScalarAffineFunction`/`ScalarQuadraticFunction`): accumulate
   `c` (affine, note the `−cᵀp` sign), `Q` triplets (quadratic; symmetrize into
   full-symmetric blocks per the 002 storage convention), `objconstant`; negate
   for `MAX_SENSE`.
3. Walk `Zeros` constraints → append equality rows to `B` triplets and `g`
   (negate coeffs, add constants); record each range in `eq_rows`.
4. Walk each cone constraint type → allocate a slack vertex `sₖ ∈ Cₖ`; append
   `fₖ(x) − sₖ = 0` rows; push `Cₖ` onto `K`; record slack column range in
   `cone_blocks` and the MOI set in `moi_cones`. Apply the PSD/Exp coordinate
   transform to the slack rows here.
5. Free-variable vertices: one **1-dim `CofreeCone` vertex per scalar variable**.
   Ensure every vertex (free and slack) gets a diagonal `Q` arc (seed with
   `allocblockdiag`-style structural zeros; off-diagonal `Q` blocks from the
   quadratic objective are added on top).
6. `IPMProblem(Q, B, c, g, K)`.

`optimize!(dest)` = `dest.result = solve(prob, dest.settings)`.

### 3. Supported sets / objective

```julia
MOI.supports_constraint(::Optimizer{T}, ::Type{MOI.VectorAffineFunction{T}},
    ::Type{<:Union{MOI.Zeros, MOI.Nonnegatives, MOI.SecondOrderCone,
        MOI.Scaled{MOI.PositiveSemidefiniteConeTriangle},
        MOI.ExponentialCone, MOI.PowerCone{T}}}) = true
MOI.supports(::Optimizer{T}, ::MOI.ObjectiveFunction{MOI.ScalarQuadraticFunction{T}}) = true
# also ScalarAffineFunction objective, VariableIndex; ObjectiveSense.
```

### 4. Result getters + status map

`TerminationStatus`, `PrimalStatus`, `DualStatus`, `VariablePrimal` (`result.p`),
`ConstraintDual` for `Zeros` (`result.y[eq_rows[k]]`), `ConstraintDual`/
`ConstraintPrimal` for cones (`untransform(moi_cones[k], result.d[...] / slack)`),
`ObjectiveValue` (`sense·(result.pobj) + objconstant`), `DualObjectiveValue`,
`SolveTimeSec`, `BarrierIterations` (`niter`), `RawStatusString`.

Status map: `OPTIMAL→OPTIMAL`, `NEAR_OPTIMAL→ALMOST_OPTIMAL`,
`PRIMAL_INFEASIBLE→INFEASIBLE`, `NEAR_PRIMAL_INFEASIBLE→ALMOST_INFEASIBLE`,
`DUAL_INFEASIBLE→DUAL_INFEASIBLE`, `NEAR_DUAL_INFEASIBLE→ALMOST_DUAL_INFEASIBLE`,
`ITERATION_LIMIT→ITERATION_LIMIT`, `STALLED→SLOW_PROGRESS`,
`NUMERICAL_FAILURE→NUMERICAL_ERROR`, `ILL_POSED→INVALID_MODEL`,
`NEAR_ILL_POSED→OTHER_ERROR`.

### 5. Attributes

`MOI.Silent` ↔ `verbose` (0/1); `RawOptimizerAttribute(name)` ↔ the matching
`*Settings` field (`feas_tol`, `gap_tol`, `itmax`, `scale_itmax`, `step_frac`,
`near_factor`, `stall_tol`, and HSD-only `illposed_tol`/`infeas_abs`/`infeas_rel`).
No `TimeLimitSec` field exists — declare unsupported. `SolverName = "CellularSheaves.IPM"`.

## Tests to Write

```julia
# Differential vs Clarabel-through-MOI on random programs (the strong oracle).
for _ in 1:20
    model_ipm = Model(IPM.Optimizer); model_cla = Model(Clarabel.Optimizer)
    build_random_qp!(model)      # LP, QP, SOCP, PSD, Exp, Power variants
    @test objective_value(ipm) ≈ objective_value(cla) rtol=1e-6
    @test value.(x_ipm)        ≈ value.(x_cla)        rtol=1e-5
end

# Each cone in isolation, dual correctness.
@test MOI.get(opt, MOI.ConstraintDual(), soc_con) ≈ known_dual   rtol=1e-6
@test MOI.get(opt, MOI.ConstraintDual(), psd_con) ≈ known_dual   rtol=1e-6  # exercises the permutation

# PSD svec permutation: n≥3 off-diagonal placement round-trips.
@test recovered_matrix ≈ X_star   rtol=1e-6

# The convex-clustering program built via JuMP matches build_convex_clustering.
@test objective_value(Model(IPM.Optimizer) |> cc) ≈ res.pobj + 0.5*sum(abs2, A)

# MOI.Test conformance subset (supported cones only).
MOI.Test.runtests(BRIDGED, CONFIG; include=[...])
```

## Verification Checklist

- [ ] `MathOptInterface` added to `CellularSheaves` `Project.toml` `[deps]` + `[compat]`.
- [ ] `Model(IPM.Optimizer)` solves an LP, a QP (with cross-terms), and an SOCP to Clarabel agreement.
- [ ] PSD: `Scaled{PSDTriangle}` accepted; the lower↔upper svec permutation verified for n=3,4 (primal *and* dual).
- [ ] Exponential-cone `(u,v,w)→(w,v,u)` permutation verified against Clarabel/known duals.
- [ ] Power cone α round-trips.
- [ ] `ConstraintDual` correct for both `Zeros` (from `y`) and cones (from `d`); `CofreeCone` dual is 0.
- [ ] Status map covers all 11 returnable `IPMStatus` values; infeasible/unbounded models report the right `TerminationStatus`.
- [ ] `Silent`, tolerance `RawOptimizerAttribute`s take effect.
- [ ] `convex_clustering.jl` benchmark rebuilt on the wrapper (retire the duplicated `build_cc_jump`).
- [ ] A `MOI.Test.runtests` subset (supported cones) passes under a bridged/cached optimizer.

## Out of Scope

- Incremental (`add_constraint`) interface — copy_to-only; JuMP's caching layer covers it.
- `DualExponentialCone`/`DualPowerCone` native support (bridge or reformulate).
- `TimeLimitSec` (no solver time limit exists).
- Migration to a weakdep package-extension (do before registration, not now).
- Efficient free-variable block grouping / exploiting model sparsity in the vertex
  partition beyond correctness — a follow-up performance issue.
