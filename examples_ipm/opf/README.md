# AC Optimal Power Flow Benchmarks

AC-OPF as a ladder of covers over the power network graph. Each relaxation rung is a choice of cover—the sheaf base graph is the cover's nerve.

## Quick Start

```bash
# Full extended benchmark (Rung 1, Rung 2, SE, Bowtie)
julia --project=. examples_ipm/opf/opf_sheaf.jl

# With Mosek
julia --project=. examples_ipm/opf/opf_sheaf.jl --mosek
```

## The Ladder of Relaxations

| Rung | Cover | Cells | Tight When |
|------|-------|-------|------------|
| 0 | rank-1 W | Nonconvex | Ground truth (NP-hard) |
| 1 | lines | 4-dim Lorentz per line | **Radial networks** |
| 2 | lines ∪ cycles | + Hermitian W_C per cycle | **Chordal networks** |
| 3 | cliques | PSD clique blocks | Often (Lavaei-Low) |

**Key insight**: The dropped constraint is **angular holonomy** (Σ ±θ_e = 0 around cycles). Rung 1 keeps each line's 2×2 physics and drops exactly this.

## Cone Mix

The fullest cone mix in the corpus under one coboundary:
- `SecondOrderCone` (line cells)
- `CofreeCone` (bus cells, generators)
- `PositiveCone` (box slack leaves)
- `SemidefiniteCone` (cycle cells at rung 2+)
- Native quadratic Q (generator cost)

## Benchmarks

### OPF Dispatch (Rung 1 & 2)

Standard OPF with linear or quadratic generation cost.

| Config | IPM vs Mosek |
|--------|--------------|
| Rung 1 SOCP | IPM wins 2-3x |
| Rung 2 +SDP | Mosek wins 2-5x |
| Quadratic cost | IPM wins 1.1-1.3x |

### State Estimation (SE)

WLS state estimation over the SOCP feasible set. The objective is entirely **dense measurement Grams**—the corpus's first dense Q on SecondOrderCone cells.

| Config | IPM vs Mosek |
|--------|--------------|
| SE Rung 1 | **IPM wins 2-10x** |
| SE Rung 2 | IPM wins 1.3-3x |

### Bowtie (Chordal Overlap)

Two triangles sharing edge (1,2)—the simplest test where cycle cells overlap in a live off-diagonal, not just magnitudes.

## Test Battery

The file includes a comprehensive test battery (`opf_test_battery()`):

1. **Monotonicity**: v₁ ≤ v₂ ≤ v₃ for all seeds
2. **Tightness theorems**: radial ⇒ v₁=v₃, single-cycle ⇒ v₂=v₃
3. **Holonomy resolution**: defect → 0 at rung 2+
4. **Holonomy-gap correlation**: ρ ≈ 0.88 (defect tracks gap)
5. **Load sweeps**: gap/defect grow together with stress
6. **Rank-1 recovery**: physics round-trip verification

```julia
include("examples_ipm/opf/opf_sheaf.jl")
opf_test_battery()
```

## References

- Jabr, IEEE TPS 21:1000 (2006) — SOCP relaxation
- Lavaei & Low, IEEE TPS 27:92 (2012) — SDP relaxation
- Kocuk, Dey, Sun, Oper. Res. 64:1177 (2016) — Cycle strengthening
- Sojoudi & Lavaei — Tree-exactness
