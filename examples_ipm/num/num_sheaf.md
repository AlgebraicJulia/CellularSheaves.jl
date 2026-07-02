# num_sheaf.jl — network utility maximization

Kelly's framework as a sheaf on the source–link bipartite graph, and
the corpus's cleanest instance of *the dual is deployed engineering*:
dual decomposition of this sheaf — link prices adjusting to congestion,
sources reacting to path price — **is** distributed TCP congestion
control (Kelly–Maulloo–Tan 1998; Low–Lapsley 1999; Chiang–Low–
Calderbank–Doyle 2007). After LMPs, Kantorovich potentials, and the
obstacle reaction measure, the session's fourth named dual is the one
running on the machine you're reading this on.

## 1. Sheaf and the fairness–cone ladder

Cells: source rates x_s and link slacks σ_l (`PositiveCone`), one
capacity hyperedge per link (Σ_{s∋l} x_s + σ_l = c_l, capacities in
g), and utility leaves per mode. The α-fairness family (Mo–Walrand)
maps onto cone types:

| mode | α | utility | cones | protocol flavor |
|---|---|---|---|---|
| `:throughput` | 0 | w·x | LP | — (starves long paths) |
| `:proportional` | 1 | w·log x | `ExponentialCone` leaves (u,1,t), the `mle_spline` pattern verbatim | TCP Vegas |
| `:delay` | 2 | −w/x | `SecondOrderCone` leaves: q·x ≥ w ⟺ (q+x, q−x, 2√w) ∈ Q³ | TCP Reno |
| `:maxmin` | ∞ | shared floor | LP (`CofreeCone` t + slack leaves) | — |

## 2. Verification (num_oracle.py, CLARABEL)

Kelly's linear network — one long flow across L = 3 unit-capacity
links, three shorts — has closed forms at **every** rung:

```
mode           x_long     closed form           value: |sheaf − direct|
throughput     0.000000   0                     5.8e-9   (x non-unique: LP)
proportional   0.250003   1/(L+1)   = 0.25      4.7e-8   |Δx| 1.8e-5
delay          0.366016   1/(1+√L) = 0.366025   1.9e-9   |Δx| 1.2e-5
maxmin         0.500000   1/2                   5.6e-10  (x non-unique: LP)
```

**Kelly's identity from the duals**: w_s/x_s = Σ_{l∈path} λ_l to
9.3×10⁻⁵, with λ = 4/3 = (L+1)/L matching the analytic price; the α=2
version w/x² = path price holds to 1.4×10⁻³. On a random 8×6 mesh,
Jain's index runs 0.597 (throughput) < {0.683 delay, 0.681
proportional} < 0.787 (maxmin) — **measured non-monotone between α = 1
and 2**: Jain's index is not the α-fair order (the floor and ceiling
hold; the middle does not sort). H¹ = 0 (capacity block full-rank;
leaves are wires).

Note the LP-degeneracy honesty: throughput and maxmin have non-unique
optimal *rates* (ties among shorts; free division above the floor), so
the assembled-vs-direct comparison is in value — where agreement is
10⁻⁹ — and the ε block again acts as min-norm selection.

## 3. Field notes

* First runs: `kelly_ladder_demo()` (four solves, four closed forms —
  the whole file in eight lines of output), `price_demo()` (Kelly's
  identity from `res.y`; take |λ| pending the library's dual sign
  convention, as with the LMPs), `mesh_demo()`.
* The delay mode's SOC trick (hyperbolic constraint qx ≥ w as a
  rotated cone) is the standard bridge for all *rational* utilities;
  general α needs power cones, which the library lacks — the honest
  scope boundary, stated rather than worked around.
* Dual decomposition as an *algorithm* (price iteration = simulated
  TCP) is deliberately not implemented — the IPM solves the same
  problem centrally; the md's claim is about what the dual *is*, and
  the price identity is its checkable shadow.
* Structural profile: pure hyperedge coupling (no restriction maps
  beyond wires), H¹ = 0, tiny cells, three cone types across modes —
  after obstacle_option's dense-Q extreme, this is the opposite corner:
  all structure in B and the cones, none in Q.

## References

Kelly, Maulloo, Tan, J. Oper. Res. Soc. **49**:237 (1998). Low &
Lapsley, IEEE/ACM ToN **7**:861 (1999). Mo & Walrand, IEEE/ACM ToN
**8**:556 (2000). Chiang, Low, Calderbank, Doyle, Proc. IEEE **95**:255
(2007). Srikant, *The Mathematics of Internet Congestion Control*.

## Files

* `num_sheaf.jl` — standalone: four-mode builder, price reader, demos,
  JuMP reference.
* `num_oracle.py` — self-contained; all numbers above.
