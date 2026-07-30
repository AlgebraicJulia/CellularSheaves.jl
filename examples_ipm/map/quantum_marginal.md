# quantum_marginal.jl — certified ground-energy bounds from cluster RDMs

The 2-RDM row of `local_consistency_sheaf.md` §4, realized: reduced
density matrices as PSD cells, **partial trace** as the restriction
map, and a certified lower bound on the ground-state energy of an XXZ
spin chain/ring as the optimum. This is the spin-lattice relaxation of
Barthel & Hübener (PRL **108**, 200404, 2012; parallel: Baumgratz &
Plenio, NJP **14**, 023027) — the field calls it a quantum marginal /
local-consistency relaxation, and it is the corpus's first example
whose restriction maps are genuinely **multi-Kraus channels** rather
than compressions.

## 1. The problem

XXZ chain/ring of $N$ spin-½ sites, bond Hamiltonian
$h = \tfrac14(\sigma^x\!\otimes\!\sigma^x + \sigma^y\!\otimes\!\sigma^y
+ \Delta\,\sigma^z\!\otimes\!\sigma^z)$, total
$H = \sum_i h_{i,i+1}$. Ground energy $E_0 = \min_\rho \langle H,\rho\rangle$
over global states is exponentially hard; but $\langle H,\rho\rangle$
reads only the two-site marginals. Relax: keep one **cluster RDM**
$\rho_C$ per window of $\ell$ consecutive sites, ask only that
neighbouring windows agree on their $(\ell{-}1)$-site overlap, and
minimize the same energy. Every true global state's marginals are
feasible, so the optimum is a **strict lower bound on $E_0$** at any
$\ell$ — the natural counterpart to variational (upper-bound) methods,
and exactly Barthel–Hübener's construction.

**Real symmetric WLOG.** $\sigma^y\!\otimes\!\sigma^y$ is real, so the
XXZ Hamiltonian is real in the computational basis, and the real part
of any feasible Hermitian family is feasible at the same energy. The
library's real `SemidefiniteCone` therefore loses nothing here. (A
complex Hamiltonian would need the standard $2\times$-size real
embedding.)

## 2. The sheaf

| sheaf datum | quantum object |
|---|---|
| vertex $C_i$ | cluster RDM $\rho_i$ on sites $i..i{+}\ell{-}1$, svec coords, `SemidefiniteCone` |
| edge $(C_i, C_{i+1})$ | $\mathrm{Tr}_{\text{first}}\,\rho_i = \mathrm{Tr}_{\text{last}}\,\rho_{i+1}$ (the shared $(\ell{-}1)$-site window) |
| pin | $\mathrm{tr}\,\rho_1 = 1$ — **one** pin; the others are implied, since agreement transports the trace ($\mathrm{tr}\,\rho_i = \mathrm{tr}\,\mathrm{Tr}_{\text{first}}\rho_i = \mathrm{tr}\,\rho_{i+1}$) |
| objective | $c_i = \mathrm{svec}(H_{C_i})$, using $\langle H,\rho\rangle = \langle\mathrm{svec}\,H, \mathrm{svec}\,\rho\rangle$; each bond owned by exactly one cluster |
| $Q$ | $\varepsilon I$ (Uzawa's SPD requirement); the exact $\varepsilon$-optimum obeys $v_\varepsilon \le E_0 + \varepsilon\,n_{\mathrm{cl}}/2$ (evaluate at the true marginals, $\|\rho\|_F \le \mathrm{tr}\,\rho = 1$), so $E_0 \ge v_\varepsilon - \varepsilon\,n_{\mathrm{cl}}/2$ — `qm_bound` applies the correction |

Partial trace is the CPTP channel with Kraus family
$\{I\otimes\langle k|\}_{k=0,1}$ — two Kraus operators, the
multi-Kraus case §4.2 of the consistency note promised. Non-consecutive
overlaps need no edges: marginals of marginals commute, so agreement
between $C_i$ and $C_{i+2}$ is implied through $C_{i+1}$.

**H¹, measured then explained.** Chain: 0. Ring: exactly **1**, at
every $\ell$ tested. Mechanism: $\mathrm{PT}^\top_{\text{first}}(Y)$ is
the functional $\langle I\otimes Y,\cdot\rangle$ and
$\mathrm{PT}^\top_{\text{last}}(Y)$ is $\langle Y\otimes I,\cdot\rangle$;
a left-null cycle needs $I\otimes Y_i = Y_{i-1}\otimes I$ around the
ring, and $I\otimes Y = Y'\otimes I$ forces $Y = cI$ — only the trace
functional survives the shift holonomy. So on rings the IPM dual lives
on a 1-dimensional coset and Uzawa's CG picks the basepoint, exactly
the behavior catalogued in the tensor-mesh files.

## 3. The point: a theorem fails

The corpus's classical law (`local_consistency_sheaf.md` §5) is
*positive local sections that agree on overlaps glue iff the pattern is
chordal / a junction tree* — odd cycles break the orthant, chordless
4-cycles break PSD completion, trees break nothing. Quantumly that law
is **false on a path**, the simplest and most chordal topology there
is:

> Pin both two-site RDMs of a 3-site chain to the Bell state
> $\Phi^+$. Each one-site marginal is $I/2$; the single overlap
> agreement holds **exactly** (sheaf residual 0). Yet no 3-qubit state
> has both as its marginals — the independent global-extension SDP
> returns *infeasible*. Monogamy of entanglement is the obstruction,
> and deciding consistency of local density matrices is QMA-complete
> (Liu 2006).

The energy shadow of the same fact: on the 3-site Heisenberg chain the
$\ell{=}2$ relaxation puts **both** bonds in the singlet
(their shared marginals are $I/2$ — consistent!), certifying $-3/2$
where the truth is $-1$. The gluing obstruction has left the base
graph's topology and moved **into the cone**: same coboundary, same
tree, new failure mode. This is the P-vs-CP ladder of the consistency
note's §4.2, now with teeth.

## 4. Verification (quantum_marginal_oracle.py, CLARABEL)

```
partial-trace svec maps: 8.88e-16   XXZ reality check: 0.0e+00

energy ladder  (ring, certified E_sdp ≤ E_exact; per-bond in brackets)
  XXZ Δ=1.0  N=10 ring:  exact  -4.51544635 [-0.451545]
    ℓ=2: E_sdp  -7.50000001 [-0.750000]   |sheaf−native| 1.51e-08   ‖Bx−g‖ 1.51e-14  -(2+Δ)/4·N = -7.500000
    ℓ=3: E_sdp  -5.00000002 [-0.500000]   |sheaf−native| 1.68e-08   ‖Bx−g‖ 7.58e-14
    ℓ=4: E_sdp  -5.00000000 [-0.500000]   |sheaf−native| 6.91e-07   ‖Bx−g‖ 1.03e-12
  XXZ Δ=0.5  N=10 ring:  exact  -3.81903278 [-0.381903]
    ℓ=2: E_sdp  -6.25000000 [-0.625000]
    ℓ=3: E_sdp  -4.21535167 [-0.421535]
    ℓ=4: E_sdp  -4.21535169 [-0.421535]

TI single-cluster ladder (per bond; thermodynamic exact Δ=1: ¼−ln2 ≈ −0.443147):
  Δ=1.0: ℓ=2: -0.750000  ℓ=3: -0.500000  ℓ=4: -0.500000  ℓ=5: -0.467129  ℓ=6: -0.467129
  Δ=0.5: ℓ=2: -0.625000  ℓ=3: -0.421535  ℓ=4: -0.421535  ℓ=5: -0.394672  ℓ=6: -0.394672

3-site Heisenberg chain: exact E0 = -1.000000   ℓ=2 bound = -1.500000   (monogamy gap 0.5000)

monogamy witness (double Bell on a 3-site chain — a PATH, trivially chordal):
  sheaf overlap residual  : 0.00e+00   (locally consistent ✓)
  global-extension SDP    : infeasible   (no global state — gluing FAILS)

H1 = rows − rank(B):  chain 0 (ℓ=2,3);  ring 1 (ℓ=2,3)
```

Three measured facts worth keeping:

* **$\ell{=}2$ captures no physics beyond one bond.** The bound is
  $-(2{+}\Delta)/4$ per bond *exactly*: every cluster sits
  independently in its bond ground state, and the one-site overlaps —
  all $I/2$ — object to nothing. The whole relaxation is doing the
  monogamy-blind thing globally that the 3-site witness does minimally.
* **The ladder improves only at odd $\ell$** (both $\Delta$): pairs
  $\{3,4\}$ and $\{5,6\}$ give bit-identical bounds. Ring optimum
  equals the translation-invariant single-cluster optimum (cyclic
  symmetrization of any ring solution is feasible with the same
  energy), verified numerically, which is what makes the large-$\ell$
  ladder cheap to compute. The even-step stall is reported as a
  measured fact; we did not chase the mechanism.
* $\ell{=}5$ reaches $-0.4671$ per bond against the thermodynamic
  $\tfrac14 - \ln 2 \approx -0.4431$ — the slow-but-certified
  convergence Barthel–Hübener report, reproduced in the sheaf idiom.

## 5. Field notes

* Cluster growth is brutal: $\rho_C$ is $2^\ell{\times}2^\ell$, svec
  dimension $2^{\ell-1}(2^\ell{+}1)$; $\ell = 5$ means 528 coordinates
  per vertex. The literature's route to tighter bounds at fixed
  $\ell$ is *more conditions* (the fermionic P/Q/G hierarchy, symmetry
  reduction), not larger clusters — deliberately out of scope here, and
  the reason the fermionic 2-RDM flavor was not built.
* The ring's $H^1 = 1$ makes it (after the tensor 2×2 grid) the
  cheapest live probe of Uzawa on consistent rank-deficient $B$ —
  and the first with the deficiency caused by *channel holonomy*
  rather than mesh corners.
* First Julia runs: `qm_monogamy_demo()` (tiny; the −3/2 vs −1 gap is
  the whole story in one solve), then
  `qm_ladder_demo(N = 10, Δ = 1.0)` — exact diagonalization is done
  in-file by dense `eigvals` (fine to $N \approx 12$).
* A rigorous certificate independent of solver tolerance would read the
  IPM's dual $(y, d)$ as an explicit lower-bound witness; the demos
  report the $\varepsilon$-corrected primal value, which is certified
  up to the duality gap the solver leaves.

## References

Barthel & Hübener, PRL **108**, 200404 (2012). Baumgratz & Plenio,
NJP **14**, 023027 (2012). Y.-K. Liu, *Consistency of local density
matrices is QMA-complete* (2006). Mazziotti PRL **93**, 213001 (2004)
for the fermionic variational 2-RDM line (not built). Anderson bounds
for the classical ancestor of the cluster-ground-state bound.

## Files

* `quantum_marginal.jl` — standalone builder (no spline-family
  includes), demos, JuMP reference.
* `quantum_marginal_oracle.py` — self-contained; all numbers above.
