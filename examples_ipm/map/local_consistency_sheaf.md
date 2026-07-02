# Local‑consistency inference as a positive section of a stochastic cellular sheaf

This note explains the problem that `map_lp.jl` solves and why it is, *by construction*, an instance of the corpus's normal form: the constrained global section of an ordered‑vector‑space cellular sheaf, solved as a conic program

$$
\min_{\mu}\ \tfrac12\,\mu^\top Q\,\mu + c^\top\mu
\quad\text{s.t.}\quad B\mu = g,\ \ \mu\in K=\textstyle\prod_v K_v,
\qquad B=\delta\ \text{(the sheaf coboundary)}.
$$

It is the branching, topologically non‑trivial companion to the safe‑corridor note: the corridor lives on a chain, where the sheaf is tidy but decorative; here the sheaf's cohomology is the whole story.

---

## 1. The problem

A pairwise **Markov random field** on a graph $\mathbb G=(\mathcal V,\mathcal E)$ assigns a discrete variable $x_i\in\{1,\dots,c_i\}$ to each node and a potential to each node and edge. **MAP inference** seeks the assignment of least energy,

$$
\min_{x}\ \sum_{i\in\mathcal V}\theta_i(x_i)+\sum_{(i,j)\in\mathcal E}\theta_{ij}(x_i,x_j),
$$

an NP‑hard combinatorial problem in general. The standard convex handle (Schlesinger; Wainwright–Jordan) is to lift each discrete variable to a **local belief** — a distribution $\mu_i$ per node and $\mu_{ij}$ per edge — and minimise the *linear* functional $\langle\theta,\mu\rangle$ over the **local polytope** $\mathcal L(\mathbb G)$:

$$
\mu\ge 0,\qquad
\sum_{x_i}\mu_i(x_i)=1\ \ \forall i,\qquad
\sum_{x_j}\mu_{ij}(x_i,x_j)=\mu_i(x_i),\qquad
\sum_{x_i}\mu_{ij}(x_i,x_j)=\mu_j(x_j).
$$

The last two families — **marginalisation** (edge beliefs agree with node beliefs where they overlap) — are the only coupling in the problem; everything else is per‑cell (nonnegativity, normalisation, the cost). That single fact is what makes this a sheaf.

$\mathcal L(\mathbb G)$ is an **outer bound** on the true **marginal polytope** $\mathcal M(\mathbb G)$ (the convex hull of genuine global joint distributions restricted to the cells). Optimising over $\mathcal L$ instead of $\mathcal M$ is tractable but *loose*: the two coincide exactly when the region complex is a junction tree (treewidth $\le 1$ for this first‑order polytope), and on a frustrated cycle the containment is strict. The gap $\mathcal L\supsetneq\mathcal M$ is the phenomenon this note is about.

---

## 2. The sheaf

Build a cellular sheaf on the **region graph** of the MRF.

| sheaf datum | MAP‑LP object |
|---|---|
| cell (vertex) $\{i\}$ | a node belief $\mu_i$, stalk $\mathbb R^{c_i}$ |
| cell (vertex) $\{i,j\}$ | an edge belief $\mu_{ij}$, stalk $\mathbb R^{c_ic_j}$ |
| base‑edge $\{i\}\subset\{i,j\}$ | an overlap; edge stalk $\mathbb R^{c_i}$ |
| restriction from $\{i,j\}$ | **marginalise** (sum out $x_j$) |
| restriction from $\{i\}$ | the **identity** (the node belief *is* that marginal) |
| coboundary $\delta$ | $(\delta\mu)_e=(\text{marg }\mu_{ij})-\mu_i$ |
| $\ker\delta = H^0$ | overlap‑consistent belief families |
| cone $K_v$ on a stalk | nonnegativity (PositiveCone) |
| normalisation pin | a degree‑1 edge, $\mathbf 1^\top\mu_i=1$ (RHS $1$) |
| objective | $Q,\,c$ — the cost carried on the cells |

So $\delta\mu = g$ *is* the local polytope's equality block: $\delta\mu=0$ on the marginalisation edges, $=1$ on the normalisation pins. Because the marginalisation maps are surjective, the **linear** system $\delta\mu=g$ is always solvable — indeed the **uniform independent** distribution ($\mu_i=1/c_i$, $\mu_{ij}=1/(c_ic_j)$) is a strictly interior point, a Slater point exactly analogous to the corridor's feasible reference. The normalisation pins play the role of the corridor's terminal jet pins (and, equivalently, of the apex/homogeniser that broadcasts the constant $1$).

Two book‑keeping notes. Edge normalisation is implied — $\sum_{ab}\mu_{ij}=\sum_a\sum_b\mu_{ij}=\sum_a\mu_i=1$ once marginalisation holds — so only node normalisation is imposed. And no base‑edge touches more than two cells (marginalisation edges are degree‑2, pins are degree‑1), which is the cellular‑sheaf‑on‑a‑graph discipline the solver requires.

### Restriction maps are Markov kernels

Marginalisation is the *deterministic* special case of a general admissible restriction map. In an **ordered** cellular sheaf the restriction maps must be **positive** — they must carry the positive cone of a stalk into the positive cone of the edge, so that the pushforward of a positive section stays positive. For the simplex the mass‑preserving positive maps are exactly the **column‑stochastic** matrices, i.e. **Markov kernels / channels**. Marginalisation is the 0/1 kernel that lumps a joint state $(a,b)$ onto its first coordinate $a$; a general column‑stochastic $S$ is a *soft* coarse‑graining or a noisy observation channel.

A base‑edge in the general form joins cells $u,w$ with a shared edge stalk $\mathbb R^{m}$ and two column‑stochastic maps $P,R$, and demands

$$
P\,\mu_u \;=\; R\,\mu_w
$$

— the two cells agree after each is pushed through its own channel to the shared observable. This is precisely the empirical‑model reading of Abramsky–Brandenburger: each context reports a distribution, and consistency is agreement on overlaps. `map_lp.jl` validates column‑stochasticity on construction and ships both the marginalisation sheaf (`mrf_spec`) and a soft coarse‑graining sheaf (`coarse_graining_spec`).

One qualitative change from marginalisation is worth flagging. Marginalisation is surjective and holonomy‑free, so its *linear* $H^1$ vanishes and the only obstruction is positivity. A general kernel can be non‑surjective, or — for permutation/relabeling kernels — carry non‑trivial holonomy around a cycle, so with stochastic maps the "no global section" phenomenon can arise **linearly** as well, before positivity enters. The coarse‑graining example stays feasible (agreement on a shared observable is transitive around a cycle); the synchronisation flavour is where linear holonomy would bite.

---

## 3. The objective dial: linear → native quadratic (SparseMAP)

The coboundary and cones fix the feasible set; the objective is an independent axis.

**Linear ($Q=0$).** $\min\langle\theta,\mu\rangle$ over $\mathcal L$ is the MAP‑LP. Its optimum is an extreme point of $\mathcal L$, which on a frustrated cycle is a *fractional pseudomarginal*.

**Native quadratic ($Q=\rho I$).** Adding $\tfrac{\rho}{2}\lVert\mu\rVert^2$ is **SparseMAP** (Niculae–Martins): completing the square,

$$
\min_{\mu\in\mathcal L}\ \langle\theta,\mu\rangle+\tfrac{\rho}{2}\lVert\mu\rVert^2
\;=\;\tfrac{\rho}{2}\min_{\mu\in\mathcal L}\big\lVert\mu-(-\theta/\rho)\big\rVert^2,
$$

so the solution is the **Euclidean projection** of $-\theta/\rho$ onto the sheaf‑cut local polytope. It is the structured analogue of *sparsemax*: strictly convex (a **unique** optimum), differentiable in $\theta$ (usable inside end‑to‑end learning), and **sparse** — the projection lands on a low‑dimensional face, producing exact zeros in the active set. As $\rho\to0$ it recovers the LP vertex; as $\rho\to\infty$ it is pulled to the min‑norm interior point (uniform‑independent).

Crucially, **the quadratic is not lifted to a cone.** It could be written as a rotated second‑order‑cone constraint, but a solver that carries $Q$ natively should not: the KKT system's only new term is $-(Q+H)$ in each cell's diagonal block (native curvature $Q$ plus the barrier's cone scaling $H$), both cell‑block‑diagonal, so $\delta$'s bipartite sparsity is untouched. Lifting instead manufactures an auxiliary cell and a wiring edge per quadratic term and replaces an exactly‑differentiated term with a barrier approximation. The rule of thumb: **lift only the cost whose curvature the solver cannot already see.** A quadratic has constant, known curvature ($Q$) and rides in the objective; the barrier's job is the non‑smooth boundary of $K$.

---

## 4. The cone dial: one assembly, a graded family

The cone $K_v$ encodes the local convex structure; the coboundary‑agreement backbone is invariant. Swapping the cone (and the matching "shared readout" kernel) gives a family of problems in the *same* assembly.

| cell object | cone | overlap coboundary | a global section is | tight when |
|---|---|---|---|---|
| local marginal $\mu_R$ | orthant / simplex | marginalise / channel | global joint distribution | junction tree |
| pseudo‑moment $M_R$ | PSD (SDP) | shared‑block consistency | representing measure (PSD completion) | chordal |
| reduced density $\rho_R$ | PSD, trace 1 | partial trace | $N$‑representable state | chordal |
| entropic MAP | exponential | marginalise (+ $-H(\mu)$ cost) | Gibbs distribution | junction tree |

Two remarks. The PSD row is an exact twin of the simplex story: chordal‑sparsity SDP decomposition rewrites one large PSD constraint as local PSD blocks on cliques that must agree on overlaps, and "does a global PSD section exist" is **PSD matrix completion**, solvable for every consistent local specification iff the pattern is **chordal** (Grone–Johnson–Sá–Wolkowicz) — the SDP analogue of "junction tree." The *shape* of the obstruction changes with the cone, though: the anti‑ferromagnetic LP frustrates on **odd** cycles (2‑colourability), while PSD completion frustrates on **chordless** cycles of length $\ge4$ (a triangle is chordal, so the smallest PSD witness is the $4$‑cycle). And the objective dial of §3 is orthogonal to this table: any row may carry a native quadratic.

### 4.1 The PSD row in full: the moment tightening

Take $\pm1$ variables for cleanliness. Where the LP held a *distribution* on each region's configurations, the SDP holds a **pseudo‑moment matrix** $M_R$: over the region's variables form the vector of degree‑$\le1$ monomials $\phi=(1,x_i,x_j,\dots)$, and $M_R$ is the pseudo‑expectation $\tilde{\mathbb E}[\phi\phi^\top]$ — its $(1,1)$ entry is $\tilde{\mathbb E}[1]$, its degree‑1 entries the means $\tilde{\mathbb E}[x_i]$, its off‑diagonals the pair moments $\tilde{\mathbb E}[x_ix_j]$. An edge region $\{i,j\}$ carries a $3\times3$ matrix indexed by $(1,x_i,x_j)$ — exactly the quantities the pairwise energy reads — and a node region $\{i\}$ a $2\times2$ indexed by $(1,x_i)$. The stalk over a cell is thus the space of symmetric matrices of that size, and the cone $K_R$ is the **PSD cone**, replacing the orthant.

The per‑cell affine data mirror §2. The normalisation pin becomes $M_R[1,1]=1$ (the constant monomial — the same "broadcast the $1$" role the apex/homogeniser played), and a **localising** constraint pins the diagonal: $x_i^2=1$ gives $M_R[x_i,x_i]=1$ (for $\{0,1\}$ variables it would read $x_i^2=x_i$, diagonal $=$ mean).

The restriction map is **principal‑submatrix extraction**: to restrict $M_R$ to a shared variable set $S=R\cap R'$, take the principal block on the monomials in $S$, a compression $M\mapsto PMP^\top$ with $P$ a coordinate projection. A compression carries PSD to PSD — it is the positive map that makes this an ordered sheaf with the PSD order — and it is the moment‑level image of marginalisation. The coboundary keeps its shape: "the two cells agree where they overlap," now meaning $M_R$ and $M_{R'}$ carry the **same shared block**. For an edge region $\{i,j\}$ and node region $\{i\}$ that block is the $2\times2$ on $(1,x_i)$, so agreement forces both to report the same $\tilde{\mathbb E}[x_i]$ — precisely the moment image of the old marginalisation constraint.

The MAP energy is affine in these moments (node potential linear in $\tilde{\mathbb E}[x_i]$; edge potential linear in $\tilde{\mathbb E}[x_i],\tilde{\mathbb E}[x_j],\tilde{\mathbb E}[x_ix_j]$), so the objective is a linear functional $\sum_R\langle C_R,M_R\rangle$ — the conic program verbatim, cone swapped to $\prod_R\mathbb S^{k_R}_+$, with an optional native $\tfrac\rho2\lVert M\rVert_F^2$ on the independent objective dial.

The section question now has **two layers**. First, do the overlap‑agreeing local PSD blocks assemble into a single global PSD matrix (one moment matrix over all monomials whose principal blocks are the $M_R$)? That is **PSD matrix completion**, and it succeeds for every consistent local specification iff the pattern is **chordal** (Grone–Johnson–Sá–Wolkowicz) — gluing. Second, is that global PSD matrix the moment matrix of an actual $\pm1$ measure? That needs a **rank / flat‑extension** condition — exactness, the SDP counterpart of the LP's integrality. Chordality handles gluing; flatness handles realisability, and the two are kept apart exactly as "$\delta M=g$ solvable" was kept apart from "$\mathcal L=\mathcal M$." This tightening is Sherali–Adams degree $1$ (the local polytope) $\subset$ Lasserre degree $2$ (the moment block): the PSD constraint implies all the linear consistency $\mathcal L$ enforces and then couples the moments in a way marginalisation cannot — *provided* a region is large enough to hold the correlation being coupled. A clique region $\{i,j,k\}$ carrying one $4\times4$ moment matrix binds the three pairwise correlations jointly; three *edge* regions that overlap only on single nodes share only means and leave the correlations decoupled, so on a triangle they fall back to the LP (§5). Coupling strength is thus set by region size, completability by region‑graph topology. The algorithmic fact behind the sheaf, not merely an analogy, is **chordal SDP decomposition** (Fukuda–Kojima–Murota–Nakata; Vandenberghe–Andersen): one large PSD constraint on a chordal pattern *is* a family of small clique‑PSD constraints that must agree on overlaps, with the clique tree as the junction tree and running intersection as the gluing condition. Goemans–Williamson MAX‑CUT is the degenerate single‑region instance (one global PSD block, no overlaps to glue).

### 4.2 Stochastic maps become channels

Do the general column‑stochastic kernels of §2 carry over to moments? The moment embedding $\mathcal A:\mu\mapsto M(\mu)=\sum_x\mu(x)\,\Phi_x$ (with $\Phi_x=\phi(x)\phi(x)^\top$ a rank‑one PSD atom) is **linear**, so pushing a distribution through a column‑stochastic $S$ induces a unique linear $\mathcal T$ with $\mathcal T(\Phi_x)=\sum_y S(y,x)\,\psi(y)\psi(y)^\top$, and $\mathcal A_{\text{coarse}}\!\circ S=\mathcal T\circ\mathcal A_{\text{fine}}$ commutes by construction. Each image is a **convex** combination of coarse rank‑one atoms (columns of $S$ sum to $1$), hence PSD. So on **genuine** moment matrices $M=\sum_x\mu(x)\Phi_x$ with $\mu\ge0$, every stochastic kernel induces a **positive** map — functoriality of the moment embedding, automatic.

The catch is that the SDP relaxes *off* the moment cone: its variables are PSD matrices that need not be moments of any nonnegative $\mu$. For such a pseudo‑moment, $\mu=\mathcal A^{-1}(M)$ can have **negative** entries, and then $\mathcal T(M)$ is a *signed* combination of PSD atoms — it can leave the PSD cone. Naive reuse of a column‑stochastic matrix is thus positive only on real moments; extending it to pseudo‑moments reintroduces exactly the flatness condition the SDP was escaping.

What survives on *all* pseudo‑moments is the **completely positive, trace‑preserving (CPTP) maps — quantum channels** — read in the adjoint direction as unital CP. By Choi–Kraus each is $M\mapsto\sum_k A_kMA_k^\top$ with $\sum_k A_k^\top A_k=I$; trace‑preservation is the moment analogue of columns‑sum‑to‑one, keeping $M[1,1]=1$ coherent around overlaps. The dictionary against the LP kernels:

| LP kernel (on distributions) | PSD channel (on moments) |
|---|---|
| marginalise (drop a variable) | **partial trace** |
| deterministic relabel (permute states) | **unitary / isometric conjugation** $M\mapsto UMU^\top$ (single Kraus) |
| soft coarse‑graining / noisy channel | **pinching / dephasing / depolarising** (multi‑Kraus) |

This is the reframing worth keeping: the honest SDP analogue of "column‑stochastic map on the simplex" is "CPTP channel on the PSD cone," and a given $S$ lifts to a legitimate restriction map exactly when it can be realised as such a channel. The clean ones are conjugations (marginalisation and relabeling, single‑Kraus, valid on the whole cone); the mixing ones need you to *supply* a Kraus form rather than read one off.

And there is a genuinely new phenomenon with no LP shadow. Classically, positivity and complete positivity **coincide** over a commutative algebra, so on the simplex there was nothing to distinguish — every positive map on distributions is automatically CP. The moment lift is noncommutative, and there **positive‑but‑not‑CP** maps exist (the transpose is the canonical one). So "does this map induce a positive map on moments" bifurcates into *positive on the cone* (true always, on real moments) versus *completely positive / channel* (valid on pseudo‑moments, and stable under the tensoring and composition the sheaf performs around overlaps). That P‑vs‑CP gap is the extra difficulty making the quantum‑marginal / $N$‑representability version strictly harder than PSD completion, itself strictly harder than the LP — the same ladder, now in the **morphisms** rather than the objects.

Two grades of non‑triviality follow, on independent knobs. **Kraus count**: a selector (principal‑submatrix extraction) is a *single*‑Kraus channel; genuine mixing (pinching, depolarising, soft coarse‑graining) is *multi*‑Kraus. **Non‑commutativity**: a permutation/projection is "classical," while a single‑Kraus but **basis‑changing** isometry $M\mapsto UMU^\top$ that rotates rather than permutes the monomials is deeply non‑trivial in the SDP despite being one operator — and has no classical shadow, since the only mass‑preserving conjugations on distributions are permutations. If the aim is a restriction map that exercises the part of the problem the LP could never see, the basis‑changing isometries are the sharp instrument; a cycle of them whose composition is $\ne I$ manufactures a genuinely quantum **holonomy** obstruction with no positivity or classical analogue. Feasibility is then a design constraint: either push a strictly‑positive reference $M_\star$ through the chosen channels and set $g_e$ to the mismatch (the corridor's "$g=\delta$ of a feasible reference," guaranteeing a Slater point for *any* channels), or keep $g=0$ homogeneous and choose the channels compatible enough to share a common input — the failure of which *is* the holonomy witness.

---

## 5. Positive global sections, and why the relaxation is loose

Two Hodge‑theoretic freedoms sit inside the conic program: the primal indeterminacy is $H^0=\ker\delta$ (belief families that agree everywhere), the dual indeterminacy is $H_1=\ker\delta^\top$ (the cycle space), and exactness $[g]=0\in H^1$ is what makes the primal feasible and the dual coset well‑defined at once.

But — and this is the point — the **linear** cohomology is *not* where the difficulty lives. For marginalisation, $\delta\mu=g$ is always solvable over signed measures. The difficulty is entirely **positivity**: the local polytope $\mathcal L$ contains nonnegative, normalised, overlap‑consistent pseudomarginals that are *not* the restriction of any global distribution. In sheaf language: **positive local sections that agree on overlaps do not glue to a positive global section.** That is exactly the Abramsky–Brandenburger characterisation of contextuality — "no global section" — and it closes precisely when the region complex is acyclic enough (junction tree / running intersection).

So the duality note's phrase *"positive global section"* is doing real work here for the first time. In the corridor it was trivial: a chain has no cycles and its overlaps are always jointly realisable, so both the topology and the positivity gap were empty. Here "positive" and "global" genuinely fail to commute with "local + agreeing," and the failure is the reason MAP inference is hard.

### Worked witness: the frustrated triangle

Three binary variables on a triangle, with an anti‑ferromagnetic potential that costs `penalty` for agreement on every edge. The odd cycle is not 2‑colourable, so any integer assignment must agree on at least one edge:

- **Integer MAP** $=1$.
- **MAP‑LP** $=0$: put every edge belief on the anti‑diagonal $\mu_{ij}=(0,\tfrac12,\tfrac12,0)$ and every node belief at $(\tfrac12,\tfrac12)$. This satisfies $\delta\mu=g$ *exactly* (marginalisation residual $0$), yet it claims each pair disagrees with probability $1$ — impossible for a genuine joint on an odd cycle. It is a fractional vertex of $\mathcal L\setminus\mathcal M$: **locally consistent, no global section.** The gap of $1$ is the certificate.
- **SparseMAP** ($Q=\rho I$): the unique projection. For small‑to‑moderate $\rho$ it coincides with that fractional vertex — same $\mathcal L\setminus\mathcal M$ point, now with $6$ exact zeros (the agreement states) as a sparse active set; as $\rho$ grows it dissolves into the dense uniform‑independent interior.

The even $4$‑cycle, by contrast, is 2‑colourable, so $\mathcal L=\mathcal M$ there (gap $0$, integral) — the gap tracks frustration, i.e. the topology of the region complex, exactly as the junction‑tree theorem predicts. These numbers are reproduced independently, without the Julia library, by `map_lp_oracle.py`.

**The witness moves graphs.** Under the PSD cone (§4.1) the *positive‑global‑section* question — does a family of overlap‑agreeing local PSD blocks complete to one global PSD moment matrix — is governed by **chordality**, and it inverts against the LP. The frustrated triangle that broke the LP is a $3$‑cycle, a complete graph, hence **chordal**: the gluing question is trivial and there is **no section obstruction**. The chordless even $4$‑cycle — which was LP‑*tight* — is exactly where it fails: unit‑diagonal correlations $(0.8,0.8,0.8,-0.8)$ around the cycle are **pairwise valid** (each $2\times2$ block PSD) yet admit **no global PSD completion**. Positive local sections agreeing on overlaps, no positive global section — the sheaf obstruction, sitting precisely where the orthant relaxation saw nothing wrong. Odd cycles for the orthant, chordless $\ge4$‑cycles for the PSD cone.

Keep this *section* (feasibility) question apart from the relaxation's *value*, which is a separate axis. On the triangle the LP relaxes the energy to $0$; the single‑clique Shor SDP (the whole triangle as one region, a $4\times4$ moment matrix) tightens it to $\tfrac34$ — the Goemans–Williamson bound — still short of the integer $1$; and the *pairwise* moment sheaf, coupling the three edge blocks only through their shared node means, leaves it back at $0$. So value‑tightening depends on the region (clique) size, while completability depends on the region‑graph topology; the two are orthogonal, and the sheaf's cohomological content is the second one. All three SDP numbers are reproduced by `map_sdp_oracle.py`.

---

## 6. Relation to the corpus normal form

Everything above lands in the same conic program the corridor did:

- **cost on the cells** — linear $\langle\theta,\mu\rangle$, or native quadratic $\tfrac\rho2\lVert\mu\rVert^2$ (a first‑class objective, not a lifted cone);
- **cones on the cells** — orthant here, PSD/exp for the graded family;
- **coupling = pure overlap agreement** — the coboundary $\delta$ of a stochastic cellular sheaf, with marginalisation the deterministic special case of a Markov kernel;
- **inhomogeneity via pins/apex** — normalisation as degree‑1 pins, $g$ exact by construction with the uniform‑independent Slater point.

What distinguishes this instance from the corridor is topology. The corridor is a chain: $H^1=0$ automatically, no positivity gap, the sheaf tidy but not load‑bearing. The region graph of an MRF branches and carries cycles, and the local‑vs‑global gap — the whole reason the relaxation is loose — is the positive‑global‑section obstruction living on those cycles. Swapping the cone renames the global object (joint distribution → representing measure → PSD completion → $N$‑representable state) without changing the topology of why it can fail to exist.

---

## Code

- **`map_lp.jl`** — the builder in the `safe_corridor.jl` idiom: stochastic‑kernel cells (`StochEdge`/`Pin`/`SheafSpec`), `mrf_spec` (marginalisation) and `coarse_graining_spec` (general channels), `build_ipm_problem(spec; rho)` assembling `IPMProblem(Q,B,c,g,K)` with $Q=\rho I$ for SparseMAP, a JuMP reference, and `verify_mrf` for the $\mathcal L\supsetneq\mathcal M$ report.
- **`map_lp_oracle.py`** — a self‑contained NumPy/SciPy/quadprog oracle (no dependency on the Julia package) that reproduces the LP and SparseMAP optima, the integer MAP, the frustration/tightness pattern across cycles, and the feasibility of a general column‑stochastic sheaf.
- **`map_sdp.jl`** — the PSD sibling (§4.1–4.2): clique moment matrices as `SemidefiniteCone` cells stored in the library's scaled `svec`, `shared_selection` as the compression restriction, diagonal pins for $M[1,1]{=}1$ and $x_i^2{=}1$, `region_cost` expanding $\pm1$ potentials into an affine moment objective, `moment_spec(inst; cliques)` (edge cliques by default, maximal cliques for the Shor tightening), `build_sdp_problem(spec; rho)`, and `completion_spec` for the chordless‑$C_4$ obstruction.
- **`map_sdp_oracle.py`** — a cvxpy/CLARABEL oracle reproducing the two orthogonal axes: the value ladder (LP $0 \le$ pairwise $0 \le$ Shor $\tfrac34 \le$ integer $1$ on the triangle) and the completability obstruction (the $(0.8,0.8,0.8,-0.8)$ $C_4$ spec that is pairwise‑valid but has no global PSD completion).

The numerics quoted here were verified with the oracles. The Julia is written against the CellularSheaves.IPM API but has not been executed in this environment (the `IPM` submodule lives on PR #67, and no Julia toolchain was available); treat the first run as needing a real solve, and reuse the tuned `IPMSettings` from `corridor_benchmark.jl` if the default Uzawa settings converge slowly. Every index, cost, and selection routine in `map_sdp.jl` was cross‑checked by a Python mirror against native SDPs (svec isometry, sub‑block extraction, the $\tfrac34$/$0$ objective split, the completion status).
