# Reference verification ledger

Each reference in `references.bib` is listed with:
- what the paper attributes to it;
- which text was read, and in which version;
- the verdict.

Theorem numbers come from the version named. Before submission, check those
marked *arXiv numbering* against the published version.

Verdicts:
- **verified**: the claim was read in the primary text;
- **secondary**: the claim is confirmed only through a source that quotes the
  work;
- **bibliographic**: only the citation data was checked (the paper cites the work
  without attributing a specific result).

## Category theory and semantics

| Key | Claim in the paper | Source read | Verdict |
|---|---|---|---|
| `lawvere1973` | `([0,∞], ≥, +, 0)` is a symmetric monoidal closed category with truncated subtraction as hom; generalized metric spaces are categories enriched in it; no symmetry, finiteness or separation is needed; functors are distance-decreasing maps. | TAC Reprints No. 1 (2002): §1, orig. pp. 136, 138, 145, 147–149 (reprint pp. 6, 8, 15, 17–19) and the 2002 commentary | verified |
| `benabou2000` | For an ordinary functor `P : X → I`, a normalized lax functor `I^op → Dist` with fibres as values on objects and morphisms over `α` as values on morphisms. | Streicher's notes, §7, pp. 24–25 (stated informally, without theorem number or proof) | verified (Set version only; the paper claims only that) |
| `street1983` | Bénabou's observation that a category enriched in a bicategory `W`, with object set `X`, is a lax functor from the chaotic category on `X` into `W`. | TAC Reprints No. 14 (2005), §2, p. 6, and commentary p. 1 | verified |
| `joyal1996` | Zig-zag morphisms (lifting each transition from every reachable state) are the open maps for branches; bisimilarity is spans of open maps. | BRICS report RS-94-7 (1994): §2, pp. 11–14, Proposition 1, Theorem 2. The definition of P-open is unnumbered. The I&C 1996 text was not accessed, and its numbering may differ, so the paper cites no numbers. | verified |
| `tarski1955` | Fixed points of a monotone map on a complete lattice form a complete lattice. | Theorem 1, p. 286 | verified |
| `paige1987` | Stability: `B ⊆ E⁻¹(S)` or `B ∩ E⁻¹(S) = ∅`; the relational coarsest partition is the coarsest stable refinement. | pp. 977–979 (§3; uniqueness in Theorem 2) | verified |

The claim that an *enriched* Bénabou correspondence exists was **not found** in
Stubbe (TAC 14, 2005) or Street, "Fibrations in bicategories". Betti, Carboni,
Street and Walters (1983) could not be accessed. The paper therefore proves its
cost-enriched statement directly and cites Bénabou only for the Set prototype.

## Sheaves and lattices

| Key | Claim in the paper | Source read | Verdict |
|---|---|---|---|
| `ghrist2022` | Restriction maps are Galois connections; sections are computed with the lower adjoints; sheaves of lattices. | arXiv 2007.04099v2: §2.4–2.5, eq. (4)–(5), Definition 3.1, Theorem 3.1 (*arXiv numbering*) | verified |
| `ghrist2025diffusion` | Lawvere Laplacian; weighting `W : V×V → Q` with `W`-fuzzy sections in both directions (Def. 3.2); suffix points ≅ fuzzy sections for a crisp Laplacian (Thm. 3.17); shortest-path example in §4.2 with a constant sheaf on `R^op`, symmetric weights, Laplacian `min_w W(v,w) + x_w`, Dijkstra sketched via flow weights. | arXiv 2501.03890v3 LaTeX source (preprint only; numbering counted from v3) | verified |
| `hansen2019` | Kron reduction of sheaf Laplacians is not possible in general; positive result for vertex stalks of dimension at most one (§4.2). | arXiv 1808.01513v2, §4.2 and Theorem 4.3 | verified |
| `curry2014` | A cellular sheaf is a functor from the face poset to a data category. | arXiv 1303.3255v2, Definition 4.1.6 | verified |
| `desilva2016` | Reeb graphs are equivalent to constructible set-valued cosheaves on ℝ. | arXiv 1501.04147, Theorem 3.22 (*arXiv numbering*) | verified |
| `reeb1946` | Reeb graphs (cited for the construction). | zbMATH Zbl 0063.06453 | bibliographic |
| `dorfler2013` | Kron reduction is the Schur complement of a (loopy) Laplacian onto boundary nodes, and the result is again a Laplacian. | arXiv 1102.2950, §2.1 and Lemma 2.1 (*arXiv numbering*) | verified |

## Planning, control and multi-agent path finding

| Key | Claim in the paper | Source read | Verdict |
|---|---|---|---|
| `holte1996` | Distances in an abstracted space are admissible heuristics; Hierarchical A*. | AAAI-96 PDF, Introduction and §"Hierarchical Search using A*" | verified |
| `culberson1998` | Pattern databases. | Crossref | bibliographic |
| `helmert2007` | Merge-and-shrink abstraction heuristics. The paper does not use that name; it describes merging and shrinking. | ICAPS 2007 PDF | verified |
| `larsen1988` | May/must transitions of modal transition systems. | Primary text not accessed. Confirmed by Beneš, Křetínský et al., arXiv 1304.5278: "Modal transition systems (MTS) were introduced by Larsen and Thomsen [LT88]". | secondary |
| `tabuada2009` | Symbolic control by abstraction and refinement. | Crossref | bibliographic |
| `reissig2017` | Feedback refinement relations; abstract controllers refine to concrete ones. | arXiv 1503.03715, Def. V.2 and Thm. VI.3 | verified |
| `matni2024` | Layered control architectures with planning and tracking layers. | arXiv 2401.15185 (published title differs from the arXiv title; the bib uses the published one) | verified |
| `delling2011` | Overlay of shortcut cliques between the boundary vertices of each cell (the default overlay uses reduced cliques). | SEA 2011 PDF, §3 | verified |
| `li2019` | WDG heuristic: edge-weighted minimum vertex cover with `x_i + x_j ≥ Δ_ij`, pair costs solved exactly. The paper does **not** attribute CG or any matching to this work. | IJCAI-19 PDF, §§4–5 | verified |
| `sharon2015` | Conflict-based search. | Crossref | bibliographic |
| `standley2010` | Independence detection. | AAAI-10 PDF, abstract and §"Independent Subproblems" | verified |

## Path algebra and dynamic programming

Pending: being verified.

## Software

`PathPack.jl` (AlgebraicJulia, unregistered) is mentioned in the repository
documentation but not cited in the paper. Confirm its public availability before
citing it.
