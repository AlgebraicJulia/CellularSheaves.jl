# Benchmarks

Two scripts produce everything on this page, and every table regenerates verbatim:

```
julia --project=docs docs/scripts/distributed_solve_benchmarks.jl   # solve cost, memory
julia --project=docs docs/scripts/distributed_solve_comparison.jl   # makespans, ordering, hops, wall clock
```

All measurements are on real harmonic-extension systems: sheaf Laplacians
(stalk dimension 2, identity restriction maps) over each formation, targets
pinned, ``H`` the free–free Dirichlet block. Slot counts and memory are exact
and machine-independent. Timings are minimum-of-runs on one machine.

## Exactness

The workspace and distributed solves reproduce the monolithic
`CliqueTrees.Multifrontal` answer to machine precision on every topology and
size tested (max relative residual ``\sim 2\times10^{-16}``), and the real
four-process run below agrees to ``9.5\times10^{-17}``. Nothing on this page
trades accuracy for anything.

## Communication makespan: the full tables

Slots to complete one solve (or reach ``\varepsilon = 10^{-6}``, for the
iterative methods), under the model of the [comparison page](comparison.md).
**tree (ND)** is the distributed tree under a nested-dissection ordering.
The **fewest-slots** cell in each row is highlighted.

#### grid

```@raw html
<table class="bench">
<thead><tr><th>agents</th><th>tree</th><th>tree (ND)</th><th>central (base)</th><th>central (in-net)</th><th>Richardson</th><th>Chebyshev</th><th>CG (+all-reduce)</th><th>&kappa;(H)</th></tr></thead>
<tbody>
<tr><td>59</td><td>26</td><td class="win">16</td><td>118</td><td>116</td><td>2060</td><td>252</td><td>2520</td><td>74.4</td></tr>
<tr><td>139</td><td>56</td><td class="win">22</td><td>278</td><td>276</td><td>5980</td><td>428</td><td>5992</td><td>216.4</td></tr>
<tr><td>251</td><td>92</td><td class="win">26</td><td>502</td><td>504</td><td>12340</td><td>616</td><td>11088</td><td>446.6</td></tr>
<tr><td>395</td><td>146</td><td class="win">30</td><td>790</td><td>792</td><td>21360</td><td>808</td><td>17776</td><td>772.9</td></tr>
<tr><td>571</td><td>184</td><td class="win">34</td><td>1142</td><td>1148</td><td>33196</td><td>1008</td><td>26208</td><td>1201.4</td></tr>
</tbody></table>
```

#### chain

```@raw html
<table class="bench">
<thead><tr><th>agents</th><th>tree</th><th>tree (ND)</th><th>central (base)</th><th>central (in-net)</th><th>Richardson</th><th>Chebyshev</th><th>CG (+all-reduce)</th><th>&kappa;(H)</th></tr></thead>
<tbody>
<tr><td>14</td><td>22</td><td class="win">8</td><td>28</td><td>26</td><td>1252</td><td>138</td><td>2070</td><td>90.5</td></tr>
<tr><td>30</td><td>54</td><td class="win">12</td><td>60</td><td>58</td><td>5372</td><td>286</td><td>8866</td><td>388.8</td></tr>
<tr><td>62</td><td>118</td><td class="win">16</td><td>124</td><td>122</td><td>22216</td><td>582</td><td>36666</td><td>1607.9</td></tr>
<tr><td>126</td><td>246</td><td class="win">20</td><td>252</td><td>250</td><td>90302</td><td>1174</td><td>149098</td><td>6536.2</td></tr>
<tr><td>254</td><td>502</td><td class="win">24</td><td>508</td><td>506</td><td>364080</td><td>2356</td><td>600780</td><td>26353.0</td></tr>
<tr><td>510</td><td>1014</td><td class="win">28</td><td>1020</td><td>1018</td><td>1462064</td><td>4720</td><td>2411920</td><td>105827.7</td></tr>
</tbody></table>
```

#### ring

```@raw html
<table class="bench">
<thead><tr><th>agents</th><th>tree</th><th>tree (ND)</th><th>central (base)</th><th>central (in-net)</th><th>Richardson</th><th>Chebyshev</th><th>CG (+all-reduce)</th><th>&kappa;(H)</th></tr></thead>
<tbody>
<tr><td>15</td><td>24</td><td class="win">8</td><td>30</td><td>40</td><td>1426</td><td>148</td><td>2516</td><td>103.1</td></tr>
<tr><td>31</td><td>56</td><td class="win">12</td><td>62</td><td>88</td><td>5726</td><td>296</td><td>9768</td><td>414.3</td></tr>
<tr><td>63</td><td>120</td><td class="win">16</td><td>126</td><td>184</td><td>22926</td><td>592</td><td>38480</td><td>1659.4</td></tr>
<tr><td>127</td><td>248</td><td class="win">20</td><td>254</td><td>376</td><td>91730</td><td>1184</td><td>152736</td><td>6639.5</td></tr>
<tr><td>255</td><td>504</td><td class="win">24</td><td>510</td><td>760</td><td>366942</td><td>2366</td><td>608062</td><td>26560.1</td></tr>
<tr><td>511</td><td>1016</td><td class="win">28</td><td>1022</td><td>1528</td><td>1467792</td><td>4730</td><td>2426490</td><td>106242.3</td></tr>
</tbody></table>
```

#### star

The degenerate case: one hub sits in every separator, so the tree collapses to a
single chunk and ties the coordinator, and nothing beats either.

```@raw html
<table class="bench">
<thead><tr><th>agents</th><th>tree</th><th>tree (ND)</th><th>central (base)</th><th>central (in-net)</th><th>Richardson</th><th>Chebyshev</th><th>CG (+all-reduce)</th><th>&kappa;(H)</th></tr></thead>
<tbody>
<tr><td>15</td><td class="win">26</td><td class="win">26</td><td>30</td><td>28</td><td>24570</td><td>1624</td><td>8120</td><td>254.0</td></tr>
<tr><td>63</td><td class="win">122</td><td class="win">122</td><td>126</td><td>124</td><td>1753422</td><td>28830</td><td>144150</td><td>4094.0</td></tr>
<tr><td>255</td><td class="win">506</td><td class="win">506</td><td>510</td><td>508</td><td>114984022</td><td>471932</td><td>2359660</td><td>65534.0</td></tr>
</tbody></table>
```

#### random geometric

No nested-dissection construction is applied here, so the two tree columns
coincide and share the win.

```@raw html
<table class="bench">
<thead><tr><th>agents</th><th>tree</th><th>tree (ND)</th><th>central (base)</th><th>central (in-net)</th><th>Richardson</th><th>Chebyshev</th><th>CG (+all-reduce)</th><th>&kappa;(H)</th></tr></thead>
<tbody>
<tr><td>30</td><td class="win">8</td><td class="win">8</td><td>60</td><td>58</td><td>24541</td><td>2070</td><td>8910</td><td>154.4</td></tr>
<tr><td>62</td><td class="win">22</td><td class="win">22</td><td>124</td><td>122</td><td>144396</td><td>6300</td><td>28000</td><td>580.6</td></tr>
<tr><td>126</td><td class="win">22</td><td class="win">22</td><td>252</td><td>250</td><td>144320</td><td>6952</td><td>34128</td><td>474.8</td></tr>
<tr><td>254</td><td class="win">52</td><td class="win">52</td><td>508</td><td>506</td><td>415428</td><td>12844</td><td>43472</td><td>1156.4</td></tr>
</tbody></table>
```

![Four log-log panels (grid, random geometric, ring, chain) showing makespan versus fleet size for all six methods, with the distributed tree at the best ordering sitting lowest in every panel](../../assets/figures/distributed_solve/comparison.svg)

Two readings worth pulling out of the tables:

- On every formation with sublinear separators, tree-ND grows like the plots'
  flattest curve (``O(\log n)``) while both centralized variants grow linearly
  and the iterative methods grow with ``\kappa(H)``.
- ``\kappa(H)`` itself grows with formation size (it is over ``10^5`` on the
  510-agent chain), so diffusion's disadvantage *compounds* with scale. Its
  locality is not the bottleneck. Its conditioning bill is.

## The ordering study

The elimination *ordering* sets the tree's depth, and the depth sets the
makespan. Nested dissection here is three lines of Julia (recursively eliminate
interval/box middles last), with no external graph partitioner:

Smaller of each metric pair highlighted: nested dissection wins depth and slots,
the default order wins fill, which is the trade the study is about.

```@raw html
<table class="bench">
<thead><tr><th>topology</th><th>agents</th><th>depth</th><th>depth (ND)</th><th>slots</th><th>slots (ND)</th><th>fill</th><th>fill (ND)</th></tr></thead>
<tbody>
<tr><td>chain</td><td>62</td><td>60</td><td class="win">6</td><td>118</td><td class="win">16</td><td class="win">240</td><td>412</td></tr>
<tr><td>chain</td><td>254</td><td>252</td><td class="win">8</td><td>502</td><td class="win">24</td><td class="win">1008</td><td>1916</td></tr>
<tr><td>chain</td><td>510</td><td>508</td><td class="win">9</td><td>1014</td><td class="win">28</td><td class="win">2032</td><td>3948</td></tr>
<tr><td>ring</td><td>511</td><td>509</td><td class="win">10</td><td>1016</td><td class="win">28</td><td class="win">2036</td><td>3952</td></tr>
<tr><td>grid</td><td>251</td><td>45</td><td class="win">10</td><td>92</td><td class="win">26</td><td class="win">6516</td><td>8076</td></tr>
<tr><td>grid</td><td>571</td><td>93</td><td class="win">14</td><td>184</td><td class="win">34</td><td class="win">19816</td><td>24708</td></tr>
</tbody></table>
```

![Left: elimination-tree depth versus chain length for the default order and nested dissection, against n and log n guide lines. Right: the resulting makespans, with the centralized cost for scale](../../assets/figures/distributed_solve/ordering.svg)

The trade is explicit: on the 510-agent chain, nested dissection buys a
**36× makespan reduction** (1014 → 28 slots) for a **1.9× fill increase**
(2032 → 3948 stored entries). Depth tracks the ``\log_2 n`` guide line almost
exactly. For a deployment this is nearly always the right trade, since
communication slots are wall-clock and battery, and a doubled factor is
kilobytes.

## The physical routing overlay

Measured hop distributions for tree messages over the physical radio graph
(discussion and caveats on the [comparison page](comparison.md)):

| formation, ordering | messages | mean hops | 1-hop share | total transmissions |
|---|---:|---:|---:|---:|
| 20×20 grid, default order | 355 | 12.3 | 3% | 4368 |
| 20×20 grid, nested dissection | 266 | 9.8 | 18% | 2605 |
| random geometric, ``n = 254`` | 71 | 3.1 | 7% | 224 |

![Hop-count histograms for the 20×20 grid and the 256-agent random geometric graph](../../assets/figures/distributed_solve/hops.svg)

## Solve cost on one machine

![Left: per-solve time versus problem size for the CliqueTrees monolithic solve, the recursive tree solve, and the workspace tree solve. Right: one-time factorization cost versus the per-step re-solve cost](../../assets/figures/distributed_solve/solve_scaling.svg)

The preallocated [`TreeWorkspace`](@ref) solve tracks the monolithic
`CliqueTrees` solve within a few percent (the allocating recursive reference is
``\approx 5\times`` slower), and the per-step re-solve is roughly an order of
magnitude cheaper than the factorization it reuses. Every "per-step" number in
this guide is a re-solve against a cached factor. The factorization is paid
once per formation and reported separately.

## A real multi-process run

| solve | wall time | vs single process | agrees to |
|---|---:|---:|---|
| single process, workspace | 0.053 ms | 1× | n/a |
| 4 worker processes, `RemoteChannel` | 3.078 ms | 58× slower | ``9.5\times10^{-17}`` |

Local worker processes share one machine, so this measures serialization and
channel protocol overhead and verifies exactness. It is not a radio
simulation (the slot tables are the communication model). The number worth
keeping is that this is real inter-process message passing, exact to the last
bit.

## Per-agent memory

![Left: factor-slice size held by each worker for a 24×24 grid, with the mean marked. Right: centralized single-node storage versus worst-case and mean per-worker storage across sizes](../../assets/figures/distributed_solve/memory.svg)

The per-worker slices are balanced around their mean and sum *exactly* to the
centralized factor, so there is zero duplication, by the disjoint-residuals
argument of the [multifrontal page](theory_multifrontal.md). On a 24×24 grid
across twelve workers, the busiest worker holds **12%** of the full factor.
