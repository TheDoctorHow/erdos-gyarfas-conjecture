# Erdős–Gyárfás Conjecture — Research Brief

**Compiled 2026-07-26.** Statement: every finite graph with minimum degree ≥ 3 contains a cycle
whose length is a power of 2.

---

## 0. The single most important reframe

**The problem is no longer "the conjecture." It is "the bounded-degree case."**

Liu & Montgomery ([arXiv:2010.15802](https://arxiv.org/abs/2010.15802), *JAMS* 2023, Cor. 1.3)
proved: there is an absolute constant `d₀` such that every graph of **average** degree ≥ `d₀`
contains a cycle of length `2^k`. Their machinery (sublinear expanders) gives, for any sequence
`σ` with `σ_{i+1} ≤ exp(σ_i^{1/10})`, a cycle of length in the sequence — powers of 2 have
*enormous* slack relative to that growth condition.

Three consequences that most of the literature does not state explicitly:

1. **The strong Erdős–Gyárfás conjecture is dead.** Erdős and Gyárfás believed the answer was
   *negative*, and conjectured that for every `r` there is a graph of min degree `r` with no
   `2^k` cycle. That is false. So a counterexample can only live at **min degree 3** (or, at most,
   min degree < `d₀`). Anyone hunting the bounty is hunting in a much smaller room than the
   original problem statement suggests.

2. **A counterexample is sparse and bounded-degenerate.** *(Derivation — I did not find this
   stated anywhere, but it is immediate.)* Let `G` be a counterexample. Every subgraph `H ⊆ G`
   also has no `2^k` cycle, so by Liu–Montgomery `H` has average degree `< d₀`. Hence `G` has
   **degeneracy < d₀** and `|E(G)| < d₀·n/2`. Combined with Carr's minimal-counterexample results
   (§2), a minimal counterexample is a *sparse, predominantly cubic, degree-critical* graph.

3. **Sublinear expanders are the wrong tool for what's left.** All the recent heavy machinery
   (Liu–Montgomery, Gao–Huo–Liu–Ma, Letzter's survey [arXiv:2401.10865](https://arxiv.org/abs/2401.10865))
   needs growing average degree. At degree 3 it gives literally nothing. **Making `d₀` explicit is a
   legitimate paper but does not approach 3** — the authors deliberately do not compute it
   ("we do not calculate the function `d₀(ε₁,ε₂)`", §2.1), and the tower of constants in sublinear
   expander arguments is astronomically far from 3.

**So: the open problem is a bounded-degree, essentially cubic, essentially finite-flavoured
problem.** The right toolkit is structural reduction + discharging + exhaustive computation +
construction search — not extremal/probabilistic machinery.

---

## 1. Exact frontier (what is actually proven)

| Result | Who | Where |
|---|---|---|
| min degree ≥ `d₀` (absolute const., non-explicit) | Liu, Montgomery 2020/23 | [2010.15802](https://arxiv.org/abs/2010.15802) |
| 3-connected **cubic planar** | Heckman, Krakovski 2013 | [EJC 20(2)#P7](https://www.combinatorics.org/ojs/index.php/eljc/article/view/v20i2p7) — discharging, partly computer-based |
| planar claw-free | Daniel, Shauger 2001 | — |
| claw-free (cubic: cycle of length `2^k` or `3·2^k`) | Nowbandegani et al. 2014 | — |
| `K_{1,m}`-free with min deg ≥ m+1 or max deg ≥ 2m−1 | Shauger 1998 | — |
| Cayley graphs on generalized quaternion / dihedral / semidihedral groups, groups of order `p³` | 2017 | [Aequat. Math.](https://link.springer.com/article/10.1007/s00010-017-0518-3) |
| `P₈`-free | Gao, Shan 2022 | [2109.01277](https://arxiv.org/pdf/2109.01277) |
| `P₁₀`-free | Hu, Shen 2024 | [2308.05675](https://arxiv.org/pdf/2308.05675) |
| **`P₁₃`-free** (computer-assisted) | 2024 | [2410.22842](https://arxiv.org/abs/2410.22842) |
| **diameter 2** (gives `C₄` or `C₈`) | Carr 2025 | [2508.19302](https://arxiv.org/abs/2508.19302) |
| minimal counterexample: ≥ 4/7 vertices cubic | Carr 2026 | [2605.22844](https://arxiv.org/abs/2605.22844) |
| counterexample ≥ 17 vertices; **cubic ≥ 30**; bipartite ≥ 30 | Markström **2004** | *Extremal graphs for some problems on cycles in graphs* |
| Four 24-vertex cubic graphs whose **only** `2^k` cycle is `C₁₆` | Markström 2004 | one is planar — the "Markström graph" |

Weaker general results: Verstraëte 2005 (avg degree ≥ 10 ⟹ cycles from a set of size `O(n^0.99)`);
Sudakov–Verstraëte 2008 (`|C(G)| ≥ Ω(d^⌊(g−1)/2⌋)`); Gao–Huo–Liu–Ma [IMRN 2022](https://arxiv.org/abs/1904.08126)
(min degree ≥ k+1 ⟹ all even lengths mod k — Thomassen's conjecture); Bai–Grzesik
[2511.03085](https://arxiv.org/abs/2511.03085) (stability, exceptions `K_{k+1}`, `K_{k,n−k}`, Petersen).

---

## 2. What has gone overlooked — ranked openings

### ⭐ 2.1 The computational front is **22 years stale**

Every "counterexample has ≥ N vertices" bound in the literature traces to Markström 2004.
Nothing has been redone since. In that window:

- **Generators got 4–30× faster.** Brinkmann–Goedgebeur–McKay (`snarkhunter`, [dmtcs.episciences.org/551](https://dmtcs.episciences.org/551))
  and Brinkmann–Goedgebeur [*Generation of cubic graphs and snarks with large girth*](https://arxiv.org/abs/1512.05944),
  JGT 2017, with native girth ≥ 4 / ≥ 5 restriction.
- **Hardware got ~10³× cheaper per core-hour.**
- **The search space is far smaller than "all cubic graphs".** A counterexample has no `C₄`.
  For the girth-≥5 slice, connected cubic counts (OEIS A014372) are:

  | n | 24 | 26 | 28 | 30 | 32 |
  |---|---|---|---|---|---|
  | count | 1.6×10⁶ | 3.1×10⁷ | 6.6×10⁸ | **1.46×10¹⁰** | 3.46×10¹¹ |

  `n = 30` is a laptop-weekend today. `n = 32` is a small cluster. (`C₄`-free-with-triangles is a
  strictly smaller companion class, handled by triangle-contraction.)

- **Better still: don't filter, constrain.** `C₄` and `C₈` are *local* — they live inside a
  radius-4 ball, which in a cubic graph has ≤ 46 vertices. So forbidden-cycle constraints can be
  propagated **during** orderly generation rather than tested afterwards. This is exactly the
  design of the `P₁₃`-free paper's backtracking engine, just aimed at a different class. Expect
  to reach `n ≈ 36–40` for cubic.

**The unclaimed prize here is not the bound — it's the catalogue.** Markström's four 24-vertex
near-misses are the *only* published data points about what near-counterexamples look like. A
complete catalogue of cubic graphs whose only `2^k` cycle is `C₁₆` (resp. only `C₃₂`) up to
`n ≈ 36` is the raw material for either (a) an infinite construction, or (b) the reducible-
configuration list for a discharging proof. Nobody has it.

### ⭐ 2.2 Nobody has run EG against the existing censuses

These catalogues already exist, are downloadable, and have apparently never been tested against
the conjecture. This is *days* of compute, not months:

- **Cubic vertex-transitive graphs up to 1280 vertices** (Potočnik–Spiga–Verret census).
- **Cubic arc-transitive** (Foster census, extended by Conder to 10 000 vertices).
- **Snarks up to 36 vertices** (Brinkmann–Goedgebeur–Hägglund–Markström) — cyclically 4-connected,
  girth ≥ 5, exactly the "few short cycles" regime.
- **Cubic graphs of girth ≥ 7, 8, 9**, cages, generalized Petersen graphs.
- **All cubic Cayley graphs on groups of order ≤ 2000** — the 2017 paper did four families by hand;
  a systematic sweep is trivial in GAP/Magma and has not been done.
- [House of Graphs 2.0](https://houseofgraphs.org) as the delivery mechanism.

Either it finds a counterexample (in which case: done), or it produces the single strongest piece
of evidence for the conjecture in 20 years. Both outcomes are publishable. Note the natural
intuition — "a counterexample hides among high-girth cubic graphs" — is *testable* here, and I
suspect it fails, because high-girth cubic graphs have notoriously **full** cycle spectra.

### ⭐ 2.3 The cycle-spectrum-gap literature has never been connected to EG

This is the sharpest quantitative framing available and it sits in a disjoint literature.

Merker, [*Gaps in the cycle spectrum of 3-connected cubic planar graphs*](https://arxiv.org/abs/1905.09101)
(JCTB 2020): for every `k`, every sufficiently large 3-connected cubic planar graph has a cycle
of length in `[k, 2k+9]`. And — crucially — for every even `k ≥ 4` there is an **infinite family**
of 3-connected cubic planar graphs with **no** cycle of length in `[k, 2k+1]`.
See also [2009.00423](https://arxiv.org/abs/2009.00423) (counterexamples to a Merker conjecture)
and [2009.02503](https://arxiv.org/abs/2009.02503) (tight gaps).

Two things follow, and both are worth internalising **before** spending months:

- **The obvious approach is provably dead.** "Show the cycle spectrum contains an interval of
  multiplicative ratio 2" would give EG instantly (any integer interval `[a, 2a]` contains a power
  of 2). Merker's families have multiplicative spectrum gaps of ratio `> 2`. So no such theorem
  exists — not even for cubic *planar* graphs.
- **But it tells you the true difficulty is razor-thin.** The spectrum can dodge `[k, 2k+1]` and
  no more (planar case: `2k+9`). EG asks whether a spectrum that is that dense can nonetheless
  dodge the *entire* geometric sequence `{4, 8, 16, 32, …}` simultaneously. **The non-planar
  analogue of Merker's theorem — is there `f(k) = O(k)` such that every 3-connected cubic graph
  has a cycle length in `[k, f(k)]`? — appears to be open, and is the right measuring stick.**

### 2.4 The induced-path ladder is at its natural end (stop climbing it)

`P₈ → P₁₀ → P₁₃`. Two facts kill the extrapolation:

- **`P₁₂`-free only yields "`C₄` or `C₈`". `P₁₃` is the first level where `C₁₆` is genuinely
  required.** Beyond that every level must reason about `C₁₆`, `C₃₂`, … — a different order of
  combinatorial explosion.
- The authors' own numbers: `k=13` took 11h56m sequential / 17m17s with Cilk; **`k=14` died on
  memory exhaustion**. And they state plainly: *"The technique that we used is not helpful for
  proving the conjecture for H-free graphs when H is not a path."*

**The reusable asset is the algorithm, not the theorem.** Their `explore` backtracking (sound
certificate that a hereditary class contains no counterexample; [code on GitHub](https://arxiv.org/abs/2410.22842))
retargets to other classes with modest work — and **nobody has done this**:

- **diameter 3** (diameter 2 is done — by hand, in a separate paper)
- **treewidth ≤ 3, ≤ 4**; pathwidth ≤ k
- `K_{1,4}`-free; `(C₄, C₈)`-free ∩ bounded something
- bounded genus

Also: the "`C₄` or `C₈`" strengthening that half these papers actually prove is **false in
general** — Markström's 24-vertex graphs refute it. Worth knowing before you try to prove it.

### 2.5 Bounded treewidth / bounded genus: structurally natural, entirely untouched

Heckman–Krakovski did *planar* cubic 3-connected by discharging. There is no result for:
bounded genus (torus!), `K₅`-minor-free, `K_{3,3}`-minor-free, bounded treewidth, bounded pathwidth.

Speculative but concrete: for treewidth ≤ w there is a DP over a tree decomposition tracking
partial-path lengths, so *deciding* EG for a given graph is polynomial. **Question worth an
afternoon:** can the DP state space be quotiented (Myhill–Nerode style) to a finite automaton, so
that a single finite computation settles *all* graphs of treewidth ≤ w with min degree 3?
If yes, that is a genuine theorem-generating machine. I found no attempt at this.

### 2.6 The minimal-counterexample program is one paper old

Carr's [2605.22844](https://arxiv.org/abs/2605.22844) derives its 4/7 bound from a two-line degree
count on top of an observation of Markström's. The obvious continuations are unwritten:

- Girth: a minimal counterexample has no `C₄`, so **girth ∉ {4, 8, 16, …}**, i.e. girth ∈
  {3, 5, 6, 7, 9, …}. Can triangles be excluded outright? (Triangle contraction in a cubic graph
  shifts cycle lengths by 0/1/2 — needs care but is standard.)
- Connectivity: is a minimal counterexample 3-connected? cyclically 4-connected? Each such lemma
  multiplies the reach of §2.1 by orders of magnitude.
- Reducible configurations for 2-cuts and 3-cuts.

This is the Four-Colour-Theorem pipeline — *reducible configurations + discharging + verified
computer check* — and for EG it has barely been started. Combined with §2.1 it is, in my view,
the most plausible route to an actual proof for cubic graphs.

### 2.7 Speculative: the spectral-moment angle

For a cubic graph, `tr(A^ℓ)` decomposes into closed-walk contributions from cycles of length ≤ ℓ
plus tree-like walks. Forbidding `C₄`, `C₈`, `C₁₆` therefore imposes **three exact equations on
the spectral moments** `tr(A⁴)`, `tr(A⁸)`, `tr(A¹⁶)` of a 3-regular graph, where `λ₁ = 3` is
pinned and the spectrum lies in `[−3, 3]`. Whether a positivity/moment obstruction (Alon–Boppana
flavour) can rule this out, I don't know — but it is cheap to test numerically on Markström's
24-vertex near-misses first, and I found no trace of anyone trying it.

---

## 3. Known obstructions — do not spend months rediscovering these

| Approach | Status |
|---|---|
| Expander / large-average-degree machinery | Gives nothing at degree 3. Structurally, not just quantitatively. |
| Make Liu–Montgomery's `d₀` explicit | Legitimate paper; does **not** approach 3. |
| Prove "`C₄` or `C₈`" | **False.** Markström's four 24-vertex cubic graphs. |
| "Spectrum contains an interval `[a, 2a]`" | **False**, even for 3-connected cubic planar (Merker). |
| Find a counterexample of min degree ≥ 4 | **Impossible** for degree ≥ `d₀`; and the strong conjecture is refuted. |
| Push the `P_k`-free ladder by brute force | `k=14` already OOMs; needs `C₁₆` reasoning, not more CPU. |
| Minor-closed / topological-minor arguments | "No `2^k` cycle" is closed under **no** natural graph operation. Suppressing degree-2 vertices changes lengths. This is the structural core of the difficulty. |

---

## 4. Newer methods that genuinely expedite this

**(a) Modern exhaustive generation.** `snarkhunter` / `genreg` / `plantri` / `nauty`, with
in-generation constraint propagation on the local `C₄`/`C₈` conditions. This is the single
highest-value-per-hour item on the list.

**(b) SAT / CP-SAT with symmetry breaking.** The Pythagorean-triples / Schur-number-5 / Keller
toolkit — cube-and-conquer, DRAT/LRAT proof certificates — did not exist in usable form in 2004.
Encoding "cubic on `n` vertices, no cycle of length 4/8/16/32" is natural; the interesting
engineering problem is cycle-length constraints (length-indexed encoding vs. lazy clause
generation with incremental blocking). **No one appears to have attempted a SAT attack on EG.**
The payoff is a machine-checkable certificate for the "no counterexample below `n`" bounds, which
the 2004 results lack.

**(c) ML-guided construction search — the best fit I can see.**
[PatternBoost](https://arxiv.org/abs/2411.00566) (Charton–Ellenberg–Wagner–Williamson) is built
for exactly this shape of problem and has already disproved a 30-year-old conjecture; see also
[AlphaEvolve / *Mathematical exploration and discovery at scale*](https://arxiv.org/abs/2511.02864)
and Wagner's RL constructions. EG has an unusually clean objective:

> minimize `#{cycles of length 2^k}` over cubic graphs on `n` vertices; **zero = counterexample**.

And the landscape is known to have good local minima — Markström's 24-vertex graphs score 1.
The local-search phase (edge swaps / 2-flips preserving cubicity) is trivial to implement, the
verifier is exact and fast, and the objective is dense rather than all-or-nothing. This is
arguably the best-matched open Erdős problem for current AI-for-math tooling, and it is untried.

**(d) Formalisation.** The erdosproblems.com Lean effort plus verified search certificates would
harden the computational results — the `P₁₃`-free proof and Heckman–Krakovski both rest on
unverified code.

---

## 5. Concrete first moves (ranked by value ÷ effort)

1. **Census sweep** (days). Pull the cubic vertex-transitive, arc-transitive, snark, and
   high-girth catalogues; test EG on each. Cheap, decisive either way.
2. **Rebuild Markström's search with modern tooling** (weeks). Target: cubic bound `30 → 36+`,
   general bound `17 → 24+` (using Carr's independence lemma to prune high-degree vertices).
   Emit the **full near-counterexample catalogue**, not just the bound.
3. **PatternBoost / evolutionary run** on the objective above, seeded with Markström's four
   24-vertex graphs (weeks, GPU-light).
4. **Retarget the `P_k` backtracking engine** to diameter 3 and to treewidth ≤ 3 (weeks).
   Two likely-quick papers.
5. **Push the minimal-counterexample structure theory** — girth, connectivity, 2-/3-cut
   reductions — in tandem with (2), since each lemma shrinks the search (ongoing, low cost,
   high leverage).

**Highest-probability path to actually closing it:** (2) + (5) converging into a
reducible-configurations-plus-discharging proof for cubic graphs, extending Heckman–Krakovski off
the plane. **Highest-probability path to the bounty:** (1) + (3).

---

## 6. Open questions I could not find answered anywhere

- Is there `f(k) = O(k)` with every 3-connected **cubic** (non-planar) graph having a cycle of
  length in `[k, f(k)]`? (Merker's theorem, off the plane.)
- Is EG known for **any** class of cubic graphs of girth ≥ 5 other than by exhaustive search?
- Has EG ever been checked on the cubic vertex-transitive census?
- Is a minimal counterexample 3-connected? Triangle-free?
- Does the treewidth-`w` DP admit a finite Myhill–Nerode quotient?
