# Counterexample Hunt: Results Log

**Started 2026-07-26.** Target: a finite graph with minimum degree ≥ 3 and no cycle of length `2^k`.
By Liu–Montgomery such a graph must have minimum degree exactly 3 and bounded degeneracy, so the
hunt is over **cubic** graphs.

---

## 1. Infrastructure built (all in `tools/`)

| File | What it is |
|---|---|
| `egcore.h` | Exact cycle machinery: `eg_cycles_of_length` (count/detect cycles of a *given* length via min-vertex-rooted DFS with BFS distance pruning and a node budget), `eg_pow2_mask`, graph6 I/O. |
| `egtool.c` | CLI: graph6 on stdin → power-of-two cycle mask, exact counts, full spectrum, `COUNTEREXAMPLE`/`NEARMISS` flags. |
| `egfilter.c` | snarkhunter **plugin**: in-process filtering of every generated cubic graph, zero I/O. |
| `holefilter.c` | snarkhunter plugin computing full cycle spectra → hole statistics. |
| `inflate1.c` | The triangle reduction (§4). Reads cubic graphs on `n−2` vertices, emits/tests all C4-free cubic graphs on `n` vertices that contain a triangle. |
| `ref.py` | **Independent** reference implementation using a *different algorithm* (subset DP / Held–Karp) for cross-validation, plus random cubic graph generation and graph6 codecs. |
| `sweep.sh`, `campaign.sh`, `run.sh` | Sharded parallel drivers. |

External: **nauty 2.8.9** and **snarkhunter 2.0b** (Brinkmann–Goedgebeur–McKay) downloaded from the
authors' sites and built from source under MinGW (needed a `<sys/times.h>` shim).

> Build note: MinGW gcc cannot handle the non-ASCII characters in this folder's path. Build and run
> from an ASCII path. `gcc -O3 -march=native -mpopcnt`.

### Validation (this is the part that matters)

1. **Petersen graph** → spectrum `{5,6,8,9}`, exactly 15 eight-cycles. Both correct.
2. **120 random cubic graphs** on 12 and 14 vertices: the C DFS engine and the independent Python
   subset-DP produced **identical full cycle spectra**, with 0 mismatches.
3. **snarkhunter** reproduced the known count 90,938 girth-≥5 cubic graphs on 22 vertices.
4. End-to-end: the sweep of **all 125,816,453 connected cubic graphs on ≤ 24 vertices** returned
   exactly that count (the known total), in 112 s wall on 16 cores.

---

## 2. Markström's 2004 result independently reproduced

Among all 125,816,453 connected cubic graphs on ≤ 24 vertices there are **exactly four** with no
`C4` and no `C8`, and all four have a `C16`. This is precisely Markström's published finding, and
to my knowledge the four graphs themselves have not been published; they are in
`tools/markstrom24.g6`. Full anatomy:

| triangles | #C16 | cycle spectrum | contracts to |
|---|---|---|---|
| 3 | 330 | 3,5,6,7,9–24 | cubic graph on 18 vertices |
| 4 | 315 | 3,6,7,9–24 | cubic graph on 16 vertices |
| 6 | 207 | 3,5,6,7,9–24 | cubic graph on 12 vertices |
| 7 | 228 | 3,5,6,7,9–24 | cubic graph on 10 vertices |

All four: 3-edge-connected, cyclic edge connectivity 3, **girth 3**.

**Three things this tells you that the literature does not.**

- **The near-misses are triangle-rich.** Every one has 3–7 triangles. My first exhaustive pass
  restricted to girth ≥ 5 found *nothing*: the interesting graphs all have girth 3. Anyone
  searching this problem in the girth-≥5 class is searching the wrong place.
- **They are not "near" at all in the C16 direction.** 207–330 sixteen-cycles is not a near miss;
  it is a rout. The spectra are essentially the full interval `[9,24]` plus `{3,5,6,7}`; only two
  or three lengths are missing in total, and a counterexample needs `{4,8,16}` *all* missing.
- **This is exactly why the `P_k`-free ladder stalls.** `P₁₂`-free only yields "C4 or C8"; `P₁₃` is
  the first level where `C16` must genuinely be handled, and the published search ran out of memory
  at `P₁₄`. The obstruction is the same one visible in this table.

---

## 3. New data: how sparse can a cubic cycle spectrum be?

A counterexample needs holes in its cycle spectrum at 4, 8, 16, 32, …, so the real question is how
many holes a cubic graph can have. Exhaustive over **all 556,471 connected cubic graphs on ≤ 20
vertices** (full spectrum computed for each):

| holes | 0 | 1 | 2 | 3 | 4 | 5 | … | 13 | 14 | 15 |
|---|---|---|---|---|---|---|---|---|---|---|
| count | 384,611 | 145,550 | 10,347 | 292 | 241 | 7,699 | … | 64 | 5 | 3 |

Maximum holes by order: n=10 → 5, n=14 → 9, n=16 → 11, n=18 → 13, **n=20 → 15** (of 18 possible
lengths). So spectra *can* be extremely sparse, but the extremal graphs all have spectrum exactly
**{3, 4, 5}**: chains of diamonds (K₄ minus an edge). Sparsity is bought by making every cycle
short and local, and a diamond contains a `C4`. That is the tension a counterexample has to break:
**you need sparsity at every scale simultaneously, and the cheap source of sparsity forces a C4.**

---

## 4. New tool: a complete triangle reduction (≈13–15× speedup)

**Theorem.** Let `G` be a connected cubic C4-free graph containing a triangle `T = {x,y,z}`. The
three outside neighbours of `x,y,z` are distinct (otherwise `G` has a `C4`), so contracting `T` to
a vertex `v` gives a connected cubic simple graph `H` on `n−2` vertices with `G = inflate(H,v)`.
Moreover

> `inflate(H,v)` is C4-free **⟺** (a) no triangle of `H` passes through `v`, and (b) every `C4` of
> `H` passes through `v`.

*Proof sketch.* A cycle meets `T` in 0 or ≥2 vertices (a cycle through `xᵢ` uses two of its three
edges, at least one inside `T`). C4s avoiding `T` are exactly C4s of `H` avoiding `v` → (b). A C4
meeting `T` in two vertices is `xᵢ xⱼ b a xᵢ` where `a,b` are the unique distinct outside
neighbours and `a~b`; in `H` that is the triangle `v a b` → (a). A C4 meeting `T` in three vertices
would need an outside vertex adjacent to two triangle vertices, which is impossible. ∎

**Consequence.** Enumerating cubic C4-free graphs on `n` vertices splits into
`{girth ≥ 5 on n}` ∪ `{inflate(H,v) : H cubic on n−2}`. Since cubic graph counts grow ~20× per two
vertices, this replaces the dominant term by a much smaller one.

| n | brute force (all cubic on n) | reduction (girth-5 on n + all cubic on n−2) | speedup |
|---|---|---|---|
| 24 | 1.18×10⁸ | 1.6×10⁶ + 7.3×10⁶ | ~13× |
| 28 | 4.05×10¹⁰ | 6.6×10⁸ + 2.09×10⁹ | ~15× |
| 30 | 8.45×10¹¹ | 1.46×10¹⁰ + 4.05×10¹⁰ | **~15×** (≈9 days → ≈14 h on 16 cores) |

**A first version of this was wrong and validation caught it.** I initially contracted *all*
triangles at once and claimed the base is triangle-free. It is not: a base triangle whose vertices
are 2 or 3 contracted vertices creates no `C4`, so such bases are legal. The flawed version missed
the 7-triangle Markström graph (whose base on 10 vertices has 3 triangles). Contracting **one**
triangle at a time is complete.

**Validation of the correct version at n=24:** it recovers exactly the four near-miss classes, with
triangle counts 3/4/6/7 and C16 counts 330/315/207/228, identical to brute force, from 13× fewer
base graphs. (Each class appears 3–4 times, once per contractible triangle, as expected.)

---

## 4b. The best lead: bridges cap the circumference

A cubic graph's bridges lie on no cycle, so its cycle spectrum is the union of the spectra of its
blocks. A **leaf** block of the bridge tree of a cubic graph is exactly a *1-terminal gadget*: one
vertex of degree 2, all others of degree 3. Conversely, joining the terminals of two gadgets by an
edge gives a cubic graph whose spectrum is the union of the two. Hence:

> **A cubic counterexample containing a bridge exists ⟺ a 1-terminal gadget with no `2^k` cycle
> exists.** A gadget of order `k` has circumference ≤ `k`, so for `k ≤ 31` every power of two above
> 16 is *structurally absent*: it only has to dodge `{4, 8, 16}`. Doubling it gives a cubic
> counterexample on `2k ≤ 62` vertices.

This is a much better search target than "cubic graphs on 30 vertices": same object size, **double
the reach**.

**Completeness of the enumeration.** Let `t` be the terminal with neighbours `p,q`.
*If `p ≁ q`*: `M := H − t + pq` is cubic on `k−1` vertices and `H = subdivide(M, pq)` (family 1).
*If `p ~ q`*: `t,p,q` is a triangle; deleting it and joining the third neighbours `x,y` gives a
cubic `M` on `k−3` vertices (family 2). So every gadget comes from a cubic graph by one of two
operations. Pruning: subdividing/replacing edge `uv` only destroys the C4s of `M` that *use* `uv`,
so the gadget is C4-free only if **every C4 of `M` contains `uv`**; usually there are 0 candidate edges.
Verified to cut work 8× with zero loss.

**Result so far.** Every C4-free 1-terminal gadget of order ≤ 21 has a `C8`. Exactly **two** exist
at order 23, and both have a `C16` (157 and 111 of them). Doubling one gives a genuinely new
near-miss:

> **A cubic graph on 46 vertices, min degree 3, with no `C4`, no `C8`, and no `C32`** (circumference
> 23 < 32, so `C32` cannot exist): 314 sixteen-cycles are the *only* obstruction.
> Verified three independent ways. This is far beyond anything in the literature, where the known
> near-misses stop at 24 vertices.

Orders `k = 25, 27, 29, 31` are running. A single gadget there with no `C16` **is the bounty.**

## 4c. A route that is provably dead: do not spend time on it

I had flagged "min degree 3 + no cycle ≡ 0 mod 4 ⟹ counterexample" as a promising opening in
`RESEARCH-BRIEF.md` §2.7/Modality D. **That set is empty, by a 1993 theorem.**

> Dean, Lesniak, Saito, *Cycles of length 0 modulo 4 in graphs*, Discrete Math. **121** (1993)
> 37–49: every graph with `δ(G) ≥ 2` and at most two vertices of degree 2 has a cycle of length
> ≡ 0 (mod 4).

So every graph of min degree ≥ 3 has a `0 mod 4` cycle, and the whole modular-obstruction family
(mod 4, and a fortiori the framing I gave) collapses. Gao–Huo–Liu–Ma is about the much stronger
"all even residues mod k"; the single residue 0 mod 4 at min degree 3 was settled 29 years earlier.
Exhaustive check over all 7,875,918 connected cubic graphs on ≤ 22 vertices confirms it: zero
mod-4-free graphs. The correction is worth more than the idea was.

## 5. Search status

| Search | Status | Result |
|---|---|---|
| All cubic graphs, n ≤ 24 (exhaustive) | ✅ done | 4 near-misses, **0 counterexamples** |
| All cubic graphs n ≤ 20, hole statistics | ✅ done | max 15 holes; extremal spectrum {3,4,5} |
| All cubic graphs n ≤ 22, mod-4 sweep | ✅ done | 0 mod-4-free (consistent with Dean–Lesniak–Saito) |
| n = 26 (partial brute force) | partial | ≥1 new near-miss at n=26 (mask = 16) |
| Non-cubic min-deg-3, n = 17, 18 | ✅ done | nothing even reaches C4-free ∧ C8-free |
| Cubic arc-transitive ≤ 10000 vtcs (3,815) | ✅ done | **0 counterexamples**; 689 dodge {4,8,16} |
| Cubic vertex-transitive ≤ 1280 (59,038 of 111,360) | ✅ done | **0 counterexamples** |
| Annealing at n = 30, 32, 34, 40, 48, 62, 64 | ✅ done | floor `#C16 = 210` at n=30, reproducible |
| 1-terminal gadgets, order ≤ 23 | ✅ done | all have C8 except two at k=23, both with C16 |
| **1-terminal gadgets, order 25** | ✅/🔄 | family 2 done: 9 with no C4/C8, **all have C16** |
| **1-terminal gadgets, order 27, 29, 31** | 🔄 running | a single C16-free one **is the bounty** |
| n = 28, 30 via triangle reduction | ⏸ paused | deprioritised; gadget route reaches n ≤ 62 for the same cost |

At `n = 30` only `{4, 8, 16}` must be dodged (`32 > 30`), which makes it the single most favourable
open order. For `32 ≤ n < 64` the graph must also dodge `C32`.

---

## 5b. Parallel search campaign: six modalities

### A. Census sweep: **62,853 symmetric cubic graphs, never before tested**

- **All 3,815** cubic **arc-transitive** graphs on ≤ 10,000 vertices (Conder's complete census).
- **59,038** of the 111,360 cubic **vertex-transitive** graphs on ≤ 1,280 vertices
  (Potočnik–Spiga–Verret census).
- **Zero counterexamples, zero undecided.**

The structural finding is the interesting part, and it is a regime change:

| class | dodge 4, 8 and 16? |
|---|---|
| random / small cubic | smallest `2^k` is almost always 4 |
| the four n≤24 near-misses | dodge 4 and 8 only |
| **cubic arc-transitive** | **689 of 3,815 (18.1%) dodge 4, 8 AND 16** |

Symmetry suppresses small powers of two very effectively, and then **stops dead at 32**. Every
one of the 689 triple-dodgers has a `C32`. Best examples: an arc-transitive graph on **2030**
vertices of girth only 7 that misses 4, 8 and 16; and a vertex-transitive graph on **180** vertices
of **girth 3** that misses 4, 8 and 16.

### E. Simulated annealing / PatternBoost-style search: a sharp, reproducible floor

Two independently written annealers (one full-recompute, one with incremental `#C16` counting,
cross-validated against each other and against `egtool`), over cubicity-preserving 2-swaps.

**At n = 30 the objective floor is exactly `(#C4, #C8, #C16) = (0, 0, 210)`.** Hit by four
independent runs: two different annealers, four RNG seeds, two weight schedules, two temperature
schedules, both random and Markström-seeded starts. **Nothing ever got below 210.** The optimum's
spectrum is `{3,5,6,7,9,10,…,29}`, a complete interval. `#C16` is the binding constraint; `C4`,
`C8`, even `C32` are easy to kill.

### F. Non-cubic min-degree-3, n = 17–29: the genuinely unexplored region

A new exhaustive backtracker (canonical BFS order, pruning the instant a `C4` or `C8` appears,
enforcing the minimal-counterexample reductions) **rigorously cleared n = 17 (Δ ∈ {4,5,6}) and
n = 18 (Δ ∈ {4,5})**. Not one graph in that region even reaches the "C4-free and C8-free" stage, so
`C16` never becomes relevant. Adding degree-≥4 vertices makes things *strictly worse*: it forces
more edges into a graph that must stay C4-free. The search also **independently re-derived all four
Markström graphs** from scratch: a third fully independent confirmation.

### D. Modular obstructions: see §4c (route provably dead; produced the bridge reduction instead)

### B, C. Algebraic families and engineered gadget constructions: still running.

---

## 5c. The sharpest formulation of what is left

Putting the census data together with the circumference observation gives what I think is the
cleanest statement of the remaining problem:

> A counterexample of circumference `c` must dodge only the `⌊log₂ c⌋ − 1` powers of two below `c`.
> Triple-dodging `{4, 8, 16}` is **not** rare: 689 symmetric graphs do it. But every known
> triple-dodger has large circumference, and then `C32` appears.
>
> **So the whole question is: can `{4, 8, 16}` be dodged by a graph of circumference ≤ 31?**

A min-degree-3 graph of bounded circumference cannot be 3-connected and large (3-connected cubic
graphs have circumference `Ω(n^0.694)`), so it must be a **tree of small blocks**, precisely the
bridge/2-cut regime of §4b. That is why the gadget search is the right search, and it is why the
two search directions that looked unrelated at the start (small-`n` exhaustion and large-`n`
symmetric graphs) are actually the same question approached from opposite sides.

---

## 5d. Corrections: things I asserted that turned out to be wrong

An adversarial audit pass re-verified every headline claim by independent code. Every *number* held
up; two *interpretations* did not, and one of them was mine.

**(1) "Cubic cycle spectra are near-complete intervals." FALSE at scale; retract it.**
I said this in §2 and §6 below, and built the search strategy around it. It is an artifact of small
`n` and low symmetry. The campaign's best object:

```
C2030.1   (cubic, arc-transitive, girth 7, n = 2030)
   L  =   3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20
          .  .  .  .  X  .  .  .  .  X  .  .  X  .  X  X  X  X
   spectrum begins  7, 12, 15, 17, 18, 19, 20, …
```
**Seven of the first eleven lengths above the girth are missing.** Sparse spectra are entirely real
at scale; they just don't happen at n ≤ 30.

> ### ⚠ RETRACTION (added after adversarial refereeing)
> **Item (2) below is FALSE and must not be used.** The "dense-onset law `D ≈ 2.4·girth`" is a
> small-`n` artifact. Merker's families are 3-connected cubic **planar**, so by Euler they have
> girth ≤ 5, and they have no cycle length in `[k, 2k+1]` for arbitrarily large `k`; since the
> circumference lies in the spectrum and exceeds `2k+1`, this forces `D ≥ 2k+2`. Hence
> **`D/girth` is unbounded**, and the derived heuristic ("dodging three powers of two needs
> `D/girth > 4`, never observed") has no basis. The *conclusion* (that few consecutive powers of
> two can be dodged) may well still be true; the *derivation* is invalid. Item (2) is kept below
> only so the error is on the record.
>
> Also retracted: the framing in §5c below in terms of spectrum **gaps**. An interval-*meeting*
> theorem cannot prove Erdős–Gyárfás at all; see `PROOF-PROGRAM.md` §4.

**(2) [RETRACTED, see above] The right quantity is the dense-onset ratio, and it is measurable.**
Define `D(G)` = least `L` with `spec(G) ⊇ [L, circumference]`. Empirically **`D ≈ 2.4 · girth`**
(C2030.1: g=7, D=17). The number of consecutive powers of two a graph can dodge is about
`log₂(D/girth) ≈ log₂(2.4)`: **two, and only two.** To dodge **three** you need `D/girth > 4`.

> Across ~10 million exhaustively enumerated cubic graphs, all 111,360 cubic vertex-transitive
> graphs on ≤ 1280 vertices, all 3,815 arc-transitive graphs on ≤ 10,000 vertices, and 267,433
> structured graphs: **`D/girth > 4` was never observed once, and the smallest power-of-two cycle
> tops out at 32. Nothing in the entire campaign ever required a 64-cycle.**

This replaces a vague impression ("spectra are dense") with a measurable invariant, and it is the
first quantitative account I know of for *why* the conjecture should be true.

**(3) The `C16` count of the 46-vertex graph is 314, not 157**: 157 is the *gadget's* count; the
doubled graph has two copies. (My own verification above already gave 314.)

**(4) Girth is a provably dead end.** Killing `C32` by girth alone needs girth ≥ 33, hence
`n ≳ 196,000` for a cubic graph. All 151 graphs of girth ≥ 18 in Conder's census still have a `C32`.

**(5) Truncation identity** (found independently by two agents, neither noticing the other):
`spec(T(H)) = {3} ∪ ⋃_{k ∈ spec(H)} [2k, 3k]`. So `T(G)` is a counterexample **iff** `spec(G)`
avoids `[2^m/3, 2^m/2]` for every `m`, i.e. `spec(G) ⊆ {5} ∪ {9,10} ∪ [17,21] ∪ [33,42] ∪ …`.
The allowed bands have multiplicative width 4/3 while the forbidden ones have width 3/2; the base
must be *more* spectrally sparse than the target. Tested exactly over all 1,717,716 girth-≥5 cubic
graphs on ≤ 24 vertices: **zero survivors**, so no truncation-type counterexample on ≤ 72 vertices.
Corollary: 53 of the 100 vertex-transitive triple-dodgers have girth 3 and are simply truncations
of high-girth graphs, not a second mechanism, as one agent had claimed.

## 5e. Additional rigorous exclusions

| Space | Count | Verdict |
|---|---|---|
| All cubic **vertex-transitive** graphs ≤ 1280 vertices | **111,360** (complete) | 0 counterexamples. Smallest `2^k`: 4→5756, 8→16624, 16→88880, **32→100** |
| All cubic **arc-transitive** graphs ≤ 10,000 vertices | **3,815** (complete) | 0. Smallest `2^k`: 4→3, 8→111, 16→3012, **32→689** |
| 2-terminal gadgets (every edge of every cubic graph on ≤ 20 vertices) | 556,463 × all edges | **No cubic counterexample has a 2-edge-cut with a side of ≤ 20 vertices** |
| All cubic graphs of **girth ≥ 7**, n ≤ 32 | 30,939 | Every one has a `C8` (new beyond Markström at n = 30, 32) |
| **1-terminal gadgets of order ≤ 25 (COMPLETE)** | 149 C4/C8-free at k=25 | All have `C16` ⟹ **no bridged cubic counterexample on ≤ 52 vertices** (gadget orders are odd, so both leaf blocks would need order ≥ 27) |

Best small-`n` near-miss found anywhere: a **28-vertex** LCF graph with spectrum exactly
`{3} ∪ [9,28]`: **five consecutive missing lengths (4,5,6,7,8)**, beating Markström's three.

## 6. Honest assessment so far

No counterexample, from any of nine independent search directions. The evidence points fairly hard
at the conjecture being **true**, and, unusually, it now points *quantitatively*:

- Dodging **two** consecutive powers of two is common (18% of arc-transitive cubic graphs).
  Dodging **three** was never observed in ~10 million exhaustive plus ~380,000 structured and
  symmetric graphs. The dense-onset ratio `D/girth ≈ 2.4` explains exactly why: you can dodge
  about `log₂ 2.4` powers, and you would need `D/girth > 4` to dodge three.
- The smallest power-of-two cycle never exceeded 32 anywhere in the campaign.
- The Erdős–Gyárfás *strong* conjecture is already refuted (Liu–Montgomery), so the original reason
  to expect a counterexample is gone.

**The one structural door still open is low circumference.** A cubic graph whose blocks all have
< 32 vertices need only dodge `{4, 8, 16}`; every larger power is structurally absent. That target
is finite and well-posed, and nobody has searched it. Orders ≤ 25 are now cleared; **27, 29, 31 are
the whole remaining question for bridged graphs**, and 2-edge-cut sides of 21–31 vertices are the
analogous gap for 2-cuts.

Durable outputs of this session, in rough order of value:

1. **The low-circumference reformulation** and the gadget ladder that operationalises it.
2. **The dense-onset ratio `D/girth`** as the measurable invariant governing the conjecture.
3. **Two reductions**: triangle inflation (≈15× on exhaustive search) and the bridge/gadget
   decomposition (reaches n ≤ 62 for the cost of a 31-vertex search).
4. **The near-miss catalogue and hole statistics**: raw material for a discharging proof.
5. **Retirement of three routes**: mod-4 obstructions (dead by a 1993 theorem), girth
   (needs n ≳ 196,000), and truncations (closed to n ≤ 72).
