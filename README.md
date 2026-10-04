# The Erdős–Gyárfás conjecture: a proof-programme archive

> **The Erdős–Gyárfás conjecture is open, and nothing in this repository proves
> it.** This is a research log: a long, self-auditing record of attempts, partial
> results, computations, refutations and corrections. It has not been refereed.

**The conjecture.** Every finite graph with minimum degree at least 3 contains a
cycle whose length is a power of two. By a theorem of Liu and Montgomery
([arXiv:2010.15802](https://arxiv.org/abs/2010.15802), *JAMS* 2023) it holds
whenever the average degree is large enough, so a counterexample would have to be
sparse and essentially cubic; the programme therefore works on cubic graphs.
[`RESEARCH-BRIEF.md`](RESEARCH-BRIEF.md) surveys the literature.

## What is here

| Path | What it is |
|---|---|
| [`PROOF-PROGRAM.md`](PROOF-PROGRAM.md) | The archive: 449 numbered sections (§1–§598), about 3 MB. Only things that are proved, computed, refuted or withdrawn, with the reasons. |
| [`INDEX.md`](INDEX.md) | Section number → line number → title. Start here: GitHub will not render the archive in the browser, so use the raw view or a clone. |
| [`RESEARCH-BRIEF.md`](RESEARCH-BRIEF.md) | Literature review and the reframing to the bounded-degree (cubic) case. |
| [`RESULTS.md`](RESULTS.md) | Log of the exhaustive counterexample hunt over cubic graphs, and the infrastructure built for it. |
| [`tools/`](tools) | Search and verification code: C filters for snarkhunter (Brinkmann–Goedgebeur–McKay, not included), exact cycle-spectrum machinery, an independent reference implementation, and ~180 analysis scripts. |
| [`tools/ladder/`](tools/ladder) | Scripts from the later part of the programme, up to §598, including the closed-ladder work (see "Provenance" below). |

## How to read the archive

In the later sections, each opens with a `STATUS:` line and closes with a ledger
whose headings say exactly what kind of claim was made: `PROVED`, `COMPUTED`,
`VERIFIED`, `REFUTED`, `CORRECTED`, `WITHDRAWN`, `IDENTIFIED`, `NOT PROVED`, and
restates that the conjecture is **unproved**.

**Errors are recorded loudly.** The archive's standing rule (cited throughout as
"374.1") is that a mistake is written up as a `CORRECTED` or `WITHDRAWN` entry
with the evidence, not quietly fixed. Recent examples: §594 shows that the
"margin" of §592 was computed by a search that could not reach half the
configurations, and that the quantity itself measured the wrong thing; §597.0
corrects a sign error in §596.2. The companion rule ("388.1") is to check a
published statement verbatim before using it. A claim labelled `PROVED` means
proved on the archive's own ledger, not independently checked.

## Where the work stands (through §598)

The recent sections study **closed ladders**: a Hamiltonian cycle (the *host*)
plus a perfect matching of *chords*, i.e. exactly the Hamiltonian cubic graphs.
"Girth-legal" means girth ≥ 5.

- **Verified exhaustively:** every girth-legal closed ladder on at most 32
  vertices has a cycle whose length is the largest power of two not exceeding its
  number of vertices (§591, §593).
- **The deficiency-2 telescope (§§595–598).** Find a cycle through all but two
  vertices, re-present the graph with that cycle as host, and repeat; the descent
  stops at the power of two. The step is an induced subgraph, so girth is
  inherited. It reached the power of two on every ladder tried and failed on none
  of 550 adversarial searches, which leaves two open pillars: **(T1)** every
  closed ladder has a cycle of length exactly 2N−2, and **(T2)** the descent tree
  cannot be starved at every branch at once. §596 checks the literature (verbatim)
  and finds (T1) in neither of the closest published results.
- **Refuted along the way:** an O(1) bound on the shortest power-of-two cycle
  (§576) and its proposed replacement L_min ≤ 2P(g) (§579); a bounded forcing
  number (§592); induction by chord contraction (§590); and the non-backtracking
  trace frame at the large power of two (§589). §597 shows that the number of
  chords (T1) needs is not bounded by 3 (eleven 18-vertex ladders need four), and
  §598 identifies a family measured to need exactly four at every size tested.
- The archive also develops a reduction chain from the conjecture to the
  statements it calls (LQ3)/(K4) and J7, on its own ledger with recorded
  provisos (§§431, 434, 464). Read those sections for the exact claims.

## Reproducing the computations

Python 3 (the scripts were run on 3.14, Windows 11) and a C compiler for the
`.c` files. **nauty and snarkhunter are not included**; build them from their
authors' sites as described in [`RESULTS.md`](RESULTS.md) §1.

Each of these regenerates a number recorded in the archive:

```sh
cd tools/ladder
python depth2.py 12 14        # exact margins 3, 3; 118 and 1055 search nodes   (§594)
python pushL.py 12 14 18      # no girth-legal ladder lacks the power-of-two cycle (§591)
python descend2.py U 20 22    # the telescope reaches 16 on U(20,5) and U(22,5)   (§595)
python t1fast.py 12 14 16     # no girth-legal ladder lacks a (2N-2)-cycle        (§595)
```

## Provenance

- **Authorship.** The programme was carried out by the repository owner in an
  extended collaboration with Claude (Anthropic's AI model) through Claude Code.
  The owner set the goals; Claude wrote most of the code and drafted the sections.
  Treat every mathematical claim as unreviewed by outside experts.
- **`tools/ladder/` was recovered.** Its scripts were originally written to a
  temporary folder that was later deleted. They were rebuilt by replaying the
  recorded file writes and edits from the session transcripts. As a check, 102 of
  the 338 files recovered this way were identical (up to line endings) to copies
  that had already been saved in `tools/`, and the deterministic numbers above
  reproduce exactly. Hard-coded temporary paths were rewritten to "this script's
  own folder".
- **These are research notebooks, not a library.** Many scripts are one-off
  probes, some are superseded by later ones, and a few no longer run (for example
  `bic3.py` raises a `TypeError` in a call to a sibling script). They are kept
  because they are the record of what was computed.

## License

Code (everything under `tools/`, plus `a4an.py` and `cnt04.c`): MIT, see
[`LICENSE`](LICENSE). Written material (the `.md` files): CC BY 4.0, see
[`LICENSE-TEXT.md`](LICENSE-TEXT.md), which also notes that short quotations
from published papers remain their authors'.
