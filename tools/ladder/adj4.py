# -*- coding: utf-8 -*-
"""EXHAUSTIVE search for a counterexample to (ADJ).

Build the matching left to right.  Prune on
  (i) hard case: span >= 4, span != 2^M - 1, overlapping pairs at
      L1 >= 3 and L1 != 2^M - 2;
  (ii) (ADJ)-FAILURE: no fully determined admissible pair (Q_p, D)
      may have 4 + L1(Q_p, D) a power of two.
A pair is fully determined as soon as mu(p), mu(p+1) and D are placed,
and it survives into every completion, so the prune is sound.

If the search completes a matching, (ADJ) IS FALSE.
If it exhausts, no hard-case counterexample of that size exists.
"""

import sys

POW = set(1 << i for i in range(1, 24))

def run(n, verbose=True):
    mu = [0] * (n + 2)
    chords = []
    nodes = [0]
    found = []

    def overlap(a, b, c, d):
        return not (b < c or d < a)

    def ok_pair(c1, c2):
        a, b = c1
        c, d = c2
        if not overlap(a, b, c, d):
            return True
        v = abs(a - c) + abs(b - d)
        return v >= 3 and (v + 2) not in POW

    def ok_span(a, b):
        s = b - a
        return s >= 4 and (s + 1) not in POW

    def adj_bad(newc):
        """Is some fully determined admissible pair a power of two?
        Only pairs involving the newly placed chord can be new."""
        for p in range(1, n):
            if mu[p] == 0 or mu[p + 1] == 0:
                continue
            cp = (min(p, mu[p]), max(p, mu[p]))
            cq = (min(p + 1, mu[p + 1]), max(p + 1, mu[p + 1]))
            if cp == cq:
                continue
            Q = (min(mu[p], mu[p + 1]), max(mu[p], mu[p + 1]))
            beta = (mu[p] < p) + (mu[p + 1] < p)
            # only pairs touching newc are new, but the pair (p, D)
            # is new also when newc completes M(p) or M(p+1)
            fresh = (newc == cp or newc == cq)
            for D in chords:
                if D == cp or D == cq:
                    continue
                if not (fresh or D == newc):
                    continue
                if not overlap(D[0], D[1], Q[0], Q[1]):
                    continue
                if (beta + (D[0] < p) + (D[1] < p)) % 2:
                    continue
                if 4 + abs(Q[0] - D[0]) + abs(Q[1] - D[1]) in POW:
                    return True
        return False

    def rec():
        nodes[0] += 1
        p = 1
        while p <= n and mu[p]:
            p += 1
        if p > n:
            found.append(list(chords))
            return True
        for q in range(p + 4, n + 1):
            if mu[q]:
                continue
            if not ok_span(p, q):
                continue
            c = (p, q)
            if any(not ok_pair(c, d) for d in chords):
                continue
            mu[p] = q
            mu[q] = p
            chords.append(c)
            if not adj_bad(c):
                if rec():
                    return True
            chords.pop()
            mu[p] = 0
            mu[q] = 0
        return False

    hit = rec()
    if verbose:
        if hit:
            print("2N =", n, " *** (ADJ) COUNTEREXAMPLE ***",
                  found[0], " nodes", nodes[0])
        else:
            print("2N =", n, " no counterexample; nodes", nodes[0])
    return hit

def main():
    sys.setrecursionlimit(30000)
    for n in range(10, 41, 2):
        if run(n):
            break

main()
