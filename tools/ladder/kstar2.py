# -*- coding: utf-8 -*-
"""516: k* in the CLOSED ladder, with a monotone prune.

k*(G) = the least number of chords used by a power-of-two cycle of the
closed ladder G = C(1..2N) + perfect matching.  k* = 0 means the host
cycle itself is a power of two, so 2N is taken to be a NON-power of
two throughout - which is the faithful case, since a J7 rung's host
cycle has length 2^t m with m odd >= 3.

The apparatus of 495-515 lives entirely at k* <= 3:
  k* = 1 : a one-chord cycle          (C + chord  = theta)
  k* = 2 : a two-chord cycle          (C + crossing pair = K_4)
  k* = 3 : a three-chord cycle        (C + crossing triple = K_{3,3})
A ladder with k* >= 4 is one on which every through-C template the
archive owns is silent.

PRUNE.  "Some cycle with at most j chords has power-of-two length" is
MONOTONE in the chord set: adding chords never deletes a cycle.  So
the DFS over matchings may abandon a branch as soon as the placed
chords give a power-of-two cycle with <= 2 chords.  Cycles are found
by brute-force enumeration of the small subgraph, never by a formula.
"""

import sys
from itertools import combinations


def cyc_lengths(chords, n, want):
    """all cycle lengths of C(1..n) + the given chords, as a dict
    length -> min chords used.  Brute force DFS, no formulas."""
    adj = dict((v, []) for v in range(1, n + 1))
    for i in range(1, n + 1):
        j = i % n + 1
        adj[i].append(j)
        adj[j].append(i)
    cs = set()
    for a, b in chords:
        adj[a].append(b)
        adj[b].append(a)
        cs.add((a, b) if a < b else (b, a))
    best = {}
    for s in range(1, n + 1):
        stack = [(s, 1, 0, 1 << s)]
        while stack:
            v, L, nc, msk = stack.pop()
            for w in adj[v]:
                e = (v, w) if v < w else (w, v)
                inc = 1 if e in cs else 0
                if w == s:
                    if L >= 3:
                        c = nc + inc
                        if L not in best or best[L] > c:
                            best[L] = c
                elif w > s and not (msk >> w & 1):
                    stack.append((w, L + 1, nc + inc, msk | (1 << w)))
    if want is None:
        return best
    return best


def small_hit(chords, n, P2, jmax):
    """does some cycle using at most jmax chords have 2-power length?
    enumerated over subsets of size <= jmax, brute force."""
    for r in range(0, jmax + 1):
        for S in combinations(chords, r):
            for L, c in cyc_lengths(list(S), n, None).items():
                if c == r and L in P2:
                    return True
    return False


def gen(n, P2):
    """matchings of 1..n, closed girth >= 5, pruned on k* <= 2."""
    mate = [0] * (n + 1)
    chords = []
    out = []

    def cdist(x, y):
        d = abs(x - y) % n
        return min(d, n - d)

    def rec():
        a = 0
        for i in range(1, n + 1):
            if not mate[i]:
                a = i
                break
        if not a:
            out.append(list(chords))
            return
        for b in range(a + 1, n + 1):
            if mate[b]:
                continue
            if cdist(a, b) < 4:
                continue
            bad = False
            for c, d in chords:
                for x, y in ((c, d), (d, c)):
                    if cdist(a, x) + cdist(b, y) + 2 <= 4:
                        bad = True
            if bad:
                continue
            # 2-power test for the NEW chord alone and paired
            newc = (a, b)
            if any(L in P2 and cc == 1
                   for L, cc in cyc_lengths([newc], n, None).items()):
                continue
            hit = False
            for old in chords:
                for L, cc in cyc_lengths([newc, old], n, None).items():
                    if cc == 2 and L in P2:
                        hit = True
            if hit:
                continue
            mate[a] = b
            mate[b] = a
            chords.append(newc)
            rec()
            chords.pop()
            mate[a] = 0
            mate[b] = 0

    rec()
    return out


def sweep(n, show=4):
    assert n & (n - 1), "2N must not itself be a power of two"
    P2 = set()
    v = 4
    while v <= 2 * n:
        P2.add(v)
        v *= 2
    surv = gen(n, P2)
    hist = {}
    wit = []
    for ch in surv:
        best = cyc_lengths(ch, n, None)
        k = min([c for L, c in best.items() if L in P2] or [99])
        hist[k] = hist.get(k, 0) + 1
        if k >= 4:
            wit.append((k, ch))
    print("2N=%2d  k*>=3 survivors=%6d   k* histogram: %s"
          % (n, len(surv), sorted(hist.items())))
    for k, ch in wit[:show]:
        print("   k*=%d : %s" % (k, ch))
    sys.stdout.flush()
    return surv, hist, wit


if __name__ == '__main__':
    for nn in (12, 14, 18, 20):
        sweep(nn)
