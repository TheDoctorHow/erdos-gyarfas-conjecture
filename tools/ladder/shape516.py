# -*- coding: utf-8 -*-
"""516: for the k* >= 3 survivors, WHICH SHAPE of triple delivers the
power of two?  C + three chords has cycle space of dimension 4 and
nine edges after contraction, so it is a cubic graph on six vertices:
K_{3,3} exactly when the triple is pairwise crossing, otherwise the
prism or a multigraph, neither of which is a template.

  crossing  -> C + triple IS a K_{3,3}-subdivision (515.4), and the
               archive's fifteen-length ledger reaches the cycle.
  other     -> the cycle is a genuine cycle of the ladder, but it is
               in NO template through C, so every through-C statement
               the archive owns is silent on it.

This decides whether the layered approaches cover the k* = 3 world.
"""

import sys
from itertools import combinations
from kstar2 import cyc_lengths, sweep


def shape(t):
    (a1, b1), (a2, b2), (a3, b3) = sorted(t)
    if a1 < a2 < a3 < b1 < b2 < b3:
        return 'crossing'
    cr = 0
    for (x1, y1), (x2, y2) in combinations(sorted(t), 2):
        if x1 < x2 < y1 < y2:
            cr += 1
    return 'cross%d' % cr


def analyse(n):
    surv, hist, wit = sweep(n, show=0)
    P2 = set()
    v = 4
    while v <= 2 * n:
        P2.add(v)
        v *= 2
    tally = {}
    only_noncross = []
    for ch in surv:
        shapes = set()
        for T in combinations(ch, 3):
            got = cyc_lengths(list(T), n, None)
            for L, c in got.items():
                if c == 3 and L in P2:
                    shapes.add(shape(T))
        for s in shapes:
            tally[s] = tally.get(s, 0) + 1
        if 'crossing' not in shapes:
            only_noncross.append((ch, sorted(shapes)))
    print("2N=%2d  survivors=%d  shapes delivering the 2-power: %s"
          % (n, len(surv), sorted(tally.items())))
    print("       survivors with NO crossing triple delivering it: %d"
          % len(only_noncross))
    for ch, s in only_noncross[:5]:
        print("   ", s, ch)
    sys.stdout.flush()


if __name__ == '__main__':
    analyse(18)
