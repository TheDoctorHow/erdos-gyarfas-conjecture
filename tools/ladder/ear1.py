# -*- coding: utf-8 -*-
"""(EAR): does EVERY cycle of a girth-legal ladder have two internally
disjoint ears with interleaving attachments?  Where does it fail?"""

import sys, random
from itertools import combinations
from k4lib import (ov, hard_pair, hard_span, build, cycle_edges,
                   cyc_order, all_ears, seven_from)

POW = set(1 << i for i in range(1, 24))

def soft_pair(c, d):
    if not ov(c, d):
        return True
    return abs(c[0] - d[0]) + abs(c[1] - d[1]) >= 3

def gen_soft(n, cap, rng):
    out, ch = [], []
    used = [False] * (n + 2)
    stop = []

    def rec():
        if len(out) >= cap:
            stop.append(1)
            return
        p = 1
        while p <= n and used[p]:
            p += 1
        if p > n:
            out.append(list(ch))
            return
        used[p] = True
        cand = list(range(p + 4, n + 1))
        rng.shuffle(cand)
        for q in cand:
            if stop:
                break
            if used[q]:
                continue
            c = (p, q)
            if c[1] - c[0] < 4:
                continue
            if any(not soft_pair(c, d) for d in ch):
                continue
            used[q] = True
            ch.append(c)
            rec()
            ch.pop()
            used[q] = False
        used[p] = False

    rec()
    return out

def has_two_ears(adj, order, E):
    ears = all_ears(adj, order, E)
    for e1, e2 in combinations(ears, 2):
        if seven_from(order, e1, e2) or seven_from(order, e2, e1):
            return True
    return False

def main():
    sys.setrecursionlimit(30000)
    rng = random.Random(20260811)
    from collections import Counter
    for n in (14, 16, 18):
        lads = gen_soft(n, 60, rng)
        tot = Counter()
        bad = Counter()
        examples = {}
        endsonZ = Counter()
        for ch in lads:
            adj = build(ch, n)
            for r in range(1, len(ch) + 1):
                for T in combinations(ch, r):
                    E = cycle_edges(list(T), n)
                    if E is None:
                        continue
                    L = len(E)
                    order = cyc_order(E)
                    tot[L] += 1
                    if not has_two_ears(adj, order, E):
                        bad[L] += 1
                        k = (1 in order) + (n in order)
                        endsonZ[k] += 1
                        examples.setdefault(L, (ch, sorted(order)))
        print("2N =", n)
        print("   cycles by length :",
              dict(sorted(tot.items())))
        print("   (EAR) FAILURES   :",
              dict(sorted(bad.items())))
        print("   failures by #ends(1 or 2N) on Z:", dict(endsonZ))
        for L in sorted(examples)[:3]:
            ch, o = examples[L]
            print("      len", L, "cycle", o, " in", ch)

main()
