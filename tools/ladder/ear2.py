# -*- coding: utf-8 -*-
"""Large verification of (EAR), plus the structural statistics the
proof needs: how often is the winning ear pair of the SIMPLE kinds -
(i) two chords of Z that cross, (ii) one ear left of min(Z) and one
right of max(Z), (iii) anything else?"""

import sys, random
from itertools import combinations
from k4lib import (ov, build, cycle_edges, cyc_order, all_ears,
                   seven_from)

def soft_pair(c, d):
    if not ov(c, d):
        return True
    return abs(c[0] - d[0]) + abs(c[1] - d[1]) >= 3

def enum_all(n):
    out, ch = [], []
    used = [False] * (n + 2)

    def rec():
        p = 1
        while p <= n and used[p]:
            p += 1
        if p > n:
            out.append(list(ch))
            return
        used[p] = True
        for q in range(p + 4, n + 1):
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

def sample(n, cap, rng):
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

def classify(adj, order, E, n):
    """returns (found, kind) for the FIRST interleaving disjoint pair,
    preferring kind 1 then 2 then 3"""
    ears = all_ears(adj, order, E)
    mn, mx = min(order), max(order)
    best = None
    for e1, e2 in combinations(ears, 2):
        if not (seven_from(order, e1, e2) or seven_from(order, e2,
                                                        e1)):
            continue
        len1, len2 = len(e1[2]), len(e2[2])
        if len1 == 0 and len2 == 0:
            return True, 1                      # two chords of Z
        i1 = set(e1[2])
        i2 = set(e2[2])
        if (i1 and max(i1) < mn and i2 and min(i2) > mx) or \
           (i2 and max(i2) < mn and i1 and min(i1) > mx):
            best = best or 2                    # one left, one right
        else:
            best = best or 3
    return (best is not None), (best or 0)

def run(lads, n, tag):
    from collections import Counter
    bad = 0
    kinds = Counter()
    tot = 0
    for ch in lads:
        adj = build(ch, n)
        for r in range(1, len(ch) + 1):
            for T in combinations(ch, r):
                E = cycle_edges(list(T), n)
                if E is None:
                    continue
                tot += 1
                order = cyc_order(E)
                ok, k = classify(adj, order, E, n)
                kinds[k] += 1
                if not ok:
                    bad += 1
                    if bad < 4:
                        print("   FAIL", sorted(order), "in", ch)
    print(tag, "2N =", n, " ladders", len(lads), " cycles", tot,
          " (EAR) failures:", bad, " winning kind:", dict(kinds))

def main():
    sys.setrecursionlimit(30000)
    rng = random.Random(777)
    for n in (12, 14):
        run(enum_all(n), n, "EXHAUSTIVE")
    for n in (16, 18, 20, 22):
        run(sample(n, 40, rng), n, "sampled")

main()
