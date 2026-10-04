# -*- coding: utf-8 -*-
"""(EAR*) for ONE-CHORD cycles, case by case.

Z = chord (a,b) + host [a,b].  In the CLOSED ladder the outside
O = [1,a-1] + closing edge + [b+1,2N] is a single path P_out, oriented
from the b-end (b+1) round to the a-end (a-1).  Write X for the
interior positions whose chord leaves [a,b] (the CROSSERS) and f(p)
for the position of mu(p) along P_out.

PROVED CASES (to be confirmed by machine):
 (i)   two inside chords cross                        -> interleaved
 (ii)  an inside chord strictly contains a crosser    -> interleaved
       (with the ear {a, crosser})
 (iii) no inside chords, and f has an ASCENT (q < p in X with
       f(q) < f(p))                                   -> interleaved
       (ears {a,p} and {q,b}, disjoint segments)
 RESIDUAL: no inside chords and f strictly decreasing.

We check that (i)-(iii) really do produce an interleaving disjoint
pair, count the residual, and test (EAR*) there.
"""

import sys, random
from itertools import combinations
from k4lib import ov, cycle_edges, cyc_order, all_ears, seven_from

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

def closed_adj(ch, n):
    adj = {v: [] for v in range(1, n + 1)}
    for j in range(1, n):
        adj[j].append(j + 1)
        adj[j + 1].append(j)
    for a, b in ch:
        adj[a].append(b)
        adj[b].append(a)
    adj[1].append(n)
    adj[n].append(1)
    return adj

def ok_cycle(adj, order, E):
    ears = all_ears(adj, order, E)
    for e1, e2 in combinations(ears, 2):
        if seven_from(order, e1, e2) or seven_from(order, e2, e1):
            return True
    return False

def classify_one_chord(ch, n, a, b):
    mu = {}
    for x, y in ch:
        mu[x] = y
        mu[y] = x
    inside = []
    X = []
    for p in range(a + 1, b):
        m = mu[p]
        if a < m < b:
            if p < m:
                inside.append((p, m))
        else:
            X.append(p)
    for (p, q), (r, s) in combinations(inside, 2):
        if p < r < q < s or r < p < s < q:
            return 'i'
    for (p, q) in inside:
        if any(p < x < q for x in X):
            return 'ii'
    if inside:
        return 'ii?'          # inside chord with no crosser inside
    # f = position along P_out, oriented b+1 ... 2N, 1 ... a-1
    def f(x):
        m = mu[x]
        return (m - (b + 1)) if m > b else (n - b - 1 + m)
    for q, p in combinations(X, 2):     # q < p
        if f(q) < f(p):
            return 'iii'
    return 'residual'

def main():
    sys.setrecursionlimit(30000)
    rng = random.Random(31415)
    from collections import Counter
    for n, lads in [(12, enum_all(12)), (14, enum_all(14)),
                    (16, sample(16, 60, rng)),
                    (18, sample(18, 40, rng)),
                    (20, sample(20, 30, rng))]:
        kinds = Counter()
        bad = 0
        resfail = 0
        for ch in lads:
            adj = closed_adj(ch, n)
            for (a, b) in ch:
                E = cycle_edges([(a, b)], n)
                if E is None:
                    continue
                k = classify_one_chord(ch, n, a, b)
                kinds[k] += 1
                good = ok_cycle(adj, cyc_order(E), E)
                if not good:
                    bad += 1
                    print("   (EAR*) FAILS", (a, b), "in", ch)
                if k == 'residual' and not good:
                    resfail += 1
        print("2N =", n, " ladders", len(lads),
              " one-chord cycles by case:", dict(kinds),
              " (EAR*) failures:", bad)

main()
