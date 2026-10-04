# -*- coding: utf-8 -*-
"""(XC*) - the unified one-jump-each criterion.

Two outside chords with ord endpoints (alpha,beta) and (gamma,delta)
- each usable in either orientation - such that

    gamma < min(u, alpha) ,   max(u, alpha) < delta <= v < beta .

Then
    A = {a,p} : a-end ... beta -chord- alpha ... u -> p
    B = {b,q} : b-end ... gamma -chord- delta ... v -> q
occupy  [min(alpha,u), max(alpha,u)] + [beta, max]  and
        [0, gamma] + [delta, v] ,
which are pairwise disjoint, so the ears are vertex-disjoint, and
a < q < p < b interleaves.  (XC) of 512.0 is the sub-case alpha >= u;
the sliver of 512.1 is the sub-case alpha <= u.
"""

import sys, random
from itertools import combinations
from k4lib import ov

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

def check(ch, n, a, b):
    """returns None if not residual, else ('XC','XC2','XCstar-no')"""
    mu = {}
    for x, y in ch:
        mu[x] = y
        mu[y] = x

    def ordf(m):
        return (m - (b + 1)) if m > b else (n - b - 1 + m)

    X = [p for p in range(a + 1, b) if not (a < mu[p] < b)]
    if len(X) != b - a - 1:
        return None
    f = {p: ordf(mu[p]) for p in X}
    for q, p in combinations(X, 2):
        if f[q] < f[p]:
            return None
    oc = []
    for (x, y) in ch:
        if not (a <= x <= b) and not (a <= y <= b):
            oc.append((ordf(x), ordf(y)))
    ors = []
    for (x, y) in oc:
        ors.append((x, y))
        ors.append((y, x))
    best = None
    for q, p in combinations(X, 2):
        u, v = f[p], f[q]
        for (al, be) in ors:
            for (ga, de) in ors:
                if {al, be} == {ga, de}:
                    continue
                if ga < min(u, al) and max(u, al) < de <= v < be:
                    tag = 'XC' if al >= u else 'XC2'
                    if best != 'XC':
                        best = tag
                    if best == 'XC':
                        return 'XC'
    return best or 'XCstar-no'

def main():
    sys.setrecursionlimit(30000)
    rng = random.Random(161803)
    from collections import Counter
    for n, lads in [(14, enum_all(14)), (16, sample(16, 120, rng)),
                    (18, sample(18, 80, rng)),
                    (20, sample(20, 60, rng)),
                    (22, sample(22, 40, rng)),
                    (24, sample(24, 30, rng))]:
        c = Counter()
        for ch in lads:
            for (a, b) in ch:
                r = check(ch, n, a, b)
                if r:
                    c[r] += 1
        tot = sum(c.values())
        print("2N =", n, " residual one-chord cycles:", tot,
              " settled by (XC*):", tot - c['XCstar-no'],
              " breakdown:", dict(c))

main()
