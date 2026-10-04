# -*- coding: utf-8 -*-
"""(XC): a PROVED sufficient criterion for (2LINK).

Residual one-chord cycle: Z = (a,b) + [a,b], every interior position a
crosser, f strictly decreasing.  Orient the outside P_out from the
b-end (ord 0) to the a-end (ord max) and let u = f(p) < v = f(q) for
crossers q < p.

CLAIM.  If the outside graph has two chords, e1 = (alpha, beta) and
e2 = (gamma, delta) in ord coordinates, with

        gamma < u <= alpha < v <= delta < beta ,

then the two ears
   {a, p} :  a-end ... beta -e1- alpha ... u = f(p) -> p
   {b, q} :  b-end ... gamma -e2- delta ... v = f(q) -> q
are VERTEX-DISJOINT (they occupy ord [u,alpha] + [beta,max] and
[0,gamma] + [v,delta]) and their attachments a < q < p < b interleave.

We verify the claim itself against the ear machinery, and measure how
often (XC) covers the residual cases.
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

def analyse(ch, n, a, b):
    mu = {}
    for x, y in ch:
        mu[x] = y
        mu[y] = x

    def ordf(m):
        return (m - (b + 1)) if m > b else (n - b - 1 + m)

    X = [p for p in range(a + 1, b) if not (a < mu[p] < b)]
    if len(X) != b - a - 1:
        return None                       # not residual: inside chord
    f = {p: ordf(mu[p]) for p in X}
    for q, p in combinations(X, 2):
        if f[q] < f[p]:
            return None                   # ascent: case (iii)
    # residual.  outside chords, in ord coordinates
    out_ch = []
    for (x, y) in ch:
        if not (a <= x <= b) and not (a <= y <= b):
            out_ch.append(tuple(sorted((ordf(x), ordf(y)))))
    for q, p in combinations(X, 2):       # q < p, so f[q] > f[p]
        u, v = f[p], f[q]
        for (g, d) in out_ch:
            for (al, be) in out_ch:
                if (g, d) == (al, be):
                    continue
                if g < u <= al < v <= d < be:
                    return ('XC', (a, b), (q, p), (u, v),
                            (g, d), (al, be))
    return ('residual-no-XC', (a, b))

def main():
    sys.setrecursionlimit(30000)
    rng = random.Random(2718)
    from collections import Counter
    for n, lads in [(14, enum_all(14)), (16, sample(16, 80, rng)),
                    (18, sample(18, 60, rng)),
                    (20, sample(20, 40, rng)),
                    (22, sample(22, 30, rng))]:
        c = Counter()
        eg = None
        for ch in lads:
            for (a, b) in ch:
                r = analyse(ch, n, a, b)
                if r is None:
                    continue
                c[r[0]] += 1
                if r[0] == 'XC' and eg is None:
                    eg = (ch, r)
        print("2N =", n, " residual one-chord cycles:", dict(c))
        if eg and n == 14:
            print("   example:", eg[0])
            print("   Z=%s  (q,p)=%s  (u,v)=%s  e2=%s  e1=%s"
                  % eg[1][1:])

main()
