# -*- coding: utf-8 -*-
"""Final audit for the (TRI) section.

(1) MIXED-CROSSER FORMULA.  For any chord A=(a,b), any X crossing A on
    the left (y < a < p < b) and any X' crossing on the right
    (a < p' < b < y'), the triple {A,X,X'} is ALWAYS a single cycle of
    length  3 + (a-y) + |p-p'| + (y'-b).   Connectivity free.
(2) Does the mixed-crosser family alone settle (TRI)?
(3) Does the ADJACENT-BLOCK family alone settle (TRI)?
(4) Maximum internal gap of the full three-chord spectrum.
"""

import sys, random
from itertools import combinations

POW = set()
x = 1
while x < 1 << 22:
    POW.add(x)
    x <<= 1

def overlap(c, d):
    (a, b), (p, q) = c, d
    return not (b < p or q < a)

def l1(c, d):
    return abs(c[0] - d[0]) + abs(c[1] - d[1])

def hard_ok(c, d):
    if not overlap(c, d):
        return True
    v = l1(c, d)
    return v >= 3 and (v + 2) not in POW

def span_hard(c):
    s = c[1] - c[0]
    return s >= 4 and (s + 1) not in POW

def gen_hard(n, cap, rng):
    out, chords = [], []
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
            out.append(list(chords))
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
            if not span_hard(c):
                continue
            if any(not hard_ok(c, d) for d in chords):
                continue
            used[q] = True
            chords.append(c)
            rec()
            chords.pop()
            used[q] = False
        used[p] = False

    rec()
    return out

def completion(T, n):
    ends = sorted([p for c in T for p in c])
    edges = [tuple(c) for c in T]
    cnt, idx = 0, 0
    for j in range(1, n):
        while idx < len(ends) and ends[idx] <= j:
            idx += 1
            cnt += 1
        if cnt % 2 == 1:
            edges.append((j, j + 1))
    deg = {}
    for u, v in edges:
        deg[u] = deg.get(u, 0) + 1
        deg[v] = deg.get(v, 0) + 1
    par = {v: v for v in deg}

    def find(v):
        while par[v] != v:
            par[v] = par[par[v]]
            v = par[v]
        return v

    for u, v in edges:
        ru, rv = find(u), find(v)
        if ru != rv:
            par[ru] = rv
    return len(edges), len(set(find(v) for v in deg)) <= 1

def mixed(ch):
    """weights of all one-chord + left-crosser + right-crosser triples,
    plus the (formula, connectivity) audit data."""
    W, audit = set(), []
    for A in ch:
        a, b = A
        Ls = [c for c in ch if c[0] < a < c[1] < b]
        Rs = [c for c in ch if a < c[0] < b < c[1]]
        for X in Ls:
            for Xp in Rs:
                y, p = X
                pp, yp = Xp
                w = 3 + (a - y) + abs(p - pp) + (yp - b)
                W.add(w)
                audit.append(((A, X, Xp), w))
    return W, audit

def adjacent(ch, n):
    """weights of triples containing an adjacent block {p, p+1}."""
    mu = {}
    for a, b in ch:
        mu[a] = b
        mu[b] = a
    S = set(map(tuple, ch))
    W = set()
    for p in range(1, n):
        c = tuple(sorted((p, mu[p])))
        cp = tuple(sorted((p + 1, mu[p + 1])))
        if c == cp:
            continue
        for D in ch:
            D = tuple(D)
            if D == c or D == cp:
                continue
            q = sorted([p, p + 1, mu[p], mu[p + 1], D[0], D[1]])
            if len(set(q)) != 6:
                continue
            blocks = [(q[0], q[1]), (q[2], q[3]), (q[4], q[5])]
            if (p, p + 1) not in blocks:
                continue
            if any(bl in S for bl in blocks):
                continue
            W.add(3 + sum(bl[1] - bl[0] for bl in blocks))
    return W

def tri_all(ch):
    S = set(map(tuple, ch))
    W = set()
    for T in combinations(sorted(ch), 3):
        q = sorted([p for c in T for p in c])
        bl = [(q[0], q[1]), (q[2], q[3]), (q[4], q[5])]
        if any(b in S for b in bl):
            continue
        W.add(3 + sum(b[1] - b[0] for b in bl))
    return W

def main():
    sys.setrecursionlimit(20000)
    rng = random.Random(777001)
    ok = bad = 0
    for n in (16, 18, 20):
        for ch in gen_hard(n, 60, rng):
            _, audit = mixed(ch)
            for T, w in audit:
                ww, conn = completion(list(T), n)
                if conn and ww == w:
                    ok += 1
                else:
                    bad += 1
    print("mixed-crosser formula + connectivity:", ok, "ok,", bad, "bad")
    print()
    for n in (18, 22, 26, 30, 36, 42, 50):
        lads = gen_hard(n, 200, rng)
        fm = fa = ft = 0
        mg = 0
        for ch in lads:
            Wm, _ = mixed(ch)
            Wa = adjacent(ch, n)
            Wt = tri_all(ch)
            if not any(w in POW for w in Wm):
                fm += 1
            if not any(w in POW for w in Wa):
                fa += 1
            if not any(w in POW for w in Wt):
                ft += 1
            s = sorted(Wt)
            if len(s) > 1:
                mg = max(mg, max(s[i + 1] - s[i]
                                 for i in range(len(s) - 1)))
        print("2N =", n, "ladders", len(lads),
              " failures - mixed:", fm, " adjacent-block:", fa,
              " all triples:", ft, " max spectrum gap:", mg)

main()
