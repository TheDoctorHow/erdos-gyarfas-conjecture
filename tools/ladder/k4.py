# -*- coding: utf-8 -*-
"""THE K4 TEMPLATE.

Three pairwise crossing chords a1<a2<a3<b1<b2<b3 of a ladder form a
SUBDIVISION OF K4 with branch vertices a2, a3, b1, b2 and branch
lengths
    p = a3-a2   (a2 a3)      u = 1           (a2 b2, chord c2)
    q = b1-a3   (a3 b1)      v = 1+b3-b2     (a3 b2, via c3)
    r = b2-b1   (b1 b2)      w = 1+a2-a1     (a2 b1, via c1)
K4 has exactly seven cycles: four triangles S - star(x) and three
quadrilaterals S - M_i (M_i the three perfect matchings).  So the same
six numbers give SEVEN cycle lengths, where 501.4's XXX formula used
one - and that one is the quadrilateral S - (p+r).

Everything is checked against the completion machinery.
"""

import sys, random
from itertools import combinations

POW = set(1 << i for i in range(1, 24))

def ov(c, d):
    return not (c[1] < d[0] or d[1] < c[0])

def soft_pair(c, d):
    if not ov(c, d):
        return True
    return abs(c[0] - d[0]) + abs(c[1] - d[1]) >= 3

def hard_pair(c, d):
    if not ov(c, d):
        return True
    v = abs(c[0] - d[0]) + abs(c[1] - d[1])
    return v >= 3 and (v + 2) not in POW

def hard_span(c):
    s = c[1] - c[0]
    return s >= 4 and (s + 1) not in POW

def gen(n, cap, rng, hard=True):
    out, ch = [], []
    used = [False] * (n + 2)
    stop = []
    pair = hard_pair if hard else soft_pair

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
            if hard and not hard_span(c):
                continue
            if not hard and c[1] - c[0] < 4:
                continue
            if any(not pair(c, d) for d in ch):
                continue
            used[q] = True
            ch.append(c)
            rec()
            ch.pop()
            used[q] = False
        used[p] = False

    rec()
    return out

def completion(T, n):
    ends = sorted([x for c in T for x in c])
    edges = [tuple(c) for c in T]
    cnt, idx = 0, 0
    for j in range(1, n):
        while idx < len(ends) and ends[idx] <= j:
            idx += 1
            cnt += 1
        if cnt % 2 == 1:
            edges.append((j, j + 1))
    deg = {}
    for a, b in edges:
        deg[a] = deg.get(a, 0) + 1
        deg[b] = deg.get(b, 0) + 1
    par = {v: v for v in deg}

    def find(v):
        while par[v] != v:
            par[v] = par[par[v]]
            v = par[v]
        return v

    for a, b in edges:
        ra, rb = find(a), find(b)
        if ra != rb:
            par[ra] = rb
    return len(edges), len(set(find(v) for v in deg)) <= 1

def seven(c1, c2, c3):
    """the seven K4 cycle lengths of a pairwise crossing triple"""
    (a1, b1), (a2, b2), (a3, b3) = c1, c2, c3
    p = a3 - a2
    q = b1 - a3
    r = b2 - b1
    u = 1
    v = 1 + b3 - b2
    w = 1 + a2 - a1
    S = p + q + r + u + v + w
    stars = [p + w + u, p + q + v, w + q + r, u + v + r]
    match = [p + r, w + v, u + q]
    return [S - x for x in stars] + [S - x for x in match], S

def xxx(c1, c2, c3):
    return 3 + (c1[1] - c1[0]) - (c2[1] - c2[0]) + (c3[1] - c3[0])

def triples(ch):
    for T in combinations(sorted(ch), 3):
        (a1, b1), (a2, b2), (a3, b3) = T
        if a1 < a2 < a3 < b1 < b2 < b3:
            yield T

def main():
    sys.setrecursionlimit(20000)
    rng = random.Random(606060)

    # ---- A: the XXX formula is the quadrilateral S-(p+r)
    ok = bad = 0
    for n in (16, 18, 20):
        for ch in gen(n, 40, rng):
            for T in triples(ch):
                L, S = seven(*T)
                if L[6 - 2] == xxx(*T) or xxx(*T) in L:
                    pass
                # quadrilateral S-(p+r) is index 4
                if L[4] == xxx(*T):
                    ok += 1
                else:
                    bad += 1
    print("XXX formula = the quadrilateral S-(p+r):", ok, "ok,", bad,
          "bad")

    # ---- B: all seven are genuine cycle lengths of the ladder
    ok = bad = 0
    for n in (16, 18):
        for ch in gen(n, 25, rng):
            spec = set()
            for r_ in range(1, len(ch) + 1):
                for T in combinations(ch, r_):
                    w, c = completion(list(T), n)
                    if c:
                        spec.add(w)
            for T in triples(ch):
                L, S = seven(*T)
                for x in L:
                    if x in spec:
                        ok += 1
                    else:
                        bad += 1
    print("all seven are realised cycle lengths:", ok, "ok,", bad,
          "bad")

    # ---- C: does the K4 template alone settle (LQ3)?
    print()
    for n in (18, 22, 26, 30, 36, 44):
        lads = gen(n, 150, rng)
        f7 = f1 = notri = 0
        for ch in lads:
            got7 = got1 = False
            any_t = False
            for T in triples(ch):
                any_t = True
                L, S = seven(*T)
                if any(x in POW for x in L):
                    got7 = True
                if xxx(*T) in POW:
                    got1 = True
            if not any_t:
                notri += 1
            if not got7:
                f7 += 1
            if not got1:
                f1 += 1
        print("2N =", n, " hard ladders", len(lads),
              " no crossing triple:", notri,
              " K4-seven fails:", f7,
              " XXX-alone fails:", f1)

main()
