# -*- coding: utf-8 -*-
"""(N9) check on (PAN): the HARD CASE is defined by a condition that
mentions powers of two, so "hard-case ladders are pancyclic" has a
conclusion-mentioning subject.  The legitimate statement is (PAN) for
ALL girth-legal ladders.  Test it there, exhaustively.

Also verifies the (CIRC) theorem c >= 2g - 2 on every ladder found.
"""

import sys
from itertools import combinations

POW = set(1 << i for i in range(1, 24))

def overlap(c, d):
    (a, b), (p, q) = c, d
    return not (b < p or q < a)

def soft_pair(c, d):
    if not overlap(c, d):
        return True
    return abs(c[0] - d[0]) + abs(c[1] - d[1]) >= 3

def soft_span(c):
    return c[1] - c[0] >= 4

def enumerate_soft(n):
    out, chords = [], []
    used = [False] * (n + 2)

    def rec():
        p = 1
        while p <= n and used[p]:
            p += 1
        if p > n:
            out.append(list(chords))
            return
        used[p] = True
        for q in range(p + 4, n + 1):
            if used[q]:
                continue
            c = (p, q)
            if not soft_span(c):
                continue
            if any(not soft_pair(c, d) for d in chords):
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

def spectrum(ch, n):
    S = set()
    for r in range(1, len(ch) + 1):
        for T in combinations(ch, r):
            w, c = completion(list(T), n)
            if c:
                S.add(w)
    return sorted(S)

def is_hard(ch):
    for a, b in ch:
        if (b - a + 1) in POW:
            return False
    for c, d in combinations(ch, 2):
        if overlap(c, d):
            if (abs(c[0] - d[0]) + abs(c[1] - d[1]) + 2) in POW:
                return False
    return True

def main():
    sys.setrecursionlimit(20000)
    for n in (10, 12, 14, 16, 18):
        lads = enumerate_soft(n)
        badpan = []
        badcirc = 0
        nohard = 0
        for ch in lads:
            S = spectrum(ch, n)
            g, c = S[0], S[-1]
            if len(S) != c - g + 1:
                badpan.append((ch, S, is_hard(ch)))
            if c < 2 * g - 2:
                badcirc += 1
            if not is_hard(ch):
                nohard += 1
        print("2N =", n, " girth-legal ladders", len(lads),
              " (of which hard-case", len(lads) - nohard, ")",
              " (PAN) FAILURES:", len(badpan),
              " (CIRC) c < 2g-2:", badcirc)
        for ch, S, h in badpan[:4]:
            print("     hard?", h, ch)
            print("       spectrum", S)

main()
