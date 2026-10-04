# -*- coding: utf-8 -*-
"""(PAN) is false.  Repaired target:

  (PAN*) LOCAL PANCYCLICITY IN THE CAPTURE WINDOW.
     If H is bipartite (equivalently all spans odd) its spectrum
     contains every EVEN integer in [g, 2g-2];
     otherwise its spectrum contains EVERY integer in [g, 2g-2].

Only the window [g, 2g-2] matters, because 504.1 captures a power of
two there and 504.0 proves c >= 2g-2.  Tested exhaustively over ALL
girth-legal ladders - not just the hard case, so the (N9) test of
499.4 is satisfied.
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

def main():
    sys.setrecursionlimit(20000)
    for n in (12, 14, 16, 18):
        lads = enumerate_soft(n)
        fstar = []
        fcirc = 0
        fp2 = 0
        bip = 0
        for ch in lads:
            S = spectrum(ch, n)
            if not S:
                continue
            g, c = S[0], S[-1]
            if c < 2 * g - 2:
                fcirc += 1
            b = all((y - x) % 2 == 1 for x, y in ch)
            bip += b
            want = [t for t in range(g, 2 * g - 1)
                    if (not b) or t % 2 == 0]
            if any(t not in S for t in want):
                fstar.append((ch, S, b, want))
            if not any(t in POW for t in S):
                fp2 += 1
        print("2N =", n, " girth-legal", len(lads),
              " bipartite", bip,
              " (PAN*) failures:", len(fstar),
              " (CIRC) failures:", fcirc,
              " no power-of-two cycle:", fp2)
        for ch, S, b, want in fstar[:4]:
            print("     bipartite?", b, ch)
            print("       spectrum", S, " window wanted", want)

main()
