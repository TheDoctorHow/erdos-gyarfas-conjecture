# -*- coding: utf-8 -*-
"""Is (PAN*) over-strong?  Petersen has girth 5 and no 7-cycle, so a
ladder analogue would refute (PAN*) while leaving the captured power
of two intact.  "No cycle of length t" is monotone, so search for it.
"""

import sys
from itertools import combinations

def overlap(c, d):
    (a, b), (p, q) = c, d
    return not (b < p or q < a)

def soft_pair(c, d):
    if not overlap(c, d):
        return True
    return abs(c[0] - d[0]) + abs(c[1] - d[1]) >= 3

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

def search(n, target):
    used = [False] * (n + 2)
    chords = []
    nodes = [0]
    hit = []
    maxr = target // 2

    def makes_target(newc):
        others = [c for c in chords if c != newc]
        for r in range(0, maxr):
            for rest in combinations(others, r):
                w, conn = completion(list(rest) + [newc], n)
                if conn and w == target:
                    return True
        return False

    def rec():
        nodes[0] += 1
        p = 1
        while p <= n and used[p]:
            p += 1
        if p > n:
            hit.append(list(chords))
            return True
        used[p] = True
        for q in range(p + 4, n + 1):
            if used[q]:
                continue
            c = (p, q)
            if any(not soft_pair(c, d) for d in chords):
                continue
            used[q] = True
            chords.append(c)
            if not makes_target(c):
                if rec():
                    return True
            chords.pop()
            used[q] = False
        used[p] = False
        return False

    ok = rec()
    return (hit[0] if ok else None), nodes[0]

def main():
    sys.setrecursionlimit(30000)
    for tgt in (6, 7, 9, 10, 11, 12):
        print("=== target length", tgt)
        for n in range(12, 27, 2):
            w, nodes = search(n, tgt)
            if w is None:
                print("   2N =", n, " none; nodes", nodes)
            else:
                S = spectrum(w, n)
                g = S[0]
                win = [t for t in range(g, 2 * g - 1)]
                print("   2N =", n, " LADDER WITH NO", tgt,
                      "-CYCLE:", w)
                print("      spectrum", S)
                print("      girth", g, " window", win,
                      " window complete?",
                      all(t in S for t in win))
                break

main()
