# -*- coding: utf-8 -*-
"""Decisive test: do the HARD-CASE girth-legal ladders have a
power-of-two cycle at all?  If not, (LQ3) itself is false - which is
possible, because positions 1 and 2N have degree 2 in the ladder, so
Erdos-Gyarfas does not apply to it.

Full connected spectrum by brute force over all chord subsets."""

import sys
from itertools import combinations

POW = set()
x = 1
while x < 1 << 20:
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

def enumerate_ladders(n):
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
    conn = len(set(find(v) for v in deg)) <= 1 and len(deg) > 0
    return len(edges), conn

def spectrum(ch, n):
    S = set()
    N = len(ch)
    for r in range(1, N + 1):
        for T in combinations(ch, r):
            w, c = completion(list(T), n)
            if c:
                S.add(w)
    return S

def zigzag_pow2(ch, n):
    by_a = sorted(ch)
    found = []

    def segs(chain):
        m = len(chain)
        if m == 1:
            return [chain[0][1] - chain[0][0]]
        f = [chain[1][0] - chain[0][0]]
        for t in range(m - 2):
            f.append(chain[t + 2][0] - chain[t][1])
        f.append(chain[-1][1] - chain[-2][1])
        return f

    def rec(chain):
        w = len(chain) + sum(segs(chain))
        if w in POW:
            found.append((w, len(chain)))
        last = chain[-1]
        prev = chain[-2] if len(chain) >= 2 else None
        for c in by_a:
            if not (last[0] < c[0] < last[1] < c[1]):
                continue
            if prev is not None and not prev[1] < c[0]:
                continue
            rec(chain + [c])

    for c in by_a:
        rec([c])
    return found

def main():
    sys.setrecursionlimit(10000)
    for n in (14, 16, 18, 20):
        lads = enumerate_ladders(n)
        noz, nop2 = 0, []
        for ch in lads:
            if zigzag_pow2(ch, n):
                continue
            noz += 1
            S = spectrum(ch, n)
            hits = sorted([w for w in S if w in POW])
            if not hits:
                nop2.append((ch, sorted(S)))
        print("2N =", n, " hard ladders:", len(lads),
              " no zigzag p2:", noz,
              " NO POWER-OF-TWO CYCLE AT ALL:", len(nop2))
        for ch, S in nop2[:4]:
            print("    ", ch)
            print("       spectrum:", S)

main()
