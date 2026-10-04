# -*- coding: utf-8 -*-
"""What shape are the power-of-two cycles in the hard-case ladders
where NO zigzag works?  Classify the winning chord sets by their
overlap pattern, to find the smallest family that still has
connectivity for free."""

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

def rel(c, d):
    """'X' crossing, 'N' nested, 'D' disjoint (c left of d or not)."""
    (a, b), (p, q) = sorted([c, d])
    if b < p:
        return 'D'
    if q < b:
        return 'N'
    return 'X'

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

def is_zigzag(T):
    S = sorted(T)
    for t in range(len(S) - 1):
        if not (S[t][0] < S[t + 1][0] < S[t][1] < S[t + 1][1]):
            return False
    for t in range(len(S) - 2):
        if not S[t][1] < S[t + 2][0]:
            return False
    return True

def pattern(T):
    """Upper-triangular relation word, chords sorted by left end."""
    S = sorted(T)
    w = []
    for i in range(len(S)):
        for j in range(i + 1, len(S)):
            w.append(rel(S[i], S[j]))
    return ''.join(w)

def has_zigzag_p2(ch, n):
    by_a = sorted(ch)

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
        if len(chain) + sum(segs(chain)) in POW:
            return True
        last = chain[-1]
        prev = chain[-2] if len(chain) >= 2 else None
        for c in by_a:
            if not (last[0] < c[0] < last[1] < c[1]):
                continue
            if prev is not None and not prev[1] < c[0]:
                continue
            if rec(chain + [c]):
                return True
        return False

    return any(rec([c]) for c in by_a)

def main():
    sys.setrecursionlimit(10000)
    from collections import Counter
    shapes = Counter()
    sizes = Counter()
    examples = {}
    for n in (14, 16, 18):
        lads = enumerate_ladders(n)
        for ch in lads:
            if has_zigzag_p2(ch, n):
                continue
            best = None
            for r in range(1, len(ch) + 1):
                for T in combinations(ch, r):
                    w, c = completion(list(T), n)
                    if c and w in POW:
                        if best is None or r < best[0]:
                            best = (r, T, w)
                if best is not None:
                    break
            r, T, w = best
            sizes[r] += 1
            p = pattern(T)
            shapes[p] += 1
            examples.setdefault(p, (n, ch, T, w))
    print("minimal winning |T| histogram:", dict(sizes))
    print()
    for p, cnt in shapes.most_common(12):
        n, ch, T, w = examples[p]
        print("pattern", p, "count", cnt, " zigzag?", is_zigzag(T))
        print("   2N =", n, "T =", list(T), "weight", w)

main()
