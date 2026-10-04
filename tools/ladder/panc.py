# -*- coding: utf-8 -*-
"""Is a hard-case girth-legal ladder PANCYCLIC above its girth?
If the full connected spectrum is the whole interval [girth,
circumference] and the circumference is at least twice the girth, then
(LQ3) follows at once from the capture lemma."""

import sys, random
from itertools import combinations

POW = set(1 << i for i in range(1, 24))

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

def enumerate_all(n, cap=None, rng=None):
    out, chords = [], []
    used = [False] * (n + 2)
    stop = []

    def rec():
        if cap is not None and len(out) >= cap:
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
        if rng:
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
    rng = random.Random(1234)
    for n in (14, 16, 18, 20):
        lads = enumerate_all(n)
        bad = []
        ratios = []
        for ch in lads:
            S = spectrum(ch, n)
            g, c = S[0], S[-1]
            if len(S) != c - g + 1:
                bad.append((ch, S))
            ratios.append(c / float(g))
        ratios.sort()
        print("2N =", n, "hard ladders", len(lads),
              " NOT gapless:", len(bad),
              " circumference/girth min %.2f med %.2f"
              % (ratios[0], ratios[len(ratios) // 2]))
        for ch, S in bad[:3]:
            print("     ", ch)
            print("       spectrum", S)
    for n in (22, 26, 30):
        lads = enumerate_all(n, cap=60, rng=rng)
        bad = 0
        ratios = []
        for ch in lads:
            S = spectrum(ch, n)
            g, c = S[0], S[-1]
            if len(S) != c - g + 1:
                bad += 1
            ratios.append(c / float(g))
        ratios.sort()
        print("2N =", n, "sampled", len(lads), " NOT gapless:", bad,
              " c/g min %.2f med %.2f"
              % (ratios[0], ratios[len(ratios) // 2]))

main()
