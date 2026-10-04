# -*- coding: utf-8 -*-
"""Structure of the three-chord spectrum of a hard-case ladder.

Is it an interval?  What is its largest internal gap, and its range?
Does it contain [x, 2x]?  And which sub-families already suffice:
 (S1) triples containing an ADJACENT BLOCK {p,p+1};
 (S2) triples with two adjacent blocks;
 (S3) all triples.
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

def analyse(ch):
    S = set(map(tuple, ch))
    full, adj1, adj2 = set(), set(), set()
    for T in combinations(sorted(ch), 3):
        q = sorted([p for c in T for p in c])
        blocks = [(q[0], q[1]), (q[2], q[3]), (q[4], q[5])]
        if any(b in S for b in blocks):
            continue
        gaps = [b[1] - b[0] for b in blocks]
        w = 3 + sum(gaps)
        full.add(w)
        ones = gaps.count(1)
        if ones >= 1:
            adj1.add(w)
        if ones >= 2:
            adj2.add(w)
    return full, adj1, adj2

def gapinfo(S):
    if not S:
        return None
    s = sorted(S)
    gaps = [s[i + 1] - s[i] for i in range(len(s) - 1)]
    return s[0], s[-1], (max(gaps) if gaps else 0), len(s)

def main():
    sys.setrecursionlimit(20000)
    rng = random.Random(20260810)
    for n in (18, 22, 26, 30, 36, 42):
        lads = gen_hard(n, 300, rng)
        worst_gap = 0
        interval_all = True
        f1 = f2 = 0
        rng_lo = rng_hi = 0
        for ch in lads:
            full, adj1, adj2 = analyse(ch)
            lo, hi, mg, cnt = gapinfo(full)
            worst_gap = max(worst_gap, mg)
            if mg > 1:
                interval_all = False
            rng_lo, rng_hi = lo, hi
            if not any(w in POW for w in adj1):
                f1 += 1
            if not any(w in POW for w in adj2):
                f2 += 1
        print("2N =", n, "ladders", len(lads),
              " 3-chord spectrum: max internal gap", worst_gap,
              " sample range", (rng_lo, rng_hi),
              " interval?", interval_all)
        print("      failures with >=1 adjacent block:", f1,
              "  with >=2 adjacent blocks:", f2)

main()
