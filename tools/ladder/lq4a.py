# -*- coding: utf-8 -*-
"""(LQ4) attack, step 1.

A. Verify the ZIGZAG WEIGHT IDENTITY
       W = m + phi(c_m) - phi(c_1) - sum_{t=2}^{m-1} sigma(c_t)
         = [2 + phi_m - phi_1] - sum_{t=2}^{m-1} (sigma_t - 1)
   against the edge-count weight, for every zigzag family.

B. Exhaustively enumerate the HARD CASE: girth-legal ladders with no
   one-chord and no two-chord power-of-two cycle (span never 2^M - 1,
   non-disjoint L^1 never 2^M - 2).  For each, the zigzag spectrum,
   the powers of two it hits, and the minimal m that hits one.
"""

import sys
from itertools import combinations

POW = set()
x = 1
while x < 1 << 20:
    POW.add(x)
    x <<= 1

def is_pow2(v):
    return v in POW

def overlap(c, d):
    (a, b), (p, q) = c, d
    return not (b < p or q < a)

def l1(c, d):
    return abs(c[0] - d[0]) + abs(c[1] - d[1])

def hard_ok(c, d):
    """girth-legal AND no two-chord power of two, for a pair."""
    if not overlap(c, d):
        return True
    v = l1(c, d)
    return v >= 3 and not is_pow2(v + 2)

def soft_ok(c, d):
    if not overlap(c, d):
        return True
    return l1(c, d) >= 3

def span_hard(c):
    s = c[1] - c[0]
    return s >= 4 and not is_pow2(s + 1)

def span_soft(c):
    return c[1] - c[0] >= 4

# ------------------------------------------------------- exhaustive gen

def enumerate_ladders(n, pair_ok, span_ok, cap=None):
    out = []
    chords = []
    used = [False] * (n + 2)

    def rec():
        if cap is not None and len(out) >= cap:
            return
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
            if not span_ok(c):
                continue
            if any(not pair_ok(c, d) for d in chords):
                continue
            used[q] = True
            chords.append(c)
            rec()
            chords.pop()
            used[q] = False
        used[p] = False

    rec()
    return out

# ----------------------------------------------------------- zigzag set

def zigzag_weights(chords, n, want_min_m=True):
    """All zigzag families, by DFS extension.  Returns dict
    weight -> minimal m."""
    by_a = sorted(chords)
    res = {}
    # m = 1 is not a zigzag chain of interest on its own?  It is: the
    # single-chord cycle.  Include it.
    for c in by_a:
        w = 1 + (c[1] - c[0])
        if w not in res or res[w] > 1:
            res[w] = 1

    def extend(chain):
        m = len(chain)
        w = m + sum(f for f in segs(chain))
        if w not in res or res[w] > m:
            res[w] = m
        last, prev = chain[-1], (chain[-2] if m >= 2 else None)
        for c in by_a:
            a, b = c
            la, lb = last
            if not (la < a < lb < b):
                continue
            if prev is not None and not (prev[1] < a):
                continue
            extend(chain + [c])

    def segs(chain):
        m = len(chain)
        if m == 1:
            return [chain[0][1] - chain[0][0]]
        f = [chain[1][0] - chain[0][0]]
        for t in range(m - 2):
            f.append(chain[t + 2][0] - chain[t][1])
        f.append(chain[-1][1] - chain[-2][1])
        return f

    for c in by_a:
        extend([c])
    return res

# --------------------------------------------------------- verification

def completion_weight(T, n):
    ends = sorted([p for c in T for p in c])
    cnt, idx, tot = 0, 0, 0
    for j in range(1, n):
        while idx < len(ends) and ends[idx] <= j:
            idx += 1
            cnt += 1
        if cnt % 2 == 1:
            tot += 1
    return len(T) + tot

def identity_check(chords, n):
    """W = m + phi_m - phi_1 - sum interior sigma."""
    bad = 0
    ok = 0
    by_a = sorted(chords)

    def rec(chain):
        nonlocal bad, ok
        m = len(chain)
        if m >= 2:
            w = completion_weight(chain, n)
            phi1 = chain[0][0] + chain[0][1]
            phim = chain[-1][0] + chain[-1][1]
            ident = m + phim - phi1 - sum(c[1] - c[0]
                                          for c in chain[1:-1])
            if w != ident:
                bad += 1
            else:
                ok += 1
        last = chain[-1]
        prev = chain[-2] if m >= 2 else None
        for c in by_a:
            a, b = c
            if not (last[0] < a < last[1] < b):
                continue
            if prev is not None and not prev[1] < a:
                continue
            rec(chain + [c])

    for c in by_a:
        rec([c])
    return ok, bad

def main():
    sys.setrecursionlimit(10000)
    # ---- A: identity, over soft (merely girth-legal) ladders
    tot_ok = tot_bad = 0
    for n in (10, 12, 14):
        lads = enumerate_ladders(n, soft_ok, span_soft, cap=400)
        for ch in lads:
            o, b = identity_check(ch, n)
            tot_ok += o
            tot_bad += b
    print("A. zigzag weight identity:", tot_ok, "ok,", tot_bad, "bad")

    # ---- B: the hard case
    for n in range(10, 23, 2):
        lads = enumerate_ladders(n, hard_ok, span_hard)
        hits = 0
        minm = {}
        worst = []
        for ch in lads:
            res = zigzag_weights(ch, n)
            p2 = sorted([w for w in res if is_pow2(w)])
            if p2:
                hits += 1
                mm = min(res[w] for w in p2)
                minm[mm] = minm.get(mm, 0) + 1
            else:
                worst.append(ch)
        print("2N =", n, " hard-case ladders:", len(lads),
              " with a zigzag power of two:", hits,
              " minimal m histogram:", minm)
        for w in worst[:3]:
            print("     NO ZIGZAG POWER OF TWO:", w)

main()
