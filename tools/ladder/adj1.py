# -*- coding: utf-8 -*-
"""(ADJ), unified.

VIRTUAL-CHORD THEOREM (to verify).  For p in [1, 2N-1] let
    Q_p = ( min(mu(p), mu(p+1)),  max(mu(p), mu(p+1)) )
be the VIRTUAL CHORD pairing the partners of p and p+1.  Let D be a
chord other than M(p), M(p+1) such that
    (a) D overlaps Q_p, and
    (b) D and Q_p have the same parity in the number of endpoints
        strictly below p.
Then {M(p), M(p+1), D} is a single cycle of length  4 + L1(Q_p, D).

Also measured: the (ADJ) spectrum's range and largest internal gap,
the surplus (how many (p,D) pairs hit a power of two), and whether
the MIXED-TYPE gaps (t_p != t_{p+1}) alone suffice.
"""

import sys, random
from itertools import combinations

POW = set(1 << i for i in range(1, 22))

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

def adj_pairs(ch, n):
    """(p, D, weight, mixed) for every admissible pair, by the
    VIRTUAL-CHORD rule."""
    mu = {}
    for a, b in ch:
        mu[a] = b
        mu[b] = a
    out = []
    for p in range(1, n):
        c = tuple(sorted((p, mu[p])))
        cp = tuple(sorted((p + 1, mu[p + 1])))
        if c == cp:
            continue
        Q = (min(mu[p], mu[p + 1]), max(mu[p], mu[p + 1]))
        beta = (mu[p] < p) + (mu[p + 1] < p)
        mixed = ((mu[p] > p) != (mu[p + 1] > p + 1))
        for D in ch:
            D = tuple(D)
            if D == c or D == cp:
                continue
            if not overlap(D, Q):
                continue
            gamma = (D[0] < p) + (D[1] < p)
            if (beta + gamma) % 2:
                continue
            out.append((p, D, 4 + l1(Q, D), mixed))
    return out

def block_weight(ch, n, p, D):
    """independent recomputation via the three-chord criterion"""
    mu = {}
    for a, b in ch:
        mu[a] = b
        mu[b] = a
    S = set(map(tuple, ch))
    q = sorted([p, p + 1, mu[p], mu[p + 1], D[0], D[1]])
    if len(set(q)) != 6:
        return None
    bl = [(q[0], q[1]), (q[2], q[3]), (q[4], q[5])]
    if (p, p + 1) not in bl:
        return None
    if any(b in S for b in bl):
        return None
    return 3 + sum(b[1] - b[0] for b in bl)

def main():
    sys.setrecursionlimit(20000)
    rng = random.Random(5150)
    ok = bad = 0
    miss = 0
    for n in (16, 18, 20):
        for ch in gen_hard(n, 40, rng):
            mu = {}
            for a, b in ch:
                mu[a] = b
                mu[b] = a
            pairs = adj_pairs(ch, n)
            got = set((p, D) for p, D, _, _ in pairs)
            for p, D, w, _ in pairs:
                bw = block_weight(ch, n, p, D)
                ww, conn = completion([tuple(sorted((p, mu[p]))),
                                       tuple(sorted((p + 1,
                                                     mu[p + 1]))),
                                       D], n)
                if bw == w and conn and ww == w:
                    ok += 1
                else:
                    bad += 1
            # completeness: any (p,D) the block rule accepts that the
            # virtual-chord rule rejects?
            for p in range(1, n):
                for D in ch:
                    D = tuple(D)
                    if block_weight(ch, n, p, D) is not None:
                        if (p, D) not in got:
                            miss += 1
    print("virtual-chord theorem:", ok, "ok,", bad, "bad,",
          miss, "missed by the rule")
    print()
    for n in (18, 22, 26, 30, 36, 42, 50, 60):
        lads = gen_hard(n, 150, rng)
        f_all = f_mix = 0
        maxgap = 0
        surplus = []
        for ch in lads:
            pr = adj_pairs(ch, n)
            W = sorted(set(w for _, _, w, _ in pr))
            Wm = set(w for _, _, w, m in pr if m)
            hits = [1 for _, _, w, _ in pr if w in POW]
            surplus.append(len(hits))
            if not any(w in POW for w in W):
                f_all += 1
            if not any(w in POW for w in Wm):
                f_mix += 1
            if len(W) > 1:
                maxgap = max(maxgap, max(W[i + 1] - W[i]
                                         for i in range(len(W) - 1)))
        surplus.sort()
        print("2N =", n, "ladders", len(lads),
              " (ADJ) failures:", f_all,
              " mixed-type-only failures:", f_mix,
              " max spectrum gap:", maxgap,
              " surplus min/med:", surplus[0],
              surplus[len(surplus) // 2])

main()
