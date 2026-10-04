# -*- coding: utf-8 -*-
"""(TRI): every girth-legal ladder has a CONNECTED THREE-CHORD cycle
of power-of-two weight.

Connectivity criterion for three chords (PROVED): sort the six
endpoints q1<...<q6; the completion is a single cycle iff none of
{q1,q2}, {q3,q4}, {q5,q6} is itself a chord.  Weight is then
    3 + (q2-q1) + (q4-q3) + (q6-q5).

Tested exhaustively over hard-case ladders (no one- or two-chord power
of two) at 2N = 14..22, and by randomised search at 2N = 24..34.
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

def enumerate_ladders(n, cap=None, shuffle=False, rng=None):
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
        rng_list = list(range(p + 4, n + 1))
        if shuffle:
            rng.shuffle(rng_list)
        for q in rng_list:
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

def tri_weights(ch):
    """All connected three-chord weights, with the witnesses."""
    S = set(map(tuple, ch))
    out = {}
    for T in combinations(sorted(ch), 3):
        q = sorted([p for c in T for p in c])
        if ((q[0], q[1]) in S or (q[2], q[3]) in S
                or (q[4], q[5]) in S):
            continue
        w = 3 + (q[1] - q[0]) + (q[3] - q[2]) + (q[5] - q[4])
        out.setdefault(w, []).append(T)
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
    return len(edges), len(set(find(v) for v in deg)) <= 1

def main():
    sys.setrecursionlimit(20000)
    rng = random.Random(31337)

    # criterion cross-check on small exhaustive ladders
    bad = 0
    chk = 0
    for n in (14, 16):
        for ch in enumerate_ladders(n):
            S = set(map(tuple, ch))
            for T in combinations(sorted(ch), 3):
                q = sorted([p for c in T for p in c])
                pred = not ((q[0], q[1]) in S or (q[2], q[3]) in S
                            or (q[4], q[5]) in S)
                w, conn = completion(list(T), n)
                if pred != conn:
                    bad += 1
                else:
                    chk += 1
                if conn:
                    w2 = 3 + (q[1]-q[0]) + (q[3]-q[2]) + (q[5]-q[4])
                    if w != w2:
                        bad += 1
    print("three-chord criterion+weight:", chk, "ok,", bad, "bad")
    print()

    for n in range(14, 23, 2):
        lads = enumerate_ladders(n)
        fail = []
        for ch in lads:
            tw = tri_weights(ch)
            if not any(w in POW for w in tw):
                fail.append(ch)
        print("2N =", n, "hard ladders", len(lads),
              " (TRI) fails on", len(fail))
        for f in fail[:3]:
            print("    ", f)

    print()
    for n in range(24, 35, 2):
        lads = enumerate_ladders(n, cap=3000, shuffle=True, rng=rng)
        fail = []
        for ch in lads:
            tw = tri_weights(ch)
            if not any(w in POW for w in tw):
                fail.append(ch)
        print("2N =", n, "sampled hard ladders", len(lads),
              " (TRI) fails on", len(fail))
        for f in fail[:3]:
            print("    ", f)

main()
