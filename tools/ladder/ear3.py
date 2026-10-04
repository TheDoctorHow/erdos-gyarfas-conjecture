# -*- coding: utf-8 -*-
"""The (EAR) obstruction, and its repair.

HYPOTHESIS: (EAR) fails only for cycles through position 1 or 2N - the
two DEGREE-2 vertices, which have no leaving edge, so they contribute
no attachment.  Those two vertices are an artefact of the bare ladder:
in a J7 rung the feet attach to C, so they have degree 3.

REPAIR (EAR+): close the ladder with one extra edge between 1 and 2N
- the minimal model of the C-arc - and re-test.
"""

import sys, random
from itertools import combinations
from k4lib import ov, cycle_edges, cyc_order, all_ears, seven_from

def soft_pair(c, d):
    if not ov(c, d):
        return True
    return abs(c[0] - d[0]) + abs(c[1] - d[1]) >= 3

def enum_all(n):
    out, ch = [], []
    used = [False] * (n + 2)

    def rec():
        p = 1
        while p <= n and used[p]:
            p += 1
        if p > n:
            out.append(list(ch))
            return
        used[p] = True
        for q in range(p + 4, n + 1):
            if used[q]:
                continue
            c = (p, q)
            if c[1] - c[0] < 4:
                continue
            if any(not soft_pair(c, d) for d in ch):
                continue
            used[q] = True
            ch.append(c)
            rec()
            ch.pop()
            used[q] = False
        used[p] = False

    rec()
    return out

def sample(n, cap, rng):
    out, ch = [], []
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
            out.append(list(ch))
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
            if c[1] - c[0] < 4:
                continue
            if any(not soft_pair(c, d) for d in ch):
                continue
            used[q] = True
            ch.append(c)
            rec()
            ch.pop()
            used[q] = False
        used[p] = False

    rec()
    return out

def build_plus(ch, n, closed):
    adj = {v: [] for v in range(1, n + 1)}
    for j in range(1, n):
        adj[j].append(j + 1)
        adj[j + 1].append(j)
    for a, b in ch:
        adj[a].append(b)
        adj[b].append(a)
    if closed:
        adj[1].append(n)
        adj[n].append(1)
    return adj

def ok_cycle(adj, order, E):
    ears = all_ears(adj, order, E)
    for e1, e2 in combinations(ears, 2):
        if seven_from(order, e1, e2) or seven_from(order, e2, e1):
            return True
    return False

def audit(lads, n, tag):
    bad_open, bad_closed = [], []
    tot = 0
    ends = 0
    for ch in lads:
        A0 = build_plus(ch, n, False)
        A1 = build_plus(ch, n, True)
        for r in range(1, len(ch) + 1):
            for T in combinations(ch, r):
                E = cycle_edges(list(T), n)
                if E is None:
                    continue
                tot += 1
                order = cyc_order(E)
                if not ok_cycle(A0, order, E):
                    bad_open.append((ch, sorted(order)))
                    if 1 in order or n in order:
                        ends += 1
                    if not ok_cycle(A1, order, E):
                        bad_closed.append((ch, sorted(order)))
    print(tag, "2N =", n, " cycles", tot,
          " (EAR) failures OPEN:", len(bad_open),
          " of which contain 1 or 2N:", ends,
          " (EAR+) failures CLOSED:", len(bad_closed))
    for b in bad_open[:4]:
        print("      open-fail", b[1], " in", b[0],
              " lengths", len(b[1]))
    for b in bad_closed[:4]:
        print("      CLOSED-FAIL", b[1], " in", b[0])

def main():
    sys.setrecursionlimit(30000)
    rng = random.Random(555)
    for n in (12, 14):
        audit(enum_all(n), n, "EXHAUSTIVE")
    for n in (16, 18, 20):
        audit(sample(n, 30, rng), n, "sampled")

main()
