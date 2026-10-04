# -*- coding: utf-8 -*-
"""Exact (K4) test.

Enumerate every K4-subdivision as (cycle Z, two internally disjoint
ears whose attachments interleave on Z) - this presentation catches
all of them, since any K4-subdivision is one of its quadrilaterals
plus the two remaining branch paths.  Ears are found by explicit DFS,
so ears inside the SAME bridge are included, which the bridge test of
k4b.py missed.  Then compute all seven lengths.
"""

import sys, random
from itertools import combinations

POW = set(1 << i for i in range(1, 24))

def ov(c, d):
    return not (c[1] < d[0] or d[1] < c[0])

def hard_pair(c, d):
    if not ov(c, d):
        return True
    v = abs(c[0] - d[0]) + abs(c[1] - d[1])
    return v >= 3 and (v + 2) not in POW

def hard_span(c):
    s = c[1] - c[0]
    return s >= 4 and (s + 1) not in POW

def gen(n, cap, rng):
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
            if not hard_span(c):
                continue
            if any(not hard_pair(c, d) for d in ch):
                continue
            used[q] = True
            ch.append(c)
            rec()
            ch.pop()
            used[q] = False
        used[p] = False

    rec()
    return out

def build(ch, n):
    adj = {v: [] for v in range(1, n + 1)}
    for j in range(1, n):
        adj[j].append(j + 1)
        adj[j + 1].append(j)
    for a, b in ch:
        adj[a].append(b)
        adj[b].append(a)
    return adj

def cycle_edges(T, n):
    ends = sorted([x for c in T for x in c])
    edges = [tuple(sorted(c)) for c in T]
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
    if any(v != 2 for v in deg.values()):
        return None
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
    if len(set(find(v) for v in deg)) != 1:
        return None
    return edges

def cyc_order(edges):
    a2 = {}
    for a, b in edges:
        a2.setdefault(a, []).append(b)
        a2.setdefault(b, []).append(a)
    start = min(a2)
    order = [start]
    prev, cur = None, start
    while True:
        nxt = a2[cur][0] if a2[cur][0] != prev else a2[cur][1]
        if nxt == start:
            break
        order.append(nxt)
        prev, cur = cur, nxt
    return order

def all_ears(adj, order, cycedges):
    """every path with both ends on Z and interior off Z"""
    onZ = set(order)
    ce = set(cycedges)
    ears = []

    def dfs(start, cur, path, seen):
        for nxt in adj[cur]:
            if cur == start and tuple(sorted((cur, nxt))) in ce:
                continue
            if nxt in onZ:
                if nxt != start and (len(path) > 1 or
                                     tuple(sorted((cur, nxt)))
                                     not in ce):
                    ears.append((start, nxt, tuple(path[1:])))
                continue
            if nxt in seen:
                continue
            dfs(start, nxt, path + [nxt], seen | {nxt})

    for s in order:
        dfs(s, s, [s], {s})
    uniq = {}
    for u, v, inte in ears:
        key = (min(u, v), max(u, v), tuple(sorted(inte)))
        if key not in uniq:
            uniq[key] = (u, v, inte)
    return list(uniq.values())

def seven_from(order, e1, e2):
    """Z + two interleaving internally disjoint ears -> 7 lengths"""
    m = len(order)
    pos = {v: i for i, v in enumerate(order)}
    (u1, u2, i1) = e1
    (v1, v2, i2) = e2
    p1, p2 = sorted((pos[u1], pos[u2]))
    q1, q2 = sorted((pos[v1], pos[v2]))
    if not (p1 < q1 < p2 < q2):
        return None
    if set(i1) & set(i2):
        return None
    if len({u1, u2, v1, v2}) != 4:
        return None
    al = q1 - p1
    be = p2 - q1
    ga = q2 - p2
    de = m - q2 + p1
    pi = len(i1) + 1
    ro = len(i2) + 1
    S = al + be + ga + de + pi + ro
    stars = [al + de + pi, al + be + ro, be + ga + pi, ga + de + ro]
    match = [al + ga, be + de, pi + ro]
    return [S - x for x in stars] + [S - x for x in match]

def k4_holds(ch, n, want_all=False):
    adj = build(ch, n)
    hits = []
    for r in range(1, len(ch) + 1):
        for T in combinations(ch, r):
            E = cycle_edges(list(T), n)
            if E is None:
                continue
            order = cyc_order(E)
            ears = all_ears(adj, order, E)
            for e1, e2 in combinations(ears, 2):
                L = seven_from(order, e1, e2)
                if L is None:
                    L = seven_from(order, e2, e1)
                if L is None:
                    continue
                for x in L:
                    if x in POW:
                        if not want_all:
                            return True, x
                        hits.append(x)
    return (len(hits) > 0, hits) if want_all else (False, None)

def main():
    sys.setrecursionlimit(30000)
    rng = random.Random(9182736)
    for n in (14, 16, 18):
        lads = gen(n, 120, rng)
        fails = []
        for ch in lads:
            ok, _ = k4_holds(ch, n)
            if not ok:
                fails.append(ch)
        print("2N =", n, " hard-case ladders", len(lads),
              " (K4) FAILS on:", len(fails))
        for f in fails[:3]:
            print("     ", f)

main()
