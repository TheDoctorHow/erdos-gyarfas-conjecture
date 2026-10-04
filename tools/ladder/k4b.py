# -*- coding: utf-8 -*-
"""Is (K4) a REFORMULATION of (LQ3) or a STRENGTHENING?

A K4-subdivision may be presented as a cycle Z together with two
internally disjoint ears whose attachment pairs INTERLEAVE on Z; then
Z is one of its seven cycles (the quadrilateral S - C).  So

    Z is one of the seven of some K4-subdivision
      <==  two BRIDGES of Z have interleaving attachments.

If every power-of-two cycle has interleaving bridges, then
(LQ3) implies (K4), so (K4) is a reformulation and escapes (N16).
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
    """edge set of the completion of T, or None if disconnected"""
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

def cyclic_order(edges):
    adj = {}
    for a, b in edges:
        adj.setdefault(a, []).append(b)
        adj.setdefault(b, []).append(a)
    start = next(iter(adj))
    order = [start]
    prev, cur = None, start
    while True:
        nxt = adj[cur][0] if adj[cur][0] != prev else adj[cur][1]
        if nxt == start:
            break
        order.append(nxt)
        prev, cur = cur, nxt
    return order

def interleaving_bridges(adj, order, n):
    """Do two bridges of the cycle have interleaving attachments?"""
    pos = {v: i for i, v in enumerate(order)}
    onZ = set(order)
    seen = set()
    bridges = []
    for v in range(1, n + 1):
        if v in onZ or v in seen:
            continue
        comp, att = set(), set()
        stack = [v]
        seen.add(v)
        comp.add(v)
        while stack:
            x = stack.pop()
            for y in adj[x]:
                if y in onZ:
                    att.add(y)
                elif y not in seen:
                    seen.add(y)
                    comp.add(y)
                    stack.append(y)
        if len(att) >= 2:
            bridges.append(sorted(pos[a] for a in att))
    # chords of Z that are not cycle edges are also bridges
    cyc = set()
    for i in range(len(order)):
        cyc.add(tuple(sorted((order[i], order[(i + 1) % len(order)]))))
    for a in onZ:
        for b in adj[a]:
            if b in onZ and tuple(sorted((a, b))) not in cyc:
                if pos[a] < pos[b]:
                    bridges.append([pos[a], pos[b]])
    for i in range(len(bridges)):
        for j in range(i + 1, len(bridges)):
            Bi, Bj = bridges[i], bridges[j]
            for x1, x2 in combinations(Bi, 2):
                for y1, y2 in combinations(Bj, 2):
                    if (x1 < y1 < x2 < y2) or (y1 < x1 < y2 < x2):
                        return True
    return False

def main():
    sys.setrecursionlimit(20000)
    rng = random.Random(9182736)

    print("=== the 507.0 witness (no 8-cycle, has a 16-cycle)")
    W = [(1, 5), (2, 7), (3, 11), (4, 17), (6, 22), (8, 13), (9, 21),
         (10, 15), (12, 16), (14, 19), (18, 23), (20, 24)]
    n = 24
    adj = build(W, n)
    tot = good = 0
    for r in range(1, len(W) + 1):
        for T in combinations(W, r):
            E = cycle_edges(list(T), n)
            if E is None or len(E) not in POW:
                continue
            tot += 1
            if interleaving_bridges(adj, cyclic_order(E), n):
                good += 1
    print("   power-of-two cycles:", tot,
          "  of which lie in a K4-subdivision:", good)

    print()
    print("=== hard-case ladders: does every power-of-two cycle lie")
    print("    in a K4-subdivision, and does (K4) hold?")
    for n in (16, 18, 20):
        lads = gen(n, 120, rng)
        badcyc = 0
        totcyc = 0
        k4fail = 0
        for ch in lads:
            adj = build(ch, n)
            ok = False
            for r in range(1, len(ch) + 1):
                for T in combinations(ch, r):
                    E = cycle_edges(list(T), n)
                    if E is None or len(E) not in POW:
                        continue
                    totcyc += 1
                    if interleaving_bridges(adj, cyclic_order(E), n):
                        ok = True
                    else:
                        badcyc += 1
            if not ok:
                k4fail += 1
        print("2N =", n, " ladders", len(lads),
              " power-of-two cycles", totcyc,
              " not in any K4-subdivision:", badcyc,
              " (K4) fails on:", k4fail)

main()
