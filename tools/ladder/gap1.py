# -*- coding: utf-8 -*-
"""GAP 1: residual one-chord cycles that (XC) does not settle.
Extract them and read off what DOES settle them."""

import sys
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

def closed_adj(ch, n):
    adj = {v: [] for v in range(1, n + 1)}
    for j in range(1, n):
        adj[j].append(j + 1)
        adj[j + 1].append(j)
    for a, b in ch:
        adj[a].append(b)
        adj[b].append(a)
    adj[1].append(n)
    adj[n].append(1)
    return adj

def state(ch, n, a, b):
    mu = {}
    for x, y in ch:
        mu[x] = y
        mu[y] = x

    def ordf(m):
        return (m - (b + 1)) if m > b else (n - b - 1 + m)

    X = [p for p in range(a + 1, b) if not (a < mu[p] < b)]
    if len(X) != b - a - 1:
        return None
    f = {p: ordf(mu[p]) for p in X}
    for q, p in combinations(X, 2):
        if f[q] < f[p]:
            return None
    out_ch = []
    for (x, y) in ch:
        if not (a <= x <= b) and not (a <= y <= b):
            out_ch.append(tuple(sorted((ordf(x), ordf(y)))))
    xc = False
    for q, p in combinations(X, 2):
        u, v = f[p], f[q]
        for (g, d) in out_ch:
            for (al, be) in out_ch:
                if (g, d) != (al, be) and g < u <= al < v <= d < be:
                    xc = True
    return X, f, out_ch, xc, ordf, mu

def main():
    sys.setrecursionlimit(30000)
    n = 14
    shown = 0
    for ch in enum_all(n):
        for (a, b) in ch:
            st = state(ch, n, a, b)
            if st is None:
                continue
            X, f, out_ch, xc, ordf, mu = st
            if xc:
                continue
            shown += 1
            if shown > 6:
                continue
            adj = closed_adj(ch, n)
            E = cycle_edges([(a, b)], n)
            order = cyc_order(E)
            ears = all_ears(adj, order, E)
            wins = []
            for e1, e2 in combinations(ears, 2):
                if seven_from(order, e1, e2) or seven_from(order, e2,
                                                           e1):
                    wins.append((e1, e2))
            print("ladder", ch, " Z = chord", (a, b))
            print("   crossers", [(p, mu[p], f[p]) for p in X],
                  " outside chords(ord)", sorted(out_ch))
            print("   #winning pairs", len(wins))
            for w in wins[:3]:
                print("      {%s,%s} via ord %s   |   {%s,%s} via ord"
                      " %s" % (w[0][0], w[0][1],
                               [ordf(z) for z in w[0][2]],
                               w[1][0], w[1][1],
                               [ordf(z) for z in w[1][2]]))
    print("total residual-without-XC at 2N =", n, ":", shown)

main()
