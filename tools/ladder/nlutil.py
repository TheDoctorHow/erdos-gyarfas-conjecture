# -*- coding: utf-8 -*-
"""shared helpers for the (NONLAM) work - closed-ladder model."""

from itertools import combinations
from collections import deque


def ov(c, d):
    a, b = c
    x, y = d
    return (a < x < b < y) or (x < a < y < b)


def soft_pair(c, d):
    if not ov(c, d):
        return True
    return abs(c[0] - d[0]) + abs(c[1] - d[1]) >= 3


def enum_matchings(n, cap=None, rng=None):
    out, ch = [], []
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
            out.append(list(ch))
            return
        used[p] = True
        cand = list(range(p + 1, n + 1))
        if rng is not None:
            rng.shuffle(cand)
        for q in cand:
            if stop:
                break
            if used[q]:
                continue
            s = q - p
            if s < 4 or s > n - 4:
                continue
            c = (p, q)
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
    adj = dict((v, []) for v in range(1, n + 1))
    for j in range(1, n):
        adj[j].append(j + 1)
        adj[j + 1].append(j)
    adj[1].append(n)
    adj[n].append(1)
    for x, y in ch:
        adj[x].append(y)
        adj[y].append(x)
    return adj


def girth(adj):
    best = 10 ** 9
    for s in adj:
        dist = {s: 0}
        par = {s: None}
        dq = deque([s])
        while dq:
            v = dq.popleft()
            if 2 * dist[v] >= best:
                continue
            for w in adj[v]:
                if w == par[v]:
                    continue
                if w in dist:
                    c = dist[v] + dist[w] + 1
                    if c < best:
                        best = c
                else:
                    dist[w] = dist[v] + 1
                    par[w] = v
                    dq.append(w)
    return best


def residual(ch, n, a, b):
    mu = {}
    for x, y in ch:
        mu[x] = y
        mu[y] = x
    m = n - (b - a + 1) - 1
    if m < 1:
        return None

    def o(t):
        return t - (b + 1) if t > b else n - b - 1 + t

    X = [p for p in range(a + 1, b) if not (a < mu[p] < b)]
    if len(X) != b - a - 1:
        return None
    f = dict((p, o(mu[p])) for p in X)
    for q, p in combinations(X, 2):
        if f[q] < f[p]:
            return None
    ports = {0: b, m: a}
    for p in X:
        if f[p] in ports:
            return None
        ports[f[p]] = p
    outch = sorted(tuple(sorted((o(x), o(y)))) for x, y in ch
                   if not (a <= x <= b) and not (a <= y <= b))
    return m, f, X, ports, outch


def planar_hamiltonian(ch):
    """G = host cycle + chords is planar iff the interlacement graph
    of the chords is bipartite (classical criterion)."""
    k = len(ch)
    adj = dict((i, []) for i in range(k))
    for i, j in combinations(range(k), 2):
        if ov(ch[i], ch[j]):
            adj[i].append(j)
            adj[j].append(i)
    col = {}
    for s in range(k):
        if s in col:
            continue
        col[s] = 0
        dq = deque([s])
        while dq:
            v = dq.popleft()
            for w in adj[v]:
                if w in col:
                    if col[w] == col[v]:
                        return False
                else:
                    col[w] = 1 - col[v]
                    dq.append(w)
    return True
