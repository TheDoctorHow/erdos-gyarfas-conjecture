# -*- coding: utf-8 -*-
"""(NONLAM) - the outside chord system of a residual one-chord cycle.

CLOSED ladder model: host cycle 1-2-...-n-1 plus a perfect matching M,
girth of the CLOSED graph >= 5 (this is the faithful J7-rung model of
510.4; the open-ladder enumerations of 511-513 did NOT impose the
closed girth, so this run also re-measures the residual counts).

For Z = the cycle (chord (a,b) + host path a..b), RESIDUAL means
  - every interior vertex of Z is a crosser (mu(p) outside [a,b]);
  - f(p) = ord(mu(p)) is strictly decreasing in p.
Outside path: ord 0 = b+1, ..., ord m = a-1 (m = n-(b-a+1)-1).
Terminals (ports):  0 -> b,  m -> a,  f(p) -> p.
Outside chords: matching edges with both ends outside [a,b].

Tests
  P0  no two terminals at consecutive ord positions
  T1  the outside chord system is NON-LAMINAR (a crossing pair exists)
  T2  (GEN) a crossing pair with terminals in valid zones -> linkage
  T5  two port-straddling chords cross
  T6  the crossing graph is non-bipartite
  T4  Z-planarity: chords 2-colour D/E, D laminar, E laminar, every
      E-chord port-free.  Z-planar  <=>  (EAR*) fails for Z.
  T3  full 2-linkage brute force = (EAR*) restricted to Z
"""

import sys
import random
from itertools import combinations
from collections import deque, Counter


def ov(c, d):
    a, b = c
    x, y = d
    return (a < x < b < y) or (x < a < y < b)


def soft_pair(c, d):
    if not ov(c, d):
        return True
    return abs(c[0] - d[0]) + abs(c[1] - d[1]) >= 3


def enum_matchings(n, cap=None, rng=None):
    """perfect matchings with 4 <= span <= n-4 and soft_pair"""
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
            return ('CLASH', p)
        ports[f[p]] = p
    outch = sorted(tuple(sorted((o(x), o(y)))) for x, y in ch
                   if not (a <= x <= b) and not (a <= y <= b))
    return m, f, X, ports, outch


def zones(T, g, al, de, be):
    """T sorted terminal list; returns counts in the five zones"""
    z = [0, 0, 0, 0, 0]
    for t in T:
        if t <= g:
            z[0] += 1
        elif t < al:
            z[1] += 1
        elif t < de:
            z[2] += 1
        elif t < be:
            z[3] += 1
        else:
            z[4] += 1
    return z


def gen_criterion(T, outch):
    """(GEN): crossing chords + terminals in valid zones -> linkage."""
    for c1, c2 in combinations(outch, 2):
        for (g, de), (al, be) in ((c1, c2), (c2, c1)):
            if not (g < al < de < be):
                continue
            z = zones(T, g, al, de, be)
            for zw in (0, 1):
                for zx in (1, 2):
                    for zy in (2, 3):
                        for zz in (3, 4):
                            need = Counter([zw, zx, zy, zz])
                            if all(z[k] >= need[k] for k in need):
                                return (g, de, al, be, zw, zx, zy, zz)
    return None


def straddlers(T, outch):
    out = []
    for c in outch:
        if any(c[0] < t < c[1] for t in T):
            out.append(c)
    return out


def cross_graph(outch):
    E = []
    for i, j in combinations(range(len(outch)), 2):
        if ov(outch[i], outch[j]):
            E.append((i, j))
    return E


def zplanar(outch, T):
    """True iff chords 2-colour with all straddlers in class D."""
    k = len(outch)
    E = cross_graph(outch)
    adj = dict((i, []) for i in range(k))
    for i, j in E:
        adj[i].append(j)
        adj[j].append(i)
    st = set(i for i in range(k)
             if any(outch[i][0] < t < outch[i][1] for t in T))
    col = {}
    for s in range(k):
        if s in col:
            continue
        col[s] = 0
        dq = deque([s])
        comp = [s]
        while dq:
            v = dq.popleft()
            for w in adj[v]:
                if w in col:
                    if col[w] == col[v]:
                        return False
                else:
                    col[w] = 1 - col[v]
                    comp.append(w)
                    dq.append(w)
        cols = set(col[v] for v in comp if v in st)
        if len(cols) > 1:
            return False
    return True


def two_link(adjH, w, y, x, z):
    stack = [(w, [w], frozenset([w]))]
    while stack:
        v, pth, seen = stack.pop()
        if v == y:
            if x in seen or z in seen:
                continue
            vis = set([x])
            dq = [x]
            while dq:
                p = dq.pop()
                for u in adjH[p]:
                    if u not in seen and u not in vis:
                        vis.add(u)
                        dq.append(u)
            if z in vis:
                return True
            continue
        for u in adjH[v]:
            if u in seen:
                continue
            stack.append((u, pth + [u], seen | frozenset([u])))
    return False


def ear_star(m, outch, T):
    adjH = dict((v, []) for v in range(m + 1))
    for j in range(m):
        adjH[j].append(j + 1)
        adjH[j + 1].append(j)
    for x, y in outch:
        adjH[x].append(y)
        adjH[y].append(x)
    for w, x, y, z in combinations(T, 4):
        if two_link(adjH, w, y, x, z):
            return True
    return False


def run(n, lads, do_link):
    st = Counter()
    bad = []
    for ch in lads:
        adj = closed_adj(ch, n)
        if girth(adj) < 5:
            st['skipped-closed-girth'] += 1
            continue
        st['ladders'] += 1
        for (a, b) in ch:
            R = residual(ch, n, a, b)
            if R is None:
                continue
            if R[0] == 'CLASH':
                st['PORT-CLASH!'] += 1
                continue
            m, f, X, ports, outch = R
            st['residual'] += 1
            T = sorted(ports)
            if any(T[i + 1] - T[i] == 1 for i in range(len(T) - 1)):
                st['P0-FAIL!'] += 1
                bad.append(('P0', ch, (a, b)))
            E = cross_graph(outch)
            if E:
                st['T1-nonlaminar'] += 1
            else:
                st['T1-FAIL-laminar!'] += 1
                bad.append(('T1', ch, (a, b)))
            if gen_criterion(T, outch) is not None:
                st['T2-GEN'] += 1
            else:
                bad.append(('T2', ch, (a, b), m, outch, T))
            S = straddlers(T, outch)
            if any(ov(c, d) for c, d in combinations(S, 2)):
                st['T5-straddlers-cross'] += 1
            if zplanar(outch, T):
                st['T4-ZPLANAR!'] += 1
                bad.append(('T4', ch, (a, b)))
            else:
                st['T4-not-Zplanar'] += 1
            if do_link:
                if ear_star(m, outch, T):
                    st['T3-EARstar'] += 1
                else:
                    st['T3-FAIL!'] += 1
                    bad.append(('T3', ch, (a, b)))
    return st, bad


def main():
    sys.setrecursionlimit(40000)
    rng = random.Random(271828)
    jobs = [(14, None, True), (16, None, True), (18, 4000, True),
            (20, 1500, False), (22, 800, False), (24, 400, False)]
    keep = {}
    for n, cap, dl in jobs:
        lads = enum_matchings(n, cap, None if cap is None else rng)
        st, bad = run(n, lads, dl)
        print("=" * 62)
        print("2N =", n, " matchings enumerated:", len(lads))
        for k in sorted(st):
            print("   %-26s %d" % (k, st[k]))
        keep[n] = bad
        for tag in ('P0', 'T1', 'T4', 'T3'):
            ex = [b for b in bad if b[0] == tag]
            if ex:
                print("   !!! %s failures: %d  first: %s"
                      % (tag, len(ex), ex[0][1:3]))
        t2 = [b for b in bad if b[0] == 'T2']
        if t2:
            print("   T2 (GEN) misses:", len(t2))
            for b in t2[:4]:
                print("      ladder", b[1], " Z", b[2],
                      " m", b[3], " outch", b[4], " T", b[5])
    return keep


main()
