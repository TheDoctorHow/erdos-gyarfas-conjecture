# -*- coding: utf-8 -*-
"""516: LOCAL CAPTURE tested on the hard case - the cages.

  (LC) A cubic graph of girth g has a cycle of length equal to the
       UNIQUE power of two in [g, 2g-2].

Uniqueness: the interval has ratio 2 - 2/g < 2.  Existence: let 2^k be
the least power of two >= g; then 2^k < 2g, and 2^k != 2g-1 since
2g-1 is odd and 2^k is even for k >= 1, so 2^k <= 2g-2.

(LC) implies Erdos-Gyarfas for cubic graphs outright, and its subject
names no power of two, so it is (N9)-safe.  The (N16) history says it
should be false; test it where it is most likely to fail - on the
cages, the extremal graphs at each girth.

NAMES ARE NOT TRUSTED.  Every graph is built from an LCF string and
then its order, regularity and GIRTH are recomputed here; the
mathematical claim uses only the computed invariants.
"""

import sys


def lcf(shifts, reps):
    n = len(shifts) * reps
    s = shifts * reps
    adj = dict((i, []) for i in range(n))
    for i in range(n):
        for j in ((i + 1) % n, (i + s[i]) % n):
            if j not in adj[i]:
                adj[i].append(j)
            if i not in adj[j]:
                adj[j].append(i)
    return adj


def petersen():
    out = dict((v, []) for v in range(10))
    for i in range(5):
        for x, y in ((i, (i + 1) % 5), (i, i + 5),
                     (i + 5, (i + 2) % 5 + 5)):
            out[x].append(y)
            out[y].append(x)
    return out


def coxeter():
    """the Coxeter graph, 28 vertices, girth 7, NOT Hamiltonian.
    Standard construction: seven 'spokes' plus three heptagons
    Z_7 with steps 1, 2, 4."""
    adj = dict((v, []) for v in range(28))

    def add(u, v):
        if v not in adj[u]:
            adj[u].append(v)
            adj[v].append(u)

    for t, step in enumerate((1, 2, 4)):
        base = 7 * t
        for i in range(7):
            add(base + i, base + (i + step) % 7)
    for i in range(7):
        for t in range(3):
            add(21 + i, 7 * t + i)
    return adj


def girth(adj):
    best = 10 ** 9
    for s in adj:
        dist = {s: 0}
        par = {s: None}
        dq = [s]
        while dq:
            nxt = []
            for x in dq:
                if 2 * dist[x] >= best:
                    continue
                for y in adj[x]:
                    if y == par[x]:
                        continue
                    if y in dist:
                        g = dist[x] + dist[y] + 1
                        if g < best:
                            best = g
                    else:
                        dist[y] = dist[x] + 1
                        par[y] = x
                        nxt.append(y)
            dq = nxt
    return best


def count_len(adj, L):
    """cycles of length exactly L, depth-bounded DFS, each once."""
    tot = 0
    for s in sorted(adj):
        stack = [(s, 1, 1 << s)]
        while stack:
            v, d, msk = stack.pop()
            for w in adj[v]:
                if w == s:
                    if d == L:
                        tot += 1
                elif w > s and d < L and not (msk >> w & 1):
                    stack.append((w, d + 1, msk | (1 << w)))
    return tot // 2


def target(g):
    k = 4
    while k < g:
        k *= 2
    assert g <= k <= 2 * g - 2, (g, k)
    return k


def test(name, adj):
    n = len(adj)
    degs = sorted(set(len(adj[v]) for v in adj))
    g = girth(adj)
    t = target(g)
    c = count_len(adj, t)
    print("%-22s n=%3d deg=%s girth=%2d  target 2^k=%3d  cycles=%6d "
          "%s" % (name, n, degs, g, t, c, "OK" if c else "*** FAILS"))
    sys.stdout.flush()
    return g, t, c


if __name__ == '__main__':
    test("Petersen (3,5)-cage", petersen())
    test("dodecahedron", lcf([10, 7, 4, -4, -7, 10, -4, 7, -7, 4], 2))
    test("Heawood (3,6)-cage", lcf([5, -5], 7))
    test("Mobius-Kantor", lcf([5, -5], 8))
    test("Pappus", lcf([5, 7, -7, 7, -7, -5], 3))
    test("Desargues", lcf([5, -5, 9, -9], 5))
    test("Nauru", lcf([5, -9, 7, -7, 9, -5], 4))
    test("McGee (3,7)-cage", lcf([12, 7, -7], 8))
    test("Coxeter", coxeter())
    test("Tutte-Coxeter (3,8)", lcf([-13, -9, 7, -7, 9, 13], 5))
    test("Dyck", lcf([5, -5, 13, -13], 8))
