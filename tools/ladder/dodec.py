# -*- coding: utf-8 -*-
"""THE DECISIVE TEST.

Z (a residual one-chord cycle) has exactly ONE bridge, so in ANY planar
embedding of the closed ladder G, Z bounds a face; and a cycle bounding
a face has no two disjoint interleaving ears.  Hence

    G planar  ==>  (EAR*) FAILS for every residual one-chord cycle.

The smallest cubic planar graph of girth 5 is the dodecahedron, on 20
vertices - which is exactly why 2N <= 18 saw no planar ladder.  So the
dodecahedron is the predicted counterexample.  Test it.
"""

import sys
from itertools import combinations
from nlutil import ov, closed_adj, girth, residual, planar_hamiltonian


def lcf(shifts, reps):
    n = len(shifts) * reps
    s = shifts * reps
    ch = set()
    for i in range(n):
        j = (i + s[i]) % n
        ch.add(tuple(sorted((i + 1, j + 1))))
    return n, sorted(ch)


def two_link(adjH, w, y, x, z):
    stack = [(w, frozenset([w]))]
    while stack:
        v, seen = stack.pop()
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
            if u not in seen:
                stack.append((u, seen | frozenset([u])))
    return False


def ear_star(m, outch, T):
    adjH = dict((v, []) for v in range(m + 1))
    for j in range(m):
        adjH[j].append(j + 1)
        adjH[j + 1].append(j)
    for x, y in outch:
        adjH[x].append(y)
        adjH[y].append(x)
    wins = []
    for w, x, y, z in combinations(T, 4):
        if two_link(adjH, w, y, x, z):
            wins.append((w, x, y, z))
    return wins


def report(name, n, ch):
    adj = closed_adj(ch, n)
    deg = set(len(adj[v]) for v in adj)
    print("=" * 64)
    print(name, " n =", n, " degrees", sorted(deg),
          " girth", girth(adj), " chords", len(ch))
    print("   host cycle 1-2-...-%d-1 ; matching %s" % (n, ch))
    print("   PLANAR (interlacement bipartite):", planar_hamiltonian(ch))
    nres = 0
    for (a, b) in ch:
        R = residual(ch, n, a, b)
        if R is None:
            continue
        nres += 1
        m, f, X, ports, outch = R
        T = sorted(ports)
        wins = ear_star(m, outch, T)
        print("   Z = chord (%d,%d):  m=%d  F=%s  T=%s" %
              (a, b, m, sorted(f.values()), T))
        print("      outside chords (ord):", outch)
        S = [c for c in outch if any(c[0] < t < c[1] for t in T)]
        print("      straddlers:", S, " a crossing pair:",
              any(ov(c, d) for c, d in combinations(S, 2)))
        print("      (EAR*) interleaved linkages found:", len(wins),
              "->", "HOLDS" if wins else "*** FAILS ***")
    print("   residual one-chord cycles:", nres)


def main():
    sys.setrecursionlimit(40000)
    n, ch = lcf([10, 7, 4, -4, -7, 10, -4, 7, -7, 4], 2)
    report("DODECAHEDRON  LCF[10,7,4,-4,-7,10,-4,7,-7,4]^2", n, ch)
    # a non-planar control of the same size: the Desargues graph
    n2, ch2 = lcf([5, -5, 9, -9], 5)
    report("DESARGUES  LCF[5,-5,9,-9]^5", n2, ch2)


main()
