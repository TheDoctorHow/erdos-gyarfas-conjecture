# -*- coding: utf-8 -*-
"""523: the rank >= 2 correction above 2g, measured exactly, and
compared with the number of theta subgraphs.

452 computes c_m exactly for g <= m < 2g because tr B^m = 2 m c_m
there.  Above 2g,

    tr B^m  =  SUM_{l | m, l >= g} 2 l c_l        (rank 1, exact)
             +  R(m)                              (rank >= 2)

and 452.5 measured R(m)/2^m in [0.66, 2.53] - "of the same order as
the main term".  516.6(i) names controlling R(m) as the archive's one
live analytic problem, and 520.3 proved that in a COUNTEREXAMPLE the
rank-1 sum vanishes at m = 2^k.  So everything now turns on R(m).

THE QUESTION THIS MEASURES.  A rank >= 2 closed non-backtracking walk
of length m runs over a subgraph of cycle rank at least two.  For
m < 4g the smallest such subgraphs are thetas and dumbbells.  If
R(m) were bounded by (number of thetas of bounded size) x O(1), the
window could be widened from 2g to 4g.  Is it?

tr B^m is computed exactly from tr A^j by Ihara-Bass:
    tr B^m = SUM_i p_m(mu_i) + n  =  SUM_j [p_m]_j tr A^j  +  n ,
with p_m the Chebyshev-like polynomial of 521.  All integer
arithmetic.
"""

import sys
from itertools import combinations


def pm_coeffs(m):
    a, b = [2], [0, 1]
    if m == 0:
        return a
    for _ in range(m - 1):
        c = [0] * (len(b) + 1)
        for i, v in enumerate(b):
            c[i + 1] += v
        for i, v in enumerate(a):
            c[i] -= 2 * v
        a, b = b, c
    return b


def traces(adj, upto):
    vs = sorted(adj)
    idx = dict((v, i) for i, v in enumerate(vs))
    n = len(vs)
    A = [[0] * n for _ in range(n)]
    for v in adj:
        for w in adj[v]:
            A[idx[v]][idx[w]] += 1

    def mul(X, Y):
        Z = [[0] * n for _ in range(n)]
        for i in range(n):
            Xi, Zi = X[i], Z[i]
            for k in range(n):
                x = Xi[k]
                if x:
                    Yk = Y[k]
                    for j in range(n):
                        Zi[j] += x * Yk[j]
        return Z

    out = [None] * (upto + 1)
    P = [[1 if i == j else 0 for j in range(n)] for i in range(n)]
    out[0] = n
    for k in range(1, upto + 1):
        P = mul(P, A)
        out[k] = sum(P[i][i] for i in range(n))
    return out


def trB(adj, m, tr):
    """tr B^m = SUM_i p_m(mu_i) + (n/2)(1^m + (-1)^m); the two blocks
    of n/2 eigenvalues +1 and -1 CANCEL for odd m and give n for even
    m.  Getting this wrong adds a spurious n at odd m."""
    c = pm_coeffs(m)
    s = sum(c[j] * tr[j] for j in range(len(c)))
    return s + (len(adj) if m % 2 == 0 else 0)


def count_len(adj, L):
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


def girth(adj):
    for L in range(3, 20):
        if count_len(adj, L):
            return L
    return 99


def theta_count(adj, maxlen):
    """thetas with total edge count <= maxlen: pairs of vertices with
    three internally disjoint paths.  Counted by enumerating cycles
    up to maxlen and pairing those that meet in exactly one path."""
    from collections import defaultdict
    cyc = []
    vs = sorted(adj)
    rank = dict((v, i) for i, v in enumerate(vs))
    for s in vs:
        stack = [(s, (s,), 1 << s)]
        while stack:
            v, path, msk = stack.pop()
            for w in adj[v]:
                if w == s:
                    if len(path) >= 3 and rank[path[1]] < rank[path[-1]]:
                        cyc.append(path)
                elif (rank[w] > rank[s] and not (msk >> w & 1)
                      and len(path) < maxlen):
                    stack.append((w, path + (w,), msk | (1 << w)))
    es = []
    for p in cyc:
        k = len(p)
        es.append((set(p), set(frozenset((p[i], p[(i + 1) % k]))
                               for i in range(k))))
    tot = 0
    for (v1, e1), (v2, e2) in combinations(es, 2):
        sh = e1 & e2
        if not sh:
            continue
        # union has rank 2 and is a theta iff the shared part is a
        # single path and the union has <= maxlen edges
        if len(e1 | e2) > maxlen:
            continue
        deg = {}
        for e in sh:
            for x in e:
                deg[x] = deg.get(x, 0) + 1
        if any(d > 2 for d in deg.values()):
            continue
        tot += 1
    return tot


def report(name, adj):
    g = girth(adj)
    n = len(adj)
    hi = min(4 * g - 1, 22)
    tr = traces(adj, hi)
    print("%s  n=%d  girth=%d" % (name, n, g))
    print("    m   tr B^m         rank-1        R(m)=rank>=2     "
          "R(m)/2^m   2m c_m")
    for m in range(g, hi + 1):
        B = trB(adj, m, tr)
        r1 = 0
        for l in range(g, m + 1):
            if m % l == 0:
                r1 += 2 * l * count_len(adj, l)
        R = B - r1
        cm = count_len(adj, m)
        flag = "  <- window" if m < 2 * g else ""
        print("   %3d %12d %12d %14d %10.3f %10d%s"
              % (m, B, r1, R, float(R) / 2 ** m, 2 * m * cm, flag))
    th = theta_count(adj, min(4 * g - 1, 16))
    print("    thetas with at most %d edges: %d      "
          "R(2g)/thetas = %.2f"
          % (min(4 * g - 1, 16), th,
             float(trB(adj, 2 * g, tr) - sum(
                 2 * l * count_len(adj, l)
                 for l in range(g, 2 * g + 1) if (2 * g) % l == 0))
             / max(th, 1)))
    print()
    sys.stdout.flush()


if __name__ == '__main__':
    from cage516 import lcf, petersen
    report("Petersen", petersen())
    report("Heawood", lcf([5, -5], 7))
    report("dodecahedron", lcf([10, 7, 4, -4, -7, 10, -4, 7, -7, 4], 2))
