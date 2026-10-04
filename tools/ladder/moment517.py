# -*- coding: utf-8 -*-
"""517: the c_8 identity, and the first-moment law.

(A)  THE IDENTITY.  For a cubic graph of girth >= 5, 8 < 2g, so the
non-backtracking trace counts each 8-cycle 16 times and Ihara-Bass
gives, with p_8(mu) = mu^8 - 16 mu^6 + 80 mu^4 - 128 mu^2 + 32,

    16 c_8  =  tr A^8 - 16 tr A^6 + 80 tr A^4 - 128 tr A^2 + 33 n .

Girth >= 5 forces tr A^2 = 3n and tr A^4 = 15n, so

    16 c_8  =  tr A^8 - 16 tr A^6 + 849 n
            =  849 n - 5103 - SUM_{i>=2} mu_i^6 (16 - mu_i^2) .

Verified here against brute-force cycle counting.

(B)  THE FIRST MOMENT.  In a uniform random cubic graph the number of
k-cycles is asymptotically Poisson with mean (d-1)^k / (2k) = 2^k/2k.
At k = 8 the mean is 16 - AN ABSOLUTE CONSTANT - so c_8 = O(1) while
every term of the identity is Theta(n).  Measured here on girth-5
closed ladders.
"""

import sys
import random


def build(chords, n):
    adj = dict((v, []) for v in range(1, n + 1))
    for i in range(1, n + 1):
        j = i % n + 1
        adj[i].append(j)
        adj[j].append(i)
    for a, b in chords:
        adj[a].append(b)
        adj[b].append(a)
    return adj


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
    for L in range(3, 12):
        if count_len(adj, L):
            return L
    return 99


def traces(adj):
    """tr A^k for k = 2,4,6,8 by repeated squaring on the integer
    adjacency matrix - no eigenvalues, no floating point."""
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
            Xi = X[i]
            Zi = Z[i]
            for k in range(n):
                x = Xi[k]
                if x:
                    Yk = Y[k]
                    for j in range(n):
                        Zi[j] += x * Yk[j]
        return Z

    A2 = mul(A, A)
    A4 = mul(A2, A2)
    A6 = mul(A4, A2)
    A8 = mul(A4, A4)
    tr = lambda M: sum(M[i][i] for i in range(n))
    return tr(A2), tr(A4), tr(A6), tr(A8), n


def check(name, adj):
    g = girth(adj)
    t2, t4, t6, t8, n = traces(adj)
    c8 = count_len(adj, 8)
    lhs = 16 * c8
    rhs = t8 - 16 * t6 + 80 * t4 - 128 * t2 + 33 * n
    rhs2 = t8 - 16 * t6 + 849 * n
    print("%-20s n=%3d g=%d  trA2=%6d(3n=%5d) trA4=%7d(15n=%6d)  "
          "c_8=%4d  16c_8=%6d  identity=%6d %s  short form=%6d %s"
          % (name, n, g, t2, 3 * n, t4, 15 * n, c8, lhs, rhs,
             "OK" if rhs == lhs else "MISMATCH", rhs2,
             "OK" if rhs2 == lhs else "MISMATCH"))
    sys.stdout.flush()
    return lhs == rhs


def rand_ladder(n, tries=400):
    """a uniformly-ish random girth-5 closed ladder, or None."""
    for _ in range(tries):
        pts = list(range(1, n + 1))
        random.shuffle(pts)
        ch = [tuple(sorted((pts[2 * i], pts[2 * i + 1])))
              for i in range(n // 2)]
        adj = build(ch, n)
        if all(len(adj[v]) == 3 for v in adj) and girth(adj) >= 5:
            return ch, adj
    return None


def c8_law(n, samples):
    random.seed(20260809 + n)
    hist = {}
    tot = 0
    got = 0
    while got < samples:
        r = rand_ladder(n)
        if r is None:
            break
        ch, adj = r
        c = count_len(adj, 8)
        hist[c] = hist.get(c, 0) + 1
        tot += c
        got += 1
    print("n=%3d  samples=%5d  mean c_8 = %6.2f  (Poisson mean 16)  "
          "min=%d  histogram %s"
          % (n, got, float(tot) / max(got, 1), min(hist),
             sorted(hist.items())[:14]))
    sys.stdout.flush()


if __name__ == '__main__':
    from cage516 import lcf, petersen, coxeter
    ok = True
    ok &= check("Petersen", petersen())
    ok &= check("dodecahedron", lcf([10, 7, 4, -4, -7, 10,
                                     -4, 7, -7, 4], 2))
    ok &= check("Heawood", lcf([5, -5], 7))
    ok &= check("Mobius-Kantor", lcf([5, -5], 8))
    ok &= check("Pappus", lcf([5, 7, -7, 7, -7, -5], 3))
    ok &= check("Desargues", lcf([5, -5, 9, -9], 5))
    ok &= check("McGee", lcf([12, 7, -7], 8))
    ok &= check("Coxeter", coxeter())
    print("identity holds on all:", ok)
    print()
    for nn in (20, 26, 34, 44):
        c8_law(nn, 400)
