# -*- coding: utf-8 -*-
"""521: THE SPECTRAL CEILING.  How far the spectrum route reaches, at
every girth and every target, and why it reaches no further.

For cubic G, tr B^m = SUM_i p_m(mu_i) + n for even m, where p_m is
the Chebyshev-like polynomial with p_m(mu) = lam^m + lam'^m,
lam lam' = 2, lam + lam' = mu.  For g <= m < 2g this is 2 m c_m.

TWO FACTS ARE COMPUTED HERE.

(A) THE TREE IDENTITY.  INT p_m dKM = -1 for every even m >= 2, where
KM is the Kesten-McKay measure (= the spectral measure of the
3-regular tree).  So a graph whose spectral measure were exactly KM
would have SUM p_m = -n and hence tr B^m = 0: c_m = 0 EXACTLY, at
every m.  The tree sits on the boundary at every target, not just at
m = 8.

(B) THE BULK CAN REALLY GO DOWN.  The crude step in the trace test is
|p_m(mu)| <= 2 . 2^{m/2} on the bulk |mu| <= 2 sqrt 2, giving
2 m c_m >= 2^m + 1 + n - 2 n 2^{m/2}.  Is that lossy?  Minimise
INT p_m dnu over measures on [-3,3] matching the tree moments below
the girth.  If the minimum is close to -2 . 2^{m/2}, the crude bound
is essentially the truth and the test's reach 2^{m/2} > 2n is the
CEILING of the whole spectral route.
"""

import sys
import math
from itertools import combinations


def pm_coeffs(m):
    """coefficients of p_m as a polynomial in mu, low order first."""
    a = [2]            # p_0 = 2
    b = [0, 1]         # p_1 = mu
    if m == 0:
        return a
    if m == 1:
        return b
    for _ in range(m - 1):
        # p_k = mu p_{k-1} - 2 p_{k-2}
        c = [0] * (len(b) + 1)
        for i, v in enumerate(b):
            c[i + 1] += v
        for i, v in enumerate(a):
            c[i] -= 2 * v
        a, b = b, c
    return b


def tree_walks(d, kmax):
    W = []
    for k in range(kmax + 1):
        vec = [0] * (kmax + 2)
        vec[0] = 1
        for _ in range(k):
            nxt = [0] * (kmax + 2)
            for j, val in enumerate(vec):
                if not val:
                    continue
                out = d if j == 0 else d - 1
                nxt[j + 1] += val * out
                if j > 0:
                    nxt[j - 1] += val
            vec = nxt
        W.append(vec[0])
    return W


def part_a(mmax=20):
    W = tree_walks(3, mmax)
    print("(A) the tree identity  INT p_m dKM = -1")
    ok = True
    for m in range(2, mmax + 1, 2):
        c = pm_coeffs(m)
        val = sum(c[j] * W[j] for j in range(len(c)))
        if val != -1:
            ok = False
        print("    m=%2d   INT p_m dKM = %4d   %s"
              % (m, val, "OK" if val == -1 else "*** NOT -1"))
    print("    identity holds for every even m tested:", ok)
    sys.stdout.flush()


def evalp(c, x):
    s = 0.0
    for v in reversed(c):
        s = s * x + v
    return s


def part_b(g, m, natoms=4, grid=260):
    """min INT p_m dnu over probability measures on [-3,3] with the
    tree moments up to order g-1.  Searched over measures with at most
    `natoms` atoms on a grid, weights by exact solve where possible;
    reported against the crude bound -2 . 2^{m/2}."""
    W = tree_walks(3, m + 2)
    # even moment constraints of order 2,4,...,<= g-1  (odd ones are
    # automatic for symmetric measures, which we restrict to)
    orders = [j for j in range(2, g) if j % 2 == 0]
    xs = [-3.0 + 6.0 * i / grid for i in range(grid + 1)]
    ys = sorted(set(round(x * x, 9) for x in xs if x >= 0))
    c = pm_coeffs(m)

    def pey(y):                      # p_m is even in mu
        return evalp(c, math.sqrt(y))

    best = None
    k = len(orders) + 1              # constraints incl. mass
    for combo in combinations(ys, min(natoms, k)):
        A = [[1.0] * len(combo)]
        rhs = [1.0]
        for o in orders:
            A.append([y ** (o // 2) for y in combo])
            rhs.append(float(W[o]))
        if len(A) != len(combo):
            continue
        n = len(combo)
        M = [row[:] + [rhs[i]] for i, row in enumerate(A)]
        for i in range(n):           # Gaussian elimination
            p = max(range(i, n), key=lambda r: abs(M[r][i]))
            if abs(M[p][i]) < 1e-12:
                break
            M[i], M[p] = M[p], M[i]
            for r in range(n):
                if r == i:
                    continue
                f = M[r][i] / M[i][i]
                for cc in range(i, n + 1):
                    M[r][cc] -= f * M[i][cc]
        else:
            w = [M[i][n] / M[i][i] for i in range(n)]
            if min(w) < -1e-9:
                continue
            val = sum(w[i] * pey(combo[i]) for i in range(n))
            if best is None or val < best[0]:
                best = (val, combo, w)
    crude = -2.0 * 2 ** (m / 2.0)
    if best is None:
        print("  g=%2d m=%2d : no feasible atom set on the grid" % (g, m))
        return
    print("  g=%2d  m=%2d   LP min INT p_m = %12.1f   crude bound "
          "-2.2^{m/2} = %12.1f   ratio %.3f" % (g, m, best[0], crude,
                                                best[0] / crude))
    sys.stdout.flush()


if __name__ == '__main__':
    part_a()
    print()
    print("(B) how lossy is the crude bulk bound?")
    for g, m in ((5, 8), (6, 8), (7, 8), (9, 16), (11, 16), (13, 16)):
        part_b(g, m)
