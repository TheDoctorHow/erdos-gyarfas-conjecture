# -*- coding: utf-8 -*-
"""520: the band-3 program as a MOMENT-FEASIBILITY test, one moment
at a time, and how far up n it reaches.

519 asked "how small can the 8th moment be".  The sharper question is
FEASIBILITY: assume c_8 = 0 and ask whether the resulting moment
sequence belongs to a measure on [0,9] at all.  That is the truncated
Hausdorff moment problem, and it is a semidefinite feasibility test -
convex, exact, and it uses every moment supplied rather than only
three.

SET-UP, girth >= 7 and c_8 = 0 assumed for contradiction.  With
y = mu^2 over the n-1 non-Perron eigenvalues,

  M0 = n-1
  M1 = 3n-9            (tr A^2 - 9)
  M2 = 15n-81          (tr A^4 - 81)
  M3 = 87n-729         (tr A^6 - 729,   c_5 = c_6 = 0)
  M4 = 543n-6561       (tr A^8 - 6561,  c_8 = 0 ASSUMED)
  M5 = tr A^10 - 59049 = 3543n + 20 c_10 - 59049  >=  3543n - 59049
  M6 = tr A^12 - 3^12  >=  W_12 n - 531441

CONDITIONS.  A sequence M0..M_{2k+1} comes from a measure on [0,b]
iff the two Hankel forms

    H1 = ( M_{i+j+1} )_{i,j<=k}          (the measure y.dnu)
    H2 = ( b M_{i+j} - M_{i+j+1} )_{i,j<=k}   (the measure (b-y).dnu)

are both positive semidefinite.  H1 gives a LOWER bound on the next
moment, H2 an UPPER bound.  If the lower bound forced by H1 EXCEEDS
what the graph allows, the hypothesis c_8 = 0 is contradicted.
"""

import sys
import numpy as np


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


W = tree_walks(3, 16)
EVEN = [W[2 * j] for j in range(0, 9)]     # 1,3,15,87,543,3543,...


def moments(n, upto):
    """M_j for j = 0..upto, in y = mu^2, over the non-Perron part,
    ASSUMING c_8 = 0 and girth >= 7.  Entries beyond j = 4 are the
    graph's MINIMUM (all higher cycle counts set to zero)."""
    return [EVEN[j] * n - 9 ** j for j in range(upto + 1)]


def min_next(M):
    """least value of the next moment consistent with H1 >= 0, where
    M = [M0..M_{2k}] is given and the unknown is M_{2k+1}.  det H1 is
    linear in it."""
    k = (len(M) - 1) // 2
    # H1 is (k+1)x(k+1) with entries M_{i+j+1}; the unknown sits at
    # (k,k) = M_{2k+1}
    def build(x):
        H = np.zeros((k + 1, k + 1))
        for i in range(k + 1):
            for j in range(k + 1):
                idx = i + j + 1
                H[i, j] = x if idx == 2 * k + 1 else M[idx]
        return H
    d0 = np.linalg.det(build(0.0))
    d1 = np.linalg.det(build(1.0))
    slope = d1 - d0
    if abs(slope) < 1e-9:
        return None
    return -d0 / slope


def report():
    print("girth >= 7, hypothesis c_8 = 0.  Hankel feasibility:")
    print("the least next moment H1 forces, against the least the "
          "graph allows.")
    print()
    for upto, nxt in ((4, 5), (6, 7)):
        lab = {5: "M5 (tr A^10)", 7: "M6/M7"}[nxt if nxt == 5 else 7]
        print("  using M0..M%d, testing M%d = %s" % (upto, nxt, lab))
        print("     n     Hankel forces M%d >=      graph allows "
              "M%d >=     contradiction?" % (nxt, nxt))
        for n in (24, 26, 28, 30, 32, 34, 36, 40, 50, 80, 200):
            M = [float(x) for x in moments(n, upto)]
            lo = min_next(M)
            allow = float(EVEN[nxt] * n - 9 ** nxt)
            if lo is None:
                print("    %4d        (degenerate)" % n)
                continue
            print("    %4d   %20.1f  %20.1f      %s"
                  % (n, lo, allow, "YES" if lo > allow else "no"))
        print()
    sys.stdout.flush()


def sanity():
    print("tree even-moment sequence W_0,W_2,...,W_16 :", EVEN)
    print("Kesten-McKay check: a measure with ALL these moments "
          "exists (the tree),")
    print("so the test can only ever fire through the Perron atom, "
          "at finite n.")
    print()
    # positive control: a genuine graph must pass every test
    print("positive control - the same test on the McGee graph "
          "(n=24, girth 7, c_8=34):")
    from cage516 import lcf
    from moment517 import traces, count_len
    g = lcf([12, 7, -7], 8)
    t2, t4, t6, t8, n = traces(g)
    print("   n=%d  trA2=%d trA4=%d trA6=%d trA8=%d  c_8=%d  "
          "543n=%d" % (n, t2, t4, t6, t8, count_len(g, 8), 543 * n))
    print("   tr A^8 - 543n = %d = 16 c_8 = %d  %s"
          % (t8 - 543 * n, 16 * count_len(g, 8),
             "OK" if t8 - 543 * n == 16 * count_len(g, 8) else "BAD"))
    print()
    sys.stdout.flush()


if __name__ == '__main__':
    sanity()
    report()
