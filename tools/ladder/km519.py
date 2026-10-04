# -*- coding: utf-8 -*-
"""519: the Kesten-McKay obstruction - why NO moment inequality can
prove (E_3).

The spectral measure of the infinite 3-regular tree (Kesten-McKay for
d = 3) is a probability measure supported on [-2 sqrt 2, 2 sqrt 2]
contained in [-3, 3], whose moments are the closed-walk counts of the
tree:  W_0, W_2, W_4, W_6, W_8 = 1, 3, 15, 87, 543.

517.1 says 16 c_8 = tr A^8 - 16 tr A^6 + 849 n for cubic girth >= 5,
equivalently  16 c_8 = SUM_i mu_i^8 - 543 n - 192 c_6 .  A measure
with the tree's moments therefore has "c_8 = 0".  Hence the infimum
of SUM mu^8 over probability measures on [-3,3] with the prescribed
2nd, 4th and 6th moments is exactly 543 n, ATTAINED, and no
inequality using only those moments can conclude c_8 > 0.

Verified here: the KM moments by numerical integration against the
density, and independently as tree closed-walk counts by recursion.
"""

import sys
from fractions import Fraction


def tree_walks(d, kmax):
    """W_k = closed walks of length k from the root of the infinite
    d-regular tree, by the standard transfer recursion on distance
    from the root."""
    # state: distance from root; step out with (d-1) choices except
    # from the root where it is d; step in with 1 choice.
    W = []
    for k in range(kmax + 1):
        # vector over distances 0..k
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


def km_moments(kmax, steps=4000000):
    """moments of the Kesten-McKay density for d = 3, by Simpson.
    rho(x) = d sqrt(4(d-1) - x^2) / (2 pi (d^2 - x^2)) on
    |x| <= 2 sqrt(d-1)."""
    import math
    d = 3.0
    R = 2.0 * math.sqrt(d - 1.0)
    out = [0.0] * (kmax + 1)
    h = 2.0 * R / steps
    for i in range(steps + 1):
        x = -R + i * h
        s = 4.0 * (d - 1.0) - x * x
        if s < 0:
            s = 0.0
        rho = d * math.sqrt(s) / (2.0 * math.pi * (d * d - x * x))
        w = 1.0 if i in (0, steps) else (4.0 if i % 2 else 2.0)
        p = 1.0
        for k in range(kmax + 1):
            out[k] += w * rho * p
            p *= x
    return [v * h / 3.0 for v in out]


if __name__ == '__main__':
    W = tree_walks(3, 10)
    print("closed walks from the root of the 3-regular tree:")
    print("  W_k for k=0..10 :", W)
    print("  even ones       :", [W[k] for k in range(0, 11, 2)])
    m = km_moments(10)
    print("Kesten-McKay moments by numerical integration (d=3):")
    for k in range(0, 11, 2):
        print("   moment %2d = %12.4f   tree walk count = %d  %s"
              % (k, m[k], W[k],
                 "OK" if abs(m[k] - W[k]) < 1e-3 * max(1, W[k])
                 else "MISMATCH"))
    print()
    print("So a measure on [-3,3] with 2nd, 4th, 6th moments "
          "3, 15, 87 can have")
    print("8th moment exactly 543, i.e. c_8 = 0.  The three-moment "
          "SOS bound of 517.1")
    print("gives 519n; the true infimum is 543n; the target needs "
          "> 543n.  The gap of")
    print("24n is NOT slack to be recovered - the optimum is the "
          "threshold itself.")
    sys.stdout.flush()
