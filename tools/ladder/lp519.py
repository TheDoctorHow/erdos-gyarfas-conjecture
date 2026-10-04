# -*- coding: utf-8 -*-
"""519: the exact convex program behind (E_3), solved.

c_8 > 0 is a POSITIVITY question about the spectral measure, so the
sharpest thing any spectral argument can say is the value of

   minimise  SUM_i mu_i^8   over measures with the known moments,

a truncated moment problem: a linear program over measures, whose
dual is the sum-of-squares certificate.  It is solvable in closed
form, and the closed form settles what convex optimisation can and
cannot do here.

SET-UP, girth >= 7 so that c_5 = c_6 = 0.  The Perron eigenvalue 3 is
simple; write the remaining n-1 eigenvalues as a measure nu of total
mass n-1 on [-3,3].  Then
    SUM mu = -3, SUM mu^2 = 3n-9, SUM mu^3 = -27,
    SUM mu^4 = 15n-81, SUM mu^5 = -243, SUM mu^6 = 87n-729,
and 517.1 gives  16 c_8 = SUM_all mu^8 - 543 n, so

    c_8 > 0   <=>   SUM_{i>=2} mu^8  >  543 n - 6561 .

SYMMETRISING loses nothing at leading order, so put y = mu^2 and work
with the measure of y on [0,9] of mass N = n-1 and moments
    M1 = 3n-9,  M2 = 15n-81,  M3 = 87n-729 .
The minimum of M4 over such measures is attained at a TWO-ATOM
measure (three constraints plus mass), and the two-atom system solves
in closed form:  with mean m, variance s2 and third central moment
t3, the atoms sit at distances A, B either side of the mean with
    A - B = t3/s2 ,  A B = s2 ,
and the minimal fourth central moment is s2 (A-B)^2 + s2^2, i.e.
    mu4c_min = s2 . ( (t3/s2)^2 + s2 ) .
At the limiting moments (m, s2, t3) = (3, 6, 6) this gives 519 - so
the LP OPTIMUM IS 519 n + O(1), the 517.1 sum-of-squares bound is
TIGHT, and the 24n gap is not slack.
"""

import sys
from fractions import Fraction as F


def lp_min_raw4(N, M1, M2, M3):
    """exact minimum of the 4th raw moment of a non-negative measure
    of mass N on the line with raw moments M1, M2, M3 (no support
    constraint used; the optimum lands inside [0,9] anyway, which is
    checked)."""
    m = F(M1, N)
    m2 = F(M2, N)
    m3 = F(M3, N)
    s2 = m2 - m * m
    t3 = m3 - 3 * m * m2 + 2 * m ** 3
    # A - B = t3/s2 , A B = s2
    d = t3 / s2
    # mu4c = s2 (d^2 + s2)
    mu4c = s2 * (d * d + s2)
    mu4 = mu4c + 4 * m * t3 + 6 * m * m * s2 + m ** 4
    # atoms, for the support check
    import math
    dd = float(d)
    ss = float(s2)
    A = (dd + math.sqrt(dd * dd + 4 * ss)) / 2.0
    B = A - dd
    return mu4 * N, float(m) + A, float(m) - B


def threshold(gmin=7):
    print("girth >= %d.  exact LP minimum of SUM_{i>=2} mu^8 against "
          "the target 543n - 6561" % gmin)
    print()
    print("   n     LP min          target       proves c_8>0?   "
          "atoms y =")
    last = None
    for n in list(range(10, 60, 2)) + [80, 120, 200, 400]:
        N = n - 1
        M1 = 3 * n - 9
        M2 = 15 * n - 81
        M3 = 87 * n - 729
        val, ya, yb = lp_min_raw4(N, M1, M2, M3)
        tgt = 543 * n - 6561
        ok = val > tgt
        if last is not None and last and not ok:
            print("   ---- last n at which the LP decides: %d ----"
                  % lastn)
        if ok:
            lastn = n
        last = ok
        print("  %3d  %14.1f  %14d      %-6s     %.3f, %.3f"
              % (n, float(val), tgt, "YES" if ok else "no", ya, yb))
    sys.stdout.flush()


def limit_check():
    print()
    print("limiting case: mass 1, moments (1, 3, 15, 87) in y")
    m, s2, t3 = F(3), F(6), F(6)
    d = t3 / s2
    mu4c = s2 * (d * d + s2)
    mu4 = mu4c + 4 * m * t3 + 6 * m * m * s2 + m ** 4
    print("   minimal 4th raw moment =", mu4, " (i.e. SUM mu^8 = "
          "%s n)" % mu4)
    A = (d + (d * d + 4 * s2) ** F(1, 2)) / 2 if False else None
    # solve A-B=1, AB=6 exactly
    print("   atoms: A-B = %s, A B = %s  ->  A = 3, B = 2, "
          "y = 6 and y = 1" % (d, s2))
    print("   weights: u = AB/((A+B)A) = 6/(5*3) = 2/5 at y=6, "
          "3/5 at y=1")
    chk = F(2, 5) * 6 + F(3, 5) * 1
    chk2 = F(2, 5) * 36 + F(3, 5) * 1
    chk3 = F(2, 5) * 216 + F(3, 5) * 1
    chk4 = F(2, 5) * 1296 + F(3, 5) * 1
    print("   check moments: %s, %s, %s  (want 3, 15, 87)   "
          "4th = %s" % (chk, chk2, chk3, chk4))
    print("   so mu in {+-1, +-sqrt 6}: sqrt 6 = 2.449 < 2 sqrt 2 = "
          "2.828, inside [-3,3].")
    print("   THE LP OPTIMUM IS 519, NOT 543.  Kesten-McKay (543) is "
          "feasible but not")
    print("   extremal, and 517.1's SOS bound is exactly tight.")
    sys.stdout.flush()


if __name__ == '__main__':
    limit_check()
    print()
    threshold()
