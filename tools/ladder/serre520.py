# -*- coding: utf-8 -*-
"""520: augmenting the band-3 convex program with an EDGE MASS
constraint, and computing how much of it is needed.

519 showed the plain program is worth 519n against a target of 543n,
the extremal measure being the two atoms mu^2 = 1, 6 with weights
3/5, 2/5.  That measure has NO MASS ABOVE mu^2 = 6, i.e. none near
the Ramanujan edge 2 sqrt 2 (mu^2 = 8).  Serre-type theorems say a
finite d-regular graph must have a POSITIVE PROPORTION of its
spectrum near the edge - a statement the tree does not satisfy in
the same way, and therefore exactly the kind of input the program
lacks.

This computes, for each threshold Y and each mass fraction c, the
value of

  min  INT y^4 dnu   s.t.  INT dnu = 1, INT y = 3, INT y^2 = 15,
                           INT y^3 = 87,  nu(y >= Y) >= c,
                           supp nu contained in [0, 9],

and reports the least c at which the value exceeds 543 - i.e. the
strength of edge-mass statement that would close band 3.

METHOD.  With four equalities and one active inequality the optimum
is carried by at most five atoms.  Placing mass c at a single point
y0 >= Y and solving the residual three-moment problem in closed form
(519.1) gives the value for that y0; minimising over y0 on a fine
grid gives the program's value over this class, which is an UPPER
bound on the true optimum - so the c reported here is a LOWER bound
on the edge mass actually required.  Flagged as such.
"""

import sys


def two_atom_min(mass, M1, M2, M3):
    """closed form of 519.1: minimal 4th raw moment of a measure of
    total mass `mass` with raw moments M1, M2, M3.  Returns None if
    the data are infeasible (non-positive variance)."""
    if mass <= 0:
        return 0.0 if abs(M1) < 1e-12 else None
    m = M1 / mass
    m2 = M2 / mass
    m3 = M3 / mass
    s2 = m2 - m * m
    if s2 <= 1e-12:
        return None
    t3 = m3 - 3 * m * m2 + 2 * m ** 3
    d = t3 / s2
    mu4c = s2 * (d * d + s2)
    mu4 = mu4c + 4 * m * t3 + 6 * m * m * s2 + m ** 4
    # support check
    disc = d * d + 4 * s2
    A = (d + disc ** 0.5) / 2.0
    B = A - d
    lo, hi = m - B, m + A
    return mu4 * mass, lo, hi


def value_with_edge_mass(c, Y, grid=900):
    """min over y0 in [Y,9] of  c*y0^4 + residual"""
    best = None
    for i in range(grid + 1):
        y0 = Y + (9.0 - Y) * i / grid
        M1 = 3.0 - c * y0
        M2 = 15.0 - c * y0 ** 2
        M3 = 87.0 - c * y0 ** 3
        r = two_atom_min(1.0 - c, M1, M2, M3)
        if r is None:
            continue
        val, lo, hi = r
        if lo < -1e-9 or hi > 9 + 1e-9:
            continue
        tot = c * y0 ** 4 + val
        if best is None or tot < best[0]:
            best = (tot, y0, lo, hi)
    return best


def scan():
    print("band-3 program augmented with 'mass >= c above y = Y'")
    print("target to beat: 543.0     (plain program = 519.0)")
    print()
    for Y in (8.0, 7.5, 7.0, 6.5, 6.0):
        print("  edge threshold y >= %.1f   (|mu| >= %.4f, "
              "Ramanujan edge is 2 sqrt 2 = 2.8284, y = 8)"
              % (Y, Y ** 0.5))
        prev = None
        need = None
        for j in range(0, 61):
            c = j / 200.0
            b = value_with_edge_mass(c, Y)
            if b is None:
                continue
            tot = b[0]
            if need is None and tot > 543.0:
                need = c
            if j % 5 == 0 or (prev is not None and need == c):
                print("      c = %.3f   value = %8.2f   %s"
                      % (c, tot, "CLOSES BAND 3" if tot > 543 else ""))
            prev = tot
            if need is not None:
                break
        print("      least edge mass that closes band 3 at this "
              "threshold: %s" % ("c = %.3f" % need if need else
                                 "not reached below c = 0.30"))
        print()
    sys.stdout.flush()


if __name__ == '__main__':
    b = value_with_edge_mass(0.0, 8.0)
    print("sanity: c = 0 reproduces the plain optimum ->  %.4f  "
          "(want 519)" % b[0])
    print("        extremal atoms  y = %.4f, %.4f  (want 1 and 6)"
          % (b[2], b[3]))
    print()
    scan()
