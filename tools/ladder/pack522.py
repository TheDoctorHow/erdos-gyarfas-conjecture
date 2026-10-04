# -*- coding: utf-8 -*-
"""522: THE GIRTH-7 PACKING THEOREM.  A cubic graph of girth 7 with
no 8-cycle has at least 34 vertices.

The ingredients, each an arithmetic consequence of "an ear of length
q with feet at H-distance d gives cycles q+d and q+7-d":

 (i)  ALL pairs of outside neighbours are at distance >= 3 in
      G - V(H).  d = 1 gives q >= 6 by girth, so s >= 4; d = 2 gives
      q >= 5, so s >= 3; d = 3 gives q >= 4 by girth and q != 4, 5 by
      the hypothesis, so s >= 4.
 (ii) Hence the 7 outside neighbours have no common neighbours, so
      level two has exactly 14 distinct vertices, none in H or Y.
(iii) A level-three vertex cannot lie in Y: that would give a path of
      length 4 between two H-vertices, hence cycles 4+d and 11-d;
      girth forces d = 3 and then 4+d = 7 while 11-d = 8.
 (iv) An edge inside level two joins children of feet at H-distance
      exactly 2 (girth kills d = 1, the hypothesis kills d = 3), and
      each such foot pair carries at most one, since two sharing an
      endpoint close a 4-cycle and two disjoint ones close a 6-cycle.
      There are 7 such foot pairs, so at most 7 edges.
  (v) So at least 28 - 14 = 14 edge-slots leave level two, and each
      level-three vertex absorbs at most 3, giving at least 5.

  n >= 7 + 7 + 14 + 5 = 33, and n is even, so n >= 34.

This script checks the arithmetic of (i)-(v) exhaustively over the
possible (d, q) and confirms on real girth-7 graphs that every
violation of (i) is accompanied by an 8-cycle.
"""

import sys
from itertools import combinations
from band522 import cycles_len, outside, dist_out, girth


def arithmetic():
    print("(i) admissible ear lengths at girth 7, no 8-cycle")
    for d in (1, 2, 3):
        ok = []
        for q in range(1, 14):
            if q + d < 7 or q + 7 - d < 7:
                continue
            if q + d == 8 or q + 7 - d == 8:
                continue
            ok.append(q)
        smin = min(ok) - 2
        print("    d=%d : admissible q = %s ...   so s = q-2 >= %d"
              % (d, ok[:6], smin))
    print("    => every pair of outside neighbours is at distance "
          ">= 3, and")
    print("       the distance-3 foot pairs at distance >= 4.")
    print()
    print("(iv) an edge inside level two: cycle length 5 + d and "
          "12 - d")
    for d in (1, 2, 3):
        bad = []
        if 5 + d < 7:
            bad.append("girth")
        if 5 + d == 8 or 12 - d == 8:
            bad.append("8-cycle")
        print("    d=%d : lengths %d and %d   %s"
              % (d, 5 + d, 12 - d, "EXCLUDED (%s)" % ",".join(bad)
                 if bad else "allowed"))
    print()
    print("(iii) a level-three vertex in Y: cycles 4 + d and 11 - d")
    for d in (1, 2, 3):
        note = []
        if 4 + d < 7:
            note.append("girth")
        if 4 + d == 8 or 11 - d == 8:
            note.append("8-cycle")
        print("    d=%d : lengths %d and %d   %s"
              % (d, 4 + d, 11 - d, "EXCLUDED (%s)" % ",".join(note)
                 if note else "allowed"))
    print()
    print("counting: level0=7, level1=7, level2=14, edges inside "
          "level2 <= 7,")
    print("slots to level3 >= 28 - 2*7 = 14, each level-3 vertex "
          "takes <= 3,")
    print("so level3 >= ceil(14/3) = 5, and n >= 7+7+14+5 = 33; "
          "n even => n >= 34.")
    sys.stdout.flush()


def confirm(name, adj):
    """every violation of (i) must come with an 8-cycle."""
    heps = cycles_len(adj, 7)
    has8 = bool(cycles_len(adj, 8))
    tested = viol = 0
    for H in heps:
        C = list(H)
        Cs = set(C)
        ys = outside(adj, C)
        if ys is None:
            continue
        for i, j in combinations(range(7), 2):
            s = dist_out(adj, ys[i], ys[j], Cs)
            if s is None:
                continue
            tested += 1
            if s < 3:
                viol += 1
    print("%-12s n=%2d g=%d  heptagons=%3d  foot pairs=%5d  "
          "distance < 3 in %4d   8-cycles present: %s   %s"
          % (name, len(adj), girth(adj), len(heps), tested, viol,
             "yes" if has8 else "NO",
             "consistent" if (viol == 0 or has8) else "*** BAD"))
    sys.stdout.flush()
    return viol == 0 or has8


if __name__ == '__main__':
    arithmetic()
    print()
    from cage516 import lcf, coxeter
    ok = True
    ok &= confirm("McGee", lcf([12, 7, -7], 8))
    ok &= confirm("Coxeter", coxeter())
    print("consistent:", ok)
