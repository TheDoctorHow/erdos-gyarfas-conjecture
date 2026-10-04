# -*- coding: utf-8 -*-
"""529: does 528.2's inequality actually fire?

For each chord of a closed ladder there are three options: omit it,
or take it into the disjoint family with one of its two arcs.  A
family is admissible when the chosen arcs are pairwise disjoint arcs
of the host cycle (they can never share an endpoint, the matching
being a matching).  The shortcut length is then

     2N  -  SUM ( s_e - 1 )  =  (gaps)  +  |D| .

This enumerates all admissible families, collects the achievable
shortcut LENGTHS, and reports the longest full run of the right
parity among them - which is exactly what 528.2 needs.  It then
compares that with the run in the FULL spectrum, to see how much of
(RUN) the shortcut mechanism supplies on its own.
"""

import sys
from itertools import product


def arcs_of(a, b, n):
    """the two arcs between a and b as vertex lists (endpoints
    included), given as (length, set of interior vertices)."""
    fwd = []
    x = a
    while x != b:
        x = x % n + 1
        if x != b:
            fwd.append(x)
    bwd = []
    x = b
    while x != a:
        x = x % n + 1
        if x != a:
            bwd.append(x)
    la = (b - a) % n
    return [(la, set(fwd)), (n - la, set(bwd))]


def shortcut_lengths(chords, n):
    opts = []
    for a, b in chords:
        o = [None]
        for L, interior in arcs_of(a, b, n):
            o.append((L, interior | {a, b}))
        opts.append(o)
    out = set()
    # depth-first over chords, maintaining the used vertex set
    res = []

    def rec(i, used, tot, k):
        if i == len(opts):
            if k:
                res.append(n - tot + k)
            return
        rec(i + 1, used, tot, k)
        for choice in opts[i][1:]:
            L, S = choice
            if used & S:
                continue
            rec(i + 1, used | S, tot + L, k + 1)

    rec(0, frozenset(), 0, 0)
    return sorted(set(res))


def longest_run(vals, step):
    S = set(vals)
    best = (0, 0, 0.0)
    for u in sorted(S):
        v = u
        while v + step in S:
            v += step
        if u and v / float(u) > best[2]:
            best = (u, v, v / float(u))
    return best


def report(name, chords, n, full=None):
    L = shortcut_lengths(chords, n)
    bip = all(x % 2 == 0 for x in L) if L else False
    step = 2 if bip else 1
    u, v, r = longest_run(L, step)
    smax = n - min(L) if L else 0
    smin = n - max(L) if L else 0
    ineq = 2 * smax - smin
    print("%-16s 2N=%3d  shortcut lengths=%3d  longest run [%2d,%2d] "
          "ratio %.2f   2S_max-S_min = %3d vs 2N = %3d  %s"
          % (name, n, len(L), u, v, r, ineq, n,
             "FIRES" if r >= 2 else "does not fire"))
    print("     lengths:", L)
    sys.stdout.flush()


def dodeca_chords():
    """dodecahedron as a closed ladder: 506.4's chord list on [1,20]
    plus the closing edge, i.e. the ladder whose host cycle is the
    Hamiltonian cycle."""
    return [(1, 11), (2, 9), (3, 7), (4, 20), (5, 18), (6, 16),
            (8, 15), (10, 14), (12, 19), (13, 17)]


if __name__ == '__main__':
    report("dodecahedron", dodeca_chords(), 20)
    # the 507.0 witness, closed
    W = [(1, 5), (2, 7), (3, 11), (4, 17), (6, 22), (8, 13), (9, 21),
         (10, 15), (12, 16), (14, 19), (18, 23), (20, 24)]
    report("507.0 closed", W, 24)
    # the girth-6 witness of 524.0
    W30 = [(1, 6), (2, 10), (3, 13), (4, 17), (5, 25), (7, 19),
           (8, 29), (9, 14), (11, 23), (12, 26), (15, 28), (16, 21),
           (18, 24), (20, 30), (22, 27)]
    report("524.0 girth-6", W30, 30)
