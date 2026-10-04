# -*- coding: utf-8 -*-
"""522: the girth-6 and girth-7 analogues of the pentagram matching
lemma, and the involution behind 521.3's obstruction.

GIRTH 6.  Hexagon H = x_0..x_5, outside neighbours y_0..y_5 (distinct
by girth 6).  y_i ~ y_j closes a cycle of length 3 + d with d the
H-distance, so girth 6 forces d = 3: only ANTIPODAL adjacencies are
possible, and there are three candidate pairs.
  LEMMA.  No 8-cycle => AT MOST ONE of the three occurs.
  PROOF.  Say y_0 ~ y_3 and y_1 ~ y_4.  Then
    H1 = x_0 y_0 y_3 x_3 x_4 x_5      and
    H2 = x_1 y_1 y_4 x_4 x_5 x_0
  are hexagons meeting in exactly the path x_4 x_5 x_0, of length 2,
  so their symmetric difference is a cycle of length 6+6-4 = 8.  For
  y_0 ~ y_3 and y_2 ~ y_5 take H1 and
    H3 = x_2 y_2 y_5 x_5 x_4 x_3 ,
  which meet in exactly x_3 x_4 x_5, again a 2-path.  QED

GIRTH 7.  Heptagon x_0..x_6.  y_i ~ y_j would close a cycle of length
3 + d <= 6 < 7, so THERE ARE NO y-adjacencies at all, by girth alone.
The 8-cycle hypothesis instead constrains the outside DISTANCES: an
ear of length q with feet at H-distance d gives cycles q + d and
q + 7 - d, so
    d = 1 : q >= 6 by girth, and q != 7      (s != 5)
    d = 2 : q >= 5 by girth, and q != 6      (s != 4)
    d = 3 : q >= 4 by girth, and q != 4, 5   (s >= 4)
  LEMMA.  In a cubic graph of girth 7 with no 8-cycle, the outside
  neighbours of two heptagon vertices at distance 3 are at distance
  at least 4 in G - V(H).

THE INVOLUTION.  521.3's obstruction is exactly this: if C has the
u-adjacency u_i ~ u_j then the new pentagon C' = v_i u_i u_j v_j v_x
has, as its own forced u-adjacency, the outside neighbours of v_j and
v_i - and applying the construction to it returns C.  The map is an
INVOLUTION, and the invariant object is the theta Theta(2,3,3)
= C union C'.
"""

import sys
from itertools import combinations
from collections import deque


def cycles_len(adj, L):
    out = []
    vs = sorted(adj)
    rank = dict((v, i) for i, v in enumerate(vs))
    for s in vs:
        stack = [(s, (s,), 1 << s)]
        while stack:
            v, path, msk = stack.pop()
            for w in adj[v]:
                if w == s:
                    if len(path) == L and rank[path[1]] < rank[path[-1]]:
                        out.append(path)
                elif (rank[w] > rank[s] and not (msk >> w & 1)
                      and len(path) < L):
                    stack.append((w, path + (w,), msk | (1 << w)))
    return out


def girth(adj):
    for L in range(3, 15):
        if cycles_len(adj, L):
            return L
    return 99


def outside(adj, C):
    Cs = set(C)
    ys = []
    for v in C:
        o = [w for w in adj[v] if w not in Cs]
        if len(o) != 1:
            return None
        ys.append(o[0])
    return ys if len(set(ys)) == len(C) else None


def dist_out(adj, a, b, banned):
    if a in banned or b in banned:
        return None
    if a == b:
        return 0
    seen = set([a]) | set(banned)
    dq = deque([(a, 0)])
    while dq:
        v, d = dq.popleft()
        for w in adj[v]:
            if w == b:
                return d + 1
            if w not in seen:
                seen.add(w)
                dq.append((w, d + 1))
    return None


def check6(name, adj):
    """girth-6 lemma: at most one antipodal y-adjacency."""
    hexes = cycles_len(adj, 6)
    has8 = bool(cycles_len(adj, 8))
    viol = 0
    for H in hexes:
        C = list(H)
        ys = outside(adj, C)
        if ys is None:
            continue
        k = sum(1 for i in range(3) if ys[i + 3] in adj[ys[i]])
        if k >= 2:
            viol += 1
    print("%-16s g=%d hexagons=%4d  8-cycles=%s   hexagons with >=2 "
          "antipodal y-adjacencies=%3d   %s"
          % (name, girth(adj), len(hexes), "yes" if has8 else "NO",
             viol, "consistent" if (viol == 0 or has8) else "*** BAD"))
    sys.stdout.flush()
    return viol == 0 or has8


def check7(name, adj):
    """girth-7 lemma: distance-3 feet have outside distance >= 4."""
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
            if min((j - i) % 7, (i - j) % 7) != 3:
                continue
            s = dist_out(adj, ys[i], ys[j], Cs)
            if s is None:
                continue
            tested += 1
            if s < 4:
                viol += 1
    print("%-16s g=%d heptagons=%4d 8-cycles=%s   distance-3 foot "
          "pairs tested=%5d   outside distance < 4 in %3d   %s"
          % (name, girth(adj), len(heps), "yes" if has8 else "NO",
             tested, viol,
             "consistent" if (viol == 0 or has8) else "*** BAD"))
    sys.stdout.flush()
    return viol == 0 or has8


def involution(adj):
    """verify C -> C' -> C on every pentagon with a u-adjacency."""
    pent = cycles_len(adj, 5)
    keys = set()
    for P in pent:
        keys.add(frozenset(P))
    ok = 0
    tot = 0
    for P in pent:
        C = list(P)
        Cs = set(C)
        ys = outside(adj, C)
        if ys is None:
            continue
        for i in range(5):
            j = (i + 2) % 5
            if ys[j] not in adj[ys[i]]:
                continue
            x = [t for t in range(5)
                 if (t - i) % 5 in (1, 4) and (t - j) % 5 in (1, 4)][0]
            Cp = [C[i], ys[i], ys[j], C[j], C[x]]
            if frozenset(Cp) not in keys:
                continue
            tot += 1
            # C' order is (v_i, u_i, u_j, v_j, v_x); its outside
            # neighbours at the distance-2 pair (v_j, v_i) are the
            # C-neighbours of v_j and v_i away from v_x
            a = [t for t in range(5) if t != x and (t - j) % 5 in (1, 4)][0]
            b = [t for t in range(5) if t != x and (t - i) % 5 in (1, 4)][0]
            back = set([C[j], C[a], C[b], C[i], C[x]])
            if back == set(C):
                ok += 1
    return ok, tot


if __name__ == '__main__':
    from cage516 import lcf, petersen, coxeter
    from ear521 import rand_cubic_girth5
    good = True
    print("GIRTH 6 LEMMA")
    for nm, g in (("Heawood", lcf([5, -5], 7)),
                  ("Mobius-Kantor", lcf([5, -5], 8)),
                  ("Pappus", lcf([5, 7, -7, 7, -7, -5], 3)),
                  ("Desargues", lcf([5, -5, 9, -9], 5)),
                  ("Nauru", lcf([5, -9, 7, -7, 9, -5], 4)),
                  ("Dyck", lcf([5, -5, 13, -13], 8))):
        good &= check6(nm, g)
    print()
    print("GIRTH 7 LEMMA")
    for nm, g in (("McGee", lcf([12, 7, -7], 8)),
                  ("Coxeter", coxeter())):
        good &= check7(nm, g)
    print()
    print("THE INVOLUTION C -> C' -> C")
    for nm, g in (("Petersen", petersen()),
                  ("dodecahedron",
                   lcf([10, 7, 4, -4, -7, 10, -4, 7, -7, 4], 2))):
        ok, tot = involution(g)
        print("  %-14s u-adjacency instances=%3d   C' maps back to "
              "C in %3d of them" % (nm, tot, ok))
    for k, nn in enumerate((14, 18, 22, 26)):
        gg = rand_cubic_girth5(nn, 700 + k)
        if gg:
            ok, tot = involution(gg)
            print("  random n=%-6d u-adjacency instances=%3d   "
                  "C' maps back to C in %3d of them" % (nn, tot, ok))
    print()
    print("all lemmas consistent:", good)
