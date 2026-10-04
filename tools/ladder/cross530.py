# -*- coding: utf-8 -*-
"""530: the crossing families break the mod-d confinement.

THE ALGEBRA.  For a chord e write a_e < b_e, s_e = b_e - a_e and the
POSITION SUM phi_e = a_e + b_e, so phi_e = 2 a_e + s_e.

  CROSSING pair a_1 < a_2 < b_1 < b_2:
      L^1 = (a_2 - a_1) + (b_2 - b_1) = phi_f - phi_e ,
      cycle lengths 2 + L^1 and 2N + 2 - L^1.
  NESTED pair a_1 < a_2 < b_2 < b_1:
      L^1 = (a_2 - a_1) + (b_1 - b_2) = s_e - s_f .
  DISJOINT pair: the shortcut of 528.0.

Now suppose every saving s_e - 1 is divisible by d, i.e. s_e = 1 mod d.
Then
    shortcuts        =  2N        (mod d)
    nested pairs     =  2 + s_e - s_f  =  2      (mod d)
    crossing pairs   =  2 + phi_f - phi_e  =  2 + 2(a_f - a_e)  (mod d)

so the crossing pairs are the ONLY two-chord family whose residue is
not pinned, and it moves with the left endpoints.

THE COUNTING.  Let n_r be the number of chords with a_e = r (mod d).
Every b_e = a_e + 1 (mod d), and the endpoints are all of [1, 2N], so
for each r
      n_r + n_{r-1}  =  #{ i in [1,2N] : i = r mod d } .
When d divides 2N this forces n_{r+2} = n_r, and for ODD d that makes
all n_r equal to N/d - the left endpoints hit EVERY residue.

  THEOREM.  If d >= 3 is odd, d divides 2N, and the crossing graph of
  the chords is connected, then some crossing pair has
  a_f - a_e non-zero mod d, hence a cycle length outside 2N mod d.
    PROOF.  The left endpoints take all d residues, so the labelling
    e -> a_e mod d is non-constant; a connected graph with a
    non-constant vertex labelling has an edge whose ends differ. For
    that crossing pair 2(a_f - a_e) is non-zero mod d since d is odd.
    QED

This verifies the counting, the connectivity hypothesis, and the
residues actually realised.
"""

import sys
from math import gcd
from functools import reduce
from itertools import combinations
from verify_w36 import spectrum
import gap526 as G


def data(ch, n):
    sav = []
    for a, b in ch:
        s = min((b - a) % n, (a - b) % n)
        sav += [s - 1, n - s - 1]
    return reduce(gcd, sav)


def crossing_graph_connected(ch, n):
    k = len(ch)
    adj = dict((i, []) for i in range(k))
    for i, j in combinations(range(k), 2):
        a1, b1 = ch[i]
        a2, b2 = ch[j]
        if (a1 < a2 < b1 < b2) or (a2 < a1 < b2 < b1):
            adj[i].append(j)
            adj[j].append(i)
    seen, st = {0}, [0]
    while st:
        x = st.pop()
        for y in adj[x]:
            if y not in seen:
                seen.add(y)
                st.append(y)
    return len(seen) == k


def left_residues(ch, d):
    from collections import Counter
    return sorted(Counter(min(e) % d for e in ch).items())


def check_dodeca():
    n = 20
    ch = [(1, 11), (2, 9), (3, 7), (4, 20), (5, 18), (6, 16),
          (8, 15), (10, 14), (12, 19), (13, 17)]
    d = data(ch, n)
    sp = sorted(spectrum(ch and {} or {}) ) if False else None
    from verify_w36 import build
    S = sorted(spectrum(build(ch, n)))
    res = sorted(set(x % d for x in S))
    print("dodecahedron: d = %d,  2N mod d = %d" % (d, n % d))
    print("   left endpoints by residue mod %d: %s"
          % (d, left_residues(ch, d)))
    print("   crossing graph connected: %s"
          % crossing_graph_connected(ch, n))
    print("   spectrum %s" % S)
    print("   residues mod %d realised by the spectrum: %s"
          % (d, res))
    print("   shortcut residue alone was {%d}" % (n % d))
    sys.stdout.flush()


def census(n):
    mate = [0] * (n + 1)
    ch = []
    tot = [0]
    dbig = [0]
    nc = [0]

    def rec():
        a = 0
        for i in range(1, n + 1):
            if not mate[i]:
                a = i
                break
        if not a:
            if not G.girth_legal_closed(ch, n):
                return
            tot[0] += 1
            d = data(sorted(ch), n)
            if d >= 3:
                dbig[0] += 1
            if not crossing_graph_connected(sorted(ch), n):
                nc[0] += 1
            return
        for b in range(a + 1, n + 1):
            if mate[b]:
                continue
            mate[a] = b
            mate[b] = a
            ch.append((a, b))
            rec()
            ch.pop()
            mate[a] = mate[b] = 0
    rec()
    print("2N=%2d  girth-5 closed ladders=%5d   with gcd(savings) >= 3:"
          " %4d   with DISCONNECTED crossing graph: %4d"
          % (n, tot[0], dbig[0], nc[0]))
    sys.stdout.flush()


if __name__ == '__main__':
    check_dodeca()
    print()
    for n in (12, 14, 16):
        census(n)
