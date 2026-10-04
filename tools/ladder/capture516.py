# -*- coding: utf-8 -*-
"""516: does LOCAL CAPTURE hold for genuine cubic graphs?

504.1 (capture lemma): [g, 2g-2] always contains a power of two - and
it is unique there, since the interval has ratio below 2.  Call

   LOCAL CAPTURE:  a cubic graph of girth g has a cycle of length
                   equal to the unique power of two in [g, 2g-2].

At g = 5 that value is 8.  505/507 refuted the ladder analogue, but a
ladder is not cubic - it has two vertices of degree 2.  This settles
the cubic question directly: CLOSE the 507.0 witness (add the edge
2N-1) and test.  If the closed graph is cubic, girth 5 and 8-cycle
free, LOCAL CAPTURE is false for cubic graphs and (N16) has a witness
in the class the conjecture is actually about.
"""

import sys

W507 = [(1, 5), (2, 7), (3, 11), (4, 17), (6, 22), (8, 13), (9, 21),
        (10, 15), (12, 16), (14, 19), (18, 23), (20, 24)]


def closed(chords, n):
    adj = dict((v, []) for v in range(1, n + 1))
    for i in range(1, n + 1):
        j = i % n + 1
        adj[i].append(j)
        adj[j].append(i)
    for a, b in chords:
        adj[a].append(b)
        adj[b].append(a)
    return adj


def spectrum(adj):
    """every cycle length present, by brute-force DFS."""
    got = {}
    vs = sorted(adj)
    for s in vs:
        stack = [(s, 1, 1 << s)]
        while stack:
            v, L, msk = stack.pop()
            for w in adj[v]:
                if w == s:
                    if L >= 3:
                        got[L] = got.get(L, 0) + 1
                elif w > s and not (msk >> w & 1):
                    stack.append((w, L + 1, msk | (1 << w)))
    # each cycle found twice (two directions)
    return dict((L, c // 2) for L, c in got.items())


def girth(adj):
    return min(spectrum(adj))


def report(name, adj):
    sp = spectrum(adj)
    deg = sorted(set(len(adj[v]) for v in adj))
    ks = sorted(L for L in sp if L & (L - 1) == 0)
    print("%-22s n=%2d degrees=%s girth=%2d  lengths=%s"
          % (name, len(adj), deg, min(sp), sorted(sp)))
    print("%-22s powers of two present: %s   (8-cycles: %d)"
          % ("", ks, sp.get(8, 0)))
    sys.stdout.flush()
    return sp


if __name__ == '__main__':
    n = 24
    op = closed(W507, n)
    # open ladder = closed minus the host edge (n,1)
    op[n].remove(1)
    op[1].remove(n)
    report("507.0 OPEN ladder", op)
    report("507.0 CLOSED (cubic)", closed(W507, n))
