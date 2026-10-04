# -*- coding: utf-8 -*-
"""524: THE TWO-FACTOR TARGET - a completely unexploited object.

A closed ladder is a Hamiltonian cubic graph on an even number of
vertices, so it is 3-edge-colourable by the canonical recipe: host
edge {i,i+1} takes colour 1 if i is odd and 2 if i is even, and every
chord takes colour 3.  Verified 60/60 at 2N = 14.

That hands over TWO canonical 2-factors - colours {1,3} and {2,3} -
each consisting of the N chords together with the N host edges of one
parity.  Each is a disjoint union of cycles covering all 2N vertices,
so ITS CYCLE LENGTHS PARTITION 2N.  The archive has never looked at
them; the whole ladder programme works with chord subsets and their
completions, and these two subgraphs are not of that form (they use
alternate host edges, not a completion's arcs).

  (2F)  Does at least one of the two alternating 2-factors of a
        girth-legal closed ladder contain a cycle of power-of-two
        length?

(2F) implies Erdos-Gyarfas for closed ladders.  It names no length,
so it is (N16)-safe, and its subject names no power of two, so it is
(N9)-safe.  It is also FREE to test: two subgraphs per ladder, no
completions, no connectivity search.
"""

import sys
from itertools import combinations


def girth5(chords, n):
    def cd(x, y):
        d = abs(x - y) % n
        return min(d, n - d)
    for a, b in chords:
        if cd(a, b) < 4:
            return False
    for (a, b), (c, d) in combinations(chords, 2):
        for x, y in ((c, d), (d, c)):
            if cd(a, x) + cd(b, y) + 2 <= 4:
                return False
    return True


def parts(chords, n, keep):
    sub = dict((v, []) for v in range(1, n + 1))
    for i in range(1, n + 1):
        if i % 2 == keep:
            j = i % n + 1
            sub[i].append(j)
            sub[j].append(i)
    for a, b in chords:
        sub[a].append(b)
        sub[b].append(a)
    if any(len(sub[v]) != 2 for v in sub):
        return None
    seen, out = set(), []
    for s in sub:
        if s in seen:
            continue
        cur, prev, L = s, None, 0
        while True:
            seen.add(cur)
            nxt = [x for x in sub[cur] if x != prev]
            prev, cur = cur, nxt[0]
            L += 1
            if cur == s:
                break
        out.append(L)
    return sorted(out)


def sweep(n, cap=None):
    P2 = set()
    v = 4
    while v <= 4 * n:
        P2.add(v)
        v *= 2
    mate = [0] * (n + 1)
    ch = []
    tot = [0]
    hit = [0]
    fails = []

    def rec():
        if cap and tot[0] >= cap:
            return
        a = 0
        for i in range(1, n + 1):
            if not mate[i]:
                a = i
                break
        if not a:
            if not girth5(ch, n):
                return
            tot[0] += 1
            p0 = parts(ch, n, 0)
            p1 = parts(ch, n, 1)
            ok = any(x in P2 for x in p0) or any(x in P2 for x in p1)
            if ok:
                hit[0] += 1
            elif len(fails) < 3:
                fails.append((list(ch), p1, p0))
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
    print("2N=%2d  girth-5 closed ladders=%6d   at least one 2-factor "
          "has a power-of-two part: %6d   FAILURES: %d"
          % (n, tot[0], hit[0], tot[0] - hit[0]))
    for f in fails:
        print("     FAIL chords=%s" % (f[0],))
        print("          odd-host parts %s ; even-host parts %s"
              % (f[1], f[2]))
    sys.stdout.flush()
    return tot[0] - hit[0]


if __name__ == '__main__':
    bad = 0
    for n in (12, 14, 16, 18, 20):
        bad += sweep(n)
    print()
    print("total failures of (2F):", bad)
