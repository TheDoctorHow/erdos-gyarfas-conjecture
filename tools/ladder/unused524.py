# -*- coding: utf-8 -*-
"""524: properties of the CLOSED ladder the archive has not used.
Each claim is checked, not asserted.

(1) 3-EDGE-COLOURABILITY.  A closed ladder is a Hamiltonian cubic
    graph on an even number of vertices, so colouring the host cycle
    alternately 1, 2 and every chord 3 is proper.  NO CLOSED LADDER
    IS A SNARK.

(2) THE TWO ALTERNATING 2-FACTORS.  Colours {1,3} and {2,3} each give
    a 2-factor: N chords plus the N host edges of one parity.  Their
    cycle lengths partition 2N.  These are the two completions of the
    FULL chord set, and the archive has never looked at them.

(3) THE COHN-LEMPEL / INTERLACEMENT NULLITY FORMULA.  For a chord
    set T, the number of components of the all-arcs transition is
    1 + nullity over GF(2) of the interlacement (crossing) matrix of
    T.  501.7 set this aside because a ladder completion keeps
    ALTERNATE arcs; in the CLOSED ladder the all-arcs transition is
    exactly a 2-factor of (2), so the formula applies on the nose.

(4) THE SPAN PARITY.  Sum of spans = sum of all labels = N(2N+1),
    so the number of ODD spans is congruent to N mod 2.

(5) THE CHORD-DIAGRAM GENUS.  With c the component count of (3),
    1 - N + c = 2 - 2*gamma defines the genus of the chord diagram -
    a second genus, never compared with the embedding genus that
    516.5's fuel formula uses.
"""

import sys
from itertools import combinations


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


def gen_ladders(n, cap=200):
    out = []
    mate = [0] * (n + 1)
    ch = []

    def rec():
        if len(out) >= cap:
            return
        a = 0
        for i in range(1, n + 1):
            if not mate[i]:
                a = i
                break
        if not a:
            if girth5(ch, n):
                out.append(list(ch))
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
    return out


def check_colouring(chords, n):
    """host edge {i,i+1} gets colour 1 if i odd else 2; chords get 3."""
    col = {}
    for i in range(1, n + 1):
        j = i % n + 1
        col[frozenset((i, j))] = 1 if i % 2 else 2
    for a, b in chords:
        col[frozenset((a, b))] = 3
    adj = closed(chords, n)
    for v in adj:
        cs = [col[frozenset((v, w))] for w in adj[v]]
        if sorted(cs) != [1, 2, 3]:
            return False, None
    return True, col


def two_factor_parts(chords, n, keep):
    """cycle lengths of the 2-factor = chords + host edges {i,i+1}
    with i congruent to `keep` mod 2."""
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
    seen = set()
    parts = []
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
        parts.append(L)
    return sorted(parts)


def interlace_nullity(T):
    """GF(2) nullity of the crossing matrix of the chord set T."""
    k = len(T)
    rows = []
    for i in range(k):
        a, b = sorted(T[i])
        r = 0
        for j in range(k):
            if i == j:
                continue
            c, d = sorted(T[j])
            if a < c < b < d or c < a < d < b:
                r |= 1 << j
        rows.append(r)
    rank = 0
    piv = []
    for r in rows:
        for p in piv:
            r = min(r, r ^ p)
        if r:
            piv.append(r)
            piv.sort(reverse=True)
            rank += 1
    return k - rank


def main():
    n = 14
    lads = gen_ladders(n, cap=60)
    print("closed girth-5 ladders at 2N=%d tested: %d" % (n, len(lads)))
    okcol = 0
    parity_ok = 0
    cl_ok = cl_bad = 0
    genus = {}
    for ch in lads:
        good, col = check_colouring(ch, n)
        okcol += good
        odd = sum(1 for a, b in ch if (b - a) % 2)
        if odd % 2 == (n // 2) % 2:
            parity_ok += 1
        p0 = two_factor_parts(ch, n, 0)
        p1 = two_factor_parts(ch, n, 1)
        # Cohn-Lempel on the FULL chord set against one 2-factor
        c = len(p1)
        if c == 1 + interlace_nullity(ch):
            cl_ok += 1
        else:
            cl_bad += 1
        g = (1 - len(ch) + c - 2) // -2
        genus[g] = genus.get(g, 0) + 1
    print("(1) proper 3-edge-colouring by the canonical recipe: "
          "%d/%d   -> no closed ladder is a snark" % (okcol, len(lads)))
    print("(4) #odd spans == N (mod 2):                          "
          "%d/%d" % (parity_ok, len(lads)))
    print("(3) Cohn-Lempel  c = 1 + nullity(interlacement):      "
          "%d ok, %d bad" % (cl_ok, cl_bad))
    print("(5) chord-diagram genus distribution:", sorted(genus.items()))
    print()
    print("(2) the two alternating 2-factors, first six ladders:")
    for ch in lads[:6]:
        p0 = two_factor_parts(ch, n, 0)
        p1 = two_factor_parts(ch, n, 1)
        print("    chords %s" % (ch,))
        print("        odd-host 2-factor  %s  (sum %d)" % (p1, sum(p1)))
        print("        even-host 2-factor %s  (sum %d)" % (p0, sum(p0)))
    sys.stdout.flush()


if __name__ == '__main__':
    main()
