# -*- coding: utf-8 -*-
"""525: one check of the Decomposition Theorem, then no more machine.

CLAIM.  In the open ladder, write f_T for the prefix-parity function
of a chord set T: f_T(i) = 1 iff an odd number of endpoints of T lie
in [1, i].  The completion Z(T) uses exactly the chords of T and the
host edges at the positions where f_T = 1, so

    w(T)  =  |T|  +  wt(f_T) .

If Z(T) has components Z_1, ..., Z_k with chord sets T_1, ..., T_k,
then f_T = SUM_j f_{T_j} over GF(2), and a position carried by two of
them would cancel - so the SUPPORTS OF THE f_{T_j} ARE PAIRWISE
DISJOINT, each Z_j is the completion of its own T_j, and

    w(T)  =  SUM_j w(T_j) ,   each  w(T_j)  a CONNECTED weight.

Checked here on every chord subset of every girth-5 ladder at
2N = 12 and 14, then the argument is used on paper.
"""

import sys
from itertools import combinations


def girth5_open(chords, n):
    for a, b in chords:
        if b - a < 4:
            return False
    for (a, b), (c, d) in combinations(chords, 2):
        if not (b < c or d < a):                      # non-disjoint
            if abs(a - c) + abs(b - d) < 3:
                return False
    return True


def completion(T, n):
    """chords T plus the host edges at odd prefix parity; returns the
    edge set and the parity vector."""
    ends = sorted(v for e in T for v in e)
    f = [0] * (n + 1)
    par = 0
    j = 0
    for i in range(1, n):
        while j < len(ends) and ends[j] <= i:
            par ^= 1
            j += 1
        f[i] = par
    edges = [tuple(sorted(e)) for e in T]
    edges += [(i, i + 1) for i in range(1, n) if f[i]]
    return edges, f


def components(edges):
    nb = {}
    for a, b in edges:
        nb.setdefault(a, []).append(b)
        nb.setdefault(b, []).append(a)
    seen, comps = set(), []
    for s in nb:
        if s in seen:
            continue
        st, cur = [s], set()
        seen.add(s)
        while st:
            x = st.pop()
            cur.add(x)
            for y in nb[x]:
                if y not in seen:
                    seen.add(y)
                    st.append(y)
        comps.append(cur)
    return comps


def check(n):
    mate = [0] * (n + 1)
    ch = []
    lads = []

    def rec():
        a = 0
        for i in range(1, n + 1):
            if not mate[i]:
                a = i
                break
        if not a:
            if girth5_open(sorted(ch), n):
                lads.append(sorted(ch))
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

    tested = addit = disj = 0
    for M in lads:
        for r in range(1, len(M) + 1):
            for T in combinations(M, r):
                edges, f = completion(list(T), n)
                comps = components(edges)
                w = len(edges)
                tested += 1
                # split T by component
                tot = 0
                sup = []
                ok = True
                for C in comps:
                    Tj = [e for e in T if e[0] in C]
                    ej, fj = completion(Tj, n)
                    if set(ej) - set(edges):
                        ok = False
                    tot += len(ej)
                    sup.append(set(i for i in range(1, n) if fj[i]))
                if ok and tot == w:
                    addit += 1
                pairwise = all(not (sup[i] & sup[j])
                               for i in range(len(sup))
                               for j in range(i + 1, len(sup)))
                if pairwise:
                    disj += 1
    print("2N=%2d  girth-5 open ladders=%4d  chord subsets tested=%7d"
          % (n, len(lads), tested))
    print("       weight additive over components : %7d / %7d"
          % (addit, tested))
    print("       component parity supports disjoint: %7d / %7d"
          % (disj, tested))
    sys.stdout.flush()


if __name__ == '__main__':
    check(12)
    check(14)
