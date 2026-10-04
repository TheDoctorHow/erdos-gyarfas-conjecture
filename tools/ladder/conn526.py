# -*- coding: utf-8 -*-
"""526: the connectivity criterion in closed form.

496.1 says the components of Z(T) are the cycles of the BLOCK
MULTIGRAPH: sort T's 2k endpoints q_1 < ... < q_2k, group them into
blocks {q_1,q_2}, {q_3,q_4}, ..., and join blocks by the chords.
Every block has two endpoints and every endpoint one chord, so the
multigraph is 2-REGULAR - a disjoint union of cycles - and Z(T) is
connected iff it is a single cycle.

So the component count is the orbit count of the product of TWO
fixed-point-free involutions on T's endpoints: the block pairing and
the matching.  That is exactly the Cohn-Lempel setting, applied to
the sub-diagram T with its OWN circular order.  524.2 tested the
formula against the FULL matching with the alternate-arc transition
and it failed; the claim here is that per-subset, with the block
pairing, it holds:

    #components of Z(T)  =  1 + nullity_GF(2) ( interlacement of T )

where two chords of T interlace if their endpoints alternate in the
circular order of T's endpoints.  If true this is a closed-form
connectivity test for ALL |T|, where 501.2 gave only |T| = 3.
"""

import sys
from itertools import combinations


def girth5_open(chords, n):
    for a, b in chords:
        if b - a < 4:
            return False
    for (a, b), (c, d) in combinations(chords, 2):
        if not (b < c or d < a):
            if abs(a - c) + abs(b - d) < 3:
                return False
    return True


def ladders(n):
    mate = [0] * (n + 1)
    ch, out = [], []

    def rec():
        a = 0
        for i in range(1, n + 1):
            if not mate[i]:
                a = i
                break
        if not a:
            if girth5_open(sorted(ch), n):
                out.append(sorted(ch))
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


def block_components(T):
    """component count via the block multigraph of 496.1."""
    q = sorted(v for e in T for v in e)
    k = len(T)
    blk = {}
    for i in range(k):
        blk[q[2 * i]] = i
        blk[q[2 * i + 1]] = i
    nb = {}
    for a, b in T:
        nb.setdefault(blk[a], []).append(blk[b])
        nb.setdefault(blk[b], []).append(blk[a])
    seen, c = set(), 0
    for s in range(k):
        if s in seen:
            continue
        c += 1
        st = [s]
        seen.add(s)
        while st:
            x = st.pop()
            for y in nb.get(x, []):
                if y not in seen:
                    seen.add(y)
                    st.append(y)
    return c


def interlace_nullity(T):
    """nullity over GF(2) of the interlacement matrix of T, taken in
    the circular order of T's OWN endpoints."""
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
    piv = []
    rank = 0
    for r in rows:
        for p in piv:
            r = min(r, r ^ p)
        if r:
            piv.append(r)
            piv.sort(reverse=True)
            rank += 1
    return k - rank


def run(n):
    L = ladders(n)
    tot = ok = 0
    byk = {}
    for M in L:
        for r in range(1, len(M) + 1):
            for T in combinations(M, r):
                T = list(T)
                c = block_components(T)
                pred = 1 + interlace_nullity(T)
                tot += 1
                d = byk.setdefault(r, [0, 0])
                d[0] += 1
                if c == pred:
                    ok += 1
                    d[1] += 1
    print("2N=%2d  ladders=%4d  subsets=%7d   criterion holds "
          "%7d / %7d" % (n, len(L), tot, ok, tot))
    print("       by |T|:", " ".join("%d:%d/%d" % (k, v[1], v[0])
                                     for k, v in sorted(byk.items())))
    sys.stdout.flush()


if __name__ == '__main__':
    run(12)
    run(14)
