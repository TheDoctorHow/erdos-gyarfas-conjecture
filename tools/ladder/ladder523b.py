# -*- coding: utf-8 -*-
"""523: THE LADDER QUESTION of 449.4 pushed past mu = 16 - fast form.

Cycles of a path-plus-chords graph are enumerated by the COMPLETION
CRITERION of 496.1 rather than by path search: a set T of chords
completes to a single cycle iff, sorting its 2|T| endpoints
q_1 < ... < q_{2|T|} and grouping them into consecutive blocks
(q_1,q_2), (q_3,q_4), ..., no chord of T lies inside a block and the
multigraph whose nodes are the blocks and whose edges are the chords
is connected; the length is then

    |T|  +  SUM_i ( q_{2i} - q_{2i-1} ) .

Only subsets containing the newly placed chord need testing, and the
prunes are monotone, so a reported maximum is a proof of absence at
those parameters.
"""

import sys
from itertools import combinations


def forbidden(t, cap):
    bad = set()
    v = 4
    while v <= cap:
        bad.add(v)
        v *= 2
    q = 1 << t
    while q <= cap:
        bad.add(q)
        q += 1 << t
    return bad


def cyc_len(T):
    """length of the completion of chord set T, or None if it is not
    a single cycle."""
    q = sorted(v for e in T for v in e)
    k = len(T)
    blocks = [(q[2 * i], q[2 * i + 1]) for i in range(k)]
    idx = {}
    for bi, (lo, hi) in enumerate(blocks):
        idx[lo] = bi
        idx[hi] = bi
    for a, b in T:
        if idx[a] == idx[b]:
            return None                      # chord inside a block
    seen = set([0])
    stack = [0]
    nb = {}
    for a, b in T:
        nb.setdefault(idx[a], []).append(idx[b])
        nb.setdefault(idx[b], []).append(idx[a])
    while stack:
        x = stack.pop()
        for y in nb.get(x, []):
            if y not in seen:
                seen.add(y)
                stack.append(y)
    if len(seen) != k:
        return None
    return k + sum(hi - lo for lo, hi in blocks)


def feasible(s, mu, t):
    bad = forbidden(t, s + 2)
    mate = [0] * (s + 2)
    chords = []

    def bad_new(new):
        rest = chords
        for r in range(0, len(rest) + 1):
            for sub in combinations(rest, r):
                T = list(sub) + [new]
                L = cyc_len(T)
                if L is not None and L in bad:
                    return True
        return False

    def rec():
        a = 0
        for i in range(1, s + 1):
            if not mate[i]:
                a = i
                break
        if not a:
            return True
        if a + 4 > s:
            return False
        for b in range(a + 4, min(a + mu, s) + 1):
            if mate[b]:
                continue
            if bad_new((a, b)):
                continue
            mate[a] = b
            mate[b] = a
            chords.append((a, b))
            if rec():
                return True
            chords.pop()
            mate[a] = 0
            mate[b] = 0
        return False

    return rec()


def maxlen(mu, t, hi=28, patience=3):
    best, miss = 0, 0
    for s in range(4, hi + 1, 2):
        if feasible(s, mu, t):
            best, miss = s, 0
        else:
            miss += 1
            if miss >= patience and best:
                break
    return best


if __name__ == '__main__':
    print("LADDER QUESTION (449.4): longest path 1..s with a perfect")
    print("matching of spans in [4, mu] whose ladder avoids every")
    print("power of two and every multiple of 2^t.")
    print("449.2 measured the mu <= 16 columns; the rest is new.")
    print()
    mus = [6, 8, 10, 12, 16, 20, 24, 32, 48]
    print("  t \\ mu" + "".join("%5d" % m for m in mus))
    for t in (3, 4, 5, 6):
        row = []
        for mu in mus:
            row.append(maxlen(mu, t))
        print("    %d  " % t + "".join("%5d" % v for v in row))
        sys.stdout.flush()
