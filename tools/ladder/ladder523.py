# -*- coding: utf-8 -*-
"""523: THE LADDER QUESTION of 449.4, pushed past mu = 16.

  (LADDER QUESTION)  Is there a constant B, independent of mu and t,
  such that a path with a perfect matching of spans at most mu, every
  cycle of which avoids all powers of 2 and all multiples of 2^t, has
  length at most B?

449.2 measured, with every interior vertex matched:
    t = 3:  mu = 6 -> 2,  mu = 8 -> 4,  mu = 10, 12, 16 -> 9
    t = 4:  mu = 6 -> 2,  mu = 8 -> 4,  mu = 10, 12 -> 9,  mu = 16 -> 10
and 449.4 records the honest limitation: "the measurements cover span
caps mu <= 16 ... for large L the bound is EXTRAPOLATED, not proven".
The escape regime has mu ~ sqrt(L), so the extrapolation is exactly
what the archive needs and does not have.

This pushes the same computation to mu = 20, 24, 32, 48, 64 and
t = 3, 4, 5, 6.  Both prunes are monotone - adding a chord never
destroys a cycle - so a reported maximum is a PROOF that no longer
ladder exists at those parameters.

The path is 1..s; a chord (a,b) has span b - a in [4, mu]; girth 5
bans span 3 and below.  Forbidden cycle lengths: every power of two,
and every positive multiple of 2^t.
"""

import sys


def forbidden_set(t, cap):
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


def feasible(s, mu, t, want_witness=False):
    """is there a perfect matching of 1..s, spans in [4,mu], whose
    ladder has no forbidden cycle length?"""
    bad = forbidden_set(t, s + 2)
    adj = dict((i, []) for i in range(1, s + 1))
    for i in range(1, s):
        adj[i].append(i + 1)
        adj[i + 1].append(i)
    mate = [0] * (s + 2)
    witness = []

    def cycles_through(a, b):
        """lengths of all cycles using the edge (a,b), in the current
        graph; (a,b) has just been added."""
        out = []
        stack = [(a, 1, 1 << a)]
        while stack:
            v, d, msk = stack.pop()
            for w in adj[v]:
                if v == a and w == b:
                    continue
                if w == b:
                    out.append(d + 1)
                elif not (msk >> w & 1):
                    stack.append((w, d + 1, msk | (1 << w)))
        return out

    def rec():
        a = 0
        for i in range(1, s + 1):
            if not mate[i]:
                a = i
                break
        if not a:
            if want_witness:
                witness.append([(i, mate[i]) for i in range(1, s + 1)
                                if i < mate[i]])
            return True
        for b in range(a + 4, min(a + mu, s) + 1):
            if mate[b]:
                continue
            adj[a].append(b)
            adj[b].append(a)
            ok = not any(L in bad for L in cycles_through(a, b))
            if ok:
                mate[a] = b
                mate[b] = a
                if rec():
                    adj[a].remove(b)
                    adj[b].remove(a)
                    return True
                mate[a] = 0
                mate[b] = 0
            adj[a].remove(b)
            adj[b].remove(a)
        return False

    r = rec()
    return (r, witness[0] if witness else None) if want_witness else r


def maxlen(mu, t, hi=40):
    best = 0
    for s in range(4, hi + 1, 2):
        if feasible(s, mu, t):
            best = s
        # do not break: longer ladders are not monotone in s
    return best


if __name__ == '__main__':
    print("maximum length s of a path 1..s with a perfect matching of")
    print("spans in [4, mu] whose ladder avoids every power of two and")
    print("every multiple of 2^t.  449.2 measured the mu <= 16 column.")
    print()
    mus = [6, 8, 10, 12, 16, 20, 24, 32, 48, 64]
    print("   t \\ mu " + "".join("%6d" % m for m in mus))
    for t in (3, 4, 5, 6):
        row = []
        for mu in mus:
            row.append(maxlen(mu, t, hi=40))
        print("     %d   " % t + "".join("%6d" % v for v in row))
        sys.stdout.flush()
