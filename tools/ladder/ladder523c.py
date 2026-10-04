# -*- coding: utf-8 -*-
"""523: THE LADDER QUESTION of 449.4, at span caps the archive never
reached.  Faithful re-implementation of tools/ladder450.py.

MODEL (read from the archive's own tool, not guessed).  Positions
1, 2, 3, ... on a path; going left to right, the least undecided
position p must be matched to some q in [p+4, p+mu]; the recorded
"ladder length" is the largest prefix p-1 the construction reaches.
Forbidden cycle lengths: every multiple of 2^t and every power of
two.  A cycle through a new chord (p,q) is a p-q path plus that
chord, so its length is (path length) + 1.

449.2 reports, with every position matched:
    t = 3:  mu = 6 -> 2,  mu = 8 -> 4,  mu = 10, 12, 16 -> 9
    t = 4:  mu = 6 -> 2,  mu = 8 -> 4,  mu = 10, 12 -> 9,  mu = 16 -> 10
and 449.4 states the gap plainly: "the measurements cover span caps
mu <= 16 ... in the large-shift regime mu grows like sqrt(L), so for
large L the bound is EXTRAPOLATED, not proven."  This reruns the
archive's column as a control and then extends it.
"""

import sys
import time


def forbidden_set(t, cap):
    bad = set()
    M = 1 << t
    for v in range(M, cap + 1, M):
        bad.add(v)
    p = 2
    while p <= cap:
        bad.add(p)
        p *= 2
    return bad


def search(t, mu, smax, cap=None):
    if cap is None:
        cap = 6 * mu + 8
    bad = forbidden_set(t, cap)
    adj = dict((z, []) for z in range(1, smax + 3))
    for z in range(1, smax + 2):
        adj[z].append(z + 1)
        adj[z + 1].append(z)
    best = [0, None]
    chords = []
    matched = [False] * (smax + 3)

    def ok_new(p, q):
        """no cycle through (p,q) has a forbidden length; iterative
        depth-bounded DFS for p-q paths, cycle length = path + 1."""
        stack = [(p, 1 << p, 0)]
        while stack:
            v, msk, ln = stack.pop()
            if ln >= cap:
                continue
            for w in adj[v]:
                if v == p and w == q:
                    continue
                if w == q:
                    if ln + 2 in bad:
                        return False
                elif not (msk >> w & 1):
                    stack.append((w, msk | (1 << w), ln + 1))
        return True

    def rec(p):
        if p - 1 > best[0]:
            best[0] = p - 1
            best[1] = list(chords)
        if p > smax:
            return True
        if matched[p]:
            return rec(p + 1)
        for q in range(p + 4, min(p + mu, smax) + 1):
            if matched[q]:
                continue
            if not ok_new(p, q):
                continue
            adj[p].append(q)
            adj[q].append(p)
            matched[p] = matched[q] = True
            chords.append((p, q))
            if rec(p + 1):
                return True
            chords.pop()
            matched[p] = matched[q] = False
            adj[p].remove(q)
            adj[q].remove(p)
        return False

    unb = rec(1)
    return best[0], best[1], unb


if __name__ == '__main__':
    print("LADDER QUESTION (449.4): longest prefix of a self-chorded")
    print("run with spans in [4, mu] avoiding every multiple of 2^t")
    print("and every power of two.  smax = 60 throughout.")
    print()
    print("  t  mu    max s   unbounded?   time")
    args = sys.argv[1:]
    ts = [int(x) for x in args[0].split(',')] if args else [3, 4]
    mus = ([int(x) for x in args[1].split(',')] if len(args) > 1
           else [6, 8, 10, 12, 16])
    for t in ts:
        for mu in mus:
            t0 = time.time()
            s, wit, unb = search(t, mu, smax=60)
            print("  %d  %3d   %5s   %-10s  %.0fs"
                  % (t, mu, (">=%d" % s) if unb else str(s),
                     "YES" if unb else "no", time.time() - t0))
            sys.stdout.flush()
