# -*- coding: utf-8 -*-
"""517: full hand-checkable audit of one 2N = 20 residue witness, and
the residue census at 2N = 18, 20, 22 recomputed from scratch.

Prints, for the witness, EVERY power-of-two cycle: its vertex
sequence, the chords it uses, and the crossing pattern of those
chords.  Nothing is taken from kstar2 or shape516.
"""

import sys
from itertools import combinations

W = [(1, 5), (2, 8), (3, 19), (4, 10), (6, 12), (7, 16), (9, 18),
     (11, 15), (13, 17), (14, 20)]


def build(chords, n):
    adj = dict((v, []) for v in range(1, n + 1))
    for i in range(1, n + 1):
        j = i % n + 1
        adj[i].append(j)
        adj[j].append(i)
    for a, b in chords:
        adj[a].append(b)
        adj[b].append(a)
    return adj


def all_cycles(adj):
    out = []
    vs = sorted(adj)
    rank = dict((v, i) for i, v in enumerate(vs))
    for s in vs:
        stack = [(s, (s,), 1 << s)]
        while stack:
            v, path, msk = stack.pop()
            for w in adj[v]:
                if w == s:
                    if len(path) >= 3 and rank[path[1]] < rank[path[-1]]:
                        out.append(path)
                elif rank[w] > rank[s] and not (msk >> w & 1):
                    stack.append((w, path + (w,), msk | (1 << w)))
    return out


def girth(adj):
    return min(len(c) for c in all_cycles(adj))


def chords_used(path, cs):
    used = []
    k = len(path)
    for i in range(k):
        u, v = path[i], path[(i + 1) % k]
        e = (u, v) if u < v else (v, u)
        if e in cs:
            used.append(e)
    return used


def pattern(T):
    ch = sorted(T)
    cr = []
    for (x1, y1), (x2, y2) in combinations(ch, 2):
        cr.append(1 if x1 < x2 < y1 < y2 else 0)
    return sum(cr), cr


def audit_witness(n=20):
    cs = set(tuple(sorted(e)) for e in W)
    adj = build(W, n)
    print("witness chords:", W)
    print("n=%d  cubic=%s  girth=%d  total cycles=%d"
          % (n, all(len(adj[v]) == 3 for v in adj), girth(adj),
             len(all_cycles(adj))))
    P2 = set([4, 8, 16, 32])
    found = 0
    for path in all_cycles(adj):
        if len(path) not in P2:
            continue
        found += 1
        u = chords_used(path, cs)
        if len(u) == 3:
            k, cr = pattern(u)
            tag = ("PAIRWISE CROSSING -> C+triple is K_{3,3}" if k == 3
                   else "only %d of 3 pairs cross -> NOT a K_{3,3}" % k)
        else:
            tag = ""
        print("  length %2d  chords used %d %s  %s"
              % (len(path), len(u), u, tag))
        print("        cycle:", path)
    print("  power-of-two cycles in total:", found)
    sys.stdout.flush()


def census(n):
    """recompute the residue from scratch: over EVERY girth-5 closed
    ladder, not only the pruned survivors."""
    P2 = set()
    v = 4
    while v <= 2 * n:
        P2.add(v)
        v *= 2
    mate = [0] * (n + 1)
    chords = []
    tot = [0]
    hist = {}
    resid = []

    def cdist(x, y):
        d = abs(x - y) % n
        return min(d, n - d)

    def rec():
        a = 0
        for i in range(1, n + 1):
            if not mate[i]:
                a = i
                break
        if not a:
            tot[0] += 1
            ch = list(chords)
            cs = set(tuple(sorted(e)) for e in ch)
            adj = build(ch, n)
            if girth(adj) < 5:
                return
            best = 99
            tmpl = False
            for path in all_cycles(adj):
                if len(path) not in P2:
                    continue
                u = chords_used(path, cs)
                if len(u) < best:
                    best = len(u)
                if len(u) == 3 and pattern(u)[0] == 3:
                    tmpl = True
                if len(u) <= 2:
                    tmpl = True
            hist[best] = hist.get(best, 0) + 1
            if best >= 3 and not tmpl:
                resid.append(ch)
            return
        for b in range(a + 1, n + 1):
            if mate[b] or cdist(a, b) < 4:
                continue
            mate[a] = b
            mate[b] = a
            chords.append((a, b))
            rec()
            chords.pop()
            mate[a] = 0
            mate[b] = 0

    rec()
    print("2N=%2d  girth-5 closed ladders=%6d  k* histogram %s  "
          "RESIDUE (no through-C template) = %d"
          % (n, sum(hist.values()), sorted(hist.items()), len(resid)))
    sys.stdout.flush()
    return resid


if __name__ == '__main__':
    audit_witness()
    print()
    for nn in (18, 20):
        census(nn)
