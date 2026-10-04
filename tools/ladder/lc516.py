# -*- coding: utf-8 -*-
"""516: LOCAL CAPTURE at girth 5 - is there a cubic graph of girth 5
with NO 8-cycle?

504.1's capture lemma says [g, 2g-2] contains a power of two, and it
is unique there since the interval has ratio below 2.  At g = 5 it is
8.  507.0 refuted the located target in the OPEN ladder, but 516's
recomputation shows the closing edge restores four 8-cycles, so the
refutation used the two degree-2 vertices and does not touch the
cubic statement.

  (LC5)  every cubic graph of girth 5 has a cycle of length 8.

(LC5) implies Erdos-Gyarfas for girth-5 cubic graphs outright.  Its
subject names no power of two, so it is (N9)-safe.  Search it in the
hard case, which here means: search for the counterexample.

SEARCH.  Every HAMILTONIAN cubic graph is a closed ladder, so run the
matching DFS on the cycle 1..n with two monotone prunes - girth >= 5
and no 8-cycle - both of which only get harder as chords are added.
An 8-cycle, once present, is never destroyed, so pruning is exact.
"""

import sys


def cdist(x, y, n):
    d = abs(x - y) % n
    return min(d, n - d)


def has_short(adj, n, cap):
    """is there a cycle of length <= cap through the last-added chord
    region?  full check: BFS-based girth, plus explicit 8-cycle test.
    Returns (girth_at_most_4, has_8)."""
    g4 = False
    for s in adj:
        dist = {s: 0}
        par = {s: None}
        dq = [s]
        while dq:
            nxt = []
            for x in dq:
                if dist[x] >= 2:
                    continue
                for y in adj[x]:
                    if y == par[x]:
                        continue
                    if y in dist:
                        if dist[x] + dist[y] + 1 <= 4:
                            g4 = True
                    else:
                        dist[y] = dist[x] + 1
                        par[y] = x
                        nxt.append(y)
            dq = nxt
        if g4:
            return True, False
    return False, count_len(adj, 8) > 0


def count_len(adj, L):
    """number of cycles of length exactly L (each counted once)."""
    tot = 0
    vs = sorted(adj)
    for s in vs:
        stack = [(s, 1, 1 << s)]
        while stack:
            v, d, msk = stack.pop()
            for w in adj[v]:
                if w == s:
                    if d == L:
                        tot += 1
                elif w > s and d < L and not (msk >> w & 1):
                    stack.append((w, d + 1, msk | (1 << w)))
    return tot // 2


def build(chords, n, partial=False):
    adj = dict((v, []) for v in range(1, n + 1))
    for i in range(1, n + 1):
        j = i % n + 1
        adj[i].append(j)
        adj[j].append(i)
    for a, b in chords:
        adj[a].append(b)
        adj[b].append(a)
    return adj


def search(n, report_all=False):
    mate = [0] * (n + 1)
    chords = []
    hits = []
    nodes = [0]

    def rec():
        nodes[0] += 1
        a = 0
        for i in range(1, n + 1):
            if not mate[i]:
                a = i
                break
        if not a:
            hits.append(list(chords))
            return
        for b in range(a + 1, n + 1):
            if mate[b]:
                continue
            if cdist(a, b, n) < 4:
                continue
            bad = False
            for c, d in chords:
                for x, y in ((c, d), (d, c)):
                    if cdist(a, x, n) + cdist(b, y, n) + 2 <= 4:
                        bad = True
            if bad:
                continue
            mate[a] = b
            mate[b] = a
            chords.append((a, b))
            # exact monotone prune: an 8-cycle now is an 8-cycle later
            if count_len(build(chords, n), 8) == 0:
                rec()
            chords.pop()
            mate[a] = 0
            mate[b] = 0

    rec()
    print("n=%2d  Hamiltonian cubic girth-5 graphs with NO 8-cycle: "
          "%d   (search nodes %d)" % (n, len(hits), nodes[0]))
    for ch in hits[:3 if not report_all else len(hits)]:
        print("    ", ch)
    sys.stdout.flush()
    return hits


if __name__ == '__main__':
    for nn in (10, 12, 14, 16, 18, 20):
        search(nn)
