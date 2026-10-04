# -*- coding: utf-8 -*-
"""517: hunt the (LC) counterexample where the first moment says it is.

FIRST MOMENT.  In a uniform random cubic graph the number of k-cycles
is asymptotically Poisson with mean (d-1)^k / (2k) = 2^k / (2k).  At
k = 8 that mean is 256/16 = 16, an ABSOLUTE CONSTANT, so

     P(no 8-cycle)  ->  e^{-16}  =  1.1 x 10^{-7} ,

and conditioning on girth 5 (c_3 = c_4 = 0) does not change it, short
cycle counts being asymptotically independent.  Connected cubic
graphs of girth >= 5 number 1, 2, 9, 49, 455, 5783, 90938, 1620479,
31478584 at n = 10, 12, ..., 26.  Multiplying, the expected number of
(LC) counterexamples is below 1 up to n = 24 and is about 3.5 at
n = 26.  So the searches at n <= 22 confirmed nothing, and n = 26 is
where the object should first appear.

SEARCH.  Every Hamiltonian cubic graph is a closed ladder.  DFS over
perfect matchings of the cycle 1..n with two prunes, both monotone
because adding an edge never destroys a cycle:
   - no cycle of length <= 4;
   - no cycle of length 8.
Absence is therefore PROVED on the searched class, not merely
unobserved.  The check is incremental: a cycle created by the new
chord must use it, so it is enough to look for paths of bounded
length between its ends.
"""

import sys


def search(n, want=3, verbose=True):
    adj = dict((v, []) for v in range(1, n + 1))
    for i in range(1, n + 1):
        j = i % n + 1
        adj[i].append(j)
        adj[j].append(i)
    mate = [0] * (n + 1)
    chords = []
    hits = []
    nodes = [0]

    def path_exists(a, b, L, banned):
        """is there a path a..b with exactly L edges avoiding the
        edge set `banned` (given as a set of frozensets)?"""
        stack = [(a, 0, 1 << a)]
        while stack:
            v, d, msk = stack.pop()
            for w in adj[v]:
                e = frozenset((v, w))
                if e in banned:
                    continue
                if w == b:
                    if d + 1 == L:
                        return True
                    continue
                if d + 1 < L and not (msk >> w & 1):
                    stack.append((w, d + 1, msk | (1 << w)))
        return False

    def cdist(x, y):
        d = abs(x - y) % n
        return min(d, n - d)

    def ok(a, b):
        """after adding chord (a,b): no new cycle of length <= 4 and
        no new 8-cycle.  The caller guarantees cdist(a,b) >= 4, so
        (a,b) is not a repeated edge and the ban removes one copy.
        A new cycle uses (a,b), so it is an a-b path of length 2, 3
        (short) or 7, avoiding the chord itself."""
        ban = set([frozenset((a, b))])
        for L in (2, 3):
            if path_exists(a, b, L, ban):
                return False
        if path_exists(a, b, 7, ban):
            return False
        return True

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
            if mate[b] or cdist(a, b) < 4:
                continue
            adj[a].append(b)
            adj[b].append(a)
            if ok(a, b):
                mate[a] = b
                mate[b] = a
                chords.append((a, b))
                rec()
                chords.pop()
                mate[a] = 0
                mate[b] = 0
            adj[a].remove(b)
            adj[b].remove(a)
            if len(hits) >= want:
                return

    rec()
    print("n=%2d  Hamiltonian cubic girth-5 graphs with NO 8-cycle: "
          "%s%d   (nodes %d)"
          % (n, ">=" if len(hits) >= want else "", len(hits), nodes[0]))
    for ch in hits:
        print("    chords:", ch)
    sys.stdout.flush()
    return hits


if __name__ == '__main__':
    import sys as S
    lo = int(S.argv[1]) if len(S.argv) > 1 else 10
    hi = int(S.argv[2]) if len(S.argv) > 2 else 24
    for nn in range(lo, hi + 1, 2):
        search(nn)
