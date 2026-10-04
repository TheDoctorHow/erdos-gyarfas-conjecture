# -*- coding: utf-8 -*-
"""516: ALL cubic graphs of girth 5 with no 8-cycle - not only the
Hamiltonian ones.

Both prunes are exact because both are MONOTONE: adding edges never
destroys a cycle, so a partial graph containing a cycle of length <= 4
or a cycle of length 8 can never be repaired.  The search therefore
proves absence, not merely fails to find.

Construction order is canonical in the weak sense used throughout the
archive: always saturate the least-indexed unsaturated vertex, and
join it only to larger indices or to already-present vertices, so
every isomorphism class is reached (many times over, which costs time
but never correctness).
"""

import sys


def cycles_upto(adj, cap):
    """lengths of all cycles of length <= cap present so far."""
    got = set()
    for s in sorted(adj):
        stack = [(s, 1, 1 << s)]
        while stack:
            v, d, msk = stack.pop()
            for w in adj[v]:
                if w == s:
                    if d >= 3:
                        got.add(d)
                elif w > s and d < cap and not (msk >> w & 1):
                    stack.append((w, d + 1, msk | (1 << w)))
    return got


def bad(adj):
    g = cycles_upto(adj, 8)
    return any(L <= 4 for L in g) or (8 in g)


def search(n, verbose=False):
    adj = dict((v, []) for v in range(n))
    hits = []
    nodes = [0]

    def rec():
        nodes[0] += 1
        v = -1
        for u in range(n):
            if len(adj[u]) < 3:
                v = u
                break
        if v < 0:
            # connected?
            seen, dq = set([0]), [0]
            while dq:
                x = dq.pop()
                for y in adj[x]:
                    if y not in seen:
                        seen.add(y)
                        dq.append(y)
            if len(seen) == n:
                hits.append(dict((u, list(adj[u])) for u in adj))
            return
        cand = [u for u in range(v + 1, n)
                if len(adj[u]) < 3 and u not in adj[v]]
        # the deficiency of v must be fillable
        need = 3 - len(adj[v])
        if len(cand) < need:
            return
        from itertools import combinations
        for pick in combinations(cand, need):
            for u in pick:
                adj[v].append(u)
                adj[u].append(v)
            if not bad(adj):
                rec()
            for u in pick:
                adj[v].remove(u)
                adj[u].remove(v)

    rec()
    print("n=%2d  cubic girth->=5 graphs with NO 8-cycle: %d  "
          "(labelled hits; search nodes %d)" % (n, len(hits), nodes[0]))
    if hits and verbose:
        h = hits[0]
        print("     example:", dict((k, sorted(vv))
                                    for k, vv in sorted(h.items())))
    sys.stdout.flush()
    return hits


if __name__ == '__main__':
    for nn in (10, 12, 14, 16):
        search(nn, verbose=True)
