# -*- coding: utf-8 -*-
"""523: independent verification of the n = 36 witness, and an audit
of the off-by-one in g67.py's girth prune.

THE PRUNE BUG.  g67.py bans a-b paths of length L for L in
range(2, gmin), and a path of length L plus the chord closes a cycle
of length L+1.  So it bans cycles of length 3..gmin, i.e. it enforces
girth >= gmin+1, NOT girth >= gmin.  Every "girth >= 6" run was a
girth >= 7 run and every "girth >= 7" run was a girth >= 8 run.

Nothing here is taken on trust: the witness is rebuilt from its chord
list and its girth, regularity and complete cycle spectrum are
recomputed by brute force.
"""

import sys

W1 = [(1, 7), (2, 11), (3, 15), (4, 26), (5, 32), (6, 20), (8, 30),
      (9, 18), (10, 25), (12, 34), (13, 19), (14, 28), (16, 31),
      (17, 23), (21, 27), (22, 36), (24, 33), (29, 35)]
W2 = [(1, 7), (2, 11), (3, 15), (4, 30), (5, 18), (6, 24), (8, 33),
      (9, 19), (10, 27), (12, 21), (13, 25), (14, 32), (16, 35),
      (17, 26), (20, 31), (22, 36), (23, 29), (28, 34)]


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


def spectrum(adj):
    """complete cycle spectrum by brute-force DFS; each cycle once."""
    got = {}
    vs = sorted(adj)
    rank = dict((v, i) for i, v in enumerate(vs))
    for s in vs:
        stack = [(s, (s,), 1 << s)]
        while stack:
            v, path, msk = stack.pop()
            for w in adj[v]:
                if w == s:
                    if len(path) >= 3 and rank[path[1]] < rank[path[-1]]:
                        got[len(path)] = got.get(len(path), 0) + 1
                elif rank[w] > rank[s] and not (msk >> w & 1):
                    stack.append((w, path + (w,), msk | (1 << w)))
    return got


def bfs_girth(adj):
    """girth by BFS from every edge - a second, independent method."""
    best = 10 ** 9
    for u in adj:
        for v in adj[u]:
            if v < u:
                continue
            dist = {u: 0}
            dq = [u]
            found = None
            while dq and found is None:
                nxt = []
                for x in dq:
                    for y in adj[x]:
                        if x == u and y == v:
                            continue
                        if y == v:
                            found = dist[x] + 1
                            break
                        if y not in dist:
                            dist[y] = dist[x] + 1
                            nxt.append(y)
                    if found is not None:
                        break
                dq = nxt
            if found is not None and found + 1 < best:
                best = found + 1
    return best


def check(name, chords, n=36):
    adj = build(chords, n)
    degs = sorted(set(len(adj[v]) for v in adj))
    simple = all(len(set(adj[v])) == len(adj[v]) for v in adj)
    sp = spectrum(adj)
    g1 = min(sp)
    g2 = bfs_girth(adj)
    print("%s  n=%d  degrees=%s  simple=%s" % (name, n, degs, simple))
    print("   girth by cycle enumeration = %d ;  girth by edge-BFS "
          "= %d   %s" % (g1, g2, "AGREE" if g1 == g2 else "*** DIFFER"))
    print("   cycle counts at 3..12 :",
          [(L, sp.get(L, 0)) for L in range(3, 13)])
    ks = sorted(L for L in sp if L & (L - 1) == 0)
    print("   powers of two present :", ks)
    print("   8-cycles: %d      16-cycles: %d      32-cycles: %d"
          % (sp.get(8, 0), sp.get(16, 0), sp.get(32, 0)))
    print("   VERDICT:", end=" ")
    if g1 == 7 and sp.get(8, 0) == 0:
        print("*** (E_3) COUNTEREXAMPLE AT GIRTH 7 ***")
    elif sp.get(8, 0) == 0 and g1 >= 9:
        print("girth %d, no 8-cycle - OUTSIDE band 3, not a "
              "counterexample to (E_3)" % g1)
    elif sp.get(8, 0) == 0:
        print("girth %d with no 8-cycle" % g1)
    else:
        print("has an 8-cycle; not a witness")
    print()
    sys.stdout.flush()


if __name__ == '__main__':
    check("witness 1", W1)
    check("witness 2", W2)
