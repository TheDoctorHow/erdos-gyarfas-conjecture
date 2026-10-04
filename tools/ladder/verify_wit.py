"""Independent check of a witness: girth, and the FULL length spectrum,
computed by enumerating the cycle space over GF(2) -- no path search, no
shared code with pushL/depth2.

Graph = host cycle 0..n-1 plus the given chords.  Spanning tree = host
minus the edge (n-1,0).  Cycle space has dimension m - n + 1; every
element is a symmetric difference of fundamental cycles; an element is a
simple cycle iff every vertex it touches has degree 2 in it and the
touched vertices are connected.
"""
import sys
from itertools import combinations

def spectrum(n, chords):
    edges = [(v, (v + 1) % n) for v in range(n)] + [tuple(sorted(c))
                                                    for c in chords]
    edges = [tuple(sorted(e)) for e in edges]
    m = len(edges)
    idx = {e: i for i, e in enumerate(edges)}
    assert len(idx) == m, "duplicate edge"
    # spanning tree: host path 0-1-...-(n-1)
    tree = set(idx[tuple(sorted((v, v + 1)))] for v in range(n - 1))
    par = {0: None}
    depth = {0: 0}
    for v in range(1, n):
        par[v] = v - 1
        depth[v] = v

    def tree_path_edges(a, b):
        s = set()
        x, y = a, b
        while x != y:
            if depth[x] > depth[y]:
                s ^= {idx[tuple(sorted((x, par[x])))]}
                x = par[x]
            else:
                s ^= {idx[tuple(sorted((y, par[y])))]}
                y = par[y]
        return s

    fund = []
    for e in edges:
        if idx[e] in tree:
            continue
        s = tree_path_edges(e[0], e[1]) ^ {idx[e]}
        fund.append(frozenset(s))
    lengths = set()
    d = len(fund)
    assert d == m - n + 1, (d, m, n)
    for r in range(1, d + 1):
        for comb in combinations(range(d), r):
            s = set()
            for i in comb:
                s ^= fund[i]
            if not s:
                continue
            deg = {}
            for i in s:
                a, b = edges[i]
                deg[a] = deg.get(a, 0) + 1
                deg[b] = deg.get(b, 0) + 1
            if any(x != 2 for x in deg.values()):
                continue
            # connected?
            adj = {}
            for i in s:
                a, b = edges[i]
                adj.setdefault(a, []).append(b)
                adj.setdefault(b, []).append(a)
            start = next(iter(adj))
            seen = {start}
            st = [start]
            while st:
                u = st.pop()
                for w in adj[u]:
                    if w not in seen:
                        seen.add(w)
                        st.append(w)
            if len(seen) == len(deg):
                lengths.add(len(s))
    return sorted(lengths)

if __name__ == '__main__':
    n = int(sys.argv[1])
    L = int(sys.argv[2])
    chords = eval(sys.argv[3])
    sp = spectrum(n, chords)
    print("  n =", n, " chords =", len(chords))
    print("  cycle lengths present:", sp)
    print("  girth =", sp[0], " (legal: %s)" % (sp[0] >= 5))
    print("  contains a %d-cycle: %s" % (L, L in sp))
