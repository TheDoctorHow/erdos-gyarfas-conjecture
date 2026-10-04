# -*- coding: utf-8 -*-
"""521: THE PENTAGRAM MATCHING LEMMA.

Let G be cubic of girth 5, C = v_1..v_5 a pentagon, u_1..u_5 the five
third-neighbours (distinct by girth 5).  Girth 5 already forces
u_i ~ u_j only when v_i, v_j are at distance TWO along C, so the
possible adjacencies form the pentagram 1-3-5-2-4-1.

  LEMMA.  If G has no 8-cycle then the adjacencies among u_1..u_5
  form a MATCHING in that pentagram - no u_i has two u-neighbours -
  so at most two of the five occur.

  PROOF.  Suppose u_1 ~ u_3 and u_1 ~ u_4, the two C-distance-2 pairs
  at v_1.  Then C' = v_1 u_1 u_3 v_3 v_2 and C'' = v_1 u_1 u_4 v_4 v_5
  are pentagons, and C' meets C'' in exactly the edge v_1 u_1: their
  other vertices are v_2, v_3, u_3 and v_4, v_5, u_4, all distinct.
  Two cycles meeting in exactly one path of length 1 have a symmetric
  difference which is a cycle, here of length 5 + 5 - 2 = 8.  QED

  COROLLARY.  The Petersen graph, whose five u's form a 5-cycle in
  the pentagram, must have an 8-cycle - and does.

VERIFIED here by the contrapositive: over cubic girth-5 graphs, every
pentagon whose u-adjacencies fail to be a matching is checked to lie
in a graph with an 8-cycle, and the 8-cycle is exhibited explicitly
as the symmetric difference C' delta C''.
"""

import sys
from itertools import combinations


def all_cycles_len(adj, L):
    out = []
    vs = sorted(adj)
    rank = dict((v, i) for i, v in enumerate(vs))
    for s in vs:
        stack = [(s, (s,), 1 << s)]
        while stack:
            v, path, msk = stack.pop()
            for w in adj[v]:
                if w == s:
                    if len(path) == L and rank[path[1]] < rank[path[-1]]:
                        out.append(path)
                elif (rank[w] > rank[s] and not (msk >> w & 1)
                      and len(path) < L):
                    stack.append((w, path + (w,), msk | (1 << w)))
    return out


def is_cycle(adj, verts, edges):
    """do `edges` form a single cycle on `verts`?"""
    deg = {}
    nb = {}
    for a, b in edges:
        deg[a] = deg.get(a, 0) + 1
        deg[b] = deg.get(b, 0) + 1
        nb.setdefault(a, []).append(b)
        nb.setdefault(b, []).append(a)
    if set(deg) != set(verts) or any(d != 2 for d in deg.values()):
        return False
    s = next(iter(verts))
    seen, cur, prev = {s}, s, None
    while True:
        nxt = [x for x in nb[cur] if x != prev]
        if not nxt:
            return False
        prev, cur = cur, nxt[0]
        if cur == s:
            break
        if cur in seen:
            return False
        seen.add(cur)
    return seen == set(verts)


def check(name, adj):
    pents = all_cycles_len(adj, 5)
    has8 = len(all_cycles_len(adj, 8)) > 0
    viol = 0
    exhibited = 0
    for P in pents:
        C = list(P)
        Cs = set(C)
        us = []
        good = True
        for v in C:
            o = [w for w in adj[v] if w not in Cs]
            if len(o) != 1:
                good = False
                break
            us.append(o[0])
        if not good or len(set(us)) != 5:
            continue
        # pentagram edges: index pairs at C-distance 2
        pg = [(i, (i + 2) % 5) for i in range(5)]
        deg = dict((i, 0) for i in range(5))
        present = []
        for i, j in pg:
            if us[j] in adj[us[i]]:
                deg[i] += 1
                deg[j] += 1
                present.append((i, j))
        if max(deg.values()) < 2:
            continue
        viol += 1
        # exhibit the 8-cycle for one offending vertex
        i = max(deg, key=lambda t: deg[t])
        js = [b if a == i else a for a, b in present if i in (a, b)]
        j, k = js[0], js[1]

        def pent(i, j):
            """the pentagon v_i u_i u_j v_j v_x where v_x is the
            common C-neighbour of v_i and v_j"""
            x = [t for t in range(5)
                 if (t - i) % 5 in (1, 4) and (t - j) % 5 in (1, 4)][0]
            verts = [C[i], us[i], us[j], C[j], C[x]]
            edges = [(C[i], us[i]), (us[i], us[j]), (us[j], C[j]),
                     (C[j], C[x]), (C[x], C[i])]
            return verts, [tuple(sorted(e)) for e in edges]

        v1, e1 = pent(i, j)
        v2, e2 = pent(i, k)
        shared = set(e1) & set(e2)
        sym = sorted(set(e1) ^ set(e2))
        vv = set()
        for a, b in sym:
            vv.add(a)
            vv.add(b)
        if len(shared) == 1 and is_cycle(adj, vv, sym) and len(sym) == 8:
            exhibited += 1
    print("%-16s pentagons=%3d  8-cycles=%s  pentagons whose "
          "u-adjacencies are NOT a matching=%3d   symmetric "
          "difference is an 8-cycle in %d of them"
          % (name, len(pents), "yes" if has8 else "NO", viol,
             exhibited))
    sys.stdout.flush()
    return viol == exhibited


if __name__ == '__main__':
    from cage516 import lcf, petersen
    from ear521 import rand_cubic_girth5
    ok = True
    ok &= check("Petersen", petersen())
    ok &= check("dodecahedron",
                lcf([10, 7, 4, -4, -7, 10, -4, 7, -7, 4], 2))
    for k, nn in enumerate((14, 18, 20, 22, 24, 26, 28, 30)):
        g = rand_cubic_girth5(nn, 500 + k)
        if g:
            ok &= check("random n=%d" % nn, g)
    print("lemma's 8-cycle exhibited in every violating case:", ok)
