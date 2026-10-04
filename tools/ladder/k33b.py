# -*- coding: utf-8 -*-
"""515: exact enumeration of EVERY K_4- and K_{3,3}-subdivision of a
cubic graph, and of the cycles they cover.

CORRECTED ALGORITHM.  A subdivision need not span G.  Let O be the set
of vertices it omits.  In G - O every surviving vertex must keep degree
2 or 3, so O must leave min degree >= 2; the vertices of degree 3 in
G - O that the subdivision does not branch at must each lose exactly
one edge, and no two of those edges may share a vertex - so they form a
MATCHING M inside the degree-3 part.  Hence

    subdivisions  =  { (G - O) - M }  over O and matchings M ,

with |D3(G-O)| - 2|M| = 4 for K_4 and = 6 for K_{3,3}.  This enumerates
all of them, spanning or not, and uses no bridge or ear machinery.
"""

import sys
from itertools import combinations
from collections import deque


def edges_of(adj):
    return sorted(set(tuple(sorted((v, w))) for v in adj for w in adj[v]))


def all_cycles(adj):
    out = []
    vs = sorted(v for v in adj if adj[v])
    rank = dict((v, i) for i, v in enumerate(vs))
    for s in vs:
        stack = [(s, (s,), frozenset([s]))]
        while stack:
            v, path, seen = stack.pop()
            for w in adj[v]:
                if w == s and len(path) >= 3:
                    if rank[path[1]] < rank[path[-1]]:
                        out.append(path)
                elif w in rank and w not in seen and rank[w] > rank[s]:
                    stack.append((w, path + (w,), seen | frozenset([w])))
    return out


def ckey(path):
    """canonical key of a cycle given as a vertex sequence"""
    n = len(path)
    best = None
    for d in (1, -1):
        for i in range(n):
            r = tuple(path[(i + d * j) % n] for j in range(n))
            if best is None or r < best:
                best = r
    return best


def contract(sub, sup):
    br = [v for v in sup if len(sub[v]) == 3]
    if not br:
        return None, None
    cnt = {}
    for v in br:
        for w0 in sub[v]:
            prev, cur, L = v, w0, 1
            while len(sub[cur]) == 2:
                nxt = [x for x in sub[cur] if x != prev][0]
                prev, cur, L = cur, nxt, L + 1
            if cur == v:
                return None, None
            cnt.setdefault(tuple(sorted((v, cur))), []).append(L)
    ed = []
    for k, Ls in cnt.items():
        for L in Ls[:len(Ls) // 2]:
            ed.append((k[0], k[1], L))
    return br, ed


def is_target(br, ed, k):
    if br is None or len(br) != k:
        return False
    m = 6 if k == 4 else 9
    if len(ed) != m or len(set((a, b) for a, b, _ in ed)) != m:
        return False
    g = dict((v, []) for v in br)
    for a, b, _ in ed:
        g[a].append(b)
        g[b].append(a)
    if any(len(g[v]) != 3 for v in br):
        return False
    if k == 4:
        return True
    col = {br[0]: 0}
    dq = deque([br[0]])
    while dq:
        v = dq.popleft()
        for w in g[v]:
            if w in col:
                if col[w] == col[v]:
                    return False
            else:
                col[w] = 1 - col[v]
                dq.append(w)
    return len(col) == 6


def matchings(edges, size):
    res = []

    def rec(i, cur, used):
        if len(cur) == size:
            res.append(list(cur))
            return
        if i >= len(edges) or len(edges) - i < size - len(cur):
            return
        a, b = edges[i]
        if a not in used and b not in used:
            cur.append(edges[i])
            rec(i + 1, cur, used | set([a, b]))
            cur.pop()
        rec(i + 1, cur, used)

    rec(0, [], set())
    return res


def subdivisions(adj, k, omit_cap=None):
    V = sorted(adj)
    n = len(V)
    out = []
    cap = n - k if omit_cap is None else omit_cap
    for osz in range(0, cap + 1):
        for O in combinations(V, osz):
            Os = set(O)
            sup = [v for v in V if v not in Os]
            if len(sup) < k:
                continue
            deg = dict((v, [w for w in adj[v] if w not in Os])
                       for v in sup)
            if any(len(deg[v]) < 2 for v in sup):
                continue
            D3 = [v for v in sup if len(deg[v]) == 3]
            need = len(D3) - k
            if need < 0 or need % 2:
                continue
            D3s = set(D3)
            E3 = [e for e in edges_of(deg) if e[0] in D3s and e[1] in D3s]
            for M in matchings(E3, need // 2):
                sub = dict((v, list(deg[v])) for v in sup)
                for a, b in M:
                    sub[a].remove(b)
                    sub[b].remove(a)
                s0 = sup[0]
                seen, dq = set([s0]), [s0]
                while dq:
                    x = dq.pop()
                    for y in sub[x]:
                        if y not in seen:
                            seen.add(y)
                            dq.append(y)
                if len(seen) != len(sup):
                    continue
                br, ed = contract(sub, sup)
                if is_target(br, ed, k):
                    out.append(sub)
    return out


def report(name, adj, do4=True, do6=True, omit_cap=None):
    keys = set(ckey(c) for c in all_cycles(adj))
    line = "%-14s n=%2d cycles=%4d" % (name, len(adj), len(keys))
    for k, flag in ((4, do4), (6, do6)):
        if not flag:
            continue
        subs = subdivisions(adj, k, omit_cap)
        cov = set()
        for s in subs:
            for c in all_cycles(s):
                cov.add(ckey(c))
        miss = keys - cov
        line += " | K%s subdivs=%5d cover=%4d MISS=%3d" % (
            "4" if k == 4 else "33", len(subs), len(keys & cov),
            len(miss))
        if miss:
            ex = sorted(miss, key=len)[0]
            line += " e.g.len%d%s" % (len(ex), list(ex))
    print(line)
    sys.stdout.flush()


def petersen():
    out = dict((v, []) for v in range(10))
    for i in range(5):
        for x, y in ((i, (i + 1) % 5), (i, i + 5),
                     (i + 5, (i + 2) % 5 + 5)):
            out[x].append(y)
            out[y].append(x)
    return out


def k33graph():
    g = dict((v, []) for v in range(6))
    for a in range(3):
        for b in range(3, 6):
            g[a].append(b)
            g[b].append(a)
    return g


def k4graph():
    return dict((v, [w for w in range(4) if w != v]) for v in range(4))


if __name__ == '__main__':
    sys.setrecursionlimit(40000)
    report("K_4", k4graph())
    report("K_{3,3}", k33graph())
    report("Petersen", petersen())
