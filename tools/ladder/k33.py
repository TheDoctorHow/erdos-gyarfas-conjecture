# -*- coding: utf-8 -*-
"""515: the K_{3,3} template.

KEY ALGORITHM.  In a CUBIC graph a subgraph with all degrees 2 or 3 is
exactly G - M for a matching M, and its degree-3 vertices are exactly
the vertices M misses.  So

    K_4-subdivisions      = { G - M : |M| = (n-4)/2, contraction K_4 }
    K_{3,3}-subdivisions  = { G - M : |M| = (n-6)/2, contraction K33 }

This enumerates every one of them exactly, and is completely
independent of the bridge / ear machinery of 509-514.
"""

import sys
from itertools import combinations
from collections import deque


def edges_of(adj):
    return sorted(set(tuple(sorted((v, w))) for v in adj for w in adj[v]))


def all_cycles(adj):
    out = []
    vs = sorted(adj)
    rank = dict((v, i) for i, v in enumerate(vs))
    for s in vs:
        stack = [(s, (s,), frozenset([s]))]
        while stack:
            v, path, seen = stack.pop()
            for w in adj[v]:
                if w == s and len(path) >= 3:
                    if rank[path[1]] < rank[path[-1]]:
                        out.append(path)
                elif w not in seen and rank[w] > rank[s]:
                    stack.append((w, path + (w,), seen | frozenset([w])))
    return out


def connected(adj, vs):
    vs = [v for v in vs if adj[v]]
    if not vs:
        return False
    s = vs[0]
    seen, dq = set([s]), [s]
    while dq:
        x = dq.pop()
        for y in adj[x]:
            if y not in seen:
                seen.add(y)
                dq.append(y)
    return len(seen) == len(vs)


def contract(adj):
    """adj has degree-3 and degree-2 vertices only; return the
    multigraph on the degree-3 vertices, as a list of edges"""
    br = [v for v in adj if len(adj[v]) == 3]
    ed = []
    for v in br:
        for w0 in adj[v]:
            prev, cur = v, w0
            path = [v]
            while len(adj[cur]) == 2:
                nxt = [x for x in adj[cur] if x != prev][0]
                path.append(cur)
                prev, cur = cur, nxt
            if (v, cur) != (cur, v) or True:
                ed.append((v, cur, len(path)))
    seen = {}
    out = []
    for a, b, L in ed:
        k = tuple(sorted((a, b)))
        seen.setdefault(k, []).append(L)
    for k, Ls in seen.items():
        if k[0] == k[1]:
            return None, None
        for L in Ls[:len(Ls) // 2]:
            out.append((k[0], k[1], L))
    return br, out


def is_k33(br, ed):
    if br is None or len(br) != 6 or len(ed) != 9:
        return False
    if len(set((a, b) for a, b, _ in ed)) != 9:
        return False
    g = dict((v, []) for v in br)
    for a, b, _ in ed:
        g[a].append(b)
        g[b].append(a)
    if any(len(g[v]) != 3 for v in br):
        return False
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


def is_k4(br, ed):
    if br is None or len(br) != 4 or len(ed) != 6:
        return False
    return len(set((a, b) for a, b, _ in ed)) == 6


def subdivisions(adj, k):
    """all K_4- (k=4) or K_{3,3}- (k=6) subdivisions, as edge sets"""
    n = len(adj)
    E = edges_of(adj)
    need = (n - k) // 2
    if (n - k) % 2 or need < 0:
        return []
    out = []
    for M in combinations(E, need):
        vs = [v for e in M for v in e]
        if len(set(vs)) != 2 * need:
            continue
        rem = set(E) - set(M)
        sub = dict((v, []) for v in adj)
        for a, b in rem:
            sub[a].append(b)
            sub[b].append(a)
        if any(len(sub[v]) < 2 for v in sub):
            continue
        if not connected(sub, list(sub)):
            continue
        br, ed = contract(sub)
        if (k == 4 and is_k4(br, ed)) or (k == 6 and is_k33(br, ed)):
            out.append(sub)
    return out


def covered_cycles(adj, k):
    cov = set()
    cnt = 0
    for sub in subdivisions(adj, k):
        cnt += 1
        for c in all_cycles(sub):
            cov.add(frozenset(c))
    return cnt, cov


def report(name, adj):
    cyc = all_cycles(adj)
    keys = set(frozenset(c) for c in cyc)
    n4, c4 = covered_cycles(adj, 4)
    n6, c6 = covered_cycles(adj, 6)
    miss4 = keys - c4
    miss6 = keys - c6
    print("%-16s n=%2d cycles=%4d | K4-subdivs=%5d cover=%4d miss=%3d"
          " | K33-subdivs=%5d cover=%4d miss=%3d"
          % (name, len(adj), len(cyc), n4, len(keys & c4), len(miss4),
             n6, len(keys & c6), len(miss6)))
    if miss4:
        print("     a cycle in NO K_4-subdivision:",
              sorted(list(miss4)[0]))
    if miss6 and n6:
        ex = sorted(min(miss6, key=len))
        print("     a cycle in NO K_33-subdivision (len %d): %s"
              % (len(ex), ex))
    sys.stdout.flush()
    return keys, c4, c6


def petersen():
    out = dict((v, []) for v in range(10))
    for i in range(5):
        for x, y in ((i, (i + 1) % 5), (i, i + 5),
                     (i + 5, (i + 2) % 5 + 5)):
            out[x].append(y)
            out[y].append(x)
    return out


def lcf(shifts, reps):
    n = len(shifts) * reps
    s = shifts * reps
    adj = dict((i, []) for i in range(n))
    for i in range(n):
        for j in ((i + 1) % n, (i + s[i]) % n):
            if j not in adj[i]:
                adj[i].append(j)
            if i not in adj[j]:
                adj[j].append(i)
    return adj


def main():
    sys.setrecursionlimit(40000)
    report("K_4", dict((v, [w for w in range(4) if w != v])
                       for v in range(4)))
    K33 = dict((v, []) for v in range(6))
    for a in range(3):
        for b in range(3, 6):
            K33[a].append(b)
            K33[b].append(a)
    report("K_{3,3}", K33)
    report("Petersen", petersen())
    report("Heawood", lcf([5, -5] * 7, 1))
    report("dodecahedron",
           lcf([10, 7, 4, -4, -7, 10, -4, 7, -7, 4], 2))


main()
