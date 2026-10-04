# -*- coding: utf-8 -*-
"""516 verification: independent, graph-based recomputation of the
hard516.py columns at 2N = 14, 16, 18 (and the 2N = 18 residue
witnesses), from raw cycle enumeration with chord-usage
classification.  No weight formulas are used anywhere here.
"""

import sys
from itertools import combinations
from hard516 import gen_hard, P2


def closed_graph(chords, n):
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
    """every cycle, as a vertex tuple (rotation/reflection canonical)."""
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


def chord_usage(path, chordset, n):
    """which chords the cycle uses, and whether it uses the host
    closing edge is irrelevant here - host edges are i ~ i%n+1."""
    used = []
    k = len(path)
    for i in range(k):
        u, v = path[i], path[(i + 1) % k]
        e = (u, v) if u < v else (v, u)
        if e in chordset:
            used.append(e)
    return used


def is_crossing_triple(t):
    ch = sorted(t)
    (a1, b1), (a2, b2), (a3, b3) = ch
    return a1 < a2 < a3 < b1 < b2 < b3


def girth_of(adj):
    best = 99
    for c in all_cycles(adj):
        if len(c) < best:
            best = len(c)
    return best


def sweep(n, show_residue=False):
    ladders = []
    gen_hard(n, ladders)
    nb = nc = nd = ne = 0
    residue = []
    for ch in ladders:
        cs = set(tuple(sorted(e)) for e in ch)
        adj = closed_graph(ch, n)
        cyc = all_cycles(adj)
        hit_x = hit_xc = hit_12 = False
        tri_pw = False
        for path in cyc:
            L = len(path)
            used = chord_usage(path, cs, n)
            if L in P2:
                if len(used) <= 2:
                    hit_12 = True
                elif len(used) == 3 and is_crossing_triple(used):
                    # which completion is XXX (no closing edge n-1)?
                    wrap = False
                    k = len(path)
                    for i in range(k):
                        u, v = path[i], path[(i + 1) % k]
                        if set((u, v)) == set((n, 1)):
                            wrap = True
                    hit_xc = True
                    if not wrap:
                        hit_x = True
                if len(used) == 3:
                    tri_pw = True
        # col B: no open XXX power. col C: no crossing-triple power
        # (either completion).  col D: C and no closed <=2-chord power.
        # col E: D and closed girth >= 5.
        if not hit_x:
            nb += 1
        if not hit_xc:
            nc += 1
            if not hit_12:
                nd += 1
                if girth_of(adj) >= 5:
                    ne += 1
                    residue.append((ch, tri_pw))
    print("2N=%2d  B=%4d  C=%4d  D=%4d  E=%4d" % (n, nb, nc, nd, ne))
    if show_residue:
        for ch, tp in residue:
            print("   residue:", ch, " 3-chord p2 cycle exists:", tp)
    sys.stdout.flush()


if __name__ == '__main__':
    sweep(14)
    sweep(16)
    sweep(18, show_residue=True)
