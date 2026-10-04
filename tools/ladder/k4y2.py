# -*- coding: utf-8 -*-
"""Lean version: (1) independent re-verification that (EAR*)/(EAR+)
fail on the dodecahedron, by enumerating ACTUAL EAR PATHS in G;
(2) an explicit K_4-subdivision through the same cycle built from a Y;
(3) the Y-or-skew sweep over every cycle of every closed girth-5
ladder at 2N = 14."""

import sys
from itertools import combinations
from collections import deque
from nlutil import enum_matchings, closed_adj, girth, planar_hamiltonian


def lcf(shifts, reps):
    n = len(shifts) * reps
    s = shifts * reps
    ch = set()
    for i in range(n):
        ch.add(tuple(sorted((i + 1, (i + s[i]) % n + 1))))
    return n, sorted(ch)


def dodec_from_faces():
    """independent construction: the 12 pentagonal faces of the
    regular dodecahedron, vertices named by the standard 20-point
    labelling of the icosian graph."""
    faces = [
        (1, 2, 3, 4, 5), (1, 5, 15, 20, 6), (1, 6, 7, 8, 2),
        (2, 8, 9, 10, 3), (3, 10, 11, 12, 4), (4, 12, 13, 14, 5),
        (5, 14, 15, 20, 6)]
    # not used for the proof - kept only as a shape note
    return faces


def all_ears(adj, cyc):
    """every ear of Z: (u, v, interior tuple), interior avoids V(Z)"""
    Vz = set(cyc)
    out = []
    for u in cyc:
        for start in [x for x in adj[u] if x not in Vz]:
            stack = [(start, (start,), frozenset([start]))]
            while stack:
                v, path, seen = stack.pop()
                for x in adj[v]:
                    if x in Vz:
                        if x != u and u < x:
                            out.append((u, x, seen))
                        elif x != u:
                            out.append((x, u, seen))
                    elif x not in seen:
                        stack.append((x, path + (x,), seen | frozenset([x])))
    return out


def ear_star_direct(adj, cyc):
    """two vertex-disjoint ears with interleaving attachments?"""
    pos = dict((v, i) for i, v in enumerate(cyc))
    L = len(cyc)
    byp = {}
    for u, v, itr in all_ears(adj, cyc):
        byp.setdefault((u, v), []).append(itr)
    for q in combinations(sorted(pos, key=lambda z: pos[z]), 4):
        w, x, y, z = q
        for A, B in (((w, y), (x, z)),):
            ka = tuple(sorted(A))
            kb = tuple(sorted(B))
            for ia in byp.get(ka, []):
                for ib in byp.get(kb, []):
                    if not (ia & ib):
                        return (A, sorted(ia)), (B, sorted(ib))
    return None


def bridges_of(adj, cyc):
    Vz = set(cyc)
    L = len(cyc)
    nxt = set(frozenset((cyc[i], cyc[(i + 1) % L])) for i in range(L))
    br = []
    for u in cyc:
        for w in adj[u]:
            if w in Vz and u < w and frozenset((u, w)) not in nxt:
                br.append(('chord', frozenset([u, w]), set([u, w])))
    seen = set()
    for v in adj:
        if v in Vz or v in seen:
            continue
        comp, dq, att = set([v]), [v], set()
        seen.add(v)
        while dq:
            x = dq.pop()
            for y in adj[x]:
                if y in Vz:
                    att.add(y)
                elif y not in comp:
                    comp.add(y)
                    seen.add(y)
                    dq.append(y)
        br.append(('comp', comp, att))
    return br


def y_or_skew(adj, cyc):
    pos = dict((v, i) for i, v in enumerate(cyc))
    br = bridges_of(adj, cyc)
    for kind, body, att in br:
        if len(att) >= 3:
            return 'Y'
    for (k1, b1, a1), (k2, b2, a2) in combinations(br, 2):
        for p1, p2 in combinations(sorted(a1, key=lambda v: pos[v]), 2):
            for q1, q2 in combinations(sorted(a2,
                                              key=lambda v: pos[v]), 2):
                if len(set([p1, p2, q1, q2])) < 4:
                    continue
                p = sorted([pos[p1], pos[p2]])
                q = sorted([pos[q1], pos[q2]])
                if p[0] < q[0] < p[1] < q[1] or q[0] < p[0] < q[1] < p[1]:
                    return 'SKEW'
    return None


def bfs_path(adj, w, t, allow):
    par = {w: None}
    dq = deque([w])
    while dq:
        v = dq.popleft()
        for u in adj[v]:
            if u in par:
                continue
            if u == t:
                par[u] = v
                p, c = [], u
                while c is not None:
                    p.append(c)
                    c = par[c]
                return p[::-1]
            if u in allow:
                par[u] = v
                dq.append(u)
    return None


def explicit_Y(adj, cyc, comp, att):
    Vz = set(cyc)
    for w in sorted(comp):
        for trio in combinations(sorted(att), 3):
            used, legs, ok = set([w]), [], True
            for t in trio:
                allow = comp - used
                p = bfs_path(adj, w, t, allow)
                if p is None:
                    ok = False
                    break
                legs.append(p)
                used |= set(p[1:-1])
            if ok:
                return w, trio, legs
    return None


def verify_k4(adj, cyc, w, trio, legs):
    pos = dict((v, i) for i, v in enumerate(cyc))
    i, j, k = sorted(pos[t] for t in trio)
    arcs = [cyc[i:j + 1], cyc[j:k + 1], cyc[k:] + cyc[:i + 1]]
    parts = legs + arcs
    for p in parts:
        for x, y in zip(p, p[1:]):
            if y not in adj[x]:
                return False, "non-edge %d-%d" % (x, y)
        if len(p) < 2:
            return False, "degenerate branch path"
    branch = set([w]) | set(trio)
    ints = [set(p[1:-1]) for p in parts]
    for A in ints:
        if A & branch:
            return False, "interior hits a branch vertex"
    for A, B in combinations(ints, 2):
        if A & B:
            return False, "two branch paths share a vertex"
    ends = [frozenset([p[0], p[-1]]) for p in parts]
    if len(set(ends)) != 6:
        return False, "branch paths do not realise the six K_4 edges"
    return True, (arcs, legs)


def main():
    sys.setrecursionlimit(40000)
    n, ch = lcf([10, 7, 4, -4, -7, 10, -4, 7, -7, 4], 2)
    adj = closed_adj(ch, n)
    print("DODECAHEDRON n=%d cubic=%s girth=%d planar=%s" %
          (n, set(len(adj[v]) for v in adj) == set([3]), girth(adj),
           planar_hamiltonian(ch)))
    print("  matching:", ch)
    for (a, b) in [(3, 7), (10, 14), (13, 17)]:
        cyc = list(range(a, b + 1))
        print("-" * 62)
        print("Z =", cyc)
        es = ear_star_direct(adj, cyc)
        print("  (EAR*) by direct ear enumeration in G:",
              "HOLDS %s" % (es,) if es else "*** FAILS ***")
        br = bridges_of(adj, cyc)
        print("  bridges:", [(k, sorted(at)) for k, _, at in br])
        print("  Y-or-skew:", y_or_skew(adj, cyc))
        for kind, body, att in br:
            if kind == 'comp' and len(att) >= 3:
                Y = explicit_Y(adj, cyc, body, att)
                ok, det = verify_k4(adj, cyc, *Y)
                print("  K_4-subdivision  hub=%d feet=%s  VERIFIED=%s"
                      % (Y[0], Y[1], ok))
                if ok:
                    print("    Z arcs :", det[0])
                    print("    Y legs :", det[1])
                else:
                    print("    ", det)
                break
    print("=" * 62)
    for N in (14,):
        tot = yc = sk = bad = 0
        for m in enum_matchings(N, None):
            A = closed_adj(m, N)
            if girth(A) < 5:
                continue
            for s in range(1, N + 1):
                stack = [(s, [s], set([s]))]
                while stack:
                    v, path, seen = stack.pop()
                    for wv in A[v]:
                        if wv == s and len(path) >= 5:
                            if path[1] < path[-1]:
                                tot += 1
                                r = y_or_skew(A, list(path))
                                if r == 'Y':
                                    yc += 1
                                elif r == 'SKEW':
                                    sk += 1
                                else:
                                    bad += 1
                        elif wv not in seen and wv > s:
                            stack.append((wv, path + [wv],
                                          seen | set([wv])))
        print("2N =", N, " cycles:", tot, " Y:", yc, " skew-only:", sk,
              " NEITHER (no K_4 through Z):", bad)


main()
