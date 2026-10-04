# -*- coding: utf-8 -*-
"""Is 509.0's Presentation Lemma an "if and only if"?

509.0 claims: Z is one of the seven cycles of some K_4-subdivision IFF
Z has two internally disjoint ears with interleaving attachments.
Its proof takes Z to be a QUADRILATERAL of the K_4.  But Z can be a
TRIANGLE, and then the rest of the K_4 is a Y - three internally
disjoint paths from one vertex off Z to three vertices of Z.  A Y is
much cheaper than two skew ears: any bridge with >= 3 attachments has
one.

This script
  (1) re-verifies, by an algorithm that never uses the ord/residual
      machinery, that (EAR*) fails on the dodecahedron;
  (2) exhibits an explicit K_4-subdivision through the same cycle,
      built from a Y, and checks it edge by edge;
  (3) sweeps every cycle of every closed girth-5 ladder at 2N = 14, 16
      and asks whether Y-or-skew always holds.
"""

import sys
from itertools import combinations
from nlutil import enum_matchings, closed_adj, girth


def lcf(shifts, reps):
    n = len(shifts) * reps
    s = shifts * reps
    ch = set()
    for i in range(n):
        ch.add(tuple(sorted((i + 1, (i + s[i]) % n + 1))))
    return n, sorted(ch)


# ---------- cycles ----------

def all_cycles(adj, n):
    """every simple cycle, as a list of vertices in cyclic order"""
    out = []
    for s in range(1, n + 1):
        stack = [(s, [s], set([s]))]
        while stack:
            v, path, seen = stack.pop()
            for w in adj[v]:
                if w == s and len(path) >= 3:
                    if path[1] < path[-1]:
                        out.append(list(path))
                elif w not in seen and w > s:
                    stack.append((w, path + [w], seen | set([w])))
    return out


# ---------- bridges ----------

def bridges_of(adj, cyc):
    Vz = set(cyc)
    pos = dict((v, i) for i, v in enumerate(cyc))
    L = len(cyc)
    nxt = set()
    for i in range(L):
        nxt.add(frozenset((cyc[i], cyc[(i + 1) % L])))
    br = []
    for u in cyc:
        for w in adj[u]:
            if w in Vz and frozenset((u, w)) not in nxt and u < w:
                br.append(('chord', (u, w), set([u, w])))
    seen = set()
    for v in adj:
        if v in Vz or v in seen:
            continue
        comp, dq = set([v]), [v]
        seen.add(v)
        att = set()
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
    return pos, br


def skew(pos, A, B, L):
    for a1, a2 in combinations(sorted(A, key=lambda v: pos[v]), 2):
        for b1, b2 in combinations(sorted(B, key=lambda v: pos[v]), 2):
            if len(set([a1, a2, b1, b2])) < 4:
                continue
            p = sorted([pos[a1], pos[a2]])
            q = sorted([pos[b1], pos[b2]])
            if (p[0] < q[0] < p[1] < q[1]) or (q[0] < p[0] < q[1] < p[1]):
                return (a1, a2, b1, b2)
    return None


def y_or_skew(adj, cyc):
    pos, br = bridges_of(adj, cyc)
    for kind, body, att in br:
        if len(att) >= 3:
            return ('Y', kind, sorted(att))
    for (k1, b1, a1), (k2, b2, a2) in combinations(br, 2):
        s = skew(pos, a1, a2, len(cyc))
        if s:
            return ('SKEW', s, None)
    return None


# ---------- explicit Y ----------

def find_Y(adj, cyc, comp, att):
    """explicit w and three internally disjoint paths w->three atts"""
    Vz = set(cyc)
    for w in sorted(comp):
        for trio in combinations(sorted(att), 3):
            paths = []
            used = set([w])
            ok = True
            for t in trio:
                pr = bfs_path(adj, w, t, comp, used, Vz)
                if pr is None:
                    ok = False
                    break
                paths.append(pr)
                used |= set(pr[1:-1])
            if ok:
                return w, trio, paths
    return None


def bfs_path(adj, w, t, comp, used, Vz):
    from collections import deque
    par = {w: None}
    dq = deque([w])
    while dq:
        v = dq.popleft()
        for u in adj[v]:
            if u in par:
                continue
            if u == t:
                par[u] = v
                p = [u]
                while p[-1] is not None:
                    p.append(par[p[-1]])
                p.pop()
                return p[::-1]
            if u in Vz or u not in comp or u in used:
                continue
            par[u] = v
            dq.append(u)
    return None


def check_k4(adj, cyc, w, trio, paths):
    """verify Z + Y really is a K_4-subdivision"""
    pos = dict((v, i) for i, v in enumerate(cyc))
    L = len(cyc)
    i, j, k = sorted(pos[t] for t in trio)
    arcs = [cyc[i:j + 1], cyc[j:k + 1], cyc[k:] + cyc[:i + 1]]
    for a in arcs:
        if len(a) < 2:
            return False, "degenerate arc"
    interiors = [set(p[1:-1]) for p in paths] + [set(a[1:-1])
                                                 for a in arcs]
    for A, B in combinations(interiors, 2):
        if A & B:
            return False, "paths share an interior vertex"
    branch = set([w]) + set() if False else set([w]) | set(trio)
    for A in interiors:
        if A & branch:
            return False, "interior meets a branch vertex"
    for p in paths + arcs:
        for x, y in zip(p, p[1:]):
            if y not in adj[x]:
                return False, "non-edge %s-%s" % (x, y)
    return True, (arcs, paths)


# ---------- ears, computed directly in G ----------

def all_ears_direct(adj, cyc):
    Vz = set(cyc)
    ears = []
    for u in cyc:
        outs = [x for x in adj[u] if x not in Vz]
        for start in outs:
            stack = [(start, [u, start], set([u, start]))]
            while stack:
                v, path, seen = stack.pop()
                for x in adj[v]:
                    if x in Vz and x != u:
                        ears.append((u, x, list(path[1:])))
                    elif x not in Vz and x not in seen:
                        stack.append((x, path + [x], seen | set([x])))
    return ears


def ear_star_direct(adj, cyc):
    pos = dict((v, i) for i, v in enumerate(cyc))
    ears = all_ears_direct(adj, cyc)
    for e1, e2 in combinations(ears, 2):
        if set(e1[2]) & set(e2[2]):
            continue
        if len(set([e1[0], e1[1], e2[0], e2[1]])) < 4:
            continue
        p = sorted([pos[e1[0]], pos[e1[1]]])
        q = sorted([pos[e2[0]], pos[e2[1]]])
        if (p[0] < q[0] < p[1] < q[1]) or (q[0] < p[0] < q[1] < p[1]):
            return e1, e2
    return None


def main():
    sys.setrecursionlimit(40000)
    n, ch = lcf([10, 7, 4, -4, -7, 10, -4, 7, -7, 4], 2)
    adj = closed_adj(ch, n)
    print("DODECAHEDRON  n=%d girth=%d" % (n, girth(adj)))
    for (a, b) in [(3, 7), (10, 14), (13, 17)]:
        cyc = list(range(a, b + 1))
        print("-" * 60)
        print("Z =", cyc, " (chord (%d,%d) + host path)" % (a, b))
        es = ear_star_direct(adj, cyc)
        print("  (EAR*) direct search:",
              "HOLDS" if es else "*** FAILS ***")
        pos, br = bridges_of(adj, cyc)
        print("  bridges:", [(k, sorted(at)) for k, _, at in br])
        res = y_or_skew(adj, cyc)
        print("  Y-or-skew:", res[0] if res else "NONE")
        for kind, body, att in br:
            if kind == 'comp' and len(att) >= 3:
                Y = find_Y(adj, cyc, body, att)
                w, trio, paths = Y
                ok, det = check_k4(adj, cyc, w, trio, paths)
                print("  explicit K_4-subdivision: hub w =", w,
                      " feet", trio, " verified:", ok)
                print("    Y legs :", paths)
                print("    Z arcs :", det[0] if ok else det)
                break

    print("=" * 62)
    for N in (14, 16):
        lads = enum_matchings(N, None)
        tot = ycount = skcount = fail = ears_fail = 0
        for m in lads:
            A = closed_adj(m, N)
            if girth(A) < 5:
                continue
            for cyc in all_cycles(A, N):
                tot += 1
                r = y_or_skew(A, cyc)
                if r is None:
                    fail += 1
                elif r[0] == 'Y':
                    ycount += 1
                else:
                    skcount += 1
        print("2N =", N, " cycles:", tot, " settled by a Y:", ycount,
              " by skew bridges only:", skcount, " NEITHER:", fail)


main()
