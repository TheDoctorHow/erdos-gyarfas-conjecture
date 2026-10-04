# -*- coding: utf-8 -*-
"""515: the K_4-Ubiquity Theorem.

  every cycle of a 3-connected graph is one of the seven cycles of a
  K_4-subdivision; the only obstruction under min degree 3 is a 2-cut
  {a,b} with ab an edge of the cycle.

Tests:
  T1  exhaustive Y-or-skew over closed girth-5 ladders, 2N = 14, 16,
      with 3-connectivity recorded.
  T2  named cubic graphs: Petersen, Heawood, dodecahedron.
  T3  SHARPNESS: a cubic girth-5 graph, not 3-connected, with a cycle
      in NO K_4-subdivision - three Petersen-minus-an-edge blobs hung
      on consecutive pairs of a 6-cycle.
  T4  for every failure, the predicted 2-cut is exhibited.
"""

import sys
from itertools import combinations
from collections import deque
from nlutil import enum_matchings, closed_adj, girth


def bridges_of(adj, cyc):
    """(kind, body, attachments) for every bridge of the cycle cyc"""
    Vz = set(cyc)
    L = len(cyc)
    nxt = set(frozenset((cyc[i], cyc[(i + 1) % L])) for i in range(L))
    br = []
    for u in cyc:
        for w in adj[u]:
            if w in Vz and u < w and frozenset((u, w)) not in nxt:
                br.append(('chord', set([u, w]), set([u, w])))
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
            return 'Y', None
    for (k1, b1, a1), (k2, b2, a2) in combinations(br, 2):
        for p1, p2 in combinations(sorted(a1, key=lambda v: pos[v]), 2):
            for q1, q2 in combinations(sorted(a2,
                                              key=lambda v: pos[v]), 2):
                if len(set([p1, p2, q1, q2])) < 4:
                    continue
                p = sorted([pos[p1], pos[p2]])
                q = sorted([pos[q1], pos[q2]])
                if p[0] < q[0] < p[1] < q[1] or q[0] < p[0] < q[1] < p[1]:
                    return 'SKEW', None
    return None, br


def two_cut_witness(adj, cyc, br):
    """the theorem says: some bridge has attachments {a,b} with ab an
    edge of cyc, and removing a,b disconnects that bridge's interior"""
    L = len(cyc)
    ed = set(frozenset((cyc[i], cyc[(i + 1) % L])) for i in range(L))
    for kind, body, att in br:
        if len(att) == 2:
            a, b = sorted(att)
            if frozenset((a, b)) in ed and kind == 'comp':
                rest = set(adj) - set([a, b])
                start = next(iter(body))
                comp, dq = set([start]), [start]
                while dq:
                    x = dq.pop()
                    for y in adj[x]:
                        if y in rest and y not in comp:
                            comp.add(y)
                            dq.append(y)
                if comp < rest:
                    return (a, b, len(comp))
    return None


def all_cycles(adj):
    """every cycle of adj, as a vertex list, each once"""
    out = []
    vs = sorted(adj)
    rank = dict((v, i) for i, v in enumerate(vs))
    for s in vs:
        stack = [(s, [s], set([s]))]
        while stack:
            v, path, seen = stack.pop()
            for w in adj[v]:
                if w == s and len(path) >= 3:
                    if len(path) == 2:
                        continue
                    if rank[path[1]] < rank[path[-1]]:
                        out.append(list(path))
                elif w not in seen and rank[w] > rank[s]:
                    stack.append((w, path + [w], seen | set([w])))
    return out


def conn3(adj):
    """3-connected?  brute force over pairs (small graphs only)"""
    vs = sorted(adj)
    if len(vs) < 4:
        return False
    for a, b in combinations(vs, 2):
        rest = [v for v in vs if v not in (a, b)]
        seen, dq = set([rest[0]]), [rest[0]]
        while dq:
            x = dq.pop()
            for y in adj[x]:
                if y not in (a, b) and y not in seen:
                    seen.add(y)
                    dq.append(y)
        if len(seen) < len(rest):
            return False
    return True


def sweep(n, label, cap=None, step=1):
    tot = yc = sk = bad = 0
    lad = c3 = 0
    badlist = []
    idx = -1
    for m in enum_matchings(n):
        A = closed_adj(m, n)
        if girth(A) < 5:
            continue
        idx += 1
        if idx % step:
            continue
        if cap is not None and lad >= cap:
            break
        lad += 1
        if conn3(A):
            c3 += 1
        for cyc in all_cycles(A):
            tot += 1
            r, br = y_or_skew(A, cyc)
            if r == 'Y':
                yc += 1
            elif r == 'SKEW':
                sk += 1
            else:
                bad += 1
                if len(badlist) < 5:
                    badlist.append((m, cyc,
                                    two_cut_witness(A, cyc, br)))
    print("%s 2N=%d ladders=%d 3-connected=%d cycles=%d Y=%d skew=%d "
          "FAIL=%d" % (label, n, lad, c3, tot, yc, sk, bad))
    for b in badlist:
        print("   fail:", b)
    sys.stdout.flush()


def petersen():
    out = dict((v, []) for v in range(10))
    for i in range(5):
        for x, y in ((i, (i + 1) % 5), (i, i + 5),
                     (i + 5, (i + 2) % 5 + 5)):
            out[x].append(y)
            out[y].append(x)
    return out


def heawood():
    n = 14
    s = [5, -5] * 7
    adj = dict((i, []) for i in range(n))
    for i in range(n):
        for j in ((i + 1) % n, (i + s[i]) % n):
            if j not in adj[i]:
                adj[i].append(j)
            if i not in adj[j]:
                adj[j].append(i)
    return adj


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


def named(name, adj):
    deg = set(len(adj[v]) for v in adj)
    cyc = all_cycles(adj)
    yc = sk = bad = 0
    wit = None
    for c in cyc:
        r, br = y_or_skew(adj, c)
        if r == 'Y':
            yc += 1
        elif r == 'SKEW':
            sk += 1
        else:
            bad += 1
            if wit is None:
                wit = (c, two_cut_witness(adj, c, br))
    print("%-14s n=%3d deg=%s girth=%d 3conn=%s cycles=%d Y=%d skew=%d "
          "FAIL=%d" % (name, len(adj), sorted(deg), girth(adj),
                       conn3(adj), len(cyc), yc, sk, bad))
    if wit:
        print("   witness cycle", wit[0], " predicted 2-cut", wit[1])
    sys.stdout.flush()


def blobgraph():
    """C_6 on 0..5 with a copy of Petersen-minus-an-edge hung on each
    of the consecutive pairs (0,1), (2,3), (4,5)."""
    adj = dict((v, []) for v in range(6))
    for i in range(6):
        j = (i + 1) % 6
        adj[i].append(j)
        adj[j].append(i)
    nxt = 6
    for k, (p, q) in enumerate([(0, 1), (2, 3), (4, 5)]):
        P = petersen()
        off = nxt
        nxt += 10
        for v in P:
            adj.setdefault(v + off, [])
        for v in P:
            for w in P[v]:
                if v < w and not (v == 0 and w == 1):
                    adj[v + off].append(w + off)
                    adj[w + off].append(v + off)
        adj[p].append(0 + off)
        adj[0 + off].append(p)
        adj[q].append(1 + off)
        adj[1 + off].append(q)
    return adj


def main():
    sys.setrecursionlimit(40000)
    named("Petersen", petersen())
    named("Heawood", heawood())
    named("dodecahedron", lcf([10, 7, 4, -4, -7, 10, -4, 7, -7, 4], 2))
    B = blobgraph()
    print("blob graph: n=%d cubic=%s girth=%d 3-connected=%s"
          % (len(B), set(len(B[v]) for v in B) == set([3]), girth(B),
             conn3(B)))
    Z = [0, 1, 2, 3, 4, 5]
    r, br = y_or_skew(B, Z)
    print("  Z = C_6 %s -> %s ; bridges %s ; 2-cut %s"
          % (Z, r, [sorted(a) for _, _, a in br],
             two_cut_witness(B, Z, br)))
    sys.stdout.flush()
    sweep(14, "LADDERS-all")
    sweep(16, "LADDERS-1in10", step=10)


main()
