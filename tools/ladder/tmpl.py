# -*- coding: utf-8 -*-
"""515: (a) TEMPLATE COMPLETENESS - the only simple cubic graphs in
which every non-zero element of the cycle space is a cycle are K_4 and
K_{3,3};  (b) the crossing triple of 508.2 is a K_{3,3}-subdivision,
not merely a K_4-subdivision, and its fifteen lengths are the seven of
508.2 plus the host cycle plus seven complements."""

import sys
from itertools import combinations
from collections import deque
from nlutil import enum_matchings, closed_adj, girth
from k33b import all_cycles, ckey, contract, is_target


# ---------- (a) template completeness -------------------------------

def cubic_graphs(n):
    """every labelled simple cubic graph on 0..n-1 (connected only)"""
    adj = dict((v, []) for v in range(n))
    out = []

    def rec():
        v = -1
        for u in range(n):
            if len(adj[u]) < 3:
                v = u
                break
        if v < 0:
            s = 0
            seen, dq = set([0]), [0]
            while dq:
                x = dq.pop()
                for y in adj[x]:
                    if y not in seen:
                        seen.add(y)
                        dq.append(y)
            if len(seen) == n:
                out.append(dict((u, list(adj[u])) for u in adj))
            return
        need = 3 - len(adj[v])
        cand = [u for u in range(v + 1, n)
                if len(adj[u]) < 3 and u not in adj[v]]
        for pick in combinations(cand, need):
            for u in pick:
                adj[v].append(u)
                adj[u].append(v)
            rec()
            for u in pick:
                adj[v].remove(u)
                adj[u].remove(v)

    rec()
    return out


def all_cycle_space_elements_are_cycles(adj):
    n = len(adj)
    E = sorted(set(tuple(sorted((v, w))) for v in adj for w in adj[v]))
    idx = dict((e, i) for i, e in enumerate(E))
    par = {0: None}
    tree = set()
    dq = deque([0])
    while dq:
        v = dq.popleft()
        for w in adj[v]:
            if w not in par:
                par[w] = v
                tree.add(tuple(sorted((v, w))))
                dq.append(w)
    basis = []
    for e in E:
        if e in tree:
            continue
        a, b = e
        pa, pb = [], []
        x = a
        while x is not None:
            pa.append(x)
            x = par[x]
        x = b
        while x is not None:
            pb.append(x)
            x = par[x]
        sa, sb = set(pa), set(pb)
        lca = next(z for z in pa if z in sb)
        m = 1 << idx[e]
        for p in (pa, pb):
            for i in range(len(p) - 1):
                if p[i] == lca:
                    break
                m ^= 1 << idx[tuple(sorted((p[i], p[i + 1])))]
        basis.append(m)
    d = len(basis)
    for mask in range(1, 1 << d):
        cur = 0
        for i in range(d):
            if mask >> i & 1:
                cur ^= basis[i]
        sub = dict((v, []) for v in adj)
        for i, e in enumerate(E):
            if cur >> i & 1:
                sub[e[0]].append(e[1])
                sub[e[1]].append(e[0])
        sup = [v for v in sub if sub[v]]
        if any(len(sub[v]) != 2 for v in sup):
            return False
        s = sup[0]
        seen, dq = set([s]), [s]
        while dq:
            x = dq.pop()
            for y in sub[x]:
                if y not in seen:
                    seen.add(y)
                    dq.append(y)
        if len(seen) != len(sup):
            return False
    return True


def part_a():
    for n in (4, 6, 8, 10):
        gs = cubic_graphs(n)
        good = [g for g in gs if all_cycle_space_elements_are_cycles(g)]
        sig = set()
        for g in good:
            tri = sum(1 for c in all_cycles(g) if len(c) == 3)
            sig.add((len(all_cycles(g)), tri, girth(g)))
        print("n=%2d labelled connected cubic graphs=%7d  "
              "ALL-CYCLE-SPACE=%6d  invariants(cycles,triangles,girth)=%s"
              % (n, len(gs), len(good), sorted(sig)))
        sys.stdout.flush()
        if n == 8:
            break


# ---------- (b) the crossing triple is a K_{3,3} --------------------

def part_b(n, cap):
    lad = trip = ok_k33 = ok_len = ok_sum = 0
    for m in enum_matchings(n):
        A = closed_adj(m, n)
        if girth(A) < 5:
            continue
        lad += 1
        if lad > cap:
            break
        for c1, c2, c3 in combinations(sorted(m), 3):
            a1, b1 = c1
            a2, b2 = c2
            a3, b3 = c3
            if not (a1 < a2 < a3 < b1 < b2 < b3):
                continue
            trip += 1
            sub = dict((v, []) for v in range(1, n + 1))
            for j in range(1, n + 1):
                k = j % n + 1
                sub[j].append(k)
                sub[k].append(j)
            for x, y in (c1, c2, c3):
                sub[x].append(y)
                sub[y].append(x)
            br, ed = contract(sub, sorted(sub))
            if is_target(br, ed, 6):
                ok_k33 += 1
            got = sorted(len(c) for c in all_cycles(sub))
            s1, s2, s3 = b1 - a1, b2 - a2, b3 - a3
            L12 = abs(a1 - a2) + abs(b1 - b2)
            L13 = abs(a1 - a3) + abs(b1 - b3)
            L23 = abs(a2 - a3) + abs(b2 - b3)
            xxx = 3 + s1 - s2 + s3
            pred = sorted([n] +
                          [s + 1 for s in (s1, s2, s3)] +
                          [n + 1 - s for s in (s1, s2, s3)] +
                          [2 + L for L in (L12, L13, L23)] +
                          [n + 2 - L for L in (L12, L13, L23)] +
                          [xxx, n + 6 - xxx])
            if got == pred:
                ok_len += 1
            S = n + 3
            four = sorted(len(c) for c in all_cycles(sub)
                          if True)
            q = [x for x in got]
            if sum(q) == 8 * S:
                ok_sum += 1
    print("2N=%d ladders=%d crossing triples=%d  is-K33=%d  "
          "15-length-formula=%d  sum-of-15=8S:%d"
          % (n, min(lad, cap), trip, ok_k33, ok_len, ok_sum))
    sys.stdout.flush()


if __name__ == '__main__':
    part_a()
    for nn in (14, 16, 18):
        part_b(nn, 60)
