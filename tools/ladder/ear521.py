# -*- coding: utf-8 -*-
"""521: THE EAR DICHOTOMY AT GIRTH 5 - a genuinely combinatorial
consequence of "no 8-cycle", and its verification.

Let G be cubic of girth 5 and C = v_1..v_5 a pentagon.  Girth 5 makes
the five third-neighbours u_1..u_5 distinct.  An ear of C from u_i to
u_j - a path internally disjoint from C - of length s gives a C-ear
of length q = s + 2 with feet v_i, v_j, and hence TWO cycles, of
lengths q + d and q + 5 - d where d is the distance from v_i to v_j
along C.  If G has no 8-cycle then

    q + d != 8   and   q + 5 - d != 8 ,

which in terms of s = q - 2 reads

    d = 1 :   s not in {1, 2, 5}        (and s >= 3 by girth)
    d = 2 :   s not in {3, 4}           (so s <= 2 or s >= 5)

  THE DICHOTOMY.  In a cubic graph of girth 5 with no 8-cycle, for
  every pentagon, ADJACENT feet have their outside neighbours FAR
  APART (>= 3, and never exactly 5), while feet at distance two have
  them either CLOSE (<= 2) or FAR (>= 5) - never at distance 3 or 4.

VERIFIED here by the contrapositive, which is the only direction a
graph with 8-cycles can test: over many cubic girth-5 graphs, every
pentagon and every pair of feet, check that whenever the forbidden
distance occurs an 8-cycle really is present.
"""

import sys
from itertools import combinations
from collections import deque


def all_cycles_len(adj, L):
    out = []
    vs = sorted(adj)
    rank = dict((v, i) for i, v in enumerate(vs))
    for s in vs:
        stack = [(s, (s,), 1 << s)]
        while stack:
            v, path, msk = stack.pop()
            if len(path) > L:
                continue
            for w in adj[v]:
                if w == s:
                    if len(path) == L and rank[path[1]] < rank[path[-1]]:
                        out.append(path)
                elif (rank[w] > rank[s] and not (msk >> w & 1)
                      and len(path) < L):
                    stack.append((w, path + (w,), msk | (1 << w)))
    return out


def dist_avoiding(adj, a, b, banned):
    """shortest a-b path avoiding the vertex set `banned`"""
    if a == b:
        return 0
    seen = set([a]) | set(banned)
    if a in banned:
        return None
    dq = deque([(a, 0)])
    while dq:
        v, d = dq.popleft()
        for w in adj[v]:
            if w == b:
                return d + 1
            if w not in seen:
                seen.add(w)
                dq.append((w, d + 1))
    return None


def check(name, adj):
    n = len(adj)
    pent = all_cycles_len(adj, 5)
    has8 = len(all_cycles_len(adj, 8)) > 0
    tested = 0
    consistent = 0
    forbidden_seen = 0
    for P in pent:
        C = list(P)
        Cs = set(C)
        us = []
        ok = True
        for v in C:
            out = [w for w in adj[v] if w not in Cs]
            if len(out) != 1:
                ok = False
                break
            us.append(out[0])
        if not ok or len(set(us)) != 5:
            continue
        for i, j in combinations(range(5), 2):
            d = min((j - i) % 5, (i - j) % 5)
            s = dist_avoiding(adj, us[i], us[j], Cs)
            if s is None:
                continue
            tested += 1
            bad = (s in (1, 2, 5)) if d == 1 else (s in (3, 4))
            if bad:
                forbidden_seen += 1
                if has8:
                    consistent += 1
    print("%-16s n=%3d pentagons=%3d  8-cycles=%s  foot-pairs "
          "tested=%5d  forbidden distance seen=%4d  each accompanied "
          "by an 8-cycle: %s"
          % (name, n, len(pent), "yes" if has8 else "NO", tested,
             forbidden_seen,
             "%d/%d" % (consistent, forbidden_seen)))
    sys.stdout.flush()
    return forbidden_seen == consistent


def rand_cubic_girth5(n, seed):
    import random
    rnd = random.Random(seed)
    for _ in range(3000):
        pts = list(range(1, n + 1))
        rnd.shuffle(pts)
        adj = dict((v, []) for v in range(1, n + 1))
        for i in range(1, n + 1):
            j = i % n + 1
            adj[i].append(j)
            adj[j].append(i)
        okk = True
        for k in range(n // 2):
            a, b = pts[2 * k], pts[2 * k + 1]
            adj[a].append(b)
            adj[b].append(a)
        if any(len(adj[v]) != 3 for v in adj):
            continue
        if all_cycles_len(adj, 3) or all_cycles_len(adj, 4):
            continue
        return adj
    return None


if __name__ == '__main__':
    from cage516 import lcf, petersen
    allok = True
    allok &= check("Petersen", petersen())
    allok &= check("dodecahedron",
                   lcf([10, 7, 4, -4, -7, 10, -4, 7, -7, 4], 2))
    for k, nn in enumerate((14, 16, 18, 20, 22, 24)):
        g = rand_cubic_girth5(nn, 1000 + k)
        if g:
            allok &= check("random n=%d" % nn, g)
    print("dichotomy consistent everywhere:", allok)
