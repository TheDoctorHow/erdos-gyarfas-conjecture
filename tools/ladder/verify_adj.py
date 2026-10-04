# -*- coding: utf-8 -*-
"""INDEPENDENT verification of the (ADJ) counterexample.

Nothing here uses the virtual-chord rule.  Everything is recomputed
from the explicit edge list: completions by parity of endpoints below
each gap, connectivity by union-find, blocks by sorting.
"""

from itertools import combinations

W = [(1, 9), (2, 12), (3, 15), (4, 22), (5, 14), (6, 17), (7, 20),
     (8, 18), (10, 16), (11, 19), (13, 21)]
n = 22
POW = set(1 << i for i in range(1, 24))

# ---- 0. it is a perfect matching of [1, 2N]
pos = sorted([x for c in W for x in c])
assert pos == list(range(1, n + 1)), pos
print("perfect matching of [1,%d]: ok" % n)

# ---- 1. hard case
def overlap(c, d):
    (a, b), (p, q) = c, d
    return not (b < p or q < a)

bad = []
for a, b in W:
    s = b - a
    if s < 4 or (s + 1) in POW:
        bad.append(("span", (a, b), s))
for c, d in combinations(W, 2):
    if overlap(c, d):
        v = abs(c[0] - d[0]) + abs(c[1] - d[1])
        if v < 3 or (v + 2) in POW:
            bad.append(("pair", c, d, v))
print("hard-case violations:", bad)
print("spans:", sorted(b - a for a, b in W))

# ---- 2. completion from scratch
def completion(T):
    ends = sorted([x for c in T for x in c])
    edges = [tuple(c) for c in T]
    cnt, idx = 0, 0
    for j in range(1, n):
        while idx < len(ends) and ends[idx] <= j:
            idx += 1
            cnt += 1
        if cnt % 2 == 1:
            edges.append((j, j + 1))
    deg = {}
    for a, b in edges:
        deg[a] = deg.get(a, 0) + 1
        deg[b] = deg.get(b, 0) + 1
    assert all(v % 2 == 0 for v in deg.values())
    par = {v: v for v in deg}

    def find(v):
        while par[v] != v:
            par[v] = par[par[v]]
            v = par[v]
        return v

    for a, b in edges:
        ra, rb = find(a), find(b)
        if ra != rb:
            par[ra] = rb
    return len(edges), len(set(find(v) for v in deg)) <= 1

# ---- 3. (ADJ): any connected triple with an adjacent block and
#         power-of-two weight?
adj_hits, tri_hits = [], []
for T in combinations(W, 3):
    w, conn = completion(list(T))
    if not conn:
        continue
    q = sorted([x for c in T for x in c])
    blocks = [(q[0], q[1]), (q[2], q[3]), (q[4], q[5])]
    if w in POW:
        tri_hits.append((T, w, blocks))
        if any(b[1] - b[0] == 1 for b in blocks):
            adj_hits.append((T, w, blocks))
print()
print("(ADJ) witnesses (connected triple, adjacent block, 2^k):",
      len(adj_hits))
print("(TRI) witnesses (connected triple, 2^k):", len(tri_hits))
for t in tri_hits:
    print("    ", t)

# ---- 4. sanity: adjacent-block triples exist at all, and none is 2^k
cnt_adj = 0
weights = set()
for T in combinations(W, 3):
    w, conn = completion(list(T))
    if not conn:
        continue
    q = sorted([x for c in T for x in c])
    blocks = [(q[0], q[1]), (q[2], q[3]), (q[4], q[5])]
    if any(b[1] - b[0] == 1 for b in blocks):
        cnt_adj += 1
        weights.add(w)
print()
print("connected adjacent-block triples:", cnt_adj)
print("their weight set:", sorted(weights))

# ---- 5. (LQ3): the full spectrum
S = set()
for r in range(1, len(W) + 1):
    for T in combinations(W, r):
        w, conn = completion(list(T))
        if conn:
            S.add(w)
S = sorted(S)
print()
print("full connected spectrum:", S)
print("powers of two in it:", [x for x in S if x in POW])
print("gapless?", len(S) == S[-1] - S[0] + 1)
