"""bic300.py -- is the last tau=1 residue the same gap as 222.3's spanning lemma?

222.3: if the even path-differences span C_0, then (SC) forces the BICYCLE DIMENSION
d = dim(C n C^perp) to satisfy d >= mu - 2.  222.4 measured that d >= mu-2 never happens
once mu >= 6.  This computes mu and d for the H's of the (N3) configuration.
"""
import itertools, sys, random, collections
from collections import deque
exec(open('tools/n3search.py').read().split('def heawood()')[0].split('"""',2)[2])

def heawood():
    E = set()
    for i in range(14): E.add(frozenset((i,(i+1)%14)))
    for i in range(0,14,2): E.add(frozenset((i,(i+5)%14)))
    return E

def gf2_rank(rows, ncols):
    rows = [r for r in rows]; rank = 0; piv = []
    for col in range(ncols):
        p = None
        for i in range(rank, len(rows)):
            if (rows[i] >> col) & 1: p = i; break
        if p is None: continue
        rows[rank], rows[p] = rows[p], rows[rank]
        for i in range(len(rows)):
            if i != rank and ((rows[i] >> col) & 1): rows[i] ^= rows[rank]
        piv.append(col); rank += 1
    return rank, rows[:rank]

def cycle_space_basis(E):
    """fundamental cycles w.r.t. a spanning tree, as bitmasks over an edge ordering"""
    El = sorted(E, key=lambda s: sorted(map(str, s)))
    idx = {e: i for i, e in enumerate(El)}
    a = adj(El)
    root = sorted(a, key=str)[0]
    par = {root: None}; pare = {}; order = [root]; q = deque([root])
    tree = set()
    while q:
        u = q.popleft()
        for w in sorted(a[u], key=str):
            if w not in par:
                par[w] = u; pare[w] = frozenset((u, w)); tree.add(frozenset((u, w)))
                order.append(w); q.append(w)
    basis = []
    for e in El:
        if e in tree: continue
        u, v = tuple(e)
        # path u->root and v->root, symmetric difference
        pu = set(); z = u
        while par[z] is not None: pu ^= {pare[z]}; z = par[z]
        pv = set(); z = v
        while par[z] is not None: pv ^= {pare[z]}; z = par[z]
        cyc = (pu ^ pv) | {e}
        m = 0
        for f in cyc: m |= 1 << idx[f]
        basis.append(m)
    return basis, len(El)

def bicycle_dim(E):
    basis, m = cycle_space_basis(E)
    mu = len(basis)
    if mu == 0: return 0, 0
    # bicycle space = C n C^perp ; C^perp = cut space.  x in C is a bicycle iff
    # <x, b> = 0 for every b in the cycle-space basis (orthogonal to all of C).
    # Solve over the coordinates of C: write x = sum a_i b_i and require Gram*a = 0.
    G = []
    for i in range(mu):
        row = 0
        for j in range(mu):
            if bin(basis[i] & basis[j]).count('1') % 2: row |= 1 << j
        G.append(row)
    r, _ = gf2_rank(G, mu)
    return mu, mu - r          # nullity of the Gram matrix = bicycle dimension

stats = collections.Counter(); rows = []
E0 = heawood()
El = sorted(E0, key=lambda s: sorted(map(str, s)))
for e1, e2 in itertools.combinations(El, 2):
    if e1 & e2: continue
    H = set(E0) - {e1, e2}
    aH = adj(H)
    if not two_conn(aH) or girth(aH) < 6: continue
    cH = bip(aH)
    if cH is None: continue
    d2 = [v for v in aH if len(aH[v]) == 2]
    if len(d2) != 4: continue
    if len([v for v in d2 if cH[v] == 0]) != 2: continue
    mu, d = bicycle_dim(H)
    stats['n=%d mu=%d d=%d' % (len(aH), mu, d)] += 1
    stats['d >= mu-2 : %s' % (d >= mu-2)] += 1
for k in sorted(stats): print('  %-30s %s' % (k, stats[k]))
