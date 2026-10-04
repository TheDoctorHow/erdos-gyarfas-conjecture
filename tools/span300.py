"""span300.py -- do the x-b_2 path differences span the cycle space of H?

If yes then, H being bipartite (so every cycle is even and C_0 = C), 222.2's isotropy
would make the mod-2 intersection form vanish identically on C, forcing the bicycle
dimension d = mu.  bic300.py computed d in {0,2} against mu = 6.  So spanning + isotropy
is contradictory, and (N3) closes.
"""
import itertools, collections
from collections import deque
exec(open('tools/n3search.py').read().split('def heawood()')[0].split('"""',2)[2])
_b = open('tools/bic300.py').read().split('stats = collections.Counter()')[0]
exec(_b[_b.index('def gf2_rank'):])

def heawood():
    E = set()
    for i in range(14): E.add(frozenset((i,(i+1)%14)))
    for i in range(0,14,2): E.add(frozenset((i,(i+5)%14)))
    return E

def span_of_differences(E, s, t):
    El = sorted(E, key=lambda f: sorted(map(str, f)))
    idx = {e: i for i, e in enumerate(El)}
    a = adj(El)
    P = paths(a, s, t)
    if len(P) < 2: return None, None, len(P)
    def mask(p):
        m = 0
        for i in range(len(p)-1): m |= 1 << idx[frozenset((p[i], p[i+1]))]
        return m
    base = mask(P[0])
    rows = [mask(p) ^ base for p in P[1:]]
    r, _ = gf2_rank(rows, len(El))
    return r, len(P), None

stats = collections.Counter()
E0 = heawood()
El = sorted(E0, key=lambda f: sorted(map(str, f)))
for e1, e2 in itertools.combinations(El, 2):
    if e1 & e2: continue
    H = set(E0) - {e1, e2}
    aH = adj(H)
    if not two_conn(aH) or girth(aH) < 6: continue
    cH = bip(aH)
    if cH is None: continue
    d2 = [v for v in aH if len(aH[v]) == 2]
    if len(d2) != 4: continue
    X = [v for v in d2 if cH[v] == 0]; Y = [v for v in d2 if cH[v] == 1]
    if len(X) != 2 or len(Y) != 2: continue
    mu, d = bicycle_dim(H)
    for x in X:
        for y2 in Y:            # y2 plays the role of b_2 : DIFFERENT class from x
            r, npaths, _ = span_of_differences(H, x, y2)
            if r is None: continue
            stats['mu=%d  rank_of_differences=%d  spans=%s' % (mu, r, r == mu)] += 1
            stats['SPANS' if r == mu else 'does NOT span'] += 1
            stats['d=%d' % d] += 1
for k in sorted(stats): print('  %-46s %s' % (k, stats[k]))
