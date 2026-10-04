"""span301.py -- is (SPAN) true in general?

(SPAN)  H 2-connected, s =/= t.  The differences P xor P_0 over s-t paths span C(H).

Equivalently: the only edge sets F such that all s-t paths have the same |F n P| mod 2
are the CUTS.  Tested over ALL 2-connected graphs on <= 7 vertices and random ones above.
"""
import itertools, random, collections
from collections import deque

def adj(E):
    a = {}
    for e in E:
        u, v = tuple(e); a.setdefault(u, set()).add(v); a.setdefault(v, set()).add(u)
    return a

def conn(a, skip=None, n=None):
    vs = [v for v in a if v != skip]
    if not vs: return True
    s = vs[0]; seen = {s}; q = deque([s])
    while q:
        u = q.popleft()
        for w in a[u]:
            if w != skip and w not in seen: seen.add(w); q.append(w)
    return len(seen) == len(vs)

def two_conn(a, n):
    return len(a) == n and conn(a) and all(conn(a, v) for v in a)

def gf2_rank(rows):
    rank = 0
    for col in range(64):
        p = None
        for i in range(rank, len(rows)):
            if (rows[i] >> col) & 1: p = i; break
        if p is None: continue
        rows[rank], rows[p] = rows[p], rows[rank]
        for i in range(len(rows)):
            if i != rank and ((rows[i] >> col) & 1): rows[i] ^= rows[rank]
        rank += 1
    return rank

def paths(a, s, t):
    out = []; p = [s]; on = {s}
    def rec(u):
        for w in a[u]:
            if w == t: out.append(list(p)+[t])
            elif w not in on:
                on.add(w); p.append(w); rec(w); p.pop(); on.discard(w)
    rec(s); return out

def test(E, n, stats, label):
    a = adj(E)
    if not two_conn(a, n): return
    El = sorted(E, key=lambda f: sorted(f)); idx = {e: i for i, e in enumerate(El)}
    mu = len(El) - n + 1
    for s, t in itertools.combinations(sorted(a), 2):
        P = paths(a, s, t)
        if len(P) < 2: continue
        def mask(p):
            m = 0
            for i in range(len(p)-1): m |= 1 << idx[frozenset((p[i], p[i+1]))]
            return m
        b0 = mask(P[0])
        r = gf2_rank([mask(p) ^ b0 for p in P[1:]])
        stats['tested'] += 1
        if r == mu: stats['SPANS'] += 1
        else:
            stats['*** FAILS ***'] += 1
            if stats['*** FAILS ***'] <= 3:
                print('SPAN FAILS: %s n=%d mu=%d rank=%d s=%s t=%s  edges=%s'
                      % (label, n, mu, r, s, t, sorted(tuple(sorted(f)) for f in El)))

stats = collections.Counter()
# exhaustive: all graphs on n <= 7 vertices
for n in range(3, 7):
    allE = [frozenset(p) for p in itertools.combinations(range(n), 2)]
    for k in range(n, len(allE)+1):
        for S in itertools.combinations(allE, k):
            test(set(S), n, stats, 'exhaustive n=%d' % n)
    print('  ...done n=%d, cumulative: %s' % (n, dict(stats)))
# random larger, including subcubic
rng = random.Random(4)
for _ in range(600):
    n = rng.choice([8, 9, 10, 11, 12])
    allE = [frozenset(p) for p in itertools.combinations(range(n), 2)]
    k = rng.randint(n, min(len(allE), 2*n))
    test(set(rng.sample(allE, k)), n, stats, 'random n=%d' % n)
for k in sorted(stats): print('  %-16s %s' % (k, stats[k]))
