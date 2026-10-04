"""n3search.py -- exhaustive construction of the (N3) configuration.

299.1: B has exactly two bipartizing edges, both at a terminal b_1, so H = B - b_1 is
BIPARTITE with exactly four degree-2 vertices.  Counting degrees across the bipartition,
3(|X|-|Y|) = a-b with a+b = 4 and 3 | (a-b), so a = b = 2: EXACTLY TWO DEGREE-2 VERTICES
IN EACH CLASS.  With x,y the neighbours of b_1 (which must lie in different classes, else
B is bipartite) and b_2,b_3 the other terminals (which must lie in different classes, else
295.5 closes the case), the four are forced to be  X = {x, b_3},  Y = {y, b_2}.

So: H = a bipartite cubic graph minus two independent edges, and B = H + b_1 joined to one
X-vertex and one Y-vertex.  This script builds every such B from the Heawood graph and
from random bipartite cubic graphs, checks all the bridge hypotheses, and tests (SC).
"""
import itertools, sys, random, collections
from collections import deque

def adj(E):
    a = {}
    for e in E:
        u, v = tuple(e); a.setdefault(u, set()).add(v); a.setdefault(v, set()).add(u)
    return a

def girth(a):
    best = 10**9
    for s in a:
        d = {s: 0}; par = {s: None}; q = deque([s])
        while q:
            u = q.popleft()
            for w in a[u]:
                if w not in d: d[w] = d[u]+1; par[w] = u; q.append(w)
                elif w != par[u]: best = min(best, d[u]+d[w]+1)
    return best

def conn(a, skip=None):
    vs = [v for v in a if v != skip]
    if not vs: return True
    s = vs[0]; seen = {s}; q = deque([s])
    while q:
        u = q.popleft()
        for w in a[u]:
            if w != skip and w not in seen: seen.add(w); q.append(w)
    return len(seen) == len(vs)

def two_conn(a): return len(a) >= 3 and conn(a) and all(conn(a, v) for v in a)

def bip(a):
    col = {}
    for s in a:
        if s in col: continue
        col[s] = 0; q = deque([s])
        while q:
            u = q.popleft()
            for w in a[u]:
                if w not in col: col[w] = 1-col[u]; q.append(w)
                elif col[w] == col[u]: return None
    return col

def dist(a, s, t):
    d = {s: 0}; q = deque([s])
    while q:
        u = q.popleft()
        if u == t: return d[u]
        for w in a[u]:
            if w not in d: d[w] = d[u]+1; q.append(w)
    return 10**9

def paths(a, s, t):
    out = []; p = [s]; on = {s}
    def rec(u):
        for w in a[u]:
            if w == t: out.append(list(p)+[t])
            elif w not in on:
                on.add(w); p.append(w); rec(w); p.pop(); on.discard(w)
    rec(s); return out

def heawood():
    E = set()
    for i in range(14): E.add(frozenset((i, (i+1) % 14)))
    for i in range(0, 14, 2): E.add(frozenset((i, (i+5) % 14)))
    return E

def rand_bip_cubic(m, rng, tries=3000):
    for _ in range(tries):
        L = list(range(m)); R = list(range(m, 2*m))
        st = [v for v in L for _ in range(3)]; tg = [v for v in R for _ in range(3)]
        rng.shuffle(tg)
        E = set(); ok = True
        for u, v in zip(st, tg):
            if frozenset((u, v)) in E: ok = False; break
            E.add(frozenset((u, v)))
        if not ok or len(E) != 3*m: continue
        a = adj(E)
        if len(a) != 2*m or girth(a) < 6 or not two_conn(a): continue
        return E
    return None

def try_base(E0, label, stats, mindist):
    a0 = adj(E0); c0 = bip(a0)
    if c0 is None: return
    El = sorted(E0, key=lambda s: sorted(map(str, s)))
    for e1, e2 in itertools.combinations(El, 2):
        if e1 & e2: continue                      # must be independent
        H = set(E0) - {e1, e2}
        aH = adj(H)
        if not two_conn(aH) or girth(aH) < 6: continue
        cH = bip(aH)
        d2 = [v for v in aH if len(aH[v]) == 2]
        if len(d2) != 4: continue
        X = [v for v in d2 if cH[v] == 0]; Y = [v for v in d2 if cH[v] == 1]
        if len(X) != 2 or len(Y) != 2: continue
        for x in X:
            b3 = [v for v in X if v != x][0]
            for y in Y:
                b2 = [v for v in Y if v != y][0]
                E = set(H) | {frozenset(('B1', x)), frozenset(('B1', y))}
                a = adj(E)
                if girth(a) < 5 or not two_conn(a) or bip(a) is not None: continue
                T = ['B1', b2, b3]
                if min(dist(a, u, v) for u, v in itertools.combinations(T, 2)) < mindist:
                    stats['rejected_distance'] += 1; continue
                bipe = [f for f in E if bip(adj(E-{f})) is not None]
                if set(bipe) != {frozenset(('B1', x)), frozenset(('B1', y))}:
                    stats['rejected_bipedges_%d' % len(bipe)] += 1; continue
                stats['CONFIGURATIONS'] += 1
                # test (SC): all even connecting path lengths agree mod 4, each pair
                bad = []
                for u, v in itertools.combinations(T, 2):
                    L = [len(p)-1 for p in paths(a, u, v) if (len(p)-1) % 2 == 0]
                    cl = {l % 4 for l in L}
                    if len(cl) > 1: bad.append((u, v, sorted(cl), len(L)))
                if not bad:
                    stats['*** (SC) SATISFIED ***'] += 1
                    print('(SC) HOLDS!  n=%d  %s  terminals %s' % (len(a), label, T))
                    print('   edges:', sorted(tuple(sorted(map(str, f))) for f in E))
                else:
                    stats['SC_fails'] += 1
                    stats['first_failing_pair_%s' % ('with_b1' if 'B1' in
                          (bad[0][0], bad[0][1]) else 'b2b3')] += 1
                    stats['num_failing_pairs_%d' % len(bad)] += 1

stats = collections.Counter()
mind = int(sys.argv[1]) if len(sys.argv) > 1 else 5
try_base(heawood(), 'Heawood', stats, mind)
rng = random.Random(5)
for m in (7, 8, 9):
    for k in range(4):
        E0 = rand_bip_cubic(m, rng)
        if E0: try_base(E0, 'bip-cubic-%d' % (2*m), stats, mind)
for k in sorted(stats): print('  %-34s %s' % (k, stats[k]))
