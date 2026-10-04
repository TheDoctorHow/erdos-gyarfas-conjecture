"""cactus3.py -- the predicted ALL-THREE-PAIRS-FAIL bridge, built exactly.

From the parity identity 295.1: take three BIPARTITE cubic parts D_1,D_2,D_3, one
terminal in each, joined in a CYCLE by three single edges e_1,e_2,e_3.  Every b_i-b_j
path uses one cut edge or the other two, so for each pair the even connecting paths all
share an edge -- and all three pairs should fail (L2**).

Parts: the Heawood graph (cubic, bipartite, girth 6) minus one edge.
"""
import itertools, sys
from collections import deque

def heawood():
    E = set()
    for i in range(14): E.add(frozenset((i, (i + 1) % 14)))
    for i in range(0, 14, 2): E.add(frozenset((i, (i + 5) % 14)))
    return E

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
    s = vs[0]; seen = {s}; q = deque([s])
    while q:
        u = q.popleft()
        for w in a[u]:
            if w != skip and w not in seen: seen.add(w); q.append(w)
    return len(seen) == len(vs)

def two_conn(a): return all(conn(a, v) for v in a)

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
            if w == t: out.append(p+[t])
            elif w not in on:
                on.add(w); p.append(w); rec(w); p.pop(); on.discard(w)
    rec(s); return out

def es(p): return frozenset(frozenset((p[i], p[i+1])) for i in range(len(p)-1))

def build(cuts, subs):
    """cuts[k] = edge removed from copy k (its ends become the join points)
       subs[k] = edge of copy k to subdivide into terminal T_k"""
    E = set(); ends = []
    for k in range(3):
        H = {frozenset(('%d_%d' % (k, u), '%d_%d' % (k, v))) for u, v in map(tuple, heawood())}
        cu, cv = cuts[k]
        ce = frozenset(('%d_%d' % (k, cu), '%d_%d' % (k, cv)))
        H.discard(ce)
        su, sv = subs[k]
        se = frozenset(('%d_%d' % (k, su), '%d_%d' % (k, sv)))
        assert se in H, 'subdivided edge must survive'
        H.discard(se)
        H.add(frozenset(('%d_%d' % (k, su), 'T%d' % k)))
        H.add(frozenset(('%d_%d' % (k, sv), 'T%d' % k)))
        E |= H
        ends.append(('%d_%d' % (k, cu), '%d_%d' % (k, cv)))
    for k in range(3):
        E.add(frozenset((ends[k][1], ends[(k+1) % 3][0])))
    return E

def analyse(E, label):
    a = adj(E)
    terms = ['T0', 'T1', 'T2']
    print('--- %s ---' % label)
    print('n=%d m=%d maxdeg=%d girth=%d 2-conn=%s bipartite=%s'
          % (len(a), len(E), max(len(a[v]) for v in a), girth(a), two_conn(a),
             bip(a) is not None))
    print('degree-2 vertices:', sorted(v for v in a if len(a[v]) == 2))
    d = {(x, y): dist(a, x, y) for x, y in itertools.combinations(terms, 2)}
    print('terminal distances:', d, ' min =', min(d.values()))
    if bip(a) is not None:
        print('  BIPARTITE - not a valid bridge, skipping'); return None
    bipe = [e for e in E if bip(adj(E - {e})) is not None]
    print('bipartizing edges:', sorted(tuple(sorted(map(str, e))) for e in bipe))
    nfail = 0
    for x, y in itertools.combinations(terms, 2):
        ev = [p for p in paths(a, x, y) if (len(p)-1) % 2 == 0]
        if len(ev) < 2:
            print('  pair %s,%s : fewer than two even connecting paths (%d)' % (x, y, len(ev)))
            nfail += 1; continue
        Es = [es(p) for p in ev]
        core = set(Es[0])
        for z in Es[1:]: core &= z
        dis = any(not (Es[i] & Es[j]) for i in range(len(Es)) for j in range(i+1, len(Es)))
        nfail += (not dis)
        print('  pair %s,%s : %4d even paths, edge-disjoint pair exists = %-5s, common edges %s'
              % (x, y, len(ev), dis, sorted(tuple(sorted(map(str, e))) for e in core)))
    print('  FAILING PAIRS: %d of 3   =>  (L2**) %s'
          % (nfail, 'FAILS FOR THIS BRIDGE' if nfail == 3 else 'holds'))
    return nfail

if __name__ == '__main__':
    # cut edges and subdivided edges chosen far apart inside each Heawood copy
    for shift in range(3):
        cuts = [(0, 1)] * 3
        subs = [(7 + shift, 8 + shift)] * 3
        try:
            E = build(cuts, subs)
        except AssertionError as ex:
            print('skip', ex); continue
        analyse(E, 'three Heawood parts, cut (0,1), subdivide (%d,%d)'
                % (7 + shift, 8 + shift))
