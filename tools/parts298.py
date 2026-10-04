"""parts298.py -- can a PART of the |Bip|=3 cactus satisfy what (SC) forces on it?

296.2 says each part D is bipartite, subcubic, 2-connected, girth >= 5 (so >= 6), with
EXACTLY THREE vertices of degree 2 -- its terminal and its two cut-endpoints -- all in
ONE class.  The independence of the three segments in A_1, A_2, A_3 then forces (SC) to
require, for EVERY pair of those three special vertices, that ALL paths between them are
congruent modulo 4.  All such paths are EVEN, since the specials share a class.

This script exhaustively tests that condition on every such D it can generate.
"""
import sys, itertools, random, collections
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

def allpaths(a, s, t):
    out = []; p = [s]; on = {s}
    def rec(u):
        for w in a[u]:
            if w == t: out.append(len(p))          # length in edges
            elif w not in on:
                on.add(w); p.append(w); rec(w); p.pop(); on.discard(w)
    rec(s); return out

def mod4_rigid(a, specials):
    """all paths between each pair of specials congruent mod 4?"""
    for x, y in itertools.combinations(specials, 2):
        L = allpaths(a, x, y)
        if len({l % 4 for l in L}) > 1:
            return False, (x, y, sorted({l % 4 for l in L}), len(L))
    return True, None

# ---- candidate parts: cubic bipartite girth-6 graph minus a vertex ---------------
def heawood():
    E = set()
    for i in range(14): E.add(frozenset((i, (i+1) % 14)))
    for i in range(0, 14, 2): E.add(frozenset((i, (i+5) % 14)))
    return E

def rand_bip_cubic(m, rng, tries=4000):
    """random bipartite cubic graph on 2m vertices, girth >= 6"""
    for _ in range(tries):
        L = list(range(m)); R = list(range(m, 2*m))
        stubs = [v for v in L for _ in range(3)]
        tgt = [v for v in R for _ in range(3)]
        rng.shuffle(tgt)
        E = set(); ok = True
        for u, v in zip(stubs, tgt):
            if frozenset((u, v)) in E: ok = False; break
            E.add(frozenset((u, v)))
        if not ok or len(E) != 3*m: continue
        a = adj(E)
        if len(a) != 2*m or girth(a) < 6 or not two_conn(a): continue
        return E
    return None

def test_part(E, label, rng=None):
    a = adj(E)
    d2 = sorted(v for v in a if len(a[v]) == 2)
    c = bip(a)
    if c is None: return None
    if len(d2) != 3 or len({c[v] for v in d2}) != 1: return None
    ok, why = mod4_rigid(a, d2)
    print('  %-42s n=%2d specials=%s  mod-4 rigid: %s%s'
          % (label, len(a), d2, ok, '' if ok else '   (%s..%s hits classes %s over %d paths)' % why))
    return ok

print('PARTS THAT (SC) WOULD REQUIRE (all three special pairs rigid mod 4):')
res = []
H = heawood()
for v in range(14):
    E = {e for e in H if v not in e}
    r = test_part(E, 'Heawood - v%d' % v)
    if r is not None: res.append(r)

rng = random.Random(9)
for m in (7, 8, 9, 10, 11):
    for k in range(6):
        E = rand_bip_cubic(m, rng)
        if E is None: continue
        a = adj(E)
        v = rng.choice(sorted(a))
        E2 = {e for e in E if v not in e}
        if not two_conn(adj(E2)): continue
        r = test_part(E2, 'random bipartite cubic K_{%d,%d} - v' % (m, m))
        if r is not None: res.append(r)

print()
print('parts tested: %d   satisfying the mod-4 condition: %d' % (len(res), sum(res)))
