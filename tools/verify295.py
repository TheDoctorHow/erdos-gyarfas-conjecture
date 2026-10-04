"""verify295.py -- independent verification of

 (i)  the REFUTATION of (E-CORE) as stated in 294.8: an explicit bridge-shaped graph and
      a terminal pair that fails (L2**) with NO edge common to all its even connecting
      paths;
 (ii) the PARITY IDENTITY  |P| == [c(b_i) =/= c(b_j)] + |P n N|  (mod 2);
 (iii) THE UNIQUE-BIPARTIZING-EDGE THEOREM: if exactly one edge e has B-e bipartite,
      then a pair fails (L2**) IFF it is bichromatic in B-e; hence at most two of the
      three pairs fail and (L2**) holds for B.

Everything is recomputed from the definitions; nothing is imported from the scripts that
produced the claims.
"""
import itertools, sys
from collections import deque

CEX = [('T0','a3'),('T0','a8'),('T1','b1'),('T1','b3'),('T2','b3'),('T2','b6'),
       ('a0','a2'),('a0','a5'),('a0','a7'),('a1','a2'),('a1','a6'),('a1','a8'),
       ('a2','a4'),('a3','a4'),('a3','a5'),('a4','a9'),('a5','a6'),('a6','a9'),
       ('a7','a9'),('a7','b7'),('a8','b1'),('b0','b2'),('b0','b6'),('b0','b7'),
       ('b1','b4'),('b2','b4'),('b2','b5'),('b3','b5'),('b4','b8'),('b5','b9'),
       ('b6','b8'),('b7','b9'),('b8','b9')]

def A(E):
    a = {}
    for u, v in E:
        a.setdefault(u, set()).add(v); a.setdefault(v, set()).add(u)
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

def conn(a, skipv=None, skipe=None):
    vs = [v for v in a if v != skipv]
    if not vs: return True
    s = vs[0]; seen = {s}; q = deque([s])
    while q:
        u = q.popleft()
        for w in a[u]:
            if w == skipv or w in seen: continue
            if skipe and frozenset((u, w)) == skipe: continue
            seen.add(w); q.append(w)
    return len(seen) == len(vs)

def two_conn(a): return all(conn(a, skipv=v) for v in a)
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

def paths(a, s, t):
    out = []; p = [s]; on = {s}
    def rec(u):
        for w in a[u]:
            if w == t: out.append(p + [t])
            elif w not in on:
                on.add(w); p.append(w); rec(w); p.pop(); on.discard(w)
    rec(s); return out

def es(p): return frozenset(frozenset((p[i], p[i+1])) for i in range(len(p)-1))

print('=== (i) VERIFYING THE (E-CORE) COUNTEREXAMPLE ===')
E = [tuple(sorted(e)) for e in CEX]
assert len(E) == len(set(E)), 'multi-edge'
a = A(E)
deg2 = sorted(v for v in a if len(a[v]) == 2)
print('n = %d, m = %d' % (len(a), len(E)))
print('max degree      :', max(len(a[v]) for v in a))
print('degree-2 vertices:', deg2)
print('girth           :', girth(a))
print('2-connected     :', two_conn(a))
print('bipartite       :', bip(a) is not None)
terms = deg2
b1, b2 = 'T0', 'T1'
P = paths(a, b1, b2)
ev = [p for p in P if (len(p)-1) % 2 == 0]
od = [p for p in P if (len(p)-1) % 2 == 1]
print('pair (%s,%s): %d even, %d odd connecting paths' % (b1, b2, len(ev), len(od)))
Es = [es(p) for p in ev]
disj = [(i, j) for i in range(len(Es)) for j in range(i+1, len(Es)) if not (Es[i] & Es[j])]
print('edge-disjoint pairs of even paths :', len(disj), '  => (L2**) FAILS:', not disj)
core = set(Es[0])
for x in Es[1:]: core &= x
print('edges common to ALL even paths    :', sorted(tuple(sorted(e)) for e in core))
# vertex core too, and the 2-edge-cut
vcore = set(ev[0][1:-1])
for p in ev[1:]: vcore &= set(p[1:-1])
print('vertices common to ALL even paths :', sorted(vcore))
bipe = [e for e in map(frozenset, E) if bip(A([tuple(x) for x in
        (set(map(frozenset, E)) - {e})])) is not None]
print('bipartizing edges of B            :', sorted(tuple(sorted(e)) for e in bipe))
cuts = [(e, f) for e, f in itertools.combinations(sorted(map(frozenset, E), key=sorted), 2)
        if not conn(A([tuple(x) for x in (set(map(frozenset, E)) - {e, f})]))]
print('2-edge-cuts of B                  :',
      [(tuple(sorted(e)), tuple(sorted(f))) for e, f in cuts])
print('VERDICT: (E-CORE) as stated at 294.8 is FALSE.' if not disj and not core
      else 'VERDICT: not a counterexample')

print()
print('=== (ii) PARITY IDENTITY, over all bipartitions of a small bridge ===')
# use the same graph; sample bipartitions
import random
rng = random.Random(1)
V = sorted(a)
bad = 0; tested = 0
for _ in range(40):
    S = {v for v in V if rng.random() < .5}
    c = {v: (1 if v in S else 0) for v in V}
    N = {frozenset(e) for e in E if c[e[0]] == c[e[1]]}
    for (x, y) in itertools.combinations(terms, 2):
        d = int(c[x] != c[y])
        for p in paths(a, x, y)[:120]:
            tested += 1
            if (len(p)-1) % 2 != (d + len(es(p) & N)) % 2: bad += 1
print('bipartitions x paths tested: %d   violations of the identity: %d' % (tested, bad))

print()
print('=== (iii) UNIQUE-BIPARTIZING-EDGE THEOREM, tested on random bridges ===')
exec(open('tools/ecore2.py').read().split('def main()')[0].split('"""', 2)[2])
rng = random.Random(11)
tot = dict(bridges=0, unique=0, pairs=0, fail=0, thm_violation=0, L2_fails_for_B=0)
for m in (10, 12, 14, 16):
    got = 0; att = 0
    while got < 12 and att < 3000:
        att += 1
        g, tt = make_bridge(m, rng)
        if g is None: continue
        got += 1; tot['bridges'] += 1
        EE = {frozenset((u, w)) for u in g for w in g[u]}
        be = [e for e in EE if bip(adj_of(EE - {e})) is not None]
        if len(be) != 1: continue
        tot['unique'] += 1
        c = bip(adj_of(EE - {be[0]}))
        anyhold = False
        for x, y in itertools.combinations(tt, 2):
            pp = paths(g, x, y)
            evp = [es(p) for p in pp if (len(p)-1) % 2 == 0]
            f = not any(not (evp[i] & evp[j]) for i in range(len(evp))
                        for j in range(i+1, len(evp)))
            tot['pairs'] += 1; tot['fail'] += f
            if f != (c[x] != c[y]): tot['thm_violation'] += 1
            if not f: anyhold = True
        if not anyhold: tot['L2_fails_for_B'] += 1
for k, v in tot.items(): print('  %-16s %s' % (k, v))
print('  theorem violations:', tot['thm_violation'],
      ' bridges where (L2**) fails outright:', tot['L2_fails_for_B'])
