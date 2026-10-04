"""refute464.py -- THE DIVISIBILITY CLOSURE IS FALSE.

CLAIM TO REFUTE.  If delta* > 0 and some girth-cycle overlap is
>= c*, then some overlap lies in C* = {(g-d)/2 : d | delta*}.

THE ESCAPE, located by arithmetic.  delta* = 2 (g = 2^{j+1} - 2, so
g = 6, 14, 30, ...) gives C* = {g/2 - 1}, a SINGLE value, while the
overlaps at or above c* are {g/2 - 1, g/2}.  An adversary whose only
overlap is g/2 satisfies the hypothesis and misses C*.  Overlap g/2
has delta = g - 2c = 0: the chain has zero step and produces nothing.

THE CONSTRUCTION.  Overlap g/2 is realised by Theta(g/2, g/2, g/2):
two branch vertices joined by three paths of length g/2.  All three
of its cycles are girth cycles and they pairwise meet in exactly g/2.
For g = 6 that is Theta(3,3,3), on 8 vertices, whose six interior
vertices each carry one free half-edge (no two interior vertices may
be joined: it would make a triangle or a 4-cycle).

IMPLANT it into a cubic graph H of girth >= 8 - which therefore has
NO hexagons of its own - by deleting three pairwise distant
independent edges of H and matching the six resulting half-edges to
the theta's six.  The result is cubic.  Its only hexagons are the
theta's three, so

    O = {3},  g = 6,  c* = 2,  C* = {2},  max O = 3 >= c*,
    O n C* = EMPTY.

That refutes the closure.  This file builds it and checks every claim.
"""
import sys, itertools
from collections import Counter, deque

sys.setrecursionlimit(100000)


def lcf(n, seq):
    adj = {i: set() for i in range(n)}
    for i in range(n):
        adj[i].add((i + 1) % n); adj[(i + 1) % n].add(i)
    for i in range(n):
        j = (i + seq[i % len(seq)]) % n
        adj[i].add(j); adj[j].add(i)
    return [sorted(adj[i]) for i in range(n)]


def girth(adj):
    n = len(adj); best = 10 ** 9
    for s in range(n):
        dist = {s: 0}; par = {s: -1}; dq = deque([s])
        while dq:
            v = dq.popleft()
            for w in adj[v]:
                if w not in dist:
                    dist[w] = dist[v] + 1; par[w] = v; dq.append(w)
                elif w != par[v]:
                    best = min(best, dist[v] + dist[w] + 1)
    return best


def dist(adj, s):
    d = {s: 0}; dq = deque([s])
    while dq:
        v = dq.popleft()
        for w in adj[v]:
            if w not in d:
                d[w] = d[v] + 1; dq.append(w)
    return d


def cycles_of_length(adj, L):
    n = len(adj); out = set()
    for s in range(n):
        seen = {s}; path = [s]
        def dfs(v, ln):
            for w in adj[v]:
                if w == s and ln == L:
                    out.add(frozenset(
                        frozenset((path[i], path[(i + 1) % ln]))
                        for i in range(ln)))
                elif w > s and w not in seen and ln < L:
                    seen.add(w); path.append(w)
                    dfs(w, ln + 1); path.pop(); seen.discard(w)
        dfs(s, 1)
    return sorted(out, key=lambda f: sorted(map(sorted, f)))


# ---- H = Tutte-Coxeter, girth 8, 30 vertices, no hexagons at all
H = lcf(30, [-13, -9, 7, -7, 9, 13] * 5)
nH = len(H)
print('=== (0) THE HOST ===', flush=True)
print('  Tutte-Coxeter: n=%d, girth=%d, hexagons=%d'
      % (nH, girth(H), len(cycles_of_length(H, 6))), flush=True)

# pick 3 pairwise distant independent edges of H
edges = sorted({(min(u, v), max(u, v)) for u in range(nH)
                for v in H[u]})
D = {v: dist(H, v) for v in range(nH)}


def far_triple():
    """Three independent edges, pairwise endpoint distance >= 3.
    (Within one removed edge the endpoints end up at distance
    girth-1 = 7, which covers the theta's adjacent interior pairs.)"""
    for e1, e2, e3 in itertools.combinations(edges, 3):
        vs = set(e1) | set(e2) | set(e3)
        if len(vs) != 6: continue
        ok = True
        for ea, eb in itertools.combinations((e1, e2, e3), 2):
            for a in ea:
                for b in eb:
                    if D[a][b] < 3: ok = False
        if ok: return e1, e2, e3
    return None


trip = far_triple()
print('  three pairwise distant independent edges: %s' % (trip,),
      flush=True)

e1, e2, e3 = trip
adj = [set(x) for x in H]
for (a, b) in (e1, e2, e3):
    adj[a].discard(b); adj[b].discard(a)

# ---- Theta(3,3,3): branch v0,v3 ; paths v0-v1-v2-v3, v0-b-a-v3,
#                    v0-b'-a'-v3
base = nH
V0, V1, V2, V3, A, B, A2, B2 = (base + i for i in range(8))
for _ in range(8):
    adj.append(set())


def link(u, v):
    adj[u].add(v); adj[v].add(u)


for u, v in ((V0, V1), (V1, V2), (V2, V3),
             (V0, B), (B, A), (A, V3),
             (V0, B2), (B2, A2), (A2, V3)):
    link(u, v)

# ---- match the six theta half-edges to the six host half-edges.
# The theta's three ADJACENT interior pairs (V1,V2), (B,A), (B2,A2)
# each need host distance >= 4 to avoid a new short cycle; pairing
# each with the two endpoints of one REMOVED edge gives distance
# girth-1 = 7 for free.  Distinct removed edges are >= 3 apart, which
# covers the theta pairs at internal distance 2 and 3.
for (t1, t2), (h1, h2) in (((V1, V2), e1), ((B, A), e2),
                           ((B2, A2), e3)):
    link(t1, h1); link(t2, h2)

G = [sorted(x) for x in adj]
print()
print('=== (1) THE IMPLANTED GRAPH ===', flush=True)
print('  n = %d   all degrees 3: %s   connected: %s'
      % (len(G), all(len(x) == 3 for x in G),
         len(dist(G, 0)) == len(G)), flush=True)
g = girth(G)
print('  girth = %d' % g, flush=True)

hexes = cycles_of_length(G, g)
print('  girth cycles: %d' % len(hexes), flush=True)

O = Counter()
for X, Y in itertools.combinations(hexes, 2):
    if X & Y: O[len(X & Y)] += 1
print('  overlap multiset: %s' % dict(O), flush=True)


def cstar(gg):
    p = 1
    while p * 2 <= gg - 1: p *= 2
    return gg - p, 2 * p


cs, tgt = cstar(g)
ds = tgt - g
C = {(g - d) // 2 for d in range(1, ds + 1)
     if ds % d == 0 and (d - g) % 2 == 0 and 1 <= (g - d) // 2 <= g // 2}
print()
print('=== (2) THE VERDICT ===', flush=True)
print('  g=%d  target 2^{j+1}=%d  delta*=%d  c*=%d  C*=%s'
      % (g, tgt, ds, cs, sorted(C)))
Os = set(O)
print('  realised overlaps O = %s' % sorted(Os))
print('  hypothesis  max O >= c* : %s' % (max(Os) >= cs if Os else False))
print('  conclusion  O n C*      : %s' % (sorted(Os & C) or 'EMPTY'))
if Os and max(Os) >= cs and not (Os & C):
    print('  ==> DIVISIBILITY CLOSURE IS FALSE.', flush=True)
else:
    print('  ==> no refutation from this graph.', flush=True)

print()
print('=== (3) EDGE MULTIPLICITY, for the record ===', flush=True)
r = Counter()
for cyc in hexes:
    for e in cyc: r[e] += 1
print('  r_e over edges carrying a girth cycle: %s'
      % dict(Counter(r.values())))
print('  max r_e = %d   threshold 2^{c*-1} = %d   forcing fires: %s'
      % (max(r.values()) if r else 0, 2 ** (cs - 1),
         (max(r.values()) if r else 0) > 2 ** (cs - 1)), flush=True)
print('  (consistent with 459: the escape sits AT the threshold.)')

print()
print('=== (4) DOES THE GRAPH HAVE AN 8-CYCLE ANYWAY? ===', flush=True)
c8 = cycles_of_length(G, 8)
print('  8-cycles in the implanted graph: %d' % len(c8), flush=True)
print('  (the refutation is of the CLOSURE LEMMA, not of')
print('   Erdos-Gyarfas: the host supplies 8-cycles at another scale.)',
      flush=True)
