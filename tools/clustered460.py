"""clustered460.py -- THE CLUSTERED CUBIC GRAPH, ON ITS OWN TERMS.

THE UNIQUE PROPERTY.  In a CUBIC graph every even subgraph is a
disjoint union of cycles: max degree 3 leaves no room for a degree-4
vertex.  So for two girth cycles C, C' sharing at least one edge, the
symmetric difference C + C' (the F2 sum in the cycle code) has
2g - 2c edges, each component is a cycle of length >= g, and
2g - 2c < 2g forces EXACTLY ONE component.

  GIRTH-PAIR THEOREM.  In a cubic graph two girth cycles that share
  an edge intersect in exactly one path, of length c in [1, g/2], and
  C + C' is a SINGLE cycle of length 2(g - c).

  (c <= g/2 because the difference cycle has length 2g - 2c >= g;
   t = 1 because two shared arcs would leave two components.)

THE OPENING.  Let 2^j be the largest power of two <= g - 1.  Then
2^{j+1} > g - 1, so 2^j >= g/2, so

    c* = g - 2^j   lies in [1, g/2]   -- ALWAYS a legal overlap.

  c* CRITERION.  If some two girth cycles meet in exactly c* edges
  then G has a cycle of length 2(g - c*) = 2^{j+1}, and
  Erdos-Gyarfas holds for G.

This is elementary, local, and checkable.  It also explains the
truncation counterexample of 452.6 exactly: in a truncation the
triangles are PAIRWISE VERTEX-DISJOINT, so no overlap of any size is
realised and the criterion has nothing to work with.  Disjoint girth
cycles are the real content of "cycle-clustered".

This file: verifies the theorem, evaluates the criterion across the
cages, truncations, random cubic graphs and random graphs of girth
>= 5, and asks whether girth-cycle disjointness is exactly a
3-edge-cut phenomenon - which 448.5 already forbids in a minimal
counterexample.
"""
import sys, random, itertools
from collections import Counter

sys.setrecursionlimit(100000)


def lcf(n, seq):
    adj = {i: set() for i in range(n)}
    for i in range(n):
        adj[i].add((i + 1) % n); adj[(i + 1) % n].add(i)
    for i in range(n):
        j = (i + seq[i % len(seq)]) % n
        adj[i].add(j); adj[j].add(i)
    return [sorted(adj[i]) for i in range(n)]


def gp(n, k):
    adj = {i: set() for i in range(2 * n)}
    def add(a, b):
        adj[a].add(b); adj[b].add(a)
    for i in range(n):
        add(i, (i + 1) % n); add(i, n + i)
        add(n + i, n + (i + k) % n)
    return [sorted(adj[i]) for i in range(2 * n)]


def truncate(adj):
    n = len(adj); idx = {}; m = 0
    for v in range(n):
        for w in adj[v]:
            idx[(v, w)] = m; m += 1
    out = [set() for _ in range(m)]
    def add(a, b):
        out[a].add(b); out[b].add(a)
    for v in range(n):
        ports = [idx[(v, w)] for w in adj[v]]
        for a in range(len(ports)):
            for b in range(a + 1, len(ports)):
                add(ports[a], ports[b])
    for v in range(n):
        for w in adj[v]:
            if v < w:
                add(idx[(v, w)], idx[(w, v)])
    return [sorted(x) for x in out]


def girth(adj):
    from collections import deque
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


def cycles_of_length(adj, L):
    """all cycles of length exactly L, as frozensets of edges."""
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
                    dfs(w, ln + 1)
                    path.pop(); seen.discard(w)
        dfs(s, 1)
    return sorted(out, key=lambda f: sorted(map(sorted, f)))


def is_single_cycle(edges):
    """is this edge set one cycle?  (cubic => even subgraph is a
    disjoint union of cycles; check connected + all degrees 2)"""
    if not edges: return False
    deg = Counter()
    for e in edges:
        for v in e: deg[v] += 1
    if any(d != 2 for d in deg.values()): return False
    adj = {}
    for e in edges:
        a, b = tuple(e)
        adj.setdefault(a, []).append(b); adj.setdefault(b, []).append(a)
    start = next(iter(adj)); seen = {start}; st = [start]
    while st:
        v = st.pop()
        for w in adj[v]:
            if w not in seen: seen.add(w); st.append(w)
    return len(seen) == len(deg)


def cstar(g):
    p = 1
    while p * 2 <= g - 1: p *= 2
    return g - p, 2 * p


def analyse(name, adj, verbose=True):
    n = len(adj); g = girth(adj)
    G = cycles_of_length(adj, g)
    cs, target = cstar(g)
    overlaps = Counter()
    thm_ok = thm_tot = 0
    for A, B in itertools.combinations(G, 2):
        I = A & B
        if not I: continue
        c = len(I)
        overlaps[c] += 1
        D = A ^ B
        thm_tot += 1
        if is_single_cycle(D) and len(D) == 2 * (g - c) and c <= g / 2:
            thm_ok += 1
    realised = sorted(overlaps)
    hit = cs in overlaps
    if verbose:
        print('  %-26s n=%4d g=%2d  |C_g|=%5d  c*=%2d -> %3d   '
              'overlaps=%-18s  %s'
              % (name, n, g, len(G), cs, target,
                 str(realised)[:18],
                 'HIT' if hit else 'MISS'), flush=True)
    return g, len(G), thm_ok, thm_tot, hit, len(overlaps) == 0


ZOO = [('Petersen', gp(5, 2)), ('Heawood', lcf(14, [5, -5] * 7)),
       ('Moebius-Kantor', lcf(16, [5, -5] * 8)),
       ('Pappus', lcf(18, [5, 7, -7, 7, -7, -5] * 3)),
       ('Desargues', lcf(20, [5, -5, 9, -9] * 5)),
       ('dodecahedron', gp(10, 2)),
       ('Nauru', gp(12, 5)),
       ('McGee', lcf(24, [12, 7, -7] * 8)),
       ('Tutte-Coxeter', lcf(30, [-13, -9, 7, -7, 9, 13] * 5)),
       ('Foster', lcf(90, [17, -9, 37, -37, 9, -17] * 15))]

print('=== (1) THE GIRTH-PAIR THEOREM, VERIFIED ===', flush=True)
print('  every intersecting pair of girth cycles must give a SINGLE')
print('  cycle of length exactly 2(g-c), with c <= g/2.')
print()
tot_ok = tot_all = 0
for name, adj in ZOO:
    g, ng, ok, tt, hit, disj = analyse(name, adj, verbose=False)
    tot_ok += ok; tot_all += tt
    print('  %-26s pairs tested %6d   conform %6d   %s'
          % (name, tt, ok, 'OK' if ok == tt else '*** FAIL ***'),
          flush=True)
print('  TOTAL %d/%d' % (tot_ok, tot_all), flush=True)

print()
print('=== (2) THE c* CRITERION ON THE CAGES ===', flush=True)
for name, adj in ZOO:
    analyse(name, adj)

print()
print('=== (3) TRUNCATIONS: GIRTH CYCLES PAIRWISE DISJOINT ===',
      flush=True)
for base, nm in ((gp(5, 2), 'truncated Petersen'),
                 (lcf(14, [5, -5] * 7), 'truncated Heawood'),
                 (gp(4, 1), 'truncated cube'),
                 (lcf(24, [12, 7, -7] * 8), 'truncated McGee')):
    g, ng, ok, tt, hit, disj = analyse(nm, truncate(base))
    print('       -> intersecting pairs: %d   %s'
          % (tt, 'ALL GIRTH CYCLES PAIRWISE DISJOINT' if disj else ''),
          flush=True)


def random_cubic(n, rng, mingirth=3, tries=8000):
    for _ in range(tries):
        stubs = [v for v in range(n) for _ in range(3)]
        rng.shuffle(stubs)
        adj = [set() for _ in range(n)]
        ok = True
        for i in range(0, len(stubs), 2):
            a, b = stubs[i], stubs[i + 1]
            if a == b or b in adj[a]:
                ok = False; break
            adj[a].add(b); adj[b].add(a)
        if not ok: continue
        seen = {0}; st = [0]
        while st:
            v = st.pop()
            for w in adj[v]:
                if w not in seen: seen.add(w); st.append(w)
        if len(seen) != n: continue
        A = [sorted(a) for a in adj]
        if girth(A) < mingirth: continue
        return A
    return None


print()
print('=== (4) RANDOM CUBIC GRAPHS ===', flush=True)
rng = random.Random(90210)
for mg, label in ((3, 'girth >= 3'), (5, 'girth >= 5')):
    hits = miss = disjs = 0
    ok = tt = 0
    for _ in range(40):
        adj = random_cubic(rng=rng, n=rng.choice([20, 24, 28, 30]),
                           mingirth=mg)
        if adj is None: continue
        g, ng, o, t, hit, disj = analyse('', adj, verbose=False)
        ok += o; tt += t
        if hit: hits += 1
        else: miss += 1
        if disj: disjs += 1
    print('  %-12s  samples=%3d   c* HIT=%3d   MISS=%3d   '
          'girth cycles pairwise disjoint=%3d   theorem %d/%d'
          % (label, hits + miss, hits, miss, disjs, ok, tt), flush=True)

print()
print('=== (5) WHEN c* IS MISSED, IS THE GRAPH DISJOINT-GIRTHED? ==='
      , flush=True)
print('  graph                       g   |C_g|  intersecting pairs'
      '   c* hit?')
rng2 = random.Random(13579)
shown = 0
for _ in range(200):
    if shown >= 8: break
    adj = random_cubic(rng=rng2, n=rng2.choice([20, 24, 28]),
                       mingirth=3)
    if adj is None: continue
    g, ng, o, t, hit, disj = analyse('', adj, verbose=False)
    if not hit:
        print('  %-26s %3d %6d %18d   %s'
              % ('random cubic (miss)', g, ng, t,
                 'no' if not hit else 'yes'), flush=True)
        shown += 1
