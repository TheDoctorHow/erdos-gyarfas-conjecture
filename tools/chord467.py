"""chord467.py -- THE a = 1 LAYER, AND WHERE THE ENERGY MUST LIVE.

466 proved the Attachment Identity T_m = SUM_a R_a(m-2a) and showed
the exact weight decomposition

    S := 2m T_m / 2^m = SUM_a chi_a gamma_{m-2a} 4^{-a} m/(m-2a),

so the energy inequality T_m < 2^m/2m is dominated by SMALL a: the
weight decays like 4^{-a}.  Two consequences are pursued here.

(1) THE ESCAPE OF 460 IS EXPONENTIALLY CHEAP.  The delta = 0 fixed
    point is Theta(g/2, g/2, g/2): two girth cycles meeting in a path
    of length g/2.  In the attachment decomposition that is a = g/2,
    weight 4^{-g/2}.  So the configuration that defeats the girth-
    scale closure argument sits in the most heavily suppressed layer
    of the energy and provably cannot supply horn-two energy, which
    needs a = O(1).  This is the quantitative form of 461.7(B2).

(2) THE a = 1 LAYER IS THE CLOSURE FRAME.  An attachment of length 1
    is a CHORD, and a (cycle, chord) pair is exactly a theta with a
    unit branch - two cycles meeting in exactly one edge.  That is
    the object of 457's Girth-Pair Theorem and 459's edge
    multiplicity r_e, but at scale m instead of at the girth.  So
    452-456 and 457-460 are not two frames: the closure frame is the
    a = 1 layer of the counting frame.

    A chord of a z-cycle at arc-distance d yields cycles d+1 and
    z-d+1.  So the a = 1 layer WINS whenever some chord sits at
    arc-distance d = 2^j - 1.  With horn-two energy there are
    ~2^{z+1}/z chord incidences and only ~2 log2 z winning distances,
    so the adversary must hold an enormous family of arc-distances
    off a logarithmically thin set.  This file measures whether real
    graphs manage that.
"""
import sys, random, itertools
from collections import Counter, deque

sys.setrecursionlimit(100000)

exec(open('tools/energy466.py').read().split(
    "ZOO = [")[0].split('"""', 2)[2])


def cyc_of_len(adj, L):
    """simple cycles of length exactly L, as vertex tuples in order."""
    n = len(adj); out = []
    for s in range(n):
        path = [s]; seen = {s}
        def dfs(v):
            if len(path) == L:
                if s in adj[v] and path[1] < path[-1]:
                    out.append(tuple(path))
                return
            for w in adj[v]:
                if w > s and w not in seen:
                    seen.add(w); path.append(w); dfs(w)
                    path.pop(); seen.discard(w)
        dfs(s)
    return out


def chords(adj, cyc):
    """(arc-distance, both) for each chord of the cycle."""
    z = len(cyc); pos = {v: i for i, v in enumerate(cyc)}
    ring = set()
    for i in range(z):
        ring.add(frozenset((cyc[i], cyc[(i + 1) % z])))
    out = []
    for i, u in enumerate(cyc):
        for w in adj[u]:
            if w in pos and frozenset((u, w)) not in ring:
                j = pos[w]
                d = (j - i) % z
                if d > z - d: d = z - d
                if i < j: out.append(d)
    return out


def pow2_dists(z):
    """arc-distances d whose chord yields a power-of-two cycle:
    d + 1 = 2^j  or  z - d + 1 = 2^j."""
    S = set(); j = 1
    while 2 ** j <= z + 1:
        for d in (2 ** j - 1, z - 2 ** j + 1):
            if 1 <= d <= z - 1: S.add(min(d, z - d))
        j += 1
    return S


def girth_of(adj):
    return girth(adj)


print('=' * 70)
print('(A) THE 460 ESCAPE, WEIGHED')
print('=' * 70)
print('  Theta(g/2,g/2,g/2) is an attachment of length a = g/2 to a')
print('  cycle of length g.  Its weight in S is 4^{-a}.')
print()
print('  %6s %8s %14s %20s' % ('g', 'a = g/2', '4^{-a}', 'a=1 weight ratio'))
for g in (6, 10, 14, 18, 22, 30, 62):
    a = g // 2
    print('  %6d %8d %14.3e %20s'
          % (g, a, 4.0 ** (-a), '1 : %.3e' % (4.0 ** (-a) / 0.25)))
print()
print('  The configuration that defeats 457-459 at the girth scale is')
print('  suppressed by 4^{-g/2} in the energy.  It cannot be the')
print('  source of horn-two energy, which needs a = O(1).  461.7(B2)')
print('  conjectured this; the attachment weight proves it.')

print()
print('=' * 70)
print('(B) THE 38-VERTEX WITNESS OF 460, DECOMPOSED')
print('=' * 70)
H = lcf(30, [-13, -9, 7, -7, 9, 13] * 5)
D = {}
for v in range(30):
    d = {v: 0}; q = deque([v])
    while q:
        x = q.popleft()
        for y in H[x]:
            if y not in d: d[y] = d[x] + 1; q.append(y)
    D[v] = d
edges = sorted({(min(u, v), max(u, v)) for u in range(30) for v in H[u]})
trip = None
for e1, e2, e3 in itertools.combinations(edges, 3):
    vs = set(e1) | set(e2) | set(e3)
    if len(vs) != 6: continue
    if all(D[a][b] >= 3 for ea, eb in itertools.combinations((e1, e2, e3), 2)
           for a in ea for b in eb):
        trip = (e1, e2, e3); break
adjW = [set(x) for x in H]
for (a, b) in trip: adjW[a].discard(b); adjW[b].discard(a)
V0, V1, V2, V3, A, B, A2, B2 = (30 + i for i in range(8))
for _ in range(8): adjW.append(set())
def link(u, v): adjW[u].add(v); adjW[v].add(u)
for u, v in ((V0, V1), (V1, V2), (V2, V3), (V0, B), (B, A), (A, V3),
             (V0, B2), (B2, A2), (A2, V3)): link(u, v)
for (t1, t2), (h1, h2) in (((V1, V2), trip[0]), ((B, A), trip[1]),
                           ((B2, A2), trip[2])):
    link(t1, h1); link(t2, h2)
W = [sorted(x) for x in adjW]
gW = girth(W)
print('  witness: n=%d, girth=%d  (460: three hexagons, overlaps {3})'
      % (len(W), gW))
print()
print('  chords of its z-cycles, and the winning distances 2^j - 1:')
print()
print('  %4s %6s %8s %26s %14s %6s'
      % ('z', 'c_z', 'chords', 'arc-distances seen', 'winning set', 'hit'))
for z in range(gW, 13):
    cs = cyc_of_len(W, z)
    if not cs: continue
    ds = Counter()
    for c in cs:
        for d in chords(W, c): ds[d] += 1
    win = pow2_dists(z)
    hit = sorted(set(ds) & win)
    print('  %4d %6d %8d %26s %14s %6s'
          % (z, len(cs), sum(ds.values()),
             str(sorted(ds.items()))[:26], sorted(win),
             'YES' if hit else 'no'))

print()
print('=' * 70)
print('(C) THE CHORD-DISTANCE DISTRIBUTION: IS IT SPREAD?')
print('=' * 70)
print('  For each z-cycle and each chord, record the arc-distance d.')
print('  A chord at d = 2^j - 1 (or z-d+1 = 2^j) yields a')
print('  power-of-two cycle outright.  Does the adversary avoid it?')
print()
print('  %-18s %4s %4s %6s %7s %9s %9s %6s'
      % ('graph', 'n', 'z', 'c_z', 'chords', 'distinct d', 'winning', 'hit'))


def report(name, G, zs):
    n = len(G)
    for z in zs:
        cs = cyc_of_len(G, z)
        if not cs: continue
        ds = Counter()
        for c in cs:
            for d in chords(G, c): ds[d] += 1
        if not ds: continue
        win = pow2_dists(z)
        hit = sorted(set(ds) & win)
        print('  %-18s %4d %4d %6d %7d %9d %9s %6s'
              % (name, n, z, len(cs), sum(ds.values()), len(ds),
                 str(sorted(win)), 'YES' if hit else 'NO'))


report('Petersen', gp(5, 2), range(5, 10))
report('Heawood', lcf(14, [5, -5]), range(6, 13))
report('Moebius-Kantor', gp(8, 3), range(6, 13))
report('Nauru', gp(12, 5), range(6, 13))
report('Desargues', gp(10, 3), range(6, 13))
report('460-witness', W, range(6, 13))

rng = random.Random(1618033)
for n in (20, 26, 34):
    G = random_cubic(n, rng, mingirth=5)
    if G: report('random n=%d' % n, G, range(girth(G), 13))

base = gp(5, 2)
for nb in (14, 20):
    Bk = random_cubic(nb, rng, mingirth=5)
    if Bk: report('Pet+rand%d' % nb, splice(base, Bk), range(5, 13))

print()
print('=' * 70)
print('(D) THE RIGIDITY THAT HORN TWO DEMANDS')
print('=' * 70)
print('  Horn two needs T_m >= 2^m/2m.  With the a=1 layer carrying an')
print('  O(1) share, that is ~2^{z+1}/z chord incidences on z-cycles,')
print('  z = m-2.  The winning distance set has size ~2 log2 z.')
print()
print('  %5s %16s %14s %12s'
      % ('z', 'chord incidences', '|winning set|', 'ratio'))
for k in (4, 5, 6, 7, 8):
    m = 2 ** k; z = m - 2
    inc = 2.0 ** (z + 1) / z
    w = len(pow2_dists(z))
    print('  %5d %16.4g %14d %12.4g' % (z, inc, w, inc / max(w, 1)))
print()
print('  So the adversary must place ~10^{...} chord incidences and')
print('  miss a set of size ~2 log2 z EVERY time.  That is the exact')
print('  statement horn two now reduces to at its dominant layer, and')
print('  it is a rigidity statement about arc-distances, not a')
print('  counting statement about kappa or (NC).')

print()
print('=' * 70)
print('(E) THE CHORD BAND LEMMA, AND THE a=1 CRITERION STANDING ALONE')
print('=' * 70)
print('  BAND LEMMA.  A chord of a z-cycle at arc-distance d yields')
print('  cycles d+1 and z-d+1, both >= g, so')
print('        g - 1  <=  d  <=  z - g + 1.')
print('  The legal band has width z - 2g + 3.  At z = 2g it collapses')
print('  to THREE values and the outcome is pure arithmetic -- which')
print('  is exactly the c*/delta* calculus of 457-459, rederived from')
print('  the energy side.  Measured above: the cages show ONE distinct')
print('  arc-distance at every z.')
print()
print('  a=1 CRITERION.  G has a 2^j-cycle whenever some cycle carries')
print('  a chord at arc-distance exactly 2^j - 1.  Does it ever fail?')
print()


def truncate(adj):
    n = len(adj); idx = {}; c = 0
    for v in range(n):
        for w in adj[v]: idx[(v, w)] = c; c += 1
    out = [set() for _ in range(c)]
    for v in range(n):
        for a, b in itertools.combinations(adj[v], 2):
            out[idx[(v, a)]].add(idx[(v, b)]); out[idx[(v, b)]].add(idx[(v, a)])
    for v in range(n):
        for w in adj[v]:
            out[idx[(v, w)]].add(idx[(w, v)]); out[idx[(w, v)]].add(idx[(v, w)])
    return [sorted(x) for x in out]


def has_cycle(adj, L):
    return len(cyc_of_len(adj, L)) > 0


def a1_hits(adj, zmax):
    """set of j such that some z-cycle has a chord at distance 2^j - 1."""
    out = set()
    for z in range(girth(adj), zmax + 1):
        for c in cyc_of_len(adj, z):
            for d in chords(adj, c):
                for dd in (d, z - d):
                    x = dd + 1
                    if x >= 3 and x & (x - 1) == 0: out.add(x)
    return out


print('  %-20s %4s %3s %6s %-16s %-16s %6s'
      % ('graph', 'n', 'g', 'zmax', 'a=1 yields', 'actual 2^k in G', 'ok'))
zoo2 = [('Petersen', gp(5, 2)), ('Heawood', lcf(14, [5, -5])),
        ('Moebius-Kantor', gp(8, 3)), ('Pappus',
         lcf(18, [5, 7, -7, 7, -7, -5])), ('Desargues', gp(10, 3)),
        ('Nauru', gp(12, 5)), ('460-witness', W),
        ('trunc-K4', truncate(lcf(4, [2, -2]))),
        ('trunc-cube', truncate(lcf(8, [3, -3]))),
        ('trunc-Petersen', truncate(gp(5, 2)))]
rng3 = random.Random(2718281)
for nn in (18, 24, 30):
    Gr = random_cubic(nn, rng3, mingirth=5)
    if Gr: zoo2.append(('random g>=5 n=%d' % nn, Gr))
    Gr3 = random_cubic(nn, rng3, mingirth=3)
    if Gr3: zoo2.append(('random g>=3 n=%d' % nn, Gr3))
for name, G in zoo2:
    n = len(G); g = girth(G); zmax = min(13, n)
    got = a1_hits(G, zmax)
    act = {2 ** j for j in range(2, 8)
           if 2 ** j <= zmax and has_cycle(G, 2 ** j)}
    print('  %-20s %4d %3d %6d %-16s %-16s %6s'
          % (name, n, g, zmax, str(sorted(got)) or '{}', str(sorted(act)),
             'FIRES' if got else 'MISS'))
print()
print('  Every MISS row is a graph on which the a=1 layer alone cannot')
print('  decide; the "actual" column says whether the graph is in fact')
print('  fine by some other mechanism.  A row with an empty "actual"')
print('  column and a small n would be a counterexample to E-G.')
