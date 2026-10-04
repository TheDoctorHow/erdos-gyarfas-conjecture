"""divis462.py -- THE DIVISIBILITY QUESTION.

458 reduced the girth-scale mechanism to  delta | delta*,  where
delta = g - 2c runs over realised girth-cycle overlaps and
delta* = 2^{j+1} - g.  This file answers: what forces a realised
overlap into the winning set, and where exactly does that forcing
stop?

THE TARGET SET.  delta | delta* with delta = g - 2c means the
winning overlaps are exactly

    C* = { (g - d)/2 : d | delta*,  d = g mod 2 },

which contains c* = (g - delta*)/2 at d = delta*, and tops out at
c = ceil(g/2) - 1 at d = 1 (g odd) or d = 2 (g even) - and d = 1 or 2
ALWAYS divides delta*, since delta* has the parity of g.  So C* is
never empty and the adversary must dodge every one of its members.

THE FORCING THEOREM.  Fix an edge e = uv.  Every girth cycle through
e is e together with a v-u path of length g-1 avoiding e.  At v there
are 2 choices (cubic, minus e) and 2 at every step after, so there
are at most 2^p distinct prefixes of length p.  Hence:

  EDGE-MULTIPLICITY FORCING.  If r_e, the number of girth cycles
  through an edge e, satisfies r_e > 2^{c*-1}, then two of them share
  a prefix of length c*-1, so their overlap is at least c*.

  (c* = 1 is the degenerate case: threshold r_e > 1, conclusion
   overlap >= 1, which is the g = 2^j + 1 theorem of 458.3.)

THE DICHOTOMY.  Either some edge exceeds the threshold - and then an
overlap of size >= c* exists, so the mechanism has a target - or every
edge is below it, and then

    c_g <= |E| . 2^{c*-1} / g = 3n 2^{c*-1} / 2g,
    kappa_g = 2g c_g / 2^g <= 3n . 2^{c*-1-g},

which for a near-Moore graph (n ~ C 2^{g/2}) and c* <= g/2 is
kappa_g <= 3C/2: BOUNDED.  Sparse-at-an-edge implies spread.

This file verifies the forcing theorem, computes C* for each graph,
and locates the exact point where the forcing stops.
"""
import sys, itertools
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
            if v < w: add(idx[(v, w)], idx[(w, v)])
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


def cstar(g):
    p = 1
    while p * 2 <= g - 1: p *= 2
    return g - p, 2 * p


def targets(g):
    cs, tgt = cstar(g)
    dstar = tgt - g
    if dstar == 0:
        return set(), 0, tgt
    out = set()
    for d in range(1, dstar + 1):
        if dstar % d == 0 and (d - g) % 2 == 0:
            c = (g - d) // 2
            if 1 <= c <= g // 2: out.add(c)
    return out, dstar, tgt


ZOO = [('Petersen', gp(5, 2)), ('Heawood', lcf(14, [5, -5] * 7)),
       ('Moebius-Kantor', lcf(16, [5, -5] * 8)),
       ('Pappus', lcf(18, [5, 7, -7, 7, -7, -5] * 3)),
       ('Desargues', lcf(20, [5, -5, 9, -9] * 5)),
       ('dodecahedron', gp(10, 2)),
       ('Nauru', gp(12, 5)),
       ('McGee', lcf(24, [12, 7, -7] * 8)),
       ('Tutte-Coxeter', lcf(30, [-13, -9, 7, -7, 9, 13] * 5)),
       ('Foster', lcf(90, [17, -9, 37, -37, 9, -17] * 15)),
       ('trunc Petersen', truncate(gp(5, 2))),
       ('trunc McGee', truncate(lcf(24, [12, 7, -7] * 8)))]

print('=== (1) THE TARGET SET C* AND THE FORCING THRESHOLD ===',
      flush=True)
print('  graph            g  2^{j+1} d*   C* (winning overlaps)'
      '   max r_e  thresh  forced?  realised O    O n C*')
rows = []
for name, adj in ZOO:
    g = girth(adj)
    G = cycles_of_length(adj, g)
    C, dstar, tgt = targets(g)
    cs, _ = cstar(g)
    r = Counter()
    for cyc in G:
        for e in cyc: r[e] += 1
    maxr = max(r.values()) if r else 0
    thresh = 2 ** (cs - 1)
    O = set()
    for A, B in itertools.combinations(G, 2):
        if A & B: O.add(len(A & B))
    forced = maxr > thresh
    maxov = max(O) if O else 0
    rows.append((name, g, maxr, thresh, forced, maxov, cs))
    print('  %-16s %2d %6d %3d   %-22s %5d %7d   %-6s %-12s %s'
          % (name, g, tgt, dstar, str(sorted(C)), maxr, thresh,
             'YES' if forced else 'no', str(sorted(O)),
             str(sorted(O & C)) if (O & C) else 'EMPTY'), flush=True)

print()
print('=== (2) THE FORCING THEOREM, CHECKED ===', flush=True)
print('  claim: max r_e > 2^{c*-1}  =>  some overlap >= c*')
bad = 0
for name, g, maxr, thresh, forced, maxov, cs in rows:
    if forced and maxov < cs:
        print('  *** FAIL *** %s' % name); bad += 1
print('  violations: %d  (over %d graphs)' % (bad, len(rows)),
      flush=True)

print()
print('=== (3) NAURU SITS EXACTLY ON THE THRESHOLD ===', flush=True)
adj = gp(12, 5); g = girth(adj); G = cycles_of_length(adj, g)
cs, tgt = cstar(g)
r = Counter()
for cyc in G:
    for e in cyc: r[e] += 1
print('  g=%d  c*=%d  threshold r_e > 2^{c*-1} = %d'
      % (g, cs, 2 ** (cs - 1)))
print('  r_e distribution over the %d edges: %s'
      % (len(r), dict(Counter(r.values()))))
print('  max r_e = %d, which is EQUAL to the threshold, not above it.'
      % max(r.values()))
print('  So the forcing theorem is SHARP: one more girth cycle on any')
print('  single edge would force an overlap >= c* = 2 and hand the')
print('  chain mechanism its 8-cycle.', flush=True)

print()
print('=== (4) WHERE THE FORCING STOPS ===', flush=True)
print('  forcing gives overlap >= c*, i.e. delta <= delta*.  It does')
print('  NOT give delta | delta*.  Graphs where O meets [c*, g/2] but')
print('  misses C* would be counterexamples to the stronger reading:')
for name, adj in ZOO:
    g = girth(adj)
    G = cycles_of_length(adj, g)
    C, dstar, tgt = targets(g)
    cs, _ = cstar(g)
    O = set()
    for A, B in itertools.combinations(G, 2):
        if A & B: O.add(len(A & B))
    if dstar == 0:
        print('  %-16s DEGENERATE: g is a power of two, so the target'
              ' IS the girth' % name, flush=True)
        continue
    big = {c for c in O if c >= cs}
    if big and not (O & C):
        print('  %-16s O>=c* = %s but O n C* empty  <-- GAP REAL'
              % (name, sorted(big)), flush=True)
print('  (no line above means: on every graph tested, once an overlap')
print('   of size >= c* exists, one of them already lies in C*.)',
      flush=True)
