"""closure463.py -- IS THE DIVISIBILITY CLOSURE TRUE?

STATEMENT.  Let G be cubic, delta* = 2^{j+1} - g > 0, and
C* = {(g-d)/2 : d | delta*, d == g mod 2}.  If some pair of girth
cycles overlaps in >= c* edges, must some pair overlap in a member
of C*?

ARITHMETIC PRE-ANALYSIS.

 * delta* = 1  (g = 2^{j+1} - 1: g = 7, 15, 31, ...).  Then
   C* = {(g-1)/2} = {c*} and c <= floor(g/2) = (g-1)/2 = c*, so
   "overlap >= c*" IS "overlap = c*".  TRUE, trivially.

 * g odd.  d = 1 divides delta*, so C* contains floor(g/2), the
   largest overlap possible, as well as c*.

 * delta* = 2  (g = 2^{j+1} - 2: g = 6, 14, 30, ...).  C* = {g/2 - 1}
   is a SINGLE value, c* = g/2 - 1, and the permitted overlaps at or
   above c* are {g/2 - 1, g/2}.  An adversary holding EVERY overlap
   at exactly g/2 (delta = 0, the difference cycle is another girth
   cycle) would satisfy "some overlap >= c*" while missing C*.
   THIS IS THE ONLY ESCAPE, and g = 6 is its smallest instance:
   hexagons meeting in 3 edges but never in 2.

So the search is aimed: cubic, girth 6, overlaps in {1,3} but not 2.
A pair meeting in 3 spans Theta(3,3,3) - two branch vertices joined
by three paths of length 3 - whose three cycles are all hexagons
pairwise meeting in 3.  The question is whether such a theta can sit
in a girth-6 cubic graph with no hexagon pair meeting in exactly 2.
"""
import sys, random, itertools
from collections import Counter

sys.setrecursionlimit(100000)


def gp(n, k):
    if k % n == 0 or (2 * k) % n == 0: return None
    adj = {i: set() for i in range(2 * n)}
    def add(a, b):
        if a == b: return
        adj[a].add(b); adj[b].add(a)
    for i in range(n):
        add(i, (i + 1) % n); add(i, n + i)
        add(n + i, n + (i + k) % n)
    A = [sorted(adj[i]) for i in range(2 * n)]
    if any(len(x) != 3 for x in A): return None
    return A


def lcf(n, seq):
    adj = {i: set() for i in range(n)}
    for i in range(n):
        adj[i].add((i + 1) % n); adj[(i + 1) % n].add(i)
    for i in range(n):
        j = (i + seq[i % len(seq)]) % n
        if j == i: return None
        adj[i].add(j); adj[j].add(i)
    A = [sorted(adj[i]) for i in range(n)]
    if any(len(x) != 3 for x in A): return None
    return A


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
    ds = tgt - g
    if ds == 0: return set(), 0
    return ({(g - d) // 2 for d in range(1, ds + 1)
             if ds % d == 0 and (d - g) % 2 == 0
             and 1 <= (g - d) // 2 <= g // 2}, ds)


def verdict(adj):
    """returns (g, delta*, c*, C*, O, status)."""
    g = girth(adj)
    C, ds = targets(g)
    if ds == 0: return g, 0, None, C, None, 'degenerate'
    cs, _ = cstar(g)
    G = cycles_of_length(adj, g)
    O = set()
    for A, B in itertools.combinations(G, 2):
        if A & B: O.add(len(A & B))
    if not O: return g, ds, cs, C, O, 'disjoint-girthed'
    if max(O) < cs: return g, ds, cs, C, O, 'below threshold'
    if O & C: return g, ds, cs, C, O, 'closure holds'
    return g, ds, cs, C, O, '*** COUNTEREXAMPLE ***'


def random_cubic(n, rng, mingirth=3, tries=20000):
    for _ in range(tries):
        stubs = [v for v in range(n) for _ in range(3)]
        rng.shuffle(stubs)
        adj = [set() for _ in range(n)]
        ok = True
        for i in range(0, len(stubs), 2):
            a, b = stubs[i], stubs[i + 1]
            if a == b or b in adj[a]: ok = False; break
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


print('=== (1) SWEEP: generalized Petersen graphs ===', flush=True)
tally = Counter(); ces = []
for nn in range(5, 25):
    for kk in range(2, nn // 2 + 1):
        A = gp(nn, kk)
        if A is None: continue
        g, ds, cs, C, O, st = verdict(A)
        tally[st] += 1
        if 'COUNTER' in st:
            ces.append(('GP(%d,%d)' % (nn, kk), g, ds, cs, C, O))
            print('  GP(%d,%d): g=%d d*=%d c*=%d C*=%s O=%s  %s'
                  % (nn, kk, g, ds, cs, sorted(C), sorted(O), st),
                  flush=True)
print('  ', dict(tally), flush=True)

print()
print('=== (2) SWEEP: LCF circulant-style cubic graphs ===',
      flush=True)
tally2 = Counter()
rng0 = random.Random(2718)
seqs = []
for n in (14, 16, 18, 20, 22, 24, 26, 28, 30):
    for _ in range(60):
        L = rng0.choice([2, 3, 4, 6])
        s = [rng0.choice([x for x in range(-n // 2 + 1, n // 2)
                          if abs(x) > 1]) for _ in range(L)]
        seqs.append((n, s))
for n, s in seqs:
    A = lcf(n, s)
    if A is None: continue
    g, ds, cs, C, O, st = verdict(A)
    tally2[st] += 1
    if 'COUNTER' in st:
        ces.append(('LCF%d%s' % (n, s), g, ds, cs, C, O))
        print('  LCF %d %s: g=%d d*=%d c*=%d C*=%s O=%s  %s'
              % (n, s, g, ds, cs, sorted(C), sorted(O), st), flush=True)
print('  ', dict(tally2), flush=True)

print()
print('=== (3) SWEEP: random cubic graphs, girth filtered ===',
      flush=True)
rng = random.Random(31415)
for mg in (5, 6):
    tally3 = Counter()
    for _ in range(120):
        A = random_cubic(rng.choice([16, 20, 24, 28, 30]), rng,
                         mingirth=mg)
        if A is None: continue
        g, ds, cs, C, O, st = verdict(A)
        tally3[st] += 1
        if 'COUNTER' in st:
            ces.append(('random n=%d' % len(A), g, ds, cs, C, O))
            print('  random n=%d: g=%d d*=%d c*=%d C*=%s O=%s  %s'
                  % (len(A), g, ds, cs, sorted(C), sorted(O), st),
                  flush=True)
    print('  girth >= %d : %s' % (mg, dict(tally3)), flush=True)

print()
print('=== (4) THE delta* = 2 ESCAPE, HUNTED DIRECTLY ===',
      flush=True)
print('  need: girth 6, some hexagon pair meeting in 3, none in 2.')
found = 0; checked = 0
rng2 = random.Random(161803)
pool = []
for nn in range(7, 25):
    for kk in range(2, nn // 2 + 1):
        A = gp(nn, kk)
        if A is not None and girth(A) == 6: pool.append(A)
for _ in range(400):
    A = random_cubic(rng2.choice([18, 20, 24, 28, 30]), rng2,
                     mingirth=6)
    if A is not None: pool.append(A)
for A in pool:
    g = girth(A)
    if g != 6: continue
    checked += 1
    G = cycles_of_length(A, 6)
    O = set()
    for X, Y in itertools.combinations(G, 2):
        if X & Y: O.add(len(X & Y))
    if 3 in O and 2 not in O:
        found += 1
        print('  ESCAPE FOUND: n=%d  O=%s' % (len(A), sorted(O)),
              flush=True)
print('  girth-6 cubic graphs checked: %d   escapes found: %d'
      % (checked, found), flush=True)

print()
print('=== (5) VERDICT ===', flush=True)
print('  counterexamples to the Divisibility Closure across all'
      ' sweeps: %d' % len(ces), flush=True)
for c in ces[:10]:
    print('   ', c, flush=True)
