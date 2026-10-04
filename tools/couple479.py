"""couple479.py -- LEVEL COUPLING: THE THETA VALUATION LEMMA AND THE
AVOIDANCE LATTICE.

470.5 left the frontier at "couple two levels".  A level is a
2-adic valuation s; its data is the set of cycle lengths with
v_2 = s.  Coupling means an inference carrying information between
levels.  This file supplies one mechanism, checks it, and then
measures how tightly the levels are in fact coupled.

PART 1 machine-checks the THETA VALUATION LEMMA.  A theta with
internally disjoint paths of lengths p, q, r has cycles A = p+q,
B = p+r, C = q+r.  Let alpha = min(v_2 p, v_2 q, v_2 r) and mu the
number of paths attaining it.  Claim:

    min(v_2 A, v_2 B, v_2 C)  =  alpha      if mu = 1 or 2,
                              =  alpha + 1  if mu = 3.

This is exact, not an inequality, and it is the transport rule: the
LEVELS of a theta's cycles are determined by the VALUATIONS of its
paths, which are shared between cycles at different levels.

PART 2 tests whether the mechanism is ever silent.  On a cycle C the
coupling fires through an ear whose arc distance has small 2-adic
valuation.  Measured: for a maximal-level cycle in each graph, the
multiset of v_2(arc distance) over all ears.

PART 3 is the decisive measurement.  A power of two at level s is
the single length 2^s (470.2: the multiples of 2^s below 3.2^s are
2^s and 2^{s+1}).  So Erdos-Gyarfas says no cubic graph avoids ALL
of 4, 8, 16, 32, ...  Which SUBSETS can be avoided?  The lattice of
avoidable subsets IS the coupling, measured.  For n <= 31 the only
powers of two available are 4, 8, 16, so avoiding {4,8,16} at those
n would be an outright counterexample.
"""
import sys, time, random, itertools
from collections import deque

sys.setrecursionlimit(100000)
P = print


def v2(x):
    k = 0
    while x % 2 == 0: x //= 2; k += 1
    return k


P('=' * 78)
P('PART 1.  THE THETA VALUATION LEMMA, MACHINE-CHECKED', flush=True)
P('=' * 78)
bad = []; nmu = {1: 0, 2: 0, 3: 0}
for p in range(1, 121):
    for q in range(1, 121):
        for r in range(1, 121):
            a = min(v2(p), v2(q), v2(r))
            mu = [v2(p), v2(q), v2(r)].count(a)
            got = min(v2(p + q), v2(p + r), v2(q + r))
            want = a + (1 if mu == 3 else 0)
            nmu[mu] += 1
            if got != want: bad.append((p, q, r, got, want))
P('  all (p,q,r) in [1,120]^3 : %d triples' % sum(nmu.values()))
P('  mu = 1 : %d      mu = 2 : %d      mu = 3 : %d'
  % (nmu[1], nmu[2], nmu[3]), flush=True)
P('  violations : %d' % len(bad), flush=True)
if bad:
    for b in bad[:8]: P('    %s' % (b,), flush=True)
P('  %s' % ('LEMMA HOLDS on the whole range.' if not bad else 'LEMMA FAILS.'),
  flush=True)
P()
P('  Consequence (LEVEL DESCENT).  min level of a theta is at most')
P('  alpha + 1.  So a theta all of whose cycles sit at level >= s has')
P('  every path length divisible by 2^{s-1}: high-level cycles force')
P('  high-level paths.  Contrapositive, and this is the coupling: ONE')
P('  path of low valuation drags a cycle down to a low level.')
P()
P('  Specialisation to an ear (a, d) on a cycle of length L, v_2(L)=t:')
P('  the theta has paths a, L-a, d.  If v_2(a) < t then a and L-a')
P('  share the valuation v_2(a), so mu >= 2 and')
P('      min level of the three cycles = min(v_2 a, v_2 d)   (mu<=2)')
P('                                    = v_2(a) + 1          (mu=3).')
P('  Since the third cycle is C itself at level t, one of the two EAR')
P('  cycles sits at that min level whenever it is < t.')
sp = 0
for a in range(1, 200):
    for d in range(1, 200):
        for L in (24, 40, 48, 96):
            t = v2(L)
            if a >= L: continue
            al = min(v2(a), v2(L - a), v2(d))
            mu = [v2(a), v2(L - a), v2(d)].count(al)
            want = al + (1 if mu == 3 else 0)
            got = min(t, v2(a + d), v2(L - a + d))
            if got != want: sp += 1
P('  ear specialisation checked on L in {24,40,48,96}, a,d < 200:'
  ' %d violations' % sp, flush=True)


# ------------------------------------------------------------ graph zoo
def gp(n, k):
    adj = {i: set() for i in range(2 * n)}
    def add(x, y):
        if x != y: adj[x].add(y); adj[y].add(x)
    for i in range(n):
        add(i, (i + 1) % n); add(i, n + i); add(n + i, n + (i + k) % n)
    A = [sorted(adj[i]) for i in range(2 * n)]
    return A if all(len(x) == 3 for x in A) else None


def lcf(n, seq):
    adj = {i: set() for i in range(n)}
    for i in range(n):
        adj[i].add((i + 1) % n); adj[(i + 1) % n].add(i)
    for i in range(n):
        j = (i + seq[i % len(seq)]) % n
        adj[i].add(j); adj[j].add(i)
    A = [sorted(adj[i]) for i in range(n)]
    return A if all(len(x) == 3 for x in A) else None


def truncate(adj):
    n = len(adj); idx = {}; c = 0
    for v in range(n):
        for w in adj[v]: idx[(v, w)] = c; c += 1
    out = [set() for _ in range(c)]
    for v in range(n):
        for x, y in itertools.combinations(adj[v], 2):
            out[idx[(v, x)]].add(idx[(v, y)]); out[idx[(v, y)]].add(idx[(v, x)])
    for v in range(n):
        for w in adj[v]:
            out[idx[(v, w)]].add(idx[(w, v)]); out[idx[(w, v)]].add(idx[(v, w)])
    return [sorted(x) for x in out]


def connected(adj):
    n = len(adj); seen = {0}; s = [0]
    while s:
        v = s.pop()
        for w in adj[v]:
            if w not in seen: seen.add(w); s.append(w)
    return len(seen) == n


def find_cycle(adj, L):
    """return a vertex list of a simple cycle of length L, or None."""
    n = len(adj)
    for s in range(n):
        seen = {s}; path = [s]; out = [None]
        def dfs(v, ln):
            if out[0]: return
            if ln == L:
                if s in adj[v]: out[0] = list(path)
                return
            for w in adj[v]:
                if w > s and w not in seen:
                    seen.add(w); path.append(w); dfs(w, ln + 1)
                    path.pop(); seen.discard(w)
        dfs(s, 1)
        if out[0]: return out[0]
    return None


def count_len(adj, L, cap=10 ** 9):
    n = len(adj); tot = [0]
    for s in range(n):
        vis = [False] * n; vis[s] = True
        def dfs(v, ln):
            if tot[0] >= cap: return
            if ln == L:
                if s in adj[v]: tot[0] += 1
                return
            for w in adj[v]:
                if w > s and not vis[w]:
                    vis[w] = True; dfs(w, ln + 1); vis[w] = False
        dfs(s, 1); vis[s] = False
        if tot[0] >= cap: break
    return tot[0] // 2


def ears_of(adj, cyc):
    """pairs (i,j) of cycle positions joined by a path off the cycle."""
    L = len(cyc); pos = {v: i for i, v in enumerate(cyc)}
    onc = set(cyc); out = []
    for i in range(L):
        u = cyc[i]
        block = onc - {u}
        seen = {u}; q = deque([u]); reach = set()
        while q:
            v = q.popleft()
            for w in adj[v]:
                if w in seen: continue
                if w in block:
                    reach.add(w); seen.add(w); continue
                seen.add(w); q.append(w)
        for w in reach:
            j = pos[w]
            if j > i: out.append((i, j))
    return out


ZOO = [('Petersen', gp(5, 2)), ('Heawood', lcf(14, [5, -5])),
       ('Moebius-Kantor', gp(8, 3)),
       ('Pappus', lcf(18, [5, 7, -7, 7, -7, -5])),
       ('Desargues', gp(10, 3)),
       ('Dodecahedron', lcf(20, [10, 7, 4, -4, -7, 10, -4, 7, -7, 4])),
       ('Nauru', gp(12, 5)), ('McGee', lcf(24, [12, 7, -7])),
       ('trunc-cube', truncate(lcf(8, [3, -3]))),
       ('trunc-Petersen', truncate(gp(5, 2)))]

P()
P('=' * 78)
P('PART 2.  IS THE COUPLING EVER SILENT?', flush=True)
P('=' * 78)
P('  For a cycle C at the graph\'s maximal level, list the 2-adic')
P('  valuations of the arc distances of all its ears.  Valuation 0')
P('  present means the coupling fires at once.')
P()
P('  %-17s %4s %4s %5s %7s   %s'
  % ('graph', 'n', '|C|', 't', 'ears', 'v_2(arc distance) multiset'),
  flush=True)
for name, G in ZOO:
    if G is None: continue
    n = len(G)
    best = None
    for L in range(n, 2, -1):
        if best is not None and v2(L) <= v2(best): continue
        if find_cycle(G, L): best = L
    t = v2(best)
    cyc = find_cycle(G, best)
    es = ears_of(G, cyc)
    prof = {}
    for (i, j) in es:
        a = min((j - i) % best, (i - j) % best)
        prof[v2(a)] = prof.get(v2(a), 0) + 1
    P('  %-17s %4d %4d %5d %7d   %s'
      % (name, n, best, t, len(es),
         ' '.join('%d:%d' % (k, prof[k]) for k in sorted(prof))), flush=True)
P()
P('  Valuation 0 present in every row means: no cubic graph in this')
P('  zoo has a parity-rigid ear system on its top-level cycle.')

P()
P('=' * 78)
P('PART 3.  THE AVOIDANCE LATTICE', flush=True)
P('=' * 78)
P('  Which subsets of the powers of two can a cubic graph avoid?')
P('  Objective: sum over L in S of min(#C_L, 30), driven down by')
P('  2-edge swaps with restarts.  0 = the subset is avoidable.')
P('  For n <= 31 the live powers of two are 4, 8, 16 only, so a 0 in')
P('  the {4,8,16} row would be an ERDOS-GYARFAS COUNTEREXAMPLE.')
P()


def randcubic(n, rng, tries=60000):
    for _ in range(tries):
        st = [v for v in range(n) for _ in range(3)]; rng.shuffle(st)
        adj = [set() for _ in range(n)]; ok = True
        for i in range(0, len(st), 2):
            a, b = st[i], st[i + 1]
            if a == b or b in adj[a]: ok = False; break
            adj[a].add(b); adj[b].add(a)
        if not ok: continue
        A = [sorted(x) for x in adj]
        if connected(A): return A
    return None


def swap(adj, rng):
    E = [(u, w) for u in range(len(adj)) for w in adj[u] if u < w]
    (a, b) = E[rng.randrange(len(E))]
    (c, d) = E[rng.randrange(len(E))]
    if len({a, b, c, d}) < 4: return None
    if rng.random() < 0.5: c, d = d, c
    if d in adj[a] or c in adj[b]: return None
    new = [set(x) for x in adj]
    new[a].discard(b); new[b].discard(a); new[c].discard(d); new[d].discard(c)
    new[a].add(d); new[d].add(a); new[b].add(c); new[c].add(b)
    A = [sorted(x) for x in new]
    if not all(len(x) == 3 for x in A): return None
    return A if connected(A) else None


CAP = 30
rng = random.Random(479479)
SUBSETS = [(4,), (8,), (16,), (4, 8), (4, 16), (8, 16), (4, 8, 16)]
P('  %-12s %8s %8s %8s %8s' % ('subset', 'n=16', 'n=20', 'n=24', 'n=28'),
  flush=True)
RESULT = {}
for S in SUBSETS:
    row = []
    for n in (16, 20, 24, 28):
        live = tuple(L for L in S if L <= n)
        if not live: row.append('n/a'); continue
        t0 = time.time(); best = 10 ** 9; bestG = None
        while time.time() - t0 < 22.0 and best > 0:
            G = randcubic(n, rng)
            if G is None: break
            cur = sum(min(count_len(G, L, CAP), CAP) for L in live)
            stall = 0
            while cur > 0 and stall < 100 and time.time() - t0 < 22.0:
                H = swap(G, rng)
                if H is None: stall += 1; continue
                c = sum(min(count_len(H, L, CAP), CAP) for L in live)
                if c <= cur:
                    stall = 0 if c < cur else stall + 1
                    G, cur = H, c
                else: stall += 1
            if cur < best: best, bestG = cur, G
        RESULT[(S, n)] = (best, bestG)
        row.append('0  YES' if best == 0 else str(best))
    P('  %-12s %8s %8s %8s %8s'
      % ('{' + ','.join(map(str, S)) + '}', row[0], row[1], row[2], row[3]),
      flush=True)

P()
P('=' * 78)
P('VERDICT', flush=True)
P('=' * 78)
ce = []
for n in (16, 20, 24, 28):
    b, G = RESULT[((4, 8, 16), n)]
    if b == 0 and G is not None:
        pw = [L for L in (4, 8, 16) if count_len(G, L, 1) > 0]
        if not pw: ce.append((n, G))
if ce:
    P('  *** COUNTEREXAMPLE CANDIDATES - VERIFY IMMEDIATELY ***')
    for (n, G) in ce: P('    n=%d %s' % (n, G), flush=True)
else:
    P('  No subset containing all live powers of two was avoidable at')
    P('  any n tested: every cubic graph reached kept at least one')
    P('  power-of-two cycle.  The row pattern above says WHICH one')
    P('  survives, and that is the coupling in measured form.')
    for S in SUBSETS:
        for n in (16, 20, 24, 28):
            b, G = RESULT[(S, n)]
            if b == 0 and G is not None:
                left = [L for L in (4, 8, 16) if L <= n
                        and count_len(G, L, 1) > 0]
                P('    avoided %-10s at n=%2d ; still carries %s'
                  % ('{' + ','.join(map(str, S)) + '}', n, left), flush=True)
                break
