"""sieve454.py -- THE GEODESIC SIEVE.

spectral453 established, for a CUBIC graph:

    tr(T^m) = SUM_lam p_m(lam) + (n/2)(1 + (-1)^m)          (exact)
    tr(T^m) = 2m * c_m                       for m < 2g     (exact)

so below 2g the cycle counts are a pure spectral quantity.  Above 2g
the trace picks up NON-SIMPLE closed geodesics.  This file measures
exactly what they are, so the identity can be continued upward.

CLASSIFICATION.  The support of a closed non-backtracking tailless
walk is a connected subgraph of cycle rank r >= 1.

  r = 1 : the support is a single cycle of length l, traversed j
          times, m = j*l.  Marked-start count 2l per unoriented
          cycle, so the total rank-1 contribution is
                SUM_{j | m, l = m/j >= g}  2 (m/j) c_{m/j}.
          The j = 1 term is the 2m c_m we want; j >= 2 needs m >= 2g.

  r = 2 : the support is a THETA (two degree-3 branch vertices, three
          internally disjoint paths a,b,c) or a DUMBBELL (two
          vertex-disjoint cycles joined by a path).  In a cubic graph
          both branch vertices of a theta are saturated, so a
          traversal uses every edge exactly twice and has length
          2(a+b+c); a dumbbell traversal has length |C1|+|C2|+2*(bar
          length).  Shortest possible: 2*(3g/2) = 3g for a theta with
          all three cycles of length g... measured below.

  r >= 3: longer still.

THE POINT.  The correction terms are all counts of SHORT CYCLES
(pairs of them, joined).  In the archive's minimal counterexample
short cycles are exactly what is scarce -- the spectrum caps, the
island bounds.  So the adversary faces a dichotomy that has nothing
to do with ear geometry:

    many short cycles  ->  the correction terms are large, but they
                           are themselves cycles, and their lengths
                           are constrained by the same congruences;
    few short cycles   ->  the sieve converges and the 2^m main term
                           of tr(T^m) forces c_{2^k} > 0.

This file verifies the classification exactly by enumerating every
closed non-backtracking tailless walk.
"""
import sys
sys.setrecursionlimit(100000)


def gp(n, k):
    adj = {i: [] for i in range(2 * n)}
    def add(a, b):
        adj[a].append(b); adj[b].append(a)
    for i in range(n):
        add(i, (i + 1) % n); add(i, n + i)
        add(n + i, n + (i + k) % n)
    return [sorted(set(adj[i])) for i in range(2 * n)]


def k4():
    return [[1, 2, 3], [0, 2, 3], [0, 1, 3], [0, 1, 2]]


def k33():
    adj = [[] for _ in range(6)]
    for i in range(3):
        for j in range(3, 6):
            adj[i].append(j); adj[j].append(i)
    return adj


def heawood():
    adj = {i: [] for i in range(14)}
    def add(a, b):
        adj[a].append(b); adj[b].append(a)
    for i in range(14):
        add(i, (i + 1) % 14)
    for i in range(0, 14, 2):
        add(i, (i + 5) % 14)
    return [sorted(set(adj[i])) for i in range(14)]


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


def cycle_counts(adj, cap):
    n = len(adj); cnt = [0] * (cap + 1)
    for s in range(n):
        seen = {s}
        def dfs(v, ln):
            for w in adj[v]:
                if w == s:
                    if ln >= 3: cnt[ln] += 1
                elif w > s and w not in seen and ln < cap:
                    seen.add(w); dfs(w, ln + 1); seen.discard(w)
        dfs(s, 1)
    return [c // 2 for c in cnt]


def matmul(X, Y):
    Yt = list(zip(*Y))
    return [[sum(a * b for a, b in zip(r, c)) for c in Yt] for r in X]


def spectral_traces(adj, M):
    n = len(adj)
    A = [[0] * n for _ in range(n)]
    for v in range(n):
        for w in adj[v]:
            A[v][w] = 1
    Pm2 = [[2 if i == j else 0 for j in range(n)] for i in range(n)]
    Pm1 = [r[:] for r in A]
    ex = n // 2
    out = [None, sum(Pm1[i][i] for i in range(n))]
    for m in range(2, M + 1):
        Cur = matmul(A, Pm1)
        for i in range(n):
            for j in range(n):
                Cur[i][j] -= 2 * Pm2[i][j]
        Pm2, Pm1 = Pm1, Cur
        out.append(sum(Pm1[i][i] for i in range(n)) + ex * (1 + (-1) ** m))
    return out


def walk_census(adj, m):
    """Enumerate every closed non-backtracking TAILLESS walk of length
    m with a marked start; classify by cycle rank of its support.
    Returns (total, {rank: count})."""
    n = len(adj)
    total = 0
    byrank = {}
    for v0 in range(n):
        for v1 in adj[v0]:
            seq = [v0, v1]
            edges = {frozenset((v0, v1))}
            counts = {frozenset((v0, v1)): 1}

            def rec(depth):
                nonlocal total
                v = seq[-1]; u = seq[-2]
                if depth == m:
                    return
                for w in adj[v]:
                    if w == u:
                        continue
                    e = frozenset((v, w))
                    if depth + 1 == m:
                        if w == v0 and seq[1] != v:
                            counts[e] = counts.get(e, 0) + 1
                            verts = set()
                            for f in counts:
                                verts |= set(f)
                            r = len(counts) - len(verts) + 1
                            byrank[r] = byrank.get(r, 0) + 1
                            total += 1
                            counts[e] -= 1
                            if counts[e] == 0: del counts[e]
                        continue
                    seq.append(w)
                    counts[e] = counts.get(e, 0) + 1
                    rec(depth + 1)
                    counts[e] -= 1
                    if counts[e] == 0: del counts[e]
                    seq.pop()

            rec(1)
    return total, byrank


GRAPHS = [('K4', k4()), ('K3,3', k33()), ('prism GP(3,1)', gp(3, 1)),
          ('cube GP(4,1)', gp(4, 1)), ('Petersen GP(5,2)', gp(5, 2)),
          ('Heawood', heawood()),
          ('Moebius-Kantor GP(8,3)', gp(8, 3))]

print('=== (1) EXACT CENSUS OF CLOSED GEODESICS BY SUPPORT RANK ===',
      flush=True)
print('  graph                 g   m   tr(T^m)  rank1   rank2  rank3+'
      '   2m*c_m  rank1 predicted')
for name, adj in GRAPHS:
    g = girth(adj)
    M = min(3 * g - 1, 13)
    st = spectral_traces(adj, M)
    cc = cycle_counts(adj, M)
    for m in range(g, M + 1):
        tot, br = walk_census(adj, m)
        r1 = br.get(1, 0); r2 = br.get(2, 0)
        r3 = sum(v for k, v in br.items() if k >= 3)
        pred1 = 0
        for j in range(1, m + 1):
            if m % j == 0:
                l = m // j
                if l < len(cc): pred1 += 2 * l * cc[l]
        flag = 'OK' if (tot == st[m] and r1 == pred1) else 'FAIL'
        print('  %-20s %3d %3d %9d %7d %7d %7d %8d %10d  %s'
              % (name, g, m, st[m], r1, r2, r3, 2 * m * cc[m],
                 pred1, flag), flush=True)

print()
print('=== (2) WHERE THE FIRST NON-SIMPLE GEODESIC APPEARS ===',
      flush=True)
for name, adj in GRAPHS:
    g = girth(adj)
    first = None
    for m in range(3, min(3 * g, 14)):
        tot, br = walk_census(adj, m)
        if sum(v for k, v in br.items() if k >= 2) > 0 or \
           (br.get(1, 0) and any(m % j == 0 and j >= 2 and
                                 m // j >= g for j in range(2, m + 1))):
            first = m; break
    print('  %-20s g=%d  2g=%d  first non-simple geodesic at m=%s'
          % (name, g, 2 * g, first), flush=True)

print()
print('=== (3) SIZE OF THE CORRECTION vs THE MAIN TERM 2^m ===',
      flush=True)
print('  graph                 m   2m*c_m    corrections   2^m'
      '       corr/2^m')
for name, adj in GRAPHS:
    g = girth(adj)
    M = min(3 * g - 1, 13)
    cc = cycle_counts(adj, M)
    for m in range(2 * g, M + 1):
        tot, br = walk_census(adj, m)
        corr = tot - 2 * m * cc[m]
        print('  %-20s %3d %8d %13d %9d %11.4f'
              % (name, m, 2 * m * cc[m], corr, 2 ** m,
                 corr / 2.0 ** m), flush=True)
