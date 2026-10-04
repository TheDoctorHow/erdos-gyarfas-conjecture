"""gf456.py -- J7 IN GENERATING-FUNCTION LANGUAGE: what casts, what
does not, and one existence theorem that falls out for free.

THREE EXACT CASTINGS.

(GF1) GLOBAL, the Ihara zeta.  1/zeta_G(u) = (1-u^2)^{n/2}
      det(I - uA + 2u^2 I), and

          SUM_m N_m u^m  =  -u d/du log zeta_G(u),   N_m = tr(T^m),
          N_m = SUM_{l | m} l P_l,   P_l = # prime geodesic classes.

      Moebius inversion gives P_l, and by 452.2  P_l = 2 c_l for
      l < 2g.  So the cycle counting function IS the prime counting
      function of an explicit rational function.  Verified in (A).

(GF2) ARITHMETIC, the roots-of-unity filter.  The condition
      "length == 0 mod 2^t" is a finite Fourier condition, and the
      filter is available:

          #{Z even : 2^k divides |Z|} = 2^{-k} SUM_{j<2^k} E_G(w^j),
          E_G(y) = SUM_{Z even} y^{|Z|},   w = exp(2 pi i / 2^k).

      This is exactly 443's Counting Identity, written in the
      standard notation.  Because the condition is a congruence, the
      GF frame is not merely descriptive here - AX'S THEOREM applies
      and yields an unconditional existence result.  See (B).

(GF3) LOCAL, the ear conflict graph.  By the One-Step Theorem
      (446.2) the lengths reachable in one step from the shallow
      cycle C are  L + SUM_{i in S} delta_i  over independent sets S
      of the conflict graph H on the ears.  So

          J7  <=>  the weighted independence polynomial
                   Z_H(u) = SUM_{S indep} u^{SUM_{i in S} delta_i}
                   misses an explicit finite set of exponents.

      This is an exact casting of J7.  It does NOT bypass 451.2: it
      repackages those five statements as "the structure of H".  The
      only way to bypass them would be to need just a DEGREE bound on
      H rather than its structure (Turan + zero-sum).  Part (C)
      computes whether that is possible with the archive's own
      parameters.  It is not.

WHAT PROVABLY DOES NOT CAST.  "length == 2^k exactly" is a lacunary
support condition.  Roots-of-unity filters exist for arithmetic
progressions, not for {2^k}; the natural detector SUM_k u^{2^k} has
the unit circle as a natural boundary (Hadamard gap theorem), so it
is not an algebraic operator on a rational zeta.  This is 447.3's
scale-invariance no-go in analytic clothing: {2^k} is fixed by
doubling, and doubling is not a group operation on Z that Fourier
analysis sees.
"""
import sys
from itertools import combinations
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


def k4():
    return [[1, 2, 3], [0, 2, 3], [0, 1, 3], [0, 1, 2]]


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


def matmul(X, Y):
    Yt = list(zip(*Y))
    return [[sum(a * b for a, b in zip(r, c)) for c in Yt] for r in X]


def traces(adj, M):
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


def mobius(k):
    r = 1; d = 2; kk = k
    while d * d <= kk:
        if kk % d == 0:
            kk //= d
            if kk % d == 0: return 0
            r = -r
        d += 1
    if kk > 1: r = -r
    return r


GRAPHS = [('K4', k4()), ('Petersen', gp(5, 2)),
          ('Heawood', lcf(14, [5, -5] * 7)),
          ('Moebius-Kantor', lcf(16, [5, -5] * 8)),
          ('Pappus', lcf(18, [5, 7, -7, 7, -7, -5] * 3)),
          ('Desargues', lcf(20, [5, -5, 9, -9] * 5)),
          ('McGee', lcf(24, [12, 7, -7] * 8))]

print('=== (A) GF1: PRIME COUNTING FROM THE ZETA ===', flush=True)
print('  N_m = SUM_{l|m} l P_l ;  Moebius inversion ;  P_l = 2 c_l'
      ' below 2g')
print('  graph            g   l   P_l    2*c_l   equal?')
ok = True
for name, adj in GRAPHS:
    g = girth(adj); M = min(2 * g + 1, 13)
    N = traces(adj, M); cc = cycle_counts(adj, M)
    for l in range(3, M + 1):
        s = 0
        for d in range(1, l + 1):
            if l % d == 0:
                s += mobius(d) * N[l // d]
        assert s % l == 0, 'Moebius inversion not integral'
        P = s // l
        eq = (P == 2 * cc[l])
        if l < 2 * g and not eq: ok = False
        print('  %-16s %2d %3d %6d %8d   %s%s'
              % (name, g, l, P, 2 * cc[l], 'yes' if eq else 'NO',
                 '' if l < 2 * g else '   (l >= 2g, gap expected)'),
              flush=True)
print('  below 2g every prime is a simple cycle: %s' % ok, flush=True)

print()
print('=== (B) GF2: ROOTS-OF-UNITY FILTER + AX\'S THEOREM ===',
      flush=True)
print('  Lucas: bit j of |Z| = e_{2^j}(x) over F2, so "2^k divides')
print('  |Z|" is k polynomial conditions of degrees 1,2,...,2^{k-1}.')
print('  Ax: the number of solutions is divisible by 2^b with')
print('  b = ceil((mu - (2^k - 1)) / 2^{k-1}).  Since Z = 0 is a')
print('  solution, b >= 1 forces a NONEMPTY one.')
print()
print('  graph            n   mu   2^k   #even Z with 2^k | |Z|'
      '   Ax says 2^b |  b   nonempty forced?')
for name, adj in GRAPHS:
    n = len(adj)
    edges = sorted({(min(u, v), max(u, v)) for u in range(n)
                    for v in adj[u]})
    eidx = {e: i for i, e in enumerate(edges)}
    m = len(edges); mu = m - n + 1
    # spanning tree -> fundamental cycles -> cycle space basis
    par = {0: None}; order = [0]; tree = set()
    st = [0]
    while st:
        v = st.pop()
        for w in adj[v]:
            if w not in par:
                par[w] = v; tree.add((min(v, w), max(v, w)))
                order.append(w); st.append(w)
    depth = {0: 0}
    for v in order[1:]:
        depth[v] = depth[par[v]] + 1
    basis = []
    for e in edges:
        if e in tree: continue
        u, v = e
        vec = 1 << eidx[e]
        a, b = u, v
        while a != b:
            if depth[a] < depth[b]: a, b = b, a
            vec ^= 1 << eidx[(min(a, par[a]), max(a, par[a]))]
            a = par[a]
        basis.append(vec)
    assert len(basis) == mu
    wts = [0] * (1 << mu)
    for i in range(1, 1 << mu):
        low = i & -i
        j = low.bit_length() - 1
        wts[i] = wts[i ^ low] ^ basis[j]
    from collections import Counter
    wc = Counter(bin(w).count('1') for w in wts)
    for k in range(2, 7):
        P = 1 << k
        if P > mu + 8: break
        cnt = sum(v for w, v in wc.items() if w % P == 0)
        b = -(-(mu - (P - 1)) // (P >> 1)) if mu >= P - 1 else 0
        b = max(b, 0)
        forced = (mu >= P)
        print('  %-16s %3d %4d %5d %22d %14s %3d   %s'
              % (name, n, mu, P, cnt,
                 'yes' if (b == 0 or cnt % (1 << b) == 0) else 'NO',
                 b, 'YES' if forced and cnt >= 2 else
                 ('vacuous' if not forced else 'FAIL')), flush=True)

print()
print('=== (C) GF3: CAN A DEGREE BOUND REPLACE 451.2? ===',
      flush=True)
print('  Archive parameters for the J7 configuration:')
print('    - C has L vertices, each with exactly one third edge, so')
print('      the ear system has exactly N = L/2 ears (449).')
print('    - every ear interior has s >= L/2 (450), and by the')
print('      packing result n > L^2/4, so the total interior is')
print('      n - L ~ L^2/4 vertices.')
print('    - each interior vertex has one spare edge = one cross')
print('      edge, and each cross edge is a conflict between the two')
print('      ears it joins.')
print()
print('    L    ears N   interior ~   cross edges   pairs C(N,2)'
      '   cross/pairs   Turan alpha >= N/(d+1)')
for L in (24, 48, 96, 192, 384, 4096):
    N = L // 2
    inter = L * L // 4
    cross = inter // 2
    pairs = N * (N - 1) // 2
    d = 2 * cross / N
    print('  %5d %8d %12d %13d %14d %13.3f %19.3f'
          % (L, N, inter, cross, pairs, cross / pairs, N / (d + 1)),
          flush=True)
print()
print('  => cross edges ~ PAIRS OF EARS, to within a few percent.')
print('     The conflict multigraph on N = L/2 ears carries about')
print('     one conflict per pair: it sits exactly at the density')
print('     threshold.  Consequences:')
print('       - Turan sees only the average degree d ~ N, so it')
print('         gives alpha >= N/(N+1) < 1: NOTHING.')
print('       - whether a large compatible family exists therefore')
print('         depends on how the cross edges are DISTRIBUTED, not')
print('         on how many there are.')
print('     That distribution question is statement (B) of 451.2.')
print('     So the degree-bound bypass does not escape the five')
print('     statements - it reduces to one of them.')
