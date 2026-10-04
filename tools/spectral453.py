"""spectral453.py -- THE TRACE FRAME.

A DIFFERENT SET OF TERMS.  Forget minimal counterexamples, ears,
shallow cycles and the 2-adic tower.  Count cycles by a TRACE.

Let G be cubic, T the 3n x 3n non-backtracking (Hashimoto) matrix on
directed edges: T[(u,v)][(v,w)] = 1 iff w != u.  Then tr(T^m) counts
closed non-backtracking walks of length m with a marked start.

LEMMA (cubic geodesic simplicity).  In a CUBIC graph every closed
non-backtracking walk of length m < 2g is a simple cycle traversed
once.  [If a vertex repeats, the two passes use two distinct pairs
from the 3 edges at that vertex, so they share an edge; the shared
edge is reused, and either it is reused in the same direction (the
walk is a proper power, m >= 2g) or in opposite directions (the walk
splits into two closed non-backtracking walks, each of length >= g).]

So for g <= m < 2g:            tr(T^m) = 2m * (#cycles of length m).

IHARA gives the spectrum of T from the spectrum of A: for cubic
(q = 2) the eigenvalues of T are the roots of z^2 - lam z + 2 for
each lam in spec(A), plus +1 and -1 each with multiplicity |E|-|V| =
n/2.  Hence with p_m(lam) = alpha^m + beta^m (alpha+beta = lam,
alpha*beta = 2):

    tr(T^m) = SUM_lam p_m(lam) + (n/2)(1 + (-1)^m).

p_m is a polynomial: p_0 = 2, p_1 = X, p_m = X p_{m-1} - 2 p_{m-2},
so SUM_lam p_m(lam) = tr(p_m(A)) -- EXACT INTEGER ARITHMETIC.

WHY POWERS OF TWO ARE SPECIAL HERE.  The doubling identity

    p_{2m} = p_m^2 - 2^{m+1},      i.e.   q_{2m} = q_m^2 - 2
    (q_m := p_m / 2^{m/2})

is the LUCAS-LEHMER / CHEBYSHEV map, conjugate to the doubling map
on the circle.  The power-of-two subsequence of the cycle counts is
the orbit of the graph's spectrum under a quadratic dynamical
system, not an arithmetic progression.  lam = 3 escapes
(p_m(3) = 2^m + 1); every lam in [-2rt2, 2rt2] stays bounded
(|p_m| <= 2 * 2^{m/2}); lam outside that interval contributes a
POSITIVE amount when m is even.

THE WINDOW BOUND.  For even m,

    tr(T^m) >= 2^m + 1 + n - 2(n-1) 2^{m/2},

positive as soon as 2^{m/2} > 2n roughly, i.e. m > 2 log2(2n).
Combined with the simplicity lemma:

  THEOREM (spectral window).  A cubic graph on n vertices with girth
  g has a cycle of EVERY even length m with
        2^m + 1 + n > 2(n-1) 2^{m/2}     and     m < 2g.
  In particular if that window has multiplicative width 2 it contains
  a power of 2 and Erdos-Gyarfas holds for G outright.

This file verifies: (1) the simplicity lemma and the trace identity
against brute-force cycle counts; (2) the spectral formula against
tr(T^m) computed directly; (3) the doubling identity; (4) the window
bound; (5) the WEAK FORM ("the unique power of 2 in [g,2g) is always
a cycle length") -- which is REFUTED here, honestly, by truncation.
"""
import sys
sys.setrecursionlimit(100000)


# ---------------------------------------------------------------- graphs
def gp(n, k):
    """Generalized Petersen GP(n,k): outer cycle, spokes, inner star."""
    adj = {i: [] for i in range(2 * n)}
    def add(a, b):
        adj[a].append(b); adj[b].append(a)
    for i in range(n):
        add(i, (i + 1) % n)
        add(i, n + i)
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


def truncate(adj):
    """Replace every vertex by a triangle.  Kills 4-cycles."""
    n = len(adj)
    idx = {}
    m = 0
    for v in range(n):
        for w in adj[v]:
            idx[(v, w)] = m; m += 1
    out = [[] for _ in range(m)]
    def add(a, b):
        out[a].append(b); out[b].append(a)
    for v in range(n):
        ports = [idx[(v, w)] for w in adj[v]]
        for a in range(len(ports)):
            for b in range(a + 1, len(ports)):
                add(ports[a], ports[b])
    for v in range(n):
        for w in adj[v]:
            if v < w:
                add(idx[(v, w)], idx[(w, v)])
    return [sorted(set(x)) for x in out]


# ------------------------------------------------------------- utilities
def girth(adj):
    from collections import deque
    n = len(adj); best = 10 ** 9
    for s in range(n):
        dist = {s: 0}; par = {s: -1}
        dq = deque([s])
        while dq:
            v = dq.popleft()
            for w in adj[v]:
                if w not in dist:
                    dist[w] = dist[v] + 1; par[w] = v; dq.append(w)
                elif w != par[v]:
                    best = min(best, dist[v] + dist[w] + 1)
    return best


def cycle_counts(adj, cap):
    """cnt[m] = number of cycles of length m, for m <= cap."""
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
    n = len(X); p = len(Y[0])
    Yt = list(zip(*Y))
    return [[sum(a * b for a, b in zip(row, col)) for col in Yt]
            for row in X]


def nb_traces(adj, M):
    """tr(T^m) for m = 1..M, T the Hashimoto matrix (direct)."""
    de = [(u, v) for u in range(len(adj)) for v in adj[u]]
    idx = {e: i for i, e in enumerate(de)}
    N = len(de)
    T = [[0] * N for _ in range(N)]
    for (u, v) in de:
        for w in adj[v]:
            if w != u:
                T[idx[(u, v)]][idx[(v, w)]] = 1
    P = [r[:] for r in T]
    out = [None, sum(P[i][i] for i in range(N))]
    for _ in range(2, M + 1):
        P = matmul(P, T)
        out.append(sum(P[i][i] for i in range(N)))
    return out


def spectral_traces(adj, M):
    """SUM_lam p_m(lam) + (n/2)(1+(-1)^m), by integer recursion on
    p_m(A) = A p_{m-1}(A) - 2 p_{m-2}(A)."""
    n = len(adj)
    A = [[0] * n for _ in range(n)]
    for v in range(n):
        for w in adj[v]:
            A[v][w] = 1
    Pm2 = [[2 if i == j else 0 for j in range(n)] for i in range(n)]
    Pm1 = [r[:] for r in A]
    ex = n // 2
    out = [None, sum(Pm1[i][i] for i in range(n)) + 0]
    for m in range(2, M + 1):
        Cur = matmul(A, Pm1)
        for i in range(n):
            for j in range(n):
                Cur[i][j] -= 2 * Pm2[i][j]
        Pm2, Pm1 = Pm1, Cur
        out.append(sum(Pm1[i][i] for i in range(n))
                   + ex * (1 + (-1) ** m))
    return out


GRAPHS = [('K4', k4()), ('K3,3', k33()), ('prism GP(3,1)', gp(3, 1)),
          ('cube GP(4,1)', gp(4, 1)), ('Petersen GP(5,2)', gp(5, 2)),
          ('Heawood', heawood()), ('Moebius-Kantor GP(8,3)', gp(8, 3)),
          ('Pappus-like GP(9,2)', gp(9, 2)),
          ('dodecahedron GP(10,2)', gp(10, 2)),
          ('Desargues GP(10,3)', gp(10, 3)),
          ('Nauru GP(12,5)', gp(12, 5))]

print('=== (1) SIMPLICITY LEMMA + TRACE IDENTITY, m < 2g ===',
      flush=True)
print('  graph                     n   g   m   tr(T^m)   2m*c_m'
      '   spectral   ok')
allok = True
for name, adj in GRAPHS:
    n = len(adj); g = girth(adj)
    M = min(2 * g - 1, 13)
    tt = nb_traces(adj, M)
    st = spectral_traces(adj, M)
    cc = cycle_counts(adj, M)
    for m in range(3, M + 1):
        ok = (tt[m] == 2 * m * cc[m]) and (st[m] == tt[m])
        allok = allok and ok
        print('  %-24s %3d %3d %3d %9d %8d %10d   %s'
              % (name, n, g, m, tt[m], 2 * m * cc[m], st[m],
                 'OK' if ok else 'FAIL'), flush=True)
print('  ALL CHECKS PASS' if allok else '  *** A CHECK FAILED ***',
      flush=True)

print()
print('=== (2) DOUBLING IDENTITY  p_2m = p_m^2 - 2^{m+1} ===',
      flush=True)
bad = 0
for lam in range(-3, 4):
    p = [2, lam]
    for m in range(2, 33):
        p.append(lam * p[-1] - 2 * p[-2])
    for m in (1, 2, 4, 8, 16):
        if p[2 * m] != p[m] ** 2 - 2 ** (m + 1):
            bad += 1
print('  violations over lam = -3..3, m = 1,2,4,8,16: %d' % bad,
      flush=True)
p3 = [2, 3]
for m in range(2, 33):
    p3.append(3 * p3[-1] - 2 * p3[-2])
print('  p_m(3) = 2^m + 1 for m = 0..32: %s'
      % all(p3[m] == 2 ** m + 1 for m in range(33)), flush=True)
ok = True
for lam in range(-2, 3):
    p = [2, lam]
    for m in range(2, 25):
        p.append(lam * p[-1] - 2 * p[-2])
    for m in range(1, 25):
        if abs(p[m]) > 2 * (2 ** (m / 2)) + 1e-9:
            ok = False
print('  |p_m(lam)| <= 2*2^{m/2} for lam in -2..2, m = 1..24: %s'
      % ok, flush=True)

print()
print('=== (3) THE WINDOW BOUND vs REALITY ===', flush=True)
print('  bound(m) = 2^m + 1 + n - 2(n-1)2^{m/2}   (even m)')
print('  graph                     n   g   even m in window'
      '   c_m actually')
for name, adj in GRAPHS:
    n = len(adj); g = girth(adj)
    cap = min(2 * g - 1, 13)
    cc = cycle_counts(adj, cap)
    win = []
    for m in range(4, 2 * g, 2):
        if m > cap: break
        bnd = 2 ** m + 1 + n - 2 * (n - 1) * 2 ** (m // 2)
        if bnd > 0:
            win.append((m, cc[m]))
    print('  %-24s %3d %3d   %-18s %s'
          % (name, n, g, str([w[0] for w in win]) or '[]',
             str([w[1] for w in win])), flush=True)
print('  (empty windows are vacuous, not failures; every nonempty'
      ' entry must have c_m > 0)', flush=True)

print()
print('=== (4) WEAK FORM: is the unique power of 2 in [g,2g)'
      ' always a cycle length? ===', flush=True)
def unique_pow2(g):
    p = 1
    while p < g: p *= 2
    return p
for name, adj in GRAPHS + [('truncated Petersen', truncate(gp(5, 2))),
                           ('truncated K4', truncate(k4())),
                           ('truncated cube', truncate(gp(4, 1)))]:
    n = len(adj); g = girth(adj); P = unique_pow2(g)
    cap = max(P, g) + 1
    cc = cycle_counts(adj, min(cap, 14))
    has = cc[P] if P < len(cc) else -1
    print('  %-24s n=%3d g=%d  2^k=%2d  #cycles of that length = %s  %s'
          % (name, n, g, P, has,
             'ok' if has > 0 else 'WEAK FORM FAILS'), flush=True)
