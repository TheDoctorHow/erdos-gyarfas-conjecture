"""verify455.py -- THE SPECTRAL WINDOW CRITERION, tested where it bites.

THEOREM (spectral window).  Let G be cubic on n vertices with girth
g.  For every even m with

    (W)   2^m + 1 + n > 2(n-1) 2^{m/2}      and      m < 2g,

G has a cycle of length m; indeed
    c_m >= [ 2^m + 1 + n - 2(n-1)2^{m/2} ] / 2m.

COROLLARY.  If some power of 2 satisfies (W), Erdos-Gyarfas holds for
G outright -- no counterexample assumed, no casework, decidable from
(n, g) alone.

(W) says m > ~2 log2(2n); simplicity says m < 2g; Moore says
g <= ~2 log2 n.  So the window is an interval of multiplicative width
g/(log2 n + 1), which reaches 2 only at the Moore bound.  This file
exhibits graphs where a power of 2 lands inside.

ERROR RECORDED (374.1).  The first version of this file asserted that
the Tutte-Coxeter graph has NO 10-cycle, and billed that as showing
the threshold is sharp.  FALSE: the census below finds 72 ten-cycles,
by brute force and by the trace formula independently.  (W) is
CONSERVATIVE, not tight -- it says where cycles are FORCED, not where
they start.  On Tutte-Coxeter cycles already exist at m = 8 and 10,
while (W) first fires at m = 11.
"""
import sys
sys.setrecursionlimit(100000)


def lcf(n, seq):
    """LCF notation: Hamiltonian cycle 0..n-1 plus chord i -> i+seq[i]."""
    adj = {i: set() for i in range(n)}
    for i in range(n):
        adj[i].add((i + 1) % n); adj[(i + 1) % n].add(i)
    for i in range(n):
        j = (i + seq[i % len(seq)]) % n
        adj[i].add(j); adj[j].add(i)
    return [sorted(adj[i]) for i in range(n)]


def rep(pat, k):
    return pat * k


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


def bound(n, m):
    return 2 ** m + 1 + n - 2 * (n - 1) * 2 ** (m // 2)


# ------------------------------------------------------------------ zoo
ZOO = [
    ('Heawood',        lcf(14, rep([5, -5], 7))),
    ('Moebius-Kantor', lcf(16, rep([5, -5], 8))),
    ('Pappus',         lcf(18, rep([5, 7, -7, 7, -7, -5], 3))),
    ('Desargues',      lcf(20, rep([5, -5, 9, -9], 5))),
    ('McGee',          lcf(24, rep([12, 7, -7], 8))),
    ('Tutte-Coxeter',  lcf(30, rep([-13, -9, 7, -7, 9, 13], 5))),
    ('Foster',         lcf(90, rep([17, -9, 37, -37, 9, -17], 15))),
]

print('=== (0) THE ZOO (LCF constructions, girth verified) ===',
      flush=True)
for name, adj in ZOO:
    degs = sorted(set(len(a) for a in adj))
    print('  %-16s n=%3d  degrees=%s  girth=%d'
          % (name, len(adj), degs, girth(adj)), flush=True)

print()
print('=== (1) DOES A POWER OF 2 LAND IN THE WINDOW? ===', flush=True)
print('  graph              n   g   2g   window (even m, W holds)'
      '        powers of 2 inside')
for name, adj in ZOO:
    n = len(adj); g = girth(adj)
    win = [m for m in range(4, 2 * g, 2) if bound(n, m) > 0]
    pw = [m for m in win if (m & (m - 1)) == 0]
    print('  %-16s %3d %3d %4d   %-28s %s'
          % (name, n, g, 2 * g, str(win), str(pw) if pw else 'NONE'),
          flush=True)

print()
print('=== (2) THE FORCED COUNTS, CHECKED AGAINST THE TRUTH ===',
      flush=True)
print('  graph              m   forced c_m >=   actual c_m   verdict')
for name, adj in ZOO:
    n = len(adj); g = girth(adj)
    win = [m for m in range(4, 2 * g, 2) if bound(n, m) > 0]
    if not win: continue
    M = max(win)
    st = spectral_traces(adj, M)
    for m in win:
        forced = -(-bound(n, m) // (2 * m))
        actual = st[m] // (2 * m)
        assert st[m] % (2 * m) == 0, 'trace not divisible by 2m'
        print('  %-16s %3d %13d %12d   %s'
              % (name, m, forced, actual,
                 'OK' if actual >= forced > 0 else 'FAIL'), flush=True)

print()
print('=== (3) SHARPNESS: Tutte-Coxeter at the threshold ===',
      flush=True)
name, adj = ZOO[5]
n = len(adj); g = girth(adj)
st = spectral_traces(adj, 15)
cc = cycle_counts(adj, 15)
print('  n=%d  g=%d  2g=%d' % (n, g, 2 * g))
print('   m   bound(n,m)   spectral c_m   brute-force c_m   (W)?')
for m in range(6, 16):
    sc = st[m] // (2 * m) if st[m] % (2 * m) == 0 else None
    print('  %2d %12d %14s %17d   %s'
          % (m, bound(n, m), str(sc), cc[m],
             'YES' if (bound(n, m) > 0 and m < 2 * g) else 'no'),
          flush=True)
print('  => (W) first fires at m = 11, but cycles already exist at')
print('     m = 8 (90 of them) and m = 10 (72).  So (W) is a')
print('     CONSERVATIVE forcing threshold, not the true onset.  The')
print('     gap between onset and forcing is the slack the sieve')
print('     would have to recover.', flush=True)

print()
print('=== (4) FOSTER: a forced 16-CYCLE, i.e. Erdos-Gyarfas for it'
      ' by pure spectrum ===', flush=True)
name, adj = ZOO[6]
n = len(adj); g = girth(adj)
st = spectral_traces(adj, 16)
print('  n=%d  g=%d  2g=%d   16 < 2g: %s' % (n, g, 2 * g, 16 < 2 * g))
print('  bound(90,16) = %d  ->  c_16 >= %d'
      % (bound(n, 16), -(-bound(n, 16) // 32)))
print('  tr(T^16) = %d,  divisible by 2m=32: %s,  c_16 = %d'
      % (st[16], st[16] % 32 == 0, st[16] // 32), flush=True)
cc = cycle_counts(adj, 16)
print('  brute-force c_16 = %d   %s'
      % (cc[16], 'MATCH' if cc[16] == st[16] // 32 else 'MISMATCH'),
      flush=True)
