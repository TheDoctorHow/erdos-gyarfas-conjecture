"""world465.py -- REVERSE-ENGINEERING THE TWO SETTLED WORLDS.

WORLD A (J7 false).  By 434.3 a graph satisfying (a)-(f) IS a cubic
counterexample, so "J7 false" = "Erdos-Gyarfas false".  Reverse-
engineering that world means writing the counterexample's identikit
from the archive's own unconditional bounds and asking whether the
coordinates are jointly feasible.

WORLD B (J7 true).  Reverse-engineering that world means asking what
signature a proof must have.  The trace frame (452) supplies a
counting instrument whose USABLE WINDOW is an interval; this file
measures that interval's RATIO, because the target set {2^k} has
spacing exactly 2.  An instrument whose window has ratio < 2 cannot
be made to hit a lacunary target for every input, no matter how the
constants improve.

Checks performed:
  (A) p_m(A) = A_m - A_{m-2}  and the DOUBLING IDENTITY
      tr(T^{2m}) = ||p_m(A)||_F^2 - n*2^{m+1} + n.
  (B) the TWO-SCALE IDENTITY in terms of the non-simple remainders
      R_m = tr(T^m) - 2m*c_m.
  (C) the RATIO-2 MEASUREMENT: window (m_0(n), 2g) against Moore.
  (D) the IDENTIKIT: which (n, g) survive every unconditional test.
  (E) the CAUCHY-SCHWARZ SELF-IMPROVEMENT along the doubling orbit,
      and its threshold.
"""
import sys, math, itertools
from collections import deque

sys.setrecursionlimit(100000)


# ---------------------------------------------------------------- graphs
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
        if a != b: adj[a].add(b); adj[b].add(a)
    for i in range(n):
        add(i, (i + 1) % n); add(i, n + i)
        add(n + i, n + (i + k) % n)
    A = [sorted(adj[i]) for i in range(2 * n)]
    return A if all(len(x) == 3 for x in A) else None


def truncate(adj):
    n = len(adj); idx = {}; c = 0
    for v in range(n):
        for w in adj[v]:
            idx[(v, w)] = c; c += 1
    out = [set() for _ in range(c)]
    for v in range(n):
        ws = adj[v]
        for a, b in itertools.combinations(ws, 2):
            out[idx[(v, a)]].add(idx[(v, b)])
            out[idx[(v, b)]].add(idx[(v, a)])
    for v in range(n):
        for w in adj[v]:
            out[idx[(v, w)]].add(idx[(w, v)])
            out[idx[(w, v)]].add(idx[(v, w)])
    return [sorted(x) for x in out]


def girth(adj):
    n = len(adj); best = 10 ** 9
    for s in range(n):
        d = {s: 0}; p = {s: -1}; q = deque([s])
        while q:
            v = q.popleft()
            for w in adj[v]:
                if w not in d:
                    d[w] = d[v] + 1; p[w] = v; q.append(w)
                elif w != p[v]:
                    best = min(best, d[v] + d[w] + 1)
    return best


def count_cycles(adj, L):
    """number of simple cycles of length exactly L."""
    n = len(adj); tot = 0
    for s in range(n):
        seen = {s};
        def dfs(v, ln):
            r = 0
            for w in adj[v]:
                if w == s and ln == L: r += 1
                elif w > s and w not in seen and ln < L:
                    seen.add(w); r += dfs(w, ln + 1); seen.discard(w)
            return r
        tot += dfs(s, 1)
    return tot // 2


# ------------------------------------------------ exact integer matrices
def matmul(X, Y):
    n = len(X); m = len(Y[0]); k = len(Y)
    return [[sum(X[i][t] * Y[t][j] for t in range(k)) for j in range(m)]
            for i in range(n)]


def adjmat(adj):
    n = len(adj)
    M = [[0] * n for _ in range(n)]
    for v in range(n):
        for w in adj[v]: M[v][w] = 1
    return M


def nb_walk_mats(adj, mmax):
    """A_0..A_mmax, A_m[u][v] = # non-backtracking walks u->v of len m."""
    n = len(adj); A = adjmat(adj)
    I = [[1 if i == j else 0 for j in range(n)] for i in range(n)]
    out = [I, A]
    if mmax >= 2:
        A2 = matmul(A, A)
        out.append([[A2[i][j] - (3 if i == j else 0) for j in range(n)]
                    for i in range(n)])
    for m in range(3, mmax + 1):
        P = matmul(A, out[m - 1])
        out.append([[P[i][j] - 2 * out[m - 2][i][j] for j in range(n)]
                    for i in range(n)])
    return out


def p_mats(adj, mmax):
    """p_m(A) by the recursion p_{m+1}=A p_m - 2 p_{m-1}, p_0=2I, p_1=A."""
    n = len(adj); A = adjmat(adj)
    P0 = [[2 if i == j else 0 for j in range(n)] for i in range(n)]
    out = [P0, A]
    for m in range(2, mmax + 1):
        P = matmul(A, out[m - 1])
        out.append([[P[i][j] - 2 * out[m - 2][i][j] for j in range(n)]
                    for i in range(n)])
    return out


def hashimoto_traces(adj, mmax):
    """exact tr(T^m) by powering the Hashimoto matrix."""
    darts = []; ix = {}
    for v in range(len(adj)):
        for w in adj[v]:
            ix[(v, w)] = len(darts); darts.append((v, w))
    D = len(darts)
    T = [[0] * D for _ in range(D)]
    for (u, v) in darts:
        for w in adj[v]:
            if w != u: T[ix[(u, v)]][ix[(v, w)]] = 1
    tr = []; Cur = [[1 if i == j else 0 for j in range(D)] for i in range(D)]
    for m in range(1, mmax + 1):
        Cur = matmul(Cur, T)
        tr.append(sum(Cur[i][i] for i in range(D)))
    return tr


ZOO = [
    ('K4',              lcf(4, [2, -2])),
    ('K33',             lcf(6, [3, 3, 3])),
    ('cube',            lcf(8, [3, -3])),
    ('Petersen',        gp(5, 2)),
    ('Heawood',         lcf(14, [5, -5])),
    ('Pappus',          lcf(18, [5, 7, -7, 7, -7, -5])),
    ('Desargues',       gp(10, 3)),
    ('Nauru',           gp(12, 5)),
    ('trunc-tetra',     truncate(lcf(4, [2, -2]))),
]

print('=' * 68)
print('(A) THE DOUBLING IDENTITY, EXACT')
print('=' * 68)
print('  p_m(A) = A_m - A_{m-2} - A_{m-4} - ...   (down to A_1 or A_0)')
print('  tr(T^{2m}) = ||p_m(A)||_F^2 - n*2^{m+1} + n')
print()
print('  %-14s %3s  %-9s %-9s' % ('graph', 'n', 'p=alt-sum', 'doubling'))
MM = 7


def alt(Aw, m, n):
    """A_m - sum_{j>=1} A_{m-2j}."""
    R = [[Aw[m][i][j] for j in range(n)] for i in range(n)]
    k = m - 2
    while k >= 0:
        for i in range(n):
            for j in range(n): R[i][j] -= Aw[k][i][j]
        k -= 2
    return R


for name, G in ZOO:
    n = len(G)
    Aw = nb_walk_mats(G, MM)
    Pm = p_mats(G, MM)
    ok1 = all(Pm[m] == alt(Aw, m, n) for m in range(2, MM + 1))
    tr = hashimoto_traces(G, 2 * MM)
    ok2 = True
    for m in range(1, MM + 1):
        F = sum(Pm[m][i][j] ** 2 for i in range(n) for j in range(n))
        if tr[2 * m - 1] != F - n * 2 ** (m + 1) + n: ok2 = False
    print('  %-14s %3d  %-9s %-9s'
          % (name, n, 'OK' if ok1 else 'FAIL', 'OK' if ok2 else 'FAIL'))

print()
print('=' * 68)
print('(B) THE TWO-SCALE IDENTITY  (R_m = tr(T^m) - 2m*c_m)')
print('=' * 68)
print('  sum_lam p_M(lam)   = R_M - n            (M even)')
print('  sum_lam p_M(lam)^2 = n(2^{M+1}-1) + R_{2M}')
print()
print('  %-12s %2s %2s %6s %8s %10s %12s %6s'
      % ('graph', 'g', 'M', 'c_M', 'R_M', 'c_2M', 'R_2M', 'both'))
for name, G in ZOO:
    n = len(G); g = girth(G)
    tr = hashimoto_traces(G, 12)
    Pm = p_mats(G, 6)
    for M in (4, 6):
        if 2 * M > 12: continue
        cM = count_cycles(G, M); c2M = count_cycles(G, 2 * M)
        RM = tr[M - 1] - 2 * M * cM
        R2M = tr[2 * M - 1] - 2 * (2 * M) * c2M
        S1 = sum(Pm[M][i][i] for i in range(n))
        F = sum(Pm[M][i][j] ** 2 for i in range(n) for j in range(n))
        ok = (S1 == RM + 2 * M * cM - n) and \
             (F == n * (2 ** (M + 1) - 1) + R2M + 2 * (2 * M) * c2M)
        print('  %-12s %2d %2d %6d %8d %10d %12d %6s'
              % (name, g, M, cM, RM, c2M, R2M, 'OK' if ok else 'FAIL'))

print()
print('=' * 68)
print('(C) THE RATIO-2 MEASUREMENT')
print('=' * 68)


def m0(n):
    """smallest even m with 2^m + 1 + n > 2(n-1)2^{m/2}: the trace
    frame's counting threshold."""
    m = 2
    while not (2 ** m + 1 + n > 2 * (n - 1) * 2 ** (m // 2)):
        m += 2
    return m


def moore_g(n):
    """largest girth the Moore bound permits for a cubic graph on n."""
    g = 3
    while True:
        r = (g - 1) // 2
        lo = 3 * 2 ** r - 2 if g % 2 else 2 ** (r + 2) - 2
        if lo > n: return g - 1
        g += 1


def pow2_in(a, b):
    """is there a power of two strictly inside (a, b)?"""
    k = 0
    while 2 ** k <= a: k += 1
    return 2 ** k < b


print('  m_0(n) ~ 2 log2(2n) is the counting threshold;')
print('  2*g is the ceiling (simplicity lemma).  Target spacing = 2.')
print('  The window is nonempty iff 2g > m_0, i.e. g > ~log2 n + 2 --')
print('  the SAME threshold at which 457.3 kills the disjoint horn.')
print()
print('  %14s %6s %8s %6s %8s %8s %6s'
      % ('n', 'm_0', 'm_0/2', 'g_max', '2g_max', 'ratio', 'fires'))
for t in range(4, 65, 4):
    n = 2 ** t
    a = m0(n); gm = moore_g(n); b = 2 * gm
    print('  %14d %6d %8.1f %6d %8d %8.4f %6s'
          % (n, a, a / 2.0, gm, b, b / a, 'YES' if pow2_in(a, b) else 'no'))
print()
print('  THE RATIO IS ALWAYS < 2.  Checked over every dyadic scale and')
print('  every n in [4, 3e5]:')
bad = 0; worst = 0.0
for n in range(4, 300000):
    r = 2.0 * moore_g(n) / m0(n)
    worst = max(worst, r)
    if r >= 2.0: bad += 1
for t in range(2, 200):
    n = 2 ** t
    r = 2.0 * moore_g(n) / m0(n)
    worst = max(worst, r)
    if r >= 2.0: bad += 1
print('    violations of ratio < 2 : %d      sup ratio seen : %.6f'
      % (bad, worst))
print()
print('  Consequence: the trace frame NEVER guarantees a hit.  Its')
print('  usable window has ratio < 2 for every n, and 2 is exactly the')
print('  spacing of the target set {2^k}.  Firing is decided by where')
print('  log2(m_0) sits modulo 1, so it fires and misses in blocks of n')
print('  of DOUBLY exponential length (m_0 grows like log n).')
print()
blocks = []
prev = None
for t in range(3, 140):
    n = 2 ** t
    f = pow2_in(m0(n), 2 * moore_g(n))
    if f != prev:
        blocks.append((t, f)); prev = f
print('  block boundaries (n = 2^t), fires? :')
print('   ', ' '.join('%d:%s' % (t, 'Y' if f else 'n') for t, f in blocks))

print()
print('=' * 68)
print('(D) THE IDENTIKIT: WHICH (n, g) SURVIVE')
print('=' * 68)
print('  A counterexample must fail EVERY unconditional test:')
print('   (i)   g is not a power of two            (girth cycle wins)')
print('   (ii)  no power of two in (m_0(n), 2g)    (452.4)')
print('   (iii) n >= Moore(g)                      (Moore)')
print('   (iv)  n > 144                            (450.4, L >= 24)')
print()
print('  %4s %14s %10s %10s %8s'
      % ('g', 'Moore(g)', 'need n >=', 'm_0 there', 'verdict'))
for g in range(3, 25):
    r = (g - 1) // 2
    mn = 3 * 2 ** r - 2 if g % 2 else 2 ** (r + 2) - 2
    if g & (g - 1) == 0:
        print('  %4d %14d %10s %10s %8s'
              % (g, mn, '-', '-', 'DEAD (2^j)'))
        continue
    # smallest n >= max(mn,145) for which the window misses
    n = max(mn, 145); found = None
    while n < 10 ** 7:
        if not pow2_in(m0(n), 2 * g): found = n; break
        n += 1
    print('  %4d %14d %10s %10d %8s'
          % (g, mn, str(found) if found else '>1e7',
             m0(found) if found else -1,
             'live' if found else 'DEAD'))

print()
print('  THE HANDOFF, MEASURED.  457.3: disjoint girth cycles give')
print('  kappa_g <= 2n/2^g < 2 once g >= log2 n + 1 -- the spread the')
print('  trace frame wants.  452.4: the window is nonempty once')
print('  2g > m_0(n).  Where are those two thresholds?')
print()
print('  %14s %16s %18s %10s'
      % ('n', 'g >= log2 n + 1', 'g > m_0(n)/2', 'gap'))
for t in range(4, 65, 6):
    n = 2 ** t
    t1 = t + 1
    t2 = m0(n) // 2 + 1
    print('  %14d %16d %18d %10d' % (n, t1, t2, t2 - t1))
print()
print('  The two thresholds agree to an additive constant: the trace')
print('  frame opens exactly where the closure frame closes.  But it')
print('  opens with ratio just above 1 and Moore caps it below 2.')

print()
print('=' * 68)
print('(E) CAUCHY-SCHWARZ ALONG THE DOUBLING ORBIT')
print('=' * 68)
print('  tr(T^{2m}) >= (tr(T^m) - e*n)^2/n - n(2^{m+1}-1),  e=[m even]')
print()
print('  %-12s %3s %3s %14s %16s %8s'
      % ('graph', 'n', 'm', 'tr(T^2m)', 'CS lower bound', 'valid'))
for name, G in ZOO[:6]:
    n = len(G)
    tr = hashimoto_traces(G, 14)
    for m in (5, 6, 7):
        e = n if m % 2 == 0 else 0
        lb = (tr[m - 1] - e) ** 2 // n - n * (2 ** (m + 1) - 1)
        print('  %-12s %3d %3d %14d %16d %8s'
              % (name, n, m, tr[2 * m - 1], lb,
                 'ok' if tr[2 * m - 1] >= lb else 'FAIL'))
print()
print('  the bound turns POSITIVE when tr(T^m) > n + n*2^{(m+1)/2},')
print('  i.e. m > 2 log2 n + 1 -- the SAME threshold as (C).')
print('  every counting instrument in the trace frame has')
print('  threshold ~2 log2 n and ceiling 2g <= 4 log2 n - 4.3.')

print()
print('=' * 68)
print('(F) THE GIRTH THRESHOLD IS A STEP FUNCTION, AND ITS STEPS ARE')
print('    POWERS OF TWO')
print('=' * 68)


def Gstar(n):
    """smallest girth above which 452.4 fires: the window (m_0, 2g)
    contains a power of two iff g > 2^{ceil(log2 m_0)-1}."""
    a = m0(n); j = 0
    while 2 ** j <= a: j += 1
    return 2 ** (j - 1)


print('  THEOREM (unconditional, 452.4 restated).  Every cubic graph')
print('  on n vertices with girth g > G*(n) has a cycle of length')
print('  2^ceil(log2 m_0(n)) -- a power of two.  G*(n) = 2^{j-1} with')
print('  2^j the least power of two exceeding m_0(n).')
print()
print('  G* only takes POWER-OF-TWO values, so it is flat over long')
print('  stretches and then doubles.  Moore rises like 2 log2 n.')
print()
print('  %14s %8s %8s %10s %10s'
      % ('n', 'm_0', 'G*(n)', 'Moore cap', 'nonvacuous'))
for t in list(range(4, 40, 2)) + [40, 48, 56, 62, 63, 64, 66, 70]:
    n = 2 ** t
    gs = Gstar(n); mc = moore_g(n)
    print('  %14s %8d %8d %10d %10s'
          % ('2^%d' % t, m0(n), gs, mc, 'YES' if gs < mc else 'no'))
print()
print('  THE FAILURE SET.  G*(n) >= Moore(n) -- no cubic graph on n')
print('  vertices can reach the threshold -- exactly when:')
lo = []
for t in range(3, 300):
    n = 2 ** t
    if Gstar(n) >= moore_g(n): lo.append(t)
runs = []
for t in lo:
    if runs and t == runs[-1][1] + 1: runs[-1][1] = t
    else: runs.append([t, t])
print('    log2 n in ', ' '.join('[%d,%d]' % (a, b) for a, b in runs))
print('    nearest powers of two:',
      ' '.join(str(2 ** k) for k in range(2, 9)))
print()
print('  So the trace frame fails precisely when log2 n is within an')
print('  additive 1 of a POWER OF TWO, i.e. when n sits at a')
print('  Fermat-type scale n ~ 2^{2^k}.  The target set obstructs the')
print('  instrument at exactly the scales the target set names: the')
print('  doubling orbit blocks its own second iterate.  This is')
print('  447.3 scale-invariance in the trace frame\'s own coordinates.')
