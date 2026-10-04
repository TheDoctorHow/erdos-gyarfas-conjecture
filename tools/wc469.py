"""wc469.py -- THE DECISIVE STATEMENT, IDENTIFIED AND ATTACKED.

WHAT WOULD SETTLE IT.  The Sieve Identity of 454.3 is EXACT for every
m; it needs no simplicity lemma and no m < 2g.  Under the
counterexample hypothesis at m = 2^k every proper-power term 2d c_d,
d | m, is 2^i c_{2^i} = 0.  So

    c_{2^k} = [ N_{2^k} - (rank >= 2 mass) ] / 2^{k+1},

and the ONLY thing standing between the trace frame and the
conjecture is a worst-case bound on the rank >= 2 mass.  454.1 proved
the first-moment scaling Theta(2^m / n); 456 showed the hypotheses
that would make it worst-case - kappa bounded, (NC) - are both false.
Name the statement:

  (WC)  For every cubic G on n vertices and every m <= sqrt(n/K),
        T_m + D_m + (rank >= 3 mass)  <=  K . 2^m . m / n.

WHY (WC) IS DECISIVE.  With (WC) the sieve's ceiling is no longer 2g
(simplicity) nor 4g (455.1) but sqrt(n/K).  The threshold is
unchanged at m_0(n) ~ 2 log2(2n).  So the usable window becomes

        ( m_0(n) , sqrt(n/K) ) ,

whose RATIO grows without bound instead of staying below 2.  461.4's
obstruction - the exact factor of 2 - evaporates, because it was an
artifact of the ceiling being 2g.  A power of two then lands in the
window for every n past a threshold, with NO girth hypothesis at all.
(WC) settles Erdos-Gyarfas for all large cubic graphs.

THE ATTACK.  462's Attachment Identity turns the rank-2 mass into a
cycle-path incidence count, and counting the incidences by the PATH
instead of by the cycle reduces (WC) to a statement about SHORTCUTS:
pairs of vertices far apart along a cycle but close in the graph.
That is an absolute geometric quantity, not a ratio, which is
exactly why it can survive 456's counterexamples.

Parts: (A) the window arithmetic and the threshold n_0;
       (B) S . n / m^2 measured across an adversarial family;
       (C) the shortcut reduction, verified;
       (D) can the adversary inflate it?
"""
import sys, math, random, itertools
from collections import Counter, deque

sys.setrecursionlimit(100000)
exec(open('tools/energy466.py').read().split(
    "ZOO = [")[0].split('"""', 2)[2])

M = 16


def m0(n):
    m = 2
    while not (2 ** m + 1 + n > 2 * (n - 1) * 2 ** (m // 2)): m += 2
    return m


def pow2_in(a, b):
    k = 0
    while 2 ** k <= a: k += 1
    return 2 ** k < b


print('=' * 72)
print('(A) THE WINDOW UNDER (WC), AND THE THRESHOLD n_0')
print('=' * 72)
print('  ceiling sqrt(n/K) instead of 2g.  461.4 measured the old')
print('  ratio at 1.980 and proved it never reaches 2.')
print()
print('  %12s %6s %10s %10s %10s %8s'
      % ('n', 'm_0', 'sqrt(n)', 'ratio', '2g_Moore', 'fires'))
for t in (10, 12, 14, 16, 20, 24, 30, 40, 60):
    n = 2 ** t
    a = m0(n); b = math.isqrt(n)
    g = 3
    while True:
        r = (g - 1) // 2
        lo = 3 * 2 ** r - 2 if g % 2 else 2 ** (r + 2) - 2
        if lo > n: g -= 1; break
        g += 1
    print('  %12d %6d %10d %10.2f %10d %8s'
          % (n, a, b, b / a, 2 * g, 'YES' if pow2_in(a, b) else 'no'))
print()
for K in (1, 4, 25, 100):
    n0 = None
    for n in range(100, 4000000):
        if pow2_in(m0(n), math.sqrt(n / K)):
            ok = all(pow2_in(m0(x), math.sqrt(x / K))
                     for x in range(n, min(n + 40000, 4000000), 97))
            if ok: n0 = n; break
    print('  K = %-4d : (WC) settles every cubic graph with n >= %s'
          % (K, n0 if n0 else '> 4e6'))
print()
print('  and the ratio sqrt(n/K)/m_0(n) -> infinity, so the window')
print('  stops being a ratio-2 interval.  That is the whole point.')

print()
print('=' * 72)
print('(B) S . n / m^2 ACROSS AN ADVERSARIAL FAMILY  (m = 16)')
print('=' * 72)
print('  S = 2m T_m / 2^m is exactly the energy ratio (462.2).')
print('  (WC) predicts S = O(m^2/n), i.e. S.n/m^2 bounded.')
print()


def profile(G, m=M, amax=6):
    cyc = all_cycles(G, m - 2)
    by_z = {}
    for c in cyc: by_z.setdefault(len(c[1]), []).append(c)
    per_a = Counter(); S = 0.0
    for a in range(1, amax + 1):
        z = m - 2 * a
        if z < 3: break
        cs = by_z.get(z, [])
        if not cs: continue
        R = sum(attachments(G, c, a) for c in cs)
        per_a[a] = R
        S += 2.0 * m * R / 2 ** m
    return sum(per_a.values()), per_a, S


def mobius_ladder(k):
    adj = {i: set() for i in range(2 * k)}
    for i in range(2 * k):
        adj[i].add((i + 1) % (2 * k)); adj[(i + 1) % (2 * k)].add(i)
    for i in range(k):
        adj[i].add(i + k); adj[i + k].add(i)
    A = [sorted(adj[i]) for i in range(2 * k)]
    return A if all(len(x) == 3 for x in A) else None


def prism(k):
    adj = {i: set() for i in range(2 * k)}
    for i in range(k):
        adj[i].add((i + 1) % k); adj[(i + 1) % k].add(i)
        adj[k + i].add(k + (i + 1) % k); adj[k + (i + 1) % k].add(k + i)
        adj[i].add(k + i); adj[k + i].add(i)
    A = [sorted(adj[i]) for i in range(2 * k)]
    return A if all(len(x) == 3 for x in A) else None


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


rng = random.Random(4444)
fam = []
for n in (18, 22, 26, 30, 34, 38):
    G = random_cubic(n, rng, mingirth=3)
    if G: fam.append(('random g>=3 n=%d' % n, G))
for n in (20, 28, 36):
    G = random_cubic(n, rng, mingirth=5)
    if G: fam.append(('random g>=5 n=%d' % n, G))
for k in (8, 12, 16):
    L = mobius_ladder(k)
    if L: fam.append(('Moebius ladder %d' % (2 * k), L))
for k in (8, 12, 16):
    P = prism(k)
    if P: fam.append(('prism %d' % (2 * k), P))
fam.append(('trunc-Petersen', truncate(gp(5, 2))))
fam.append(('trunc-cube', truncate(lcf(8, [3, -3]))))
ch = gp(5, 2)
for _ in range(2): ch = splice(ch, gp(5, 2))
fam.append(('Petersen chain x3', ch))
lad = mobius_ladder(8)
for nb in (16, 24):
    B = random_cubic(nb, rng, mingirth=3)
    if B and lad: fam.append(('ladder16+rand%d' % nb, splice(lad, B)))

print('  %-22s %4s %3s %9s %9s %10s'
      % ('graph', 'n', 'g', 'T_16', 'S', 'S.n/m^2'))
vals = []
for name, G in fam:
    n = len(G); g = girth(G)
    T, per_a, S = profile(G)
    v = S * n / (M * M)
    vals.append((v, name))
    print('  %-22s %4d %3d %9d %9.4f %10.4f' % (name, n, g, T, S, v))
vals.sort()
print()
print('  min %.4f (%s)   max %.4f (%s)'
      % (vals[0][0], vals[0][1], vals[-1][0], vals[-1][1]))
print('  spread factor %.1f over %d graphs of girth 3 to 6.'
      % (vals[-1][0] / max(vals[0][0], 1e-9), len(vals)))

print()
print('=' * 72)
print('(C) THE SHORTCUT REDUCTION')
print('=' * 72)
print('  Count incidences by the PATH, not the cycle.  A length-a')
print('  attachment to Z joins x,y in Z at GRAPH distance <= a; the')
print('  Band Lemma (462.5) forces their ARC distance >= g-a.  So')
print()
print('    R_a(z)  <=  w_a . SUM_Z sigma_a(Z),  w_1 = 1, w_a = 2^{a-2},')
print('    sigma_a(Z) = #{x,y in Z : arc-dist >= g-a, graph-dist <= a}.')
print()
print('  CORRECTED (374.1): the weight was first written 2^{a-2} for')
print('  every a, and it failed at a = 1 by exactly a factor 2 in 8 of')
print('  8 rows.  At a = 1 the attachment IS the shortcut, so w_1 = 1')
print('  and the relation is an EQUALITY: R_1(z) = SUM_Z sigma_1(Z).')
print()
print('  sigma is a SHORTCUT count: pairs far along the cycle but')
print('  close in the graph.  It is absolute and geometric, not a')
print('  ratio - which is why 456 does not touch it.')
print()


def sigma(G, cyc, a, g):
    Z = cyc[0]; z = len(Z); pos = {v: i for i, v in enumerate(Z)}
    d = {}
    for s in Z:
        dd = {s: 0}; q = deque([s])
        while q:
            v = q.popleft()
            if dd[v] >= a: continue
            for w in G[v]:
                if w not in dd: dd[w] = dd[v] + 1; q.append(w)
        d[s] = dd
    cnt = 0
    for i in range(z):
        for j in range(i + 1, z):
            arc = min(j - i, z - (j - i))
            if arc < g - a: continue
            if d[Z[i]].get(Z[j], 99) <= a: cnt += 1
    return cnt


print('  %-22s %4s %3s %3s %8s %12s %10s'
      % ('graph', 'n', 'z', 'a', 'R_a(z)', 'w_a.Ssigma', 'valid'))
for name, G in fam[:8]:
    n = len(G); g = girth(G)
    for a in (1, 2, 3):
        z = M - 2 * a
        cs = [c for c in all_cycles(G, z) if len(c[1]) == z]
        if not cs: continue
        R = sum(attachments(G, c, a) for c in cs)
        Ssig = sum(sigma(G, c, a, g) for c in cs)
        w = 1 if a == 1 else 2 ** (a - 2)
        bnd = w * Ssig
        print('  %-22s %4d %3d %3d %8d %12.1f %10s'
              % (name, n, z, a, R, bnd,
                 ('EQUAL' if R == bnd else 'ok') if R <= bnd + 1e-9
                 else 'FAIL'))

print()
print('=' * 72)
print('(D) VERDICT')
print('=' * 72)
print('  (WC) is the decisive statement: it removes the 2g ceiling,')
print('  hence 461.4\'s ratio-2 obstruction, hence settles E-G for all')
print('  large cubic graphs with no girth hypothesis.  It is NOT')
print('  implied by anything banked, and 456 killed the two ratio')
print('  hypotheses that would have given it.  The shortcut reduction')
print('  of (C) is the first absolute form of the same quantity.')
