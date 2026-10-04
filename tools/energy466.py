"""energy466.py -- THE ENERGY DICHOTOMY, PURSUED.

454.3 left the trace frame resting on one inequality,
    T_m < 2^m / 2m,     m = 2^k,
with T_m the number of pairs of cycles that meet in a path and whose
lengths sum to m.  455 proved the theta half under kappa < 1.41; 456
showed kappa is unbounded and (NC) false.  This file changes the
representation of T_m and asks what the inequality then needs.

THE ATTACHMENT IDENTITY.  A theta with branches (a,b,c) carries the
cycles a+b, a+c, b+c and its length-m geodesic doubles one branch.
Doubling c gives m = (a+b) + 2c.  So the pair counted by T_m is
recovered from the OPPOSITE cycle Z = a+b together with the doubled
branch P = c, and P is a path with both ends on Z and interior off
it -- an ATTACHMENT.  Hence

    T_m  =  # { (Z, P) : Z a cycle, P an attachment, |Z| + 2|P| = m }
         =  SUM_{a>=1} R_a(m - 2a),   R_a(z) = # (z-cycle, a-attachment).

T_m is therefore an INCIDENCE count between cycles and paths, not a
convolution of the cycle spectrum.  That is the precise sense in
which the energy is geometric, and it is why BSG does not apply
directly.

WHAT THE IDENTITY BUYS.  In a cubic graph a cycle Z of length z has
exactly z legs, so r_a(Z) <= z 2^{a-2}, and with c_z <= N_z / 2z the
identity gives a term-by-term bound whose weight is GEOMETRIC in a.
Writing chi_a for the average number of a-attachments per cycle of
length m-2a and gamma_z = 2z c_z / 2^z for the cycle density,

    THE LAYER DECOMPOSITION (an IDENTITY, not a bound)
        S := 2m T_m / 2^m
           = SUM_{a>=1} chi_a . gamma_{m-2a} . 4^{-a} . m/(m-2a),
    so S < 1 is EXACTLY the energy inequality T_m < 2^m/2m, split
    into layers by attachment length with explicit weight 4^{-a}.

The a = 1 term is chords of an (m-2)-cycle, weight 1/4: the whole
inequality is dominated by the question "how many chords does the
average long cycle have?", and the answer must be fewer than four.
Every larger a is suppressed by 4^{-a}.

THE (NC) TEST.  (NC) failed in 456 because a spliced Petersen block
holds all the short cycles in ten vertices.  The Chord Criterion asks
about chords of LONG cycles, which is a different quantity.  Part (E)
builds exactly 456's counterexample and measures it.

CORRECTION (374.1), RECORDED LOUDLY.  Parts (D) and (E) below run at
m = 8, and at m = 8 the test is VACUOUS: for girth >= 5 every split
x + y = 8 needs a cycle of length <= 3, so T_8 = 0 identically and
S = 0 for every graph in those tables.  The rows discriminate
nothing.  They are left in place because the zeros are correct and
because the vacuity is the point: the first scale with content at
girth 5 is m = 16.  The real test is tools/nc468.py.
"""
import sys, random, itertools
from collections import Counter, deque

sys.setrecursionlimit(100000)


# ------------------------------------------------------------- graphs
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


def girth(adj):
    n = len(adj); best = 10 ** 9
    for s in range(n):
        d = {s: 0}; p = {s: -1}; q = deque([s])
        while q:
            v = q.popleft()
            for w in adj[v]:
                if w not in d:
                    d[w] = d[v] + 1; p[w] = v; q.append(w)
                elif w != p[v]: best = min(best, d[v] + d[w] + 1)
    return best


def splice(adjA, adjB):
    """cubic connected sum: delete one edge from each, cross-join."""
    nA = len(adjA)
    A = [set(x) for x in adjA]
    B = [set(y + nA for y in x) for x in adjB]
    a0, a1 = 0, adjA[0][0]
    b0, b1 = nA, adjB[0][0] + nA
    A[a0].discard(a1); A[a1].discard(a0)
    B[b0 - nA].discard(b1); B[b1 - nA].discard(b0)
    G = A + B
    G[a0].add(b0); G[b0].add(a0)
    G[a1].add(b1); G[b1].add(a1)
    return [sorted(x) for x in G]


def random_cubic(n, rng, mingirth=3, tries=40000):
    for _ in range(tries):
        st = [v for v in range(n) for _ in range(3)]; rng.shuffle(st)
        adj = [set() for _ in range(n)]; ok = True
        for i in range(0, len(st), 2):
            a, b = st[i], st[i + 1]
            if a == b or b in adj[a]: ok = False; break
            adj[a].add(b); adj[b].add(a)
        if not ok: continue
        seen = {0}; s = [0]
        while s:
            v = s.pop()
            for w in adj[v]:
                if w not in seen: seen.add(w); s.append(w)
        if len(seen) != n: continue
        A = [sorted(x) for x in adj]
        if girth(A) < mingirth: continue
        return A
    return None


# ------------------------------------------------------------- cycles
def all_cycles(adj, Lmax):
    """every simple cycle of length <= Lmax, as (vertex tuple, edge set)."""
    n = len(adj); out = []
    for s in range(n):
        path = [s]; seen = {s}
        def dfs(v):
            for w in adj[v]:
                if w == s and len(path) >= 3:
                    if path[1] < path[-1]:
                        es = frozenset(frozenset((path[i], path[(i + 1)
                                       % len(path)]))
                                       for i in range(len(path)))
                        out.append((tuple(path), es))
                elif w > s and w not in seen and len(path) < Lmax:
                    seen.add(w); path.append(w); dfs(w)
                    path.pop(); seen.discard(w)
        dfs(s)
    return out


def single_path(edges):
    """do these edges form one nonempty simple path?"""
    if not edges: return False
    deg = Counter()
    for e in edges:
        for v in e: deg[v] += 1
    if any(d > 2 for d in deg.values()): return False
    ends = [v for v, d in deg.items() if d == 1]
    if len(ends) != 2: return False
    nb = {}
    for e in edges:
        u, v = tuple(e)
        nb.setdefault(u, []).append(v); nb.setdefault(v, []).append(u)
    seen = {ends[0]}; st = [ends[0]]
    while st:
        v = st.pop()
        for w in nb[v]:
            if w not in seen: seen.add(w); st.append(w)
    return len(seen) == len(deg)


def T_by_pairs(cycles, m):
    """T_m by its definition: pairs meeting in one path, lengths sum m."""
    t = 0
    for i in range(len(cycles)):
        for j in range(i + 1, len(cycles)):
            X = cycles[i][1]; Y = cycles[j][1]
            if len(X) + len(Y) != m: continue
            if single_path(X & Y): t += 1
    return t


def attachments(adj, cyc, a):
    """paths of length a, both ends on the cycle, interior and edges off
    it.  returns the count (unordered)."""
    Zv = set(cyc[0]); Ze = cyc[1]
    cnt = 0
    for u in cyc[0]:
        st = [(u, [u], set([u]))]
        while st:
            v, path, used = st.pop()
            if len(path) - 1 == a:
                if v in Zv and v != u and path[-2] not in Zv or \
                   (a == 1 and v in Zv and v != u):
                    cnt += 1
                continue
            for w in adj[v]:
                if frozenset((v, w)) in Ze: continue
                if w in used: continue
                if w in Zv and len(path) - 1 + 1 < a: continue
                st.append((w, path + [w], used | {w}))
    return cnt // 2


def T_by_attach(adj, cycles, m, amax=None):
    """T_m by the Attachment Identity."""
    if amax is None: amax = m // 2
    tot = 0; per_a = Counter()
    for a in range(1, amax + 1):
        z = m - 2 * a
        if z < 3: break
        for cyc in cycles:
            if len(cyc[1]) != z: continue
            k = attachments(adj, cyc, a)
            tot += k; per_a[a] += k
    return tot, per_a


ZOO = [
    ('Petersen',        gp(5, 2),                            10),
    ('Heawood',         lcf(14, [5, -5]),                    12),
    ('Moebius-Kantor',  gp(8, 3),                            12),
    ('Pappus',          lcf(18, [5, 7, -7, 7, -7, -5]),      12),
    ('Desargues',       gp(10, 3),                           12),
    ('Nauru',           gp(12, 5),                           12),
]

print('=' * 70)
print('(A) THE ATTACHMENT IDENTITY, EXACT')
print('=' * 70)
print('  T_m = #{(Z,P) : Z cycle, P attachment, |Z| + 2|P| = m}')
print()
print('  %-16s %3s %3s %8s %10s %10s %8s'
      % ('graph', 'n', 'm', 'T (pairs)', 'T (attach)', '2m T_m', 'agree'))
for name, G, m in ZOO:
    n = len(G)
    cyc = all_cycles(G, m - 2)
    tp = T_by_pairs(cyc, m)
    ta, per_a = T_by_attach(G, cyc, m)
    print('  %-16s %3d %3d %8d %10d %10d %8s'
          % (name, n, m, tp, ta, 2 * m * tp, 'OK' if tp == ta else 'FAIL'))
print()
print('  454.3 recorded T_10(Petersen)=60, T_12(Heawood)=336,')
print('  T_12(Moebius-Kantor)=216 with rank-2 census = 2m T_m.')

print()
print('=' * 70)
print('(B) WHERE THE ENERGY SITS: THE a-PROFILE')
print('=' * 70)
print('  R_a(m-2a) and its weight 4^{-a}.  a=1 is chords.')
print()
for name, G, m in ZOO:
    cyc = all_cycles(G, m - 2)
    ta, per_a = T_by_attach(G, cyc, m)
    prof = ' '.join('a=%d:%d' % (a, per_a[a])
                    for a in sorted(per_a) if per_a[a])
    print('  %-16s m=%2d  T=%5d   %s' % (name, m, ta, prof))

print()
print('=' * 70)
print('(C) THE CHORD CRITERION, EVALUATED')
print('=' * 70)
print('  chi_a = avg # a-attachments per (m-2a)-cycle')
print('  gamma_z = 2z c_z / 2^z    (cycle density vs the main term)')
print('  criterion:  S = SUM_a chi_a gamma_{m-2a} 4^{-a} m/(m-2a) < 1')
print('  S < 1  =>  T_m < 2^m/2m  =>  c_m > 0.')
print()
print('  %-16s %3s %8s %8s %10s %9s %7s'
      % ('graph', 'm', 'chi_1', 'gamma', 'S', 'T_m', '2^m/2m'))
for name, G, m in ZOO:
    cyc = all_cycles(G, m - 2)
    ta, per_a = T_by_attach(G, cyc, m)
    S = 0.0; chi1 = 0.0
    for a in sorted(per_a):
        z = m - 2 * a
        cz = sum(1 for c in cyc if len(c[1]) == z)
        if cz == 0: continue
        chi = per_a[a] / cz
        gam = 2.0 * z * cz / 2 ** z
        S += chi * gam * 4.0 ** (-a) * m / z
        if a == 1: chi1 = chi
    tgt = 2.0 ** m / (2 * m)
    print('  %-16s %3d %8.2f %8.3f %10.4f %9d %7.0f'
          % (name, m, chi1, gam, S, ta, tgt))
print()
print('  (these are Moore-extremal graphs: m ~ 2g ~ n, exactly the')
print('   regime 455.2 showed the criterion cannot hold in.  The test')
print('   that matters is m << n, below.)')

print()
print('=' * 70)
print('(D) THE REGIME THAT MATTERS: m << n, RANDOM CUBIC')
print('=' * 70)
print('  %6s %4s %4s %9s %9s %9s %9s'
      % ('n', 'g', 'm', 'c_{m-2}', 'chi_1', 'S', 'verdict'))
rng = random.Random(20260808)
for n in (20, 30, 40, 60, 80, 120):
    G = random_cubic(n, rng, mingirth=5)
    if G is None: continue
    g = girth(G); m = 8
    while m < g: m *= 2
    if m - 2 > 16: continue
    cyc = all_cycles(G, m - 2)
    ta, per_a = T_by_attach(G, cyc, m)
    S = 0.0; chi1 = 0.0; c1 = 0
    for a in sorted(per_a):
        z = m - 2 * a
        cz = sum(1 for c in cyc if len(c[1]) == z)
        if cz == 0: continue
        chi = per_a[a] / cz
        gam = 2.0 * z * cz / 2 ** z
        S += chi * gam * 4.0 ** (-a) * m / z
        if a == 1: chi1 = chi; c1 = cz
    print('  %6d %4d %4d %9d %9.3f %9.5f %9s'
          % (n, g, m, c1, chi1, S, 'PASS' if S < 1 else 'fail'))

print()
print('=' * 70)
print('(E) THE DECISIVE TEST: DOES 456\'s (NC)-BREAKER BREAK THIS?')
print('=' * 70)
print('  456.2 killed (NC) by splicing a Petersen onto random cubic')
print('  blocks: all 5-cycles stay inside ten vertices, so the')
print('  per-vertex concentration K_obs blows up like n/y.')
print('  VACUOUS HERE (374.1): m=8 forces T_8=0 at girth>=5, so every')
print('  S below is 0 for a trivial reason.  Rerun at m=16: nc468.py.')
print()
print('  %-22s %4s %4s %9s %9s %9s %8s'
      % ('graph', 'n', 'g', 'K_obs(g)', 'chi_1', 'S', 'crit'))
rng2 = random.Random(31337)
tests = []
base = gp(5, 2)
for nb in (14, 20, 26):
    B = random_cubic(nb, rng2, mingirth=5)
    if B is None: continue
    tests.append(('Petersen+rand%d' % nb, splice(base, B)))
tests.append(('Petersen+Petersen', splice(base, gp(5, 2))))
for name, G in tests:
    n = len(G); g = girth(G)
    cg = all_cycles(G, g)
    cg = [c for c in cg if len(c[1]) == g]
    hits = Counter()
    for c in cg:
        for v in c[0]: hits[v] += 1
    mean = g * len(cg) / n if n else 0
    K = (max(hits.values()) / mean) if mean > 0 else 0
    m = 8
    while m < g: m *= 2
    cyc = all_cycles(G, m - 2)
    ta, per_a = T_by_attach(G, cyc, m)
    S = 0.0; chi1 = 0.0
    for a in sorted(per_a):
        z = m - 2 * a
        cz = sum(1 for c in cyc if len(c[1]) == z)
        if cz == 0: continue
        chi = per_a[a] / cz
        gam = 2.0 * z * cz / 2 ** z
        S += chi * gam * 4.0 ** (-a) * m / z
        if a == 1: chi1 = chi
    print('  %-22s %4d %4d %9.2f %9.3f %9.5f %8s'
          % (name, n, g, K, chi1, S, 'PASS' if S < 1 else 'FAIL'))

print()
print('=' * 70)
print('(F) CAN THE ADVERSARY PUSH chi_1 ABOVE 4?')
print('=' * 70)
print('  chi_1 is the average number of chords of a long cycle.  A')
print('  z-cycle with h chords has 2h of its z legs landing back on')
print('  itself, so h <= (z-3)/2 by 3-edge-connectivity -- but that is')
print('  one cycle.  The criterion needs the AVERAGE over all of them.')
print()
print('  %-18s %4s %4s %5s %8s %8s %8s'
      % ('graph', 'n', 'z', 'c_z', 'max ch', 'avg ch', 'bound'))
for name, G in ([('Petersen', gp(5, 2)), ('Heawood', lcf(14, [5, -5])),
                 ('Nauru', gp(12, 5)),
                 ('Moebius-Kantor', gp(8, 3))] +
                [(nm, g) for nm, g in tests[:2]]):
    n = len(G)
    for z in (6, 8, 10):
        cyc = [c for c in all_cycles(G, z) if len(c[1]) == z]
        if not cyc: continue
        ch = [attachments(G, c, 1) for c in cyc]
        print('  %-18s %4d %4d %5d %8d %8.3f %8.1f'
              % (name, n, z, len(cyc), max(ch),
                 sum(ch) / len(ch), (z - 3) / 2.0))
