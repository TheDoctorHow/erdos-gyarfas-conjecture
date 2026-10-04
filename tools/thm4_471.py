"""thm4_471.py -- THEOREM IV, THE NICER FORM.

465 left Theorem IV split: same-stratum sharp at 2^t, progress
growing.  This file asks what the progress currency actually IS once
the Gram encoding is pushed one step further, and finds that it is
not a Davenport problem at all.

STEP 1 - THE CYCLE FORM.  In 465 the model was carried on the ear
lifts z_i.  Substituting z_i = C_i - x, where C_i is the composite
cycle of ear i, and writing p = #{sigma_i = +1}:

    V(sigma) = (1 - nu) x + SUM sigma_i z_i = (1 - p) x + SUM sigma_i C_i .

The coefficient of x is 1 - nu in the z-basis and 1 - p in the
C-basis: the arc-switch involution is exactly the exchange
(C_i <-> Z_i, + <-> -), where Z_i = x - C_i is the difference cycle.
Every ear offers a PAIR OF CYCLES SUMMING TO THE BASE, C_i + Z_i = x.

STEP 2 - THE OFFSET, AND HOW TO DELETE IT.  The base contributes
y_0^2 L == 2^t . [y_0 odd] (mod 2^{t+1}).  That single term is the
whole difference between the two currencies: same-stratum is the
congruence mod 2^t, where it vanishes; progress is mod 2^{t+1},
where it does not.  It is why 445.6 called Theorem IV ASYMMETRIC.

  But y_0 = 1 - nu is ZERO at nu = 1.  At exactly one flip the base
  cycle leaves the expression entirely (the Odd-Flip Lemma of 465.1),
  and the progress target becomes

      | SUM sigma_i z_i |^2  ==  0   (mod 2^{t+1}),

  a symmetric, offset-free question about the ear lifts alone.

STEP 3 - AT TWO ELEMENTS IT IS THE DIFFERENCE CRITERION.  Take
nu = 1 and |S| = 2, sigma = (+1, -1).  Then

      V = z_i - z_j = C_i - C_j = Z_j - Z_i ,

so the composite is C_i Delta C_j, and the progress target reads

      |C_i| + |C_j| - 2 |C_i cap C_j|  ==  0   (mod 2^{t+1}).

That is 457.4's DIFFERENCE CRITERION, at the composite scale rather
than at the girth.  Theorem IV's progress currency is not a zero-sum
constant: it is the trace frame's difference closure.

This file verifies steps 1-3 on real graphs, shows why the abstract
Davenport route provably cannot deliver it (the couplings are free),
and measures the difference criterion where it actually lives.
"""
import sys, itertools
import numpy as np

sys.setrecursionlimit(100000)
P = print


def gp(n, k):
    adj = {i: set() for i in range(2 * n)}
    def add(a, b):
        if a != b: adj[a].add(b); adj[b].add(a)
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
    return [sorted(adj[i]) for i in range(n)]


def cyc_len(adj, L):
    n = len(adj); out = []
    for s in range(n):
        path = [s]; seen = {s}
        def dfs(v):
            if len(path) == L:
                if s in adj[v] and path[1] < path[-1]: out.append(tuple(path))
                return
            for w in adj[v]:
                if w > s and w not in seen:
                    seen.add(w); path.append(w); dfs(w)
                    path.pop(); seen.discard(w)
        dfs(s)
    return out


def eidx(adj):
    E = sorted({(min(u, v), max(u, v)) for u in range(len(adj)) for v in adj[u]})
    return E, {e: i for i, e in enumerate(E)}


def circ(cyc, ix, m):
    v = np.zeros(m, dtype=np.int64)
    for i in range(len(cyc)):
        a, b = cyc[i], cyc[(i + 1) % len(cyc)]
        if a < b: v[ix[(a, b)]] += 1
        else: v[ix[(b, a)]] -= 1
    return v


ZOO = [('Petersen', gp(5, 2), 5), ('Moebius-Kantor', gp(8, 3), 6),
       ('Heawood', lcf(14, [5, -5]), 6), ('Pappus',
        lcf(18, [5, 7, -7, 7, -7, -5]), 6), ('Desargues', gp(10, 3), 6)]

P('=' * 72)
P('(A) THE CYCLE FORM:  V = (1-nu)x + SUM sigma_i z_i = (1-p)x + SUM sigma_i C_i',
  flush=True)
P('=' * 72)
P('  p = #{sigma_i = +1}, nu = #{sigma_i = -1}.  Also C_i + Z_i = x.')
P()
P('  %-16s %3s %3s %8s %10s %12s'
  % ('graph', 'L', 'k', 'sigmas', 'z=C form', 'C_i+Z_i=x'), flush=True)
DATA = []
for name, G, Lb in ZOO:
    E, ix = eidx(G); m = len(E)
    base = cyc_len(G, Lb)[0]
    x = circ(base, ix, m); L = int(x @ x)
    Cs = []
    for ln in range(Lb, Lb + 4):
        for cy in cyc_len(G, ln):
            if tuple(cy) == base: continue
            cv = circ(cy, ix, m)
            if int(cv @ x) <= 0: cv = -cv
            if int(cv @ x) <= 0: continue
            Cs.append(cv)
            if len(Cs) >= 5: break
        if len(Cs) >= 5: break
    k = len(Cs)
    zs = [c - x for c in Cs]
    Zs = [x - c for c in Cs]
    ok1 = ok2 = tot = 0
    for sig in itertools.product((-1, 0, 1), repeat=k):
        if not any(sig): continue
        nu = sum(1 for s in sig if s == -1); p = sum(1 for s in sig if s == 1)
        V1 = (1 - nu) * x + sum((s * z for s, z in zip(sig, zs)),
                                np.zeros(m, dtype=np.int64))
        V2 = (1 - p) * x + sum((s * c for s, c in zip(sig, Cs)),
                               np.zeros(m, dtype=np.int64))
        tot += 1
        if np.array_equal(V1, V2): ok1 += 1
    for c, z in zip(Cs, Zs):
        if np.array_equal(c + z, x): ok2 += 1
    P('  %-16s %3d %3d %8d %10s %12s'
      % (name, L, k, tot, '%d/%d' % (ok1, tot), '%d/%d' % (ok2, k)), flush=True)
    DATA.append((name, G, x, L, Cs, zs, Zs, m))

P()
P('=' * 72)
P('(B) AT nu = 1, |S| = 2 THE COMPOSITE IS C_i DELTA C_j', flush=True)
P('=' * 72)
P('  sigma = (+1,-1) gives p = 1, so V = C_i - C_j = Z_j - Z_i, and')
P('  |V|^2 = |C_i| + |C_j| - 2|C_i cap C_j| = 457.4\'s difference.')
P()
P('  %-16s %8s %10s %14s' % ('graph', 'pairs', 'V=C_i-C_j', 'len = diff fmla'),
  flush=True)
for name, G, x, L, Cs, zs, Zs, m in DATA:
    k = len(Cs); okv = okl = tot = 0
    for i in range(k):
        for j in range(k):
            if i == j: continue
            V = zs[i] - zs[j]
            tot += 1
            if np.array_equal(V, Cs[i] - Cs[j]) and \
               np.array_equal(V, Zs[j] - Zs[i]): okv += 1
            ov = int(np.sum((Cs[i] != 0) & (Cs[j] != 0)))
            li = int(np.count_nonzero(Cs[i])); lj = int(np.count_nonzero(Cs[j]))
            if int(V @ V) == li + lj - 2 * int(Cs[i] @ Cs[j]): okl += 1
    P('  %-16s %8d %10s %14s'
      % (name, tot, '%d/%d' % (okv, tot), '%d/%d' % (okl, tot)), flush=True)

P()
P('=' * 72)
P('(C) WHY THE ABSTRACT DAVENPORT ROUTE CANNOT DELIVER IT', flush=True)
P('=' * 72)
P('  At nu=1, |S|=2 the target is l_i + l_j - 2 g_ij == 0 mod 2^{t+1}.')
P('  Abstractly g_ij is FREE mod 2^t, so 2 g_ij takes every even')
P('  residue mod 2^{t+1} and the adversary just avoids one value.')
P('  Hence the pairwise constant is INFINITE in the abstract model,')
P('  at every t and every k -- verified by construction:')
P()
for t in (1, 2, 3):
    MOD = 1 << (t + 1)
    for k in (4, 8, 16, 32):
        l = [1] * k                       # all odd -> l_i + l_j even
        bad = 0
        for i in range(k):
            for j in range(i + 1, k):
                g = ((l[i] + l[j]) // 2 + 1) % (1 << t)   # dodge the one value
                if (l[i] + l[j] - 2 * g) % MOD == 0: bad += 1
        P('    t=%d k=%-3d avoider valid: %s' % (t, k, bad == 0), flush=True)
P()
P('  So Theorem IV is NOT a zero-sum constant.  What pins g_ij is')
P('  geometry: g_ij = <z_i,z_j> is the overlap of two cycles, and')
P('  462.5\'s Band Lemma bounds where that overlap can sit.')

P()
P('=' * 72)
P('(D) THE DIFFERENCE CRITERION WHERE IT ACTUALLY LIVES', flush=True)
P('=' * 72)
P('  For each graph: base cycle x, all cycles C meeting x in a path,')
P('  and the differences |C_i Delta C_j|.  How often is the progress')
P('  target 0 mod 2^{t+1} hit, with t = v_2(L)?')
P()
P('  %-16s %3s %3s %6s %8s %10s %10s'
  % ('graph', 'L', 't', 'cycles', 'pairs', 'hits 2^{t+1}', 'rate'), flush=True)
for name, G, Lb in ZOO:
    E, ix = eidx(G); m = len(E)
    base = cyc_len(G, Lb)[0]
    x = circ(base, ix, m); L = int(x @ x)
    t = 0; LL = L
    while LL % 2 == 0: LL //= 2; t += 1
    Cs = []
    for ln in range(3, min(len(G), 13) + 1):
        for cy in cyc_len(G, ln):
            if tuple(cy) == base: continue
            cv = circ(cy, ix, m)
            if int(cv @ x) < 0: cv = -cv
            if int(cv @ x) <= 0: continue
            Cs.append(cv)
    MOD = 1 << (t + 1); hits = 0; tot = 0
    for i in range(len(Cs)):
        for j in range(i + 1, len(Cs)):
            V = Cs[i] - Cs[j]
            if not np.any(V): continue
            tot += 1
            if int(V @ V) % MOD == 0: hits += 1
    P('  %-16s %3d %3d %6d %8d %10d %10s'
      % (name, L, t, len(Cs), tot, hits,
         '%.3f' % (hits / tot) if tot else '-'), flush=True)
P()
P('  t = v_2(L) is small on these graphs, so this is a sanity')
P('  measurement of the mechanism, not a test of J7 (which needs')
P('  t >= 3, i.e. L divisible by 8).')
