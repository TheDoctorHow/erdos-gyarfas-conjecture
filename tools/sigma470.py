"""sigma470.py -- THE SIGMA-PAIRED DAVENPORT PROBLEM, RE-ENCODED.

445.5 posed it; 464.3 showed Track A's five statements reduce to it.

THE RE-ENCODING.  Ear families were carried as three separate objects
(shifts delta_i, couplings c_ij, arc data d_i and a_i) with a side
condition for shallowness.  All of it is ONE GRAM MATRIX.

Write z_i for the arc-lift of ear i.  The Arc-Switch Law (445.1) says
the partner is z^- = -x - z^+, so a selection with signs is

    sigma in {-1, 0, +1}^k,   nu = #{sigma_i = -1},
    x + W = (1 - nu) x + SUM_i sigma_i z_i =: V(sigma).

By the Ear Identity the composite length is exactly |V|^2.  Put
y = (1 - nu, sigma_1, ..., sigma_k) and let G be the Gram matrix of
(x, z_1, ..., z_k).  Then

    COMPOSITE LENGTH  =  y^T G y .

So the whole model is a symmetric (k+1)x(k+1) integer matrix, and:

  * G_00 = L == 2^t (mod 2^{t+1}), since L = 2^t m with m odd;
  * off-diagonal entries appear DOUBLED, so they matter only mod 2^t
    - this is 445.1's "couplings are always even", now automatic;
  * diagonal entries l_i = |z_i|^2 matter mod 2^{t+1};
  * PROGRESS  (composite has v_2 > t) is  y^T G y == 0 mod 2^{t+1};
  * SAME-STRATUM              is  y^T G y == 0 mod 2^t;
  * SHALLOWNESS of a single ear is just the k = 1 case of the same
    condition - no side condition at all.

Four objects and a constraint collapse to one matrix and one
congruence.  The search space per new element drops from
2^{2t+1} . 2^{tk} bookkeeping tuples to 2^t . 2^{t+1} . 2^{tk} matrix
entries, and every hit test becomes one vectorised quadratic form.

THE ODD-FLIP LEMMA (the pairing's engine).  y_0 = 1 - nu vanishes at
nu = 1, so flipping exactly one element deletes x from V entirely:
V = SUM sigma_i z_i, and the composite length is a pure Gram sum of
the ear-lifts with no L and no shifts in it.  445.5's invincible
adversary was built to dodge shift arithmetic and has no answer to
this.

Parts: (A) the reduction verified against real cycle lengths;
       (B) 445.5's invincible adversary under pairing;
       (C) the constants D(t), exhaustive.
"""
import sys, itertools
import numpy as np
from collections import deque

sys.setrecursionlimit(100000)
P = print


def gp(n, k):
    adj = {i: set() for i in range(2 * n)}
    def add(a, b):
        if a != b: adj[a].add(b); adj[b].add(a)
    for i in range(n):
        add(i, (i + 1) % n); add(i, n + i)
        add(n + i, n + (i + k) % n)
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


def cycles_of_length(adj, L):
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


def edge_index(adj):
    E = sorted({(min(u, v), max(u, v)) for u in range(len(adj)) for v in adj[u]})
    return E, {e: i for i, e in enumerate(E)}


def circ(cyc, ix, m):
    v = np.zeros(m, dtype=np.int64)
    for i in range(len(cyc)):
        a, b = cyc[i], cyc[(i + 1) % len(cyc)]
        if a < b: v[ix[(a, b)]] += 1
        else: v[ix[(b, a)]] -= 1
    return v


P('=' * 72)
P('(A) THE REDUCTION: COMPOSITE LENGTH = y^T G y', flush=True)
P('=' * 72)
P('  z_i = C_i - x with C_i oriented to agree with x on shared edges.')
P('  Then |z_i|^2 = |C_i| + L - 2|C_i cap x| = the difference cycle')
P('  (457.1), and every sigma-composite length is y^T G y.')
P()
P('  %-16s %4s %3s %8s %14s %14s'
  % ('graph', 'L', 'k', 'sigmas', 'y^TGy = |V|^2', 'V a cycle => len'), flush=True)
for name, Gr, Lb, Lmax in (('Petersen', gp(5, 2), 5, 8),
                           ('Moebius-Kantor', gp(8, 3), 6, 8),
                           ('Heawood', lcf(14, [5, -5]), 6, 8)):
    E, ix = edge_index(Gr); m = len(E)
    base = cycles_of_length(Gr, Lb)[0]
    x = circ(base, ix, m); L = int(x @ x)
    zs = []
    for ln in range(Lb, Lmax + 1):
        for cyc in cycles_of_length(Gr, ln):
            if tuple(cyc) == base: continue
            cv = circ(cyc, ix, m)
            if int(cv @ x) < 0: cv = -cv
            if int(cv @ x) <= 0: continue          # must share a path
            zs.append(cv - x)
            if len(zs) >= 4: break
        if len(zs) >= 4: break
    k = len(zs)
    M = np.zeros((k + 1, m), dtype=np.int64); M[0] = x
    for i, z in enumerate(zs): M[i + 1] = z
    G = M @ M.T
    okq = oklen = tot = 0
    for sig in itertools.product((-1, 0, 1), repeat=k):
        if not any(sig): continue
        nu = sum(1 for s in sig if s == -1)
        y = np.array([1 - nu] + list(sig), dtype=np.int64)
        V = y @ M
        tot += 1
        if int(y @ G @ y) == int(V @ V): okq += 1
        if np.all(np.abs(V) <= 1):
            if int(V @ V) == int(np.count_nonzero(V)): oklen += 1
        else: oklen += 1
    P('  %-16s %4d %3d %8d %14s %14s'
      % (name, L, k, tot, '%d/%d' % (okq, tot), '%d/%d' % (oklen, tot)),
      flush=True)

P()
P('=' * 72)
P("(B) 445.5's INVINCIBLE ADVERSARY, UNDER PAIRING", flush=True)
P('=' * 72)
P('  Unpaired at t = 2: delta_i == 2, c_ij == 4 (mod 8) gives')
P('  Q = 2|S|^2 in {0,2} mod 8 - never the progress target, at every')
P('  k.  In Gram terms that is G_ii + 2G_0i == 2 and 2G_ij == 4.')
P()


def sigmas(k):
    S = np.array(list(itertools.product((-1, 0, 1), repeat=k)), dtype=np.int64)
    return S


t = 2; MOD = 1 << (t + 1); Lv = 1 << t
P('  %4s %26s %10s' % ('k', 'surviving adversaries', 'verdict'), flush=True)
for k in (1, 2, 3, 4):
    S = sigmas(k)
    nz = np.any(S != 0, axis=1)
    S = S[nz]
    nu = np.sum(S == -1, axis=1)
    Y = np.concatenate([(1 - nu).reshape(-1, 1), S], axis=1)
    surv = 0; tot = 0
    # G_0i free mod 4, G_ii determined by delta^+ = G_ii + 2 G_0i == 2
    for p in itertools.product(range(4), repeat=k):
        G = np.zeros((k + 1, k + 1), dtype=np.int64)
        G[0, 0] = Lv
        for i in range(k):
            G[0, i + 1] = G[i + 1, 0] = p[i]
            G[i + 1, i + 1] = (2 - 2 * p[i]) % MOD
        for i in range(k):
            for j in range(i + 1, k):
                G[i + 1, j + 1] = G[j + 1, i + 1] = 2      # 2G_ij = 4
        tot += 1
        q = np.einsum('sa,ab,sb->s', Y, G, Y) % MOD
        if not np.any(q == 0): surv += 1
    P('  %4d %26s %10s'
      % (k, '%d / %d' % (surv, tot), 'DEAD' if surv == 0 else 'alive'), flush=True)
P()
P('  MECHANISM.  At nu = 1 the coefficient of x is 1 - nu = 0, so the')
P('  composite is |z_i +- z_j - z_l|^2 with no L and no shift in it.')

P()
P('=' * 72)
P('(C) THE CONSTANTS  D(t)', flush=True)
P('=' * 72)
P('  D(t) = least k such that EVERY Gram instance on k ears admits a')
P('  nonzero sigma with y^TGy == 0 (mod 2^{t+1} progress / 2^t same).')
P('  Exhaustive batched DFS; a hit at depth j survives every')
P('  extension, so pruning is sound.  Elements kept in canonical')
P('  (non-decreasing) order.')
P()


def sigset(k):
    """all nonzero sigma in {-1,0,1}^k, plus nu."""
    if k == 0: return np.zeros((0, 0), dtype=np.int64)
    S = np.array(list(itertools.product((-1, 0, 1), repeat=k)), dtype=np.int64)
    return S[np.any(S != 0, axis=1)]


def constants(t, mode, maxk=9, seconds=120.0):
    """largest avoider by batched DFS.  Returns (bestk, witness, done)."""
    import time
    MOD = 1 << (t + 1); Lv = 1 << t
    M2 = MOD if mode == 'progress' else (1 << t)
    half = 1 << t
    SIG = {k: sigset(k) for k in range(0, maxk + 2)}
    COLS = {k: (np.array(list(itertools.product(range(half), repeat=k)),
                         dtype=np.int64).reshape(-1, k) if k else
                np.zeros((1, 0), dtype=np.int64)) for k in range(0, maxk + 1)}
    best = [0, None]; done = [True]; t0 = time.time()

    def rec(G, k, lastkey):
        if k > best[0]: best[0] = k; best[1] = G.copy()
        if k == maxk or time.time() - t0 > seconds:
            if k != maxk: done[0] = False
            return
        Sp = SIG[k + 1]                       # (N, k+1) child sigmas
        s = Sp[:, k]                          # last coordinate
        sg = Sp[:, :k]                        # first k
        nu = np.sum(Sp == -1, axis=1)
        y0 = 1 - nu
        Yold = np.concatenate([y0.reshape(-1, 1), sg], axis=1)   # (N, k+1)
        base = np.einsum('na,ab,nb->n', Yold, G, Yold)           # parent part
        R = COLS[k]                                              # (C, k)
        N = len(Sp)
        chunk = max(1, 4 * 10 ** 6 // max(N, 1))
        ss = (s * s).reshape(-1, 1)
        for p in range(half):
            for lo in range(0, len(R), chunk):
                Rc = R[lo:lo + chunk]
                SR = sg @ Rc.T if k else np.zeros((N, len(Rc)), dtype=np.int64)
                lead = (base.reshape(-1, 1)
                        + 2 * s.reshape(-1, 1) * ((y0 * p).reshape(-1, 1) + SR))
                for w in range(MOD):
                    if (p, w) < lastkey: continue
                    alive = ~np.any(((lead + ss * w) % M2) == 0, axis=0)
                    if not alive.any(): continue
                    for ci in np.nonzero(alive)[0]:
                        r = Rc[ci]
                        G2 = np.zeros((k + 2, k + 2), dtype=np.int64)
                        G2[:k + 1, :k + 1] = G
                        G2[0, k + 1] = G2[k + 1, 0] = p
                        if k:
                            G2[1:k + 1, k + 1] = r; G2[k + 1, 1:k + 1] = r
                        G2[k + 1, k + 1] = w
                        rec(G2, k + 1, (p, w))
                        if time.time() - t0 > seconds:
                            done[0] = False; return
    G0 = np.zeros((1, 1), dtype=np.int64); G0[0, 0] = Lv
    rec(G0, 0, (-1, -1))
    return best[0], best[1], done[0]


P('  %-4s %-10s %6s %10s %10s %8s'
  % ('t', 'target', '2^t', 'max avoid', 'D(t)', 'complete'), flush=True)
RES = {}
for tt in (1, 2, 3):
    for mode in ('same', 'progress'):
        bk, wit, done = constants(tt, mode, maxk=9, seconds=100.0)
        RES[(tt, mode)] = (bk, done)
        P('  %-4d %-10s %6d %10d %10s %8s'
          % (tt, mode, 1 << tt, bk,
             str(bk + 1) if done else '>= %d' % (bk + 1),
             'yes' if done else 'TIMEOUT'), flush=True)
        if wit is not None and bk > 0:
            # INDEPENDENT re-verification of the witness, brute force
            MODv = 1 << (tt + 1); M2v = MODv if mode == 'progress' else (1 << tt)
            kk = len(wit) - 1; hits = 0
            for sig in itertools.product((-1, 0, 1), repeat=kk):
                if not any(sig): continue
                nuv = sum(1 for z in sig if z == -1)
                yv = np.array([1 - nuv] + list(sig), dtype=np.int64)
                if int(yv @ wit @ yv) % M2v == 0: hits += 1
            P('       witness re-verified: %d hits (0 = genuine avoider) %s'
              % (hits, 'OK' if hits == 0 else '*** SEARCH BUG ***'), flush=True)
            P('       diag %s' % list(map(int, np.diag(wit))), flush=True)
P()
P('  Against the conjecture D(t) = O(2^t):', flush=True)
for (tt, mode), (bk, done) in sorted(RES.items()):
    P('    t=%d %-9s D%s%d   2^t = %d   ratio %s'
      % (tt, mode, ' = ' if done else ' >= ', bk + 1, 1 << tt,
         '%.2f' % ((bk + 1) / (1 << tt))), flush=True)
