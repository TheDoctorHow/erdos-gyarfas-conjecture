"""moment457.py -- A FIRST-MOMENT ATTACK ON THE SIEVE GAP.

452.5 measured the non-simple geodesic corrections at m >= 2g and
found corr/2^m in the range 0.66 .. 2.53 - the same order as the main
term - and concluded the sieve does not converge.  That measurement
was taken on graphs where m was comparable to n.  This file asks
whether the correction is really Theta(2^m) or whether it is
Theta(2^m / n), which would be a completely different verdict.

THE FIRST-MOMENT PREDICTION.  In the configuration model the expected
number of copies of a fixed subgraph F in a random cubic graph scales
as (d-1)^{e(F)} n^{v(F)-e(F)} with d-1 = 2.  For a CYCLE, v = e, so
the count is Theta(2^a): this is the main term, and indeed the
expected number of a-cycles is ~ 2^a / 2a, giving 2a c_a ~ 2^a.

For a THETA (two branch vertices, paths a,b,c) and for a DUMBBELL
(two disjoint cycles joined by a path), v = e - 1.  So both are
Theta(2^{e} / n).  Every rank-2 support carries a factor 1/n that the
main term does not, and a rank-r support carries n^{1-r}.

Summing the shapes: for dumbbells with a + b + 2l = m the contribution
to tr(T^m) is about 4m . 2^{a+b+l} / (n . 4ab), and
SUM_{a+b=m-2l} 1/(ab) = O(log m / m), so the total over all l >= 1 is

    rank-2 contribution  =  O( 2^m log m / n ).

PREDICTION:   corr(m) . n / 2^m   is roughly CONSTANT in n.

If true, the sieve DOES converge whenever n is large compared with
the graph's own short-cycle density, and 452's window is not capped
at 2g at all for such graphs - it extends until the corrections catch
up.  The measurements of 452.5 would then have been an artifact of
m ~ n on tiny graphs, not a real obstruction.

The honest caveat, stated up front: this is a statement about TYPICAL
cubic graphs.  A worst-case bound cannot recover the 1/n, because an
adversary may cluster its short cycles.  What the experiment settles
is which of the two scalings is the truth in the generic case, and
therefore what a counterexample would have to look like.
"""
import random, sys
import numpy as np

sys.setrecursionlimit(100000)


def random_cubic(n, rng, tries=4000):
    """Configuration model with rejection until simple and connected."""
    for _ in range(tries):
        stubs = [v for v in range(n) for _ in range(3)]
        rng.shuffle(stubs)
        adj = [set() for _ in range(n)]
        ok = True
        for i in range(0, len(stubs), 2):
            a, b = stubs[i], stubs[i + 1]
            if a == b or b in adj[a]:
                ok = False; break
            adj[a].add(b); adj[b].add(a)
        if not ok:
            continue
        seen = {0}; st = [0]
        while st:
            v = st.pop()
            for w in adj[v]:
                if w not in seen:
                    seen.add(w); st.append(w)
        if len(seen) != n:
            continue
        return [sorted(a) for a in adj]
    return None


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


def traces_np(adj, M):
    """tr(T^m), m = 1..M, via the integer recursion p_m(A)."""
    n = len(adj)
    A = np.zeros((n, n), dtype=np.int64)
    for v in range(n):
        for w in adj[v]:
            A[v, w] = 1
    Pm2 = 2 * np.eye(n, dtype=np.int64)
    Pm1 = A.copy()
    ex = n // 2
    out = [None, int(np.trace(Pm1))]
    for m in range(2, M + 1):
        Cur = A @ Pm1 - 2 * Pm2
        Pm2, Pm1 = Pm1, Cur
        out.append(int(np.trace(Pm1)) + ex * (1 + (-1) ** m))
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


def correction(adj, m):
    """corr = tr(T^m) - (all rank-1 contributions)."""
    N = traces_np(adj, m)
    cc = cycle_counts(adj, m)
    r1 = 0
    for j in range(1, m + 1):
        if m % j == 0:
            l = m // j
            r1 += 2 * l * cc[l]
    return N[m], r1, N[m] - r1, cc


rng = random.Random(20260808)

print('=== (1) DOES THE CORRECTION SCALE LIKE 2^m OR LIKE 2^m / n?'
      ' ===', flush=True)
print('  random cubic graphs, m = 12 fixed, averaged over samples')
print()
print('     n  samples   mean 2m*c_m/2^m   mean corr/2^m'
      '   mean corr*n/2^m')
m = 12
for n in (20, 30, 40, 60, 80, 120, 160, 240):
    reps = 12 if n <= 80 else 6
    a1 = []; a2 = []
    got = 0
    for _ in range(reps):
        adj = random_cubic(n, rng)
        if adj is None: continue
        Nm, r1, corr, cc = correction(adj, m)
        a1.append(2 * m * cc[m] / 2.0 ** m)
        a2.append(corr / 2.0 ** m)
        got += 1
    if not got: continue
    mc = sum(a1) / got; cr = sum(a2) / got
    print('  %5d %8d %17.4f %15.4f %17.3f'
          % (n, got, mc, cr, cr * n), flush=True)

print()
print('  If corr = Theta(2^m) the middle column stays flat and the')
print('  right column grows linearly.  If corr = Theta(2^m / n) the')
print('  middle column decays like 1/n and the RIGHT COLUMN IS FLAT.',
      flush=True)

print()
print('=== (2) THE SAME TEST AT SEVERAL m ===', flush=True)
print('     n     m   2m*c_m/2^m    corr/2^m    corr*n/2^m', flush=True)
for n in (40, 80, 160):
    for m in (10, 12, 14):
        reps = 8 if n <= 80 else 4
        a1 = []; a2 = []; got = 0
        for _ in range(reps):
            adj = random_cubic(n, rng)
            if adj is None: continue
            Nm, r1, corr, cc = correction(adj, m)
            a1.append(2 * m * cc[m] / 2.0 ** m)
            a2.append(corr / 2.0 ** m)
            got += 1
        if not got: continue
        print('  %5d %5d %12.4f %11.4f %13.3f'
              % (n, m, sum(a1) / got, sum(a2) / got,
                 (sum(a2) / got) * n), flush=True)

print()
print('=== (3) THE MAIN TERM: does 2m c_m / 2^m -> 1 ? ===',
      flush=True)
print('  the trace main term is 2^m + 1 + n; if the corrections are')
print('  o(2^m) then 2m c_m / 2^m -> 1 and cycles of length m are')
print('  forced by counting alone.')
print()
print('     n     m   2m*c_m/2^m   (target 1.000)', flush=True)
for n in (80, 160, 240, 320):
    for m in (10, 12):
        reps = 6 if n <= 160 else 3
        a1 = []; got = 0
        for _ in range(reps):
            adj = random_cubic(n, rng)
            if adj is None: continue
            _, _, _, cc = correction(adj, m)
            a1.append(2 * m * cc[m] / 2.0 ** m)
            got += 1
        if got:
            print('  %5d %5d %12.4f' % (n, m, sum(a1) / got), flush=True)
