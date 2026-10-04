"""bound458.py -- HOW FAR PAST 2g DOES THE SIEVE ACTUALLY SURVIVE?

454 reduced the trace frame to one inequality: the total non-simple
correction must be smaller than the main term, i.e.

    corr(m) = N_m - SUM_{j | m} 2(m/j) c_{m/j}   <   N_m / 2,

since then 2m c_m > N_m/2 > 0 and a cycle of length m exists.

THE DETERMINISTIC THETA BOUND.  Every term of the sieve identity is
non-negative, so c_x <= N_x / 2x EXACTLY, for every x.  A theta at m
is a pair of cycles whose lengths sum to m, so dropping the
requirement that they actually meet,

    T_m <= (1/2) SUM_{x+y=m} c_x c_y <= (1/8) SUM_{x+y=m} N_x N_y/(xy).

Writing kappa = max_x N_x / 2^x over g <= x <= m-g and using
SUM_{x=g}^{m-g} 1/(x(m-x)) = (2/m) ln((m-g)/g),

    T_m <= (kappa^2 / 4) (2^m / m) ln((m-g)/g),

so T_m < 2^m / 2m as soon as  kappa^2 ln((m-g)/g) < 2.  For
kappa = 1 that is m < (1 + e^2) g = 8.39 g - far past the 2g cap of
452.  But rank-3 supports appear before 3g (Petersen: first at
m = 12 = 2.4 g), so the bound above is not the whole story.

This file measures the truth: the FULL correction, all ranks, against
the main term, as a function of m/g.
"""
import sys
import numpy as np

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
    n = len(adj)
    A = np.zeros((n, n), dtype=object)
    for v in range(n):
        for w in adj[v]:
            A[v, w] = 1
    Pm2 = np.zeros((n, n), dtype=object)
    for i in range(n):
        Pm2[i, i] = 2
    Pm1 = A.copy()
    ex = n // 2
    out = [None, int(np.trace(Pm1))]
    for m in range(2, M + 1):
        Cur = A @ Pm1 - 2 * Pm2
        Pm2, Pm1 = Pm1, Cur
        out.append(int(np.trace(Pm1)) + ex * (1 + (-1) ** m))
    return out


def cycle_counts(adj, cap):
    """cnt[m] = # cycles of length m, all m <= cap, one DFS pass."""
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


ZOO = [('Petersen', gp(5, 2), 18),
       ('Heawood', lcf(14, [5, -5] * 7), 20),
       ('Moebius-Kantor', lcf(16, [5, -5] * 8), 20),
       ('Pappus', lcf(18, [5, 7, -7, 7, -7, -5] * 3), 19),
       ('Desargues', lcf(20, [5, -5, 9, -9] * 5), 19),
       ('McGee', lcf(24, [12, 7, -7] * 8), 20),
       ('Tutte-Coxeter', lcf(30, [-13, -9, 7, -7, 9, 13] * 5), 20)]

print('=== (1) THE FULL CORRECTION vs THE MAIN TERM ===', flush=True)
print('  survives while corr/N_m < 1/2, i.e. 2m c_m > N_m/2 > 0')
print()
print('  graph            g    m   m/g    2m c_m / N_m   corr / N_m'
      '   pair bound / (2^m/2m)')
survive = {}
for name, adj, M in ZOO:
    g = girth(adj)
    N = traces_np(adj, M)
    cc = cycle_counts(adj, M)
    last = None
    for m in range(2 * g, M + 1):
        if N[m] == 0:
            continue
        r1 = 0
        for j in range(1, m + 1):
            if m % j == 0:
                r1 += 2 * (m // j) * cc[m // j]
        corr = N[m] - r1
        main = 2 * m * cc[m] / N[m]
        pb = 0.0
        for x in range(g, m - g + 1):
            if x <= m - x:
                if x == m - x:
                    pb += cc[x] * (cc[x] - 1) / 2.0
                else:
                    pb += cc[x] * cc[m - x]
        print('  %-16s %2d %4d %5.2f %14.4f %12.4f %20.4f'
              % (name, g, m, m / g, main, corr / N[m],
                 pb / (2.0 ** m / (2 * m))), flush=True)
        if main > 0.5:
            last = m
    survive[name] = (g, last)

print()
print('=== (2) HOW FAR DID IT SURVIVE? ===', flush=True)
print('  graph            g   last m with 2m c_m > N_m/2   m/g'
      '   (452 cap was m/g = 2)')
for name, (g, last) in survive.items():
    print('  %-16s %2d %28s %5s'
          % (name, g, str(last),
             ('%.2f' % (last / g)) if last else '-'), flush=True)

print()
print('=== (3) THE WINDOW THIS BUYS ===', flush=True)
print('  a power of 2 is GUARANTEED in (m*, B g) when B g / m* >= 2,')
print('  and m* ~ 2 log2(2n).  With B = 2 (452) that needs')
print('  g >= 2 log2(2n); with larger B it needs g >= 2 log2(2n) . 2/B.')
print()
import math
print('    B     required girth as a fraction of the Moore bound'
      '   (Moore: g ~ 2 log2 n)')
for B in (2.0, 2.5, 3.0, 4.0, 8.39):
    print('  %5.2f %30.3f' % (B, 2.0 / B), flush=True)
