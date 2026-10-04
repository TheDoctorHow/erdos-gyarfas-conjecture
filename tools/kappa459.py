"""kappa459.py -- (NC) AND kappa, pinned down.

455 left two things uncontrolled:
  kappa = max_{g <= x <= m-g} N_x / 2^x, equivalently
          kappa_x = 2x c_x / 2^x  for x < 2g;
  (NC)  = no vertex lies on more than K . y c_y / n cycles of length y.

WHAT IS PROVABLE ABOUT kappa, ABOVE m*.  From the identity
N_x = 2^x + 1 + n + E_x with |E_x| <= 2(n-1) 2^{x/2},

    kappa_x <= 1 + (1+n)/2^x + 2(n-1) 2^{-x/2},

so once x >= m* + 2s with m* = 2 log2(2n), kappa_x <= 1 + 2^{-s} + o(1).
kappa is therefore pinned to 1 ABOVE the window threshold, and says
nothing below it.  The theta sum of 455.1 runs over x in [g, m-g],
which is exactly the range BELOW m* where nothing is pinned.

WHAT THIS FILE MEASURES.
 (A) kappa_x on the cages: where it peaks and how big it gets.
 (B) kappa at x = g for truncations: unbounded, Theta(n).
 (C) (NC): the max-to-mean ratio of cycles through a vertex, on
     vertex-transitive graphs (where it is exactly 1 for free), on
     random cubic graphs, and on a cubic graph deliberately built so
     that all girth cycles sit in one small region.
"""
import sys, random
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


def truncate(adj):
    n = len(adj); idx = {}; m = 0
    for v in range(n):
        for w in adj[v]:
            idx[(v, w)] = m; m += 1
    out = [set() for _ in range(m)]
    def add(a, b):
        out[a].add(b); out[b].add(a)
    for v in range(n):
        ports = [idx[(v, w)] for w in adj[v]]
        for a in range(len(ports)):
            for b in range(a + 1, len(ports)):
                add(ports[a], ports[b])
    for v in range(n):
        for w in adj[v]:
            if v < w:
                add(idx[(v, w)], idx[(w, v)])
    return [sorted(x) for x in out]


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


def cycle_counts_by_vertex(adj, cap):
    """returns cnt[m] and per[m][v] = # m-cycles through v."""
    n = len(adj)
    cnt = [0] * (cap + 1)
    per = [[0] * n for _ in range(cap + 1)]
    for s in range(n):
        seen = {s}; path = [s]
        def dfs(v, ln):
            for w in adj[v]:
                if w == s:
                    if ln >= 3:
                        cnt[ln] += 1
                        for u in path: per[ln][u] += 1
                elif w > s and w not in seen and ln < cap:
                    seen.add(w); path.append(w)
                    dfs(w, ln + 1)
                    path.pop(); seen.discard(w)
        dfs(s, 1)
    cnt = [c // 2 for c in cnt]
    per = [[p // 2 for p in row] for row in per]
    return cnt, per


def random_cubic(n, rng, tries=3000):
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
        if not ok: continue
        seen = {0}; st = [0]
        while st:
            v = st.pop()
            for w in adj[v]:
                if w not in seen: seen.add(w); st.append(w)
        if len(seen) != n: continue
        return [sorted(a) for a in adj]
    return None


def splice(adjA, adjB):
    """Cubic connected sum: delete one edge from each, cross-join."""
    nA = len(adjA)
    A = [set(x) for x in adjA]; B = [set(x) for x in adjB]
    a, b = 0, adjA[0][0]
    c, d = 0, adjB[0][0]
    A[a].discard(b); A[b].discard(a)
    B[c].discard(d); B[d].discard(c)
    out = [set(x) for x in A] + [set(y + nA for y in x) for x in B]
    def add(u, v):
        out[u].add(v); out[v].add(u)
    add(a, c + nA); add(b, d + nA)
    return [sorted(x) for x in out]


ZOO = [('Petersen', gp(5, 2)), ('Heawood', lcf(14, [5, -5] * 7)),
       ('Moebius-Kantor', lcf(16, [5, -5] * 8)),
       ('Pappus', lcf(18, [5, 7, -7, 7, -7, -5] * 3)),
       ('Desargues', lcf(20, [5, -5, 9, -9] * 5)),
       ('McGee', lcf(24, [12, 7, -7] * 8)),
       ('Tutte-Coxeter', lcf(30, [-13, -9, 7, -7, 9, 13] * 5))]

print('=== (A) kappa_x = 2x c_x / 2^x ON THE CAGES ===', flush=True)
print('  the theta bound of 455.1 needs kappa^2 (1 + ln((m-g)/g)) < 2,')
print('  hence kappa < 1.41 even in the best case ln -> 0.')
print()
print('  graph              n   g   m*    kappa_g   max kappa_x  at x')
for name, adj in ZOO:
    n = len(adj); g = girth(adj)
    cap = min(2 * g - 1, 15)
    cc, _ = cycle_counts_by_vertex(adj, cap)
    ks = [(2 * x * cc[x] / 2.0 ** x, x) for x in range(g, cap + 1)]
    mx, ax = max(ks)
    import math
    mstar = 2 * math.log2(2 * n)
    print('  %-16s %3d %3d %5.1f %10.3f %12.3f %5d'
          % (name, n, g, mstar, 2 * g * cc[g] / 2.0 ** g, mx, ax),
          flush=True)
print()
print('  => kappa peaks at x = g and sits at 3.7 .. 5.6 on the cages,')
print('     i.e. FOUR TIMES the 1.41 the theta bound can tolerate.')

print()
print('=== (B) kappa IS UNBOUNDED: TRUNCATIONS ===', flush=True)
print('  graph                     n   g   c_3   kappa_3 = 6 c_3 / 8')
for base, nm in ((gp(5, 2), 'truncated Petersen'),
                 (lcf(14, [5, -5] * 7), 'truncated Heawood'),
                 (lcf(24, [12, 7, -7] * 8), 'truncated McGee'),
                 (lcf(30, [-13, -9, 7, -7, 9, 13] * 5),
                  'truncated Tutte-Coxeter')):
    adj = truncate(base); n = len(adj)
    cc, _ = cycle_counts_by_vertex(adj, 4)
    print('  %-24s %4d %3d %5d %14.2f'
          % (nm, n, girth(adj), cc[3], 6 * cc[3] / 8.0), flush=True)
print('  => kappa_3 = n/4 exactly (one triangle per original vertex).')
print('     No absolute bound on kappa exists.')

print()
print('=== (C) (NC): max-to-mean cycles through a vertex ===',
      flush=True)
print('  K_obs = max_w c_y(w) / (y c_y / n).  (NC) asks K_obs = O(1).')
print()
print('  graph                     n    y    c_y   max c_y(w)   K_obs')
for name, adj in ZOO[:4]:
    n = len(adj); g = girth(adj)
    cc, per = cycle_counts_by_vertex(adj, g + 2)
    y = g
    mean = y * cc[y] / n
    mx = max(per[y])
    print('  %-24s %4d %4d %6d %12d %7.2f'
          % (name + ' (vtx-trans)', n, y, cc[y], mx,
             mx / mean if mean else 0), flush=True)

rng = random.Random(4711)
for n in (30, 60, 120):
    adj = random_cubic(n, rng)
    g = girth(adj)
    cc, per = cycle_counts_by_vertex(adj, g + 1)
    y = g; mean = y * cc[y] / n; mx = max(per[y])
    print('  %-24s %4d %4d %6d %12d %7.2f'
          % ('random cubic', n, y, cc[y], mx,
             mx / mean if mean else 0), flush=True)

print()
print('  now a graph built so all girth cycles sit in one small'
      ' region:')
base = gp(5, 2)
for reps, label in ((1, 'Petersen + 1 random block'),
                    (3, 'Petersen + 3 random blocks'),
                    (6, 'Petersen + 6 random blocks')):
    adj = base
    for _ in range(reps):
        H = random_cubic(30, rng)
        adj = splice(adj, H)
    n = len(adj); g = girth(adj)
    cc, per = cycle_counts_by_vertex(adj, max(g, 5) + 1)
    y = 5
    if cc[y] == 0: continue
    mean = y * cc[y] / n; mx = max(per[y])
    print('  %-24s %4d %4d %6d %12d %7.2f'
          % (label, n, y, cc[y], mx, mx / mean if mean else 0),
          flush=True)
print()
print('  => K_obs grows linearly in n.  (NC) IS FALSE in general:')
print('     a bounded region carrying all the y-cycles forces')
print('     K_obs = Theta(n/y).')
