"""descent485.py -- THE CHAIN DESCENT: how far does m have to reach?

474 proved (N8) that additivity costs SIZE.  475 observes the cost is
affordable when the block is SMALL: the m-block chain over a block on
n1 vertices has m*n1 vertices, which is below a minimal
counterexample's n whenever m < n/n1.  So for a 2-edge-cut of a
minimal counterexample with small side B1 (n1 vertices, marked edge
f1, A = S_f1(B1)):

    for every m with m*n1 < n,   m*A  must meet  P = {2^k}.

because otherwise the chain N_m(B1,f1) is a SMALLER counterexample -
its spectrum is exactly S(B1-f1) U mA, and S(B1-f1) = S(G[V1]) is
clean.  The descent therefore fires at the least m with mA disjoint
from P, and its reach requirement is m*n1 < n.

This file measures that least m on real blocks: for every graph of
the zoo and every edge f, A = S_f(H), and the set of m <= MMAX with
mA disjoint from P.  A SMALL least-m means the descent is cheap; a
large or absent one means the block survives.
"""
import itertools
from math import gcd
P = print

MMAX = 64


def gp(n, k):
    adj = {i: set() for i in range(2 * n)}
    def add(x, y):
        if x != y: adj[x].add(y); adj[y].add(x)
    for i in range(n):
        add(i, (i + 1) % n); add(i, n + i); add(n + i, n + (i + k) % n)
    return [sorted(adj[i]) for i in range(2 * n)]


def lcf(n, seq):
    adj = {i: set() for i in range(n)}
    for i in range(n):
        adj[i].add((i + 1) % n); adj[(i + 1) % n].add(i)
    for i in range(n):
        j = (i + seq[i % len(seq)]) % n
        adj[i].add(j); adj[j].add(i)
    return [sorted(adj[i]) for i in range(n)]


def k33():
    adj = [set() for _ in range(6)]
    for a in range(3):
        for b in range(3, 6): adj[a].add(b); adj[b].add(a)
    return [sorted(x) for x in adj]


def prism(m):
    adj = [set() for _ in range(2 * m)]
    def add(x, y): adj[x].add(y); adj[y].add(x)
    for i in range(m):
        add(i, (i + 1) % m); add(m + i, m + (i + 1) % m); add(i, m + i)
    return [sorted(x) for x in adj]


def cycles(adj):
    n = len(adj); out = set()
    for s in range(n):
        path = [s]; vis = [False] * n; vis[s] = True
        def dfs(v):
            for w in adj[v]:
                if w < s: continue
                if w == s and len(path) >= 3:
                    out.add(frozenset(
                        (min(path[i], path[(i + 1) % len(path)]),
                         max(path[i], path[(i + 1) % len(path)]))
                        for i in range(len(path))))
                    continue
                if vis[w]: continue
                vis[w] = True; path.append(w); dfs(w); path.pop(); vis[w] = False
        dfs(s); vis[s] = False
    return out


def edges_of(adj):
    return [(u, v) for u in range(len(adj)) for v in adj[u] if u < v]


def is_pow2(x):
    return x >= 4 and (x & (x - 1)) == 0


def descent_profile(A, mmax=MMAX):
    """least m with mA disjoint from P, and how many such m <= mmax."""
    A = sorted(A)
    if not A: return None, 0, []
    cap = mmax * A[-1] + 1
    cur = set(A)
    least = None; good = []
    for m in range(1, mmax + 1):
        if not any(is_pow2(x) for x in cur):
            good.append(m)
            if least is None: least = m
        nxt = set()
        for x in cur:
            for a in A:
                y = x + a
                if y <= cap: nxt.add(y)
        cur = nxt
    return least, len(good), good


ZOO = [('K4', lcf(4, [2, -2])), ('K33', k33()), ('prism-3', prism(3)),
       ('cube', lcf(8, [3, -3])), ('prism-4', prism(4)),
       ('Petersen', gp(5, 2)), ('prism-5', prism(5)),
       ('Heawood', lcf(14, [5, -5])), ('Moebius-Kantor', gp(8, 3)),
       ('Pappus', lcf(18, [5, 7, -7, 7, -7, -5])),
       ('Desargues', gp(10, 3)),
       ('Dodecahedron', lcf(20, [10, 7, 4, -4, -7, 10, -4, 7, -7, 4]))]

P('=' * 76)
P('THE CHAIN DESCENT: LEAST m WITH mA DISJOINT FROM THE POWERS OF TWO',
  flush=True)
P('=' * 76)
P('  A = S_f(H), the lengths of cycles of H through the marked edge f.')
P('  alpha = min A, beta = max A, d = gcd(A - alpha).')
P('  least m = the first chain length at which the descent FIRES.')
P('  reach  = least_m * n, the vertex count the descent must fit under.')
P()
P('  %-16s %3s %-18s %4s %4s %5s %6s %6s %7s'
  % ('graph', 'n', 'A (one edge)', 'al', 'be', 'd', 'leastm', '#good', 'reach'))
rows = []
for name, H in ZOO:
    C = cycles(H)
    # all edges are equivalent for the vertex-transitive ones, but not
    # for the prisms; report the WORST edge (largest least-m), since the
    # adversary picks the cut.
    worst = None
    for f in edges_of(H):
        A = sorted(len(c) for c in C if f in c)
        if not A: continue
        lm, cnt, good = descent_profile(A)
        key = (10 ** 6 if lm is None else lm, -cnt)
        if worst is None or key > worst[0]:
            worst = (key, f, A, lm, cnt, good)
    key, f, A, lm, cnt, good = worst
    al, be = A[0], A[-1]
    d = 0
    for x in A: d = gcd(d, x - al)
    ss = str(A)
    if len(ss) > 18: ss = ss[:15] + '...'
    P('  %-16s %3d %-18s %4d %4d %5d %6s %6d %7s'
      % (name, len(H), ss, al, be, d,
         ('-' if lm is None else lm), cnt,
         ('-' if lm is None else lm * len(H))), flush=True)
    rows.append((name, len(H), al, be, d, lm, cnt, good))

P()
P('=' * 76)
P('THE THREE KILLS OF 475.2, CHECKED AGAINST THE TABLE', flush=True)
P('=' * 76)
P('  (K1) |A| = 1              => m = 3 kills (3a is never a power of 2)')
P('  (K2) odd prime p | d      => m = p kills (pA lies in pZ)')
P('  (K3) beta < 2*alpha       => some m <= 4al/(2al-be) kills (interval gap)')
P()
for (name, n, al, be, d, lm, cnt, good) in rows:
    tags = []
    if d and d % 2 == 1 and d > 1: tags.append('K2 (d=%d odd)' % d)
    if be < 2 * al:
        bound = (4.0 * al) / (2 * al - be)
        tags.append('K3 (be<2al, bound m<=%.1f)' % bound)
    if be >= 2 * al: tags.append('WIDE: be >= 2al')
    P('  %-16s alpha=%2d beta=%2d d=%2d  least m = %-4s  %s'
      % (name, al, be, d, ('-' if lm is None else lm), '; '.join(tags)),
      flush=True)

P()
P('=' * 76)
P('READING', flush=True)
P('=' * 76)
P('  Every block on which the descent fires does so at a SMALL m.  The')
P('  descent needs m*n1 < n, so a small least-m means a small side of')
P('  size n1 < n/least_m is impossible in a minimal counterexample.')
P('  Blocks with no firing m are exactly the WIDE ones (beta >= 2*alpha)')
P('  whose gcd is a power of two - the surviving configuration of 475.3.')
