"""manifold482.py -- THE SIMPLICITY GAP, MEASURED IN THE JACOBIAN.

A metric graph with unit edge lengths IS a tropical curve, and its
Jacobian is the real torus R^g / Lambda where Lambda is the lattice
of integral flows and the polarisation is Q(x) = ||x||^2 = length.
That lattice is exactly 441's flow lattice.  So the manifold-side
restatement of J7 is:

    the principally polarised tropical Jacobian of a trivalent
    tropical curve carries an INDECOMPOSABLE +-1 lattice vector of
    norm a power of two, whenever it carries one of norm divisible
    by 8.

160.4 recorded the gap in this route verbatim: the theta series
counts ALL lattice vectors of a given norm, while cycles are only
the indecomposable +-1 ones.  This file measures that gap for the
first time, at three levels:

    a_N = # lattice vectors of norm N                (theta series)
    u_N = # of those with all entries in {-1,0,1}    (even subgraphs)
    s_N = # of those that are single circuits        (= 2 c_N)

and it checks the one place where the three agree: PROVED below,
min(Lambda) = girth and the minimal vectors are exactly the girth
cycles, so at N = g the gap is closed and only there.

PROOF that min(Lambda) = g.  A nonzero integral flow has support of
minimum degree >= 2, hence containing a cycle, so ||x||^2 >= |supp x|
>= g.  Equality forces every entry to be 0 or +-1 and |supp x| = g,
so supp x has min degree >= 2 with exactly g edges: a single g-cycle.
"""
import sys, itertools
from collections import deque
P = print


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


def edges_of(adj):
    return [(u, v) for u in range(len(adj)) for v in adj[u] if u < v]


def spanning_tree(adj):
    n = len(adj); seen = {0}; par = {0: None}; q = deque([0]); tre = set()
    while q:
        v = q.popleft()
        for w in adj[v]:
            if w not in seen:
                seen.add(w); par[w] = v; tre.add((min(v, w), max(v, w)))
                q.append(w)
    return tre, par


def fundamental_cycles(adj):
    """for each non-tree edge, the +-1 flow of its fundamental cycle."""
    E = edges_of(adj); idx = {e: i for i, e in enumerate(E)}
    tre, par = spanning_tree(adj)
    def path_to_root(v):
        out = []
        while par[v] is not None: out.append((par[v], v)); v = par[v]
        return out
    basis = []
    for (u, v) in E:
        if (u, v) in tre: continue
        x = [0] * len(E)
        x[idx[(u, v)]] = 1                      # oriented u -> v
        pu = path_to_root(u); pv = path_to_root(v)
        # the flow returns v -> u along the tree: v to root, then root to u.
        # path_to_root yields pairs (parent, child).
        for (a, b) in pv:                    # traversed child -> parent
            e = (min(a, b), max(a, b))
            x[idx[e]] += 1 if (b, a) == e else -1
        for (a, b) in pu:                    # traversed parent -> child
            e = (min(a, b), max(a, b))
            x[idx[e]] += 1 if (a, b) == e else -1
        basis.append(x)
    return E, basis


def is_flow(adj, E, x):
    """check conservation at every vertex for the oriented edge list E."""
    n = len(adj); bal = [0] * n
    for i, (u, v) in enumerate(E):
        bal[u] += x[i]; bal[v] -= x[i]
    return all(b == 0 for b in bal)


def circuits(adj):
    """all simple cycles as edge sets, keyed by length."""
    n = len(adj); out = {}
    for s in range(n):
        vis = [False] * n; vis[s] = True; path = [s]
        def dfs(v):
            for w in adj[v]:
                if w < s: continue
                if w == s and len(path) >= 3:
                    out.setdefault(len(path), set()).add(
                        frozenset((min(path[i], path[(i + 1) % len(path)]),
                                   max(path[i], path[(i + 1) % len(path)]))
                                  for i in range(len(path))))
                    continue
                if vis[w]: continue
                vis[w] = True; path.append(w); dfs(w)
                path.pop(); vis[w] = False
        dfs(s); vis[s] = False
    return {L: len(v) for L, v in out.items()}


ZOO = [('K4', lcf(4, [2, -2])), ('K33', k33()),
       ('cube', lcf(8, [3, -3])), ('Petersen', gp(5, 2))]
NMAX = 12

P('=' * 76)
P('THE SIMPLICITY GAP IN THE TROPICAL JACOBIAN', flush=True)
P('=' * 76)
P('  a_N = lattice vectors of norm N   (the theta series coefficient)')
P('  u_N = those with entries in {-1,0,1}   (even subgraphs, oriented)')
P('  s_N = those that are a single circuit  (= 2 c_N, what we need)')
P()
for name, G in ZOO:
    E, basis = fundamental_cycles(G)
    r = len(basis)
    # sanity: every basis vector is a flow
    assert all(is_flow(G, E, b) for b in basis), name
    B = int(NMAX ** 0.5)
    a = {}; u = {}
    rng = range(-B, B + 1)
    for co in itertools.product(rng, repeat=r):
        if all(c == 0 for c in co): continue
        x = [0] * len(E)
        for j, c in enumerate(co):
            if c:
                bj = basis[j]
                for i in range(len(E)):
                    if bj[i]: x[i] += c * bj[i]
        N = sum(t * t for t in x)
        if N == 0 or N > NMAX: continue
        a[N] = a.get(N, 0) + 1
        if all(abs(t) <= 1 for t in x): u[N] = u.get(N, 0) + 1
    cyc = circuits(G)
    g = min(cyc)
    P('  %s  (n = %d, lattice rank g = %d, girth = %d)'
      % (name, len(G), r, g), flush=True)
    P('    %4s %10s %10s %10s %12s %12s'
      % ('N', 'a_N', 'u_N', 's_N', 'a_N/s_N', 'u_N/s_N'), flush=True)
    for N in range(g, NMAX + 1):
        s = 2 * cyc.get(N, 0)
        aN = a.get(N, 0); uN = u.get(N, 0)
        if aN == 0 and s == 0: continue
        f1 = ('%.2f' % (aN / s)) if s else 'inf'
        f2 = ('%.2f' % (uN / s)) if s else 'inf'
        P('    %4d %10d %10d %10d %12s %12s'
          % (N, aN, uN, s, f1, f2), flush=True)
    mn = min(k for k in a)
    P('    min norm = %d, girth = %d : %s ; minimal vectors %d vs 2 c_g = %d'
      % (mn, g, 'AGREE' if mn == g else 'DISAGREE',
         a[mn], 2 * cyc.get(g, 0)), flush=True)
    P()

P('=' * 76)
P('READING', flush=True)
P('=' * 76)
P('  At N = girth the three counts coincide - the lattice sees the')
P('  cycle exactly.  Above the girth a_N and u_N run away from s_N,')
P('  and every lattice-theoretic invariant (theta series, Voronoi')
P('  cell, Torelli data) reads a_N.  That divergence IS 160.4\'s')
P('  recorded gap, now with numbers on it.')
