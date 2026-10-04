"""manifold483.py -- THE MINIMALITY RATE IS NOT A CONSTANT FACTOR.

161.1 chose the cycle-code translation over the Outer-space one with
this argument: among closed geodesics, cycles are EXPONENTIALLY rare
(2^{n/2} cycles against 2^n/n geodesics), whereas in the cycle code
"the same filter costs a CONSTANT FACTOR (160.1: cycles are 57-90%
of Z(G))".

The comparison is sound.  The constant-factor claim is an AGGREGATE
over all weights, and J7 needs one weight, 2^k.  This file computes
the minimality rate WEIGHT BY WEIGHT: for each w, the number of
codewords of the cycle code of weight w, and how many of them are
single circuits (minimal codewords).

If some weight has codewords but no minimal ones, then no theorem of
the form "the code has a codeword of weight 2^k" can yield a cycle,
and the aggregate rate is the wrong statistic for this application.
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


def cycle_space_basis(adj):
    """fundamental cycles as edge bitmasks over F_2."""
    E = [(u, v) for u in range(len(adj)) for v in adj[u] if u < v]
    idx = {e: i for i, e in enumerate(E)}
    seen = {0}; par = {0: None}; q = deque([0]); tre = set()
    while q:
        v = q.popleft()
        for w in adj[v]:
            if w not in seen:
                seen.add(w); par[w] = v; tre.add((min(v, w), max(v, w)))
                q.append(w)
    def up(v):
        out = []
        while par[v] is not None:
            a, b = par[v], v
            out.append((min(a, b), max(a, b))); v = a
        return out
    basis = []
    for (u, v) in E:
        if (u, v) in tre: continue
        m = 1 << idx[(u, v)]
        for e in up(u): m ^= 1 << idx[e]
        for e in up(v): m ^= 1 << idx[e]
        basis.append(m)
    return E, basis


def is_single_cycle(E, mask):
    """support has all degrees 2 and is connected."""
    es = [E[i] for i in range(len(E)) if (mask >> i) & 1]
    if not es: return False
    deg = {}
    for (u, v) in es:
        deg[u] = deg.get(u, 0) + 1; deg[v] = deg.get(v, 0) + 1
    if any(d != 2 for d in deg.values()): return False
    nb = {v: [] for v in deg}
    for (u, v) in es: nb[u].append(v); nb[v].append(u)
    s = next(iter(nb)); seen = {s}; st = [s]
    while st:
        v = st.pop()
        for w in nb[v]:
            if w not in seen: seen.add(w); st.append(w)
    return len(seen) == len(deg)


ZOO = [('K4', lcf(4, [2, -2])), ('cube', lcf(8, [3, -3])),
       ('Petersen', gp(5, 2)), ('Heawood', lcf(14, [5, -5])),
       ('Moebius-Kantor', gp(8, 3)),
       ('Pappus', lcf(18, [5, 7, -7, 7, -7, -5])),
       ('Desargues', gp(10, 3)),
       ('Dodecahedron', lcf(20, [10, 7, 4, -4, -7, 10, -4, 7, -7, 4])),
       ('McGee', lcf(24, [12, 7, -7])),
       ('trunc-Petersen', truncate(gp(5, 2)))]

P('=' * 78)
P('MINIMALITY RATE OF THE CYCLE CODE, WEIGHT BY WEIGHT', flush=True)
P('=' * 78)
P('  A_w = codewords of weight w ; M_w = those that are a single')
P('  circuit.  "dead" weights are those with A_w > 0 and M_w = 0.')
P()
DEAD = []
for name, G in ZOO:
    E, basis = cycle_space_basis(G)
    r = len(basis)
    A = {}; M = {}
    for co in range(1, 1 << r):
        m = 0; c = co; j = 0
        while c:
            if c & 1: m ^= basis[j]
            c >>= 1; j += 1
        w = bin(m).count('1')
        A[w] = A.get(w, 0) + 1
        if is_single_cycle(E, m): M[w] = M.get(w, 0) + 1
    tot = sum(A.values()); mtot = sum(M.values())
    dead = sorted(w for w in A if M.get(w, 0) == 0)
    pw = [w for w in dead if w > 0 and (w & (w - 1)) == 0]
    P('  %-16s n=%2d rank=%2d   codewords %6d   minimal %6d  (%.0f%%)'
      % (name, len(G), r, tot, mtot, 100.0 * mtot / tot), flush=True)
    P('      dead weights (A_w > 0, M_w = 0): %s'
      % (dead if dead else 'none'), flush=True)
    if pw:
        P('      *** DEAD AT A POWER OF TWO: %s ***' % pw, flush=True)
    for w in sorted(A):
        if M.get(w, 0) == 0 and A[w] > 0:
            DEAD.append((name, w, A[w]))
    P()

P('=' * 78)
P('VERDICT', flush=True)
P('=' * 78)
if DEAD:
    P('  The weight-by-weight minimality rate reaches ZERO.  Examples:')
    for (nm, w, a) in DEAD[:14]:
        P('    %-16s weight %2d : %d codewords, none minimal'
          % (nm, w, a), flush=True)
    P()
    P('  So 160.1\'s 57-90%% is an AGGREGATE and does not transfer to a')
    P('  single weight.  A theorem producing a codeword of weight 2^k')
    P('  need not produce a CYCLE of length 2^k: the filter at one')
    P('  weight can be total.  161.1\'s comparison with Outer space')
    P('  still stands - exponential versus non-uniform is still a real')
    P('  gain - but the phrase "costs a constant factor" is wrong for')
    P('  the single-weight use J7 makes of it.')
else:
    P('  No dead weight found; the constant-factor reading survives.')
