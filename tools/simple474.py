"""simple474.py -- J7 STRIPPED, AND THE DESCENT IT SUGGESTS.

THE SIMPLIFICATION.  Take the counterexample hypothesis seriously and
see where it actually enters.  G cubic, no cycle length a power of
two.  Let t = max{ v_2(|C|) : C a cycle }.  Let C be a shortest cycle
with v_2(|C|) = t, so |C| = 2^t m with m odd.  Now:

    m = 1 would make |C| = 2^t, a power of two.  So m >= 3.
    THAT IS THE ONLY PLACE THE HYPOTHESIS IS USED.

And an ear raising the valuation to t+1 contradicts maximality of t.
So J7 - "no shallow cycle at level t >= 3" - is exactly:

    J7.  There is no cubic graph in which no cycle length is a power
         of two and some cycle length is divisible by 8.

    Equivalently, in positive form:

    (J7')  If a cubic graph has a cycle of length divisible by 8,
           it has a cycle whose length is a power of two.

No islands, no shifts, no shallowness, no lattice.  One line.

THE DESCENT IT SUGGESTS.  The natural proof of (J7') is induction on
the ODD PART, not on t:

    (DESC)  In a cubic graph, the SHORTEST cycle whose length is
            divisible by 2^t is itself a power of two.

(DESC) gives (J7') immediately.  And its induction is on m, which
STRICTLY DECREASES while t STAYS FIXED - so it is untouched by
447.3's Scale-Invariance Theorem, which forbids induction on t.  The
step it needs is: given a cycle of length 2^t m with m odd >= 3,
produce a SHORTER cycle of length divisible by 2^t.  That is exactly
465.4(a)'s same-stratum-plus-shortness kill, whose congruence half
465.3 measured sharp at 2^t.

This file tests (DESC) directly.  Unlike J7 in its island form -
which is vacuous on any graph that HAS a power-of-two cycle, since
the max-valuation cycle can never have a valuation-raising ear -
(DESC) is a statement about every cubic graph and is checkable.
"""
import sys, random, itertools
from collections import deque

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
    A = [sorted(adj[i]) for i in range(n)]
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


def mobius(k):
    adj = {i: set() for i in range(2 * k)}
    for i in range(2 * k):
        adj[i].add((i + 1) % (2 * k)); adj[(i + 1) % (2 * k)].add(i)
    for i in range(k): adj[i].add(i + k); adj[i + k].add(i)
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


def girth(adj):
    n = len(adj); best = 10 ** 9
    for s in range(n):
        d = {s: 0}; p = {s: -1}; q = deque([s])
        while q:
            v = q.popleft()
            for w in adj[v]:
                if w not in d: d[w] = d[v] + 1; p[w] = v; q.append(w)
                elif w != p[v]: best = min(best, d[v] + d[w] + 1)
    return best


def has_cycle(adj, L):
    """is there a simple cycle of length exactly L?"""
    n = len(adj)
    for s in range(n):
        seen = {s}; found = [False]
        def dfs(v, ln):
            if found[0]: return
            if ln == L:
                if s in adj[v]: found[0] = True
                return
            for w in adj[v]:
                if w > s and w not in seen:
                    seen.add(w); dfs(w, ln + 1); seen.discard(w)
        dfs(s, 1)
        if found[0]: return True
    return False


def spectrum(adj, cap=None):
    n = len(adj); cap = cap or n
    return [L for L in range(3, cap + 1) if has_cycle(adj, L)]


def v2(x):
    k = 0
    while x % 2 == 0: x //= 2; k += 1
    return k


def ispow2(x): return x > 0 and (x & (x - 1)) == 0


ZOO = [('K4', lcf(4, [2, -2])), ('K33', lcf(6, [3, 3, 3])),
       ('cube', lcf(8, [3, -3])), ('Petersen', gp(5, 2)),
       ('Moebius-Kantor', gp(8, 3)), ('Heawood', lcf(14, [5, -5])),
       ('Pappus', lcf(18, [5, 7, -7, 7, -7, -5])),
       ('Desargues', gp(10, 3)), ('Nauru', gp(12, 5)),
       ('Dodecahedron', lcf(20, [10, 7, 4, -4, -7, 10, -4, 7, -7, 4])),
       ('trunc-K4', truncate(lcf(4, [2, -2]))),
       ('trunc-cube', truncate(lcf(8, [3, -3]))),
       ('trunc-Petersen', truncate(gp(5, 2))),
       ('Moebius ladder 16', mobius(8)),
       ('Moebius ladder 24', mobius(12)),
       ('prism 18', prism(9)), ('prism 24', prism(12)),
       ('GP(7,2)', gp(7, 2)), ('GP(9,2)', gp(9, 2)),
       ('GP(11,3)', gp(11, 3)), ('GP(13,5)', gp(13, 5))]
rng = random.Random(31337)
for n in (16, 20, 24, 26, 28):
    G = random_cubic(n, rng, 3)
    if G: ZOO.append(('random n=%d' % n, G))
for n in (20, 24, 28):
    G = random_cubic(n, rng, 5)
    if G: ZOO.append(('random g>=5 n=%d' % n, G))

P('=' * 76)
P('(DESC): IS THE SHORTEST CYCLE DIVISIBLE BY 2^t A POWER OF TWO?', flush=True)
P('=' * 76)
P('  For each graph and each t >= 3, find the least cycle length')
P('  divisible by 2^t.  (DESC) says it is always a power of two.')
P()
P('  %-20s %4s %4s %-26s %s'
  % ('graph', 'n', 'g', 'spectrum', 't=3 / t=4 / t=5'), flush=True)
viol = []
for name, G in ZOO:
    if G is None: continue
    n = len(G)
    sp = spectrum(G)
    cells = []
    for t in (3, 4, 5):
        M = 1 << t
        div = [L for L in sp if L % M == 0]
        if not div: cells.append('-')
        else:
            s = min(div)
            cells.append('%d%s' % (s, '' if ispow2(s) else '  <-- NOT 2^k'))
            if not ispow2(s): viol.append((name, t, s, sp))
    sps = str(sp)
    if len(sps) > 26: sps = sps[:23] + '...'
    P('  %-20s %4d %4d %-26s %s'
      % (name, n, girth(G), sps, ' / '.join(cells)), flush=True)

P()
P('=' * 76)
P('VERDICT', flush=True)
P('=' * 76)
if not viol:
    P('  (DESC) holds on every graph tested.')
else:
    P('  (DESC) FAILS on %d rows:' % len(viol))
    for (nm, t, s, sp) in viol:
        P('    %-20s t=%d: least multiple of %d is %d, spectrum %s'
          % (nm, t, 1 << t, s, sp), flush=True)
    P()
    P('  A failure of (DESC) does NOT refute J7: J7 needs the graph to')
    P('  have NO power of two anywhere, and each row above may still')
    P('  carry one at another length.  Checking that:')
    for (nm, t, s, sp) in viol:
        pw = [L for L in sp if ispow2(L)]
        P('    %-20s powers of two in the spectrum: %s'
          % (nm, pw if pw else 'NONE  *** E-G COUNTEREXAMPLE ***'), flush=True)
