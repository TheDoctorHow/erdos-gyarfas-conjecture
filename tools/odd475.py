"""odd475.py -- THE ODD PART AS THE TRANSVERSE COORDINATE.

447.3 proved that halving l |-> l/2 carries the level-t problem to
the level-(t-1) problem and carries the target P = {2^k} with it, so
induction on t is circular.  Read that as a group action.  Write
every length uniquely as

        l  =  2^t . m ,   m odd ,

so halving is the map H : (t, m) |-> (t-1, m).  It moves t and FIXES
m.  Its orbits are the sets {m . 2^s}, one per odd m, and

        P = {2^k}  IS EXACTLY THE FIBRE  m = 1 .

So the target of Erdos-Gyarfas is one fibre of the odd-part map, and
447.3's symmetry acts WITHIN fibres.  Moving along an orbit can never
change which fibre you are in - which is precisely why induction on t
cannot approach the target.  The odd part m is the unique coordinate
transverse to the symmetry, and hence the only induction parameter
that can see P at all.

447.4's integrality obstruction ("the t = 2 machinery needs shifts in
Z; halving puts them in (1/2)Z") also never fires here, because an
m-descent never halves: t, the modulus 2^{t+1}, and the shift lattice
Z are all FIXED throughout.

WHAT THIS FILE MEASURES.

  PART 1  the (t, m) grid of the cycle spectrum for a zoo of cubic
          graphs: q_t = min{ l/2^t : 2^t | l } and m_t = min{ odd part
          of l : v_2(l) = t }.
  PART 2  for which t the descent statement "q_t is a power of two"
          is true.  It is FALSE at t = 1 (Petersen: q_1 = 3).  The
          boundary is the point.
  PART 3  the base-case collapse at m = 3, checked as arithmetic, and
          how often the zoo actually sits at m = 3.
  PART 4  the base-case witness hunt: is there a cubic graph with no
          8-cycle, no 16-cycle, and a 24-cycle?  Such a graph has
          q_3 = 3 and would inhabit the base case.
"""
import sys, time, random, itertools
from collections import deque

sys.setrecursionlimit(100000)
P = print


# ---------------------------------------------------------------- graphs
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


def randcubic(n, rng, mingirth=3, tries=60000):
    for _ in range(tries):
        st = [v for v in range(n) for _ in range(3)]; rng.shuffle(st)
        adj = [set() for _ in range(n)]; ok = True
        for i in range(0, len(st), 2):
            a, b = st[i], st[i + 1]
            if a == b or b in adj[a]: ok = False; break
            adj[a].add(b); adj[b].add(a)
        if not ok: continue
        A = [sorted(x) for x in adj]
        if not connected(A): continue
        if girth(A) < mingirth: continue
        return A
    return None


def connected(adj):
    n = len(adj); seen = {0}; s = [0]
    while s:
        v = s.pop()
        for w in adj[v]:
            if w not in seen: seen.add(w); s.append(w)
    return len(seen) == n


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


# ------------------------------------------------------------- spectrum
def has_cycle(adj, L):
    """exact: is there a simple cycle of length exactly L?"""
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


def spectrum(adj, seconds=40.0):
    """all cycle lengths.  returns (set, complete?)."""
    n = len(adj); found = set(); t0 = time.time()
    state = {'stop': False}
    for s in range(n):
        vis = [False] * n; vis[s] = True
        def dfs(v, ln):
            if state['stop']: return
            if time.time() - t0 > seconds: state['stop'] = True; return
            for w in adj[v]:
                if w < s: continue
                if w == s:
                    if ln >= 3: found.add(ln)
                    continue
                if vis[w]: continue
                vis[w] = True; dfs(w, ln + 1); vis[w] = False
        dfs(s, 1); vis[s] = False
        if state['stop']: break
    return found, not state['stop']


def v2(x):
    k = 0
    while x % 2 == 0: x //= 2; k += 1
    return k


def ispow2(x): return x > 0 and (x & (x - 1)) == 0


ZOO = [('K4', lcf(4, [2, -2])), ('K33', lcf(6, [3, 3, 3])),
       ('cube', lcf(8, [3, -3])), ('Petersen', gp(5, 2)),
       ('Moebius-Kantor', gp(8, 3)), ('Heawood', lcf(14, [5, -5])),
       ('Pappus', lcf(18, [5, 7, -7, 7, -7, -5])),
       ('Desargues', gp(10, 3)),
       ('Dodecahedron', lcf(20, [10, 7, 4, -4, -7, 10, -4, 7, -7, 4])),
       ('Nauru', gp(12, 5)),
       ('McGee', lcf(24, [12, 7, -7])),
       ('trunc-K4', truncate(lcf(4, [2, -2]))),
       ('trunc-cube', truncate(lcf(8, [3, -3]))),
       ('trunc-Petersen', truncate(gp(5, 2))),
       ('Moebius ladder 16', mobius(8)),
       ('Moebius ladder 24', mobius(12)),
       ('prism 18', prism(9)), ('prism 24', prism(12)),
       ('GP(7,2)', gp(7, 2)), ('GP(9,2)', gp(9, 2)),
       ('GP(11,3)', gp(11, 3)), ('GP(13,5)', gp(13, 5))]
rng = random.Random(475475)
for nn in (16, 20, 24, 26):
    G = randcubic(nn, rng, 3)
    if G: ZOO.append(('random n=%d' % nn, G))
for nn in (20, 24, 28):
    G = randcubic(nn, rng, 5)
    if G: ZOO.append(('random g>=5 n=%d' % nn, G))
ZOO = [(nm, G) for (nm, G) in ZOO if G is not None]

P('=' * 78)
P('PART 1.  THE (t, m) GRID', flush=True)
P('=' * 78)
P('  l = 2^t . m, m odd.  q_t = min{ l / 2^t : 2^t divides l }   (the')
P('  divisibility descent target).  m_t = min{ m : v_2(l) = t exactly }')
P('  (the valuation-preserving descent target).  "." = stratum empty.')
P()
P('  %-19s %3s %2s | %-17s | %-17s'
  % ('graph', 'n', 'g', 'q_t, t=1..5', 'm_t, t=1..5'), flush=True)
DATA = []
for name, G in ZOO:
    sp, done = spectrum(G, 40.0)
    if not done:
        P('  %-19s  SPECTRUM INCOMPLETE - excluded' % name, flush=True)
        continue
    qs, ms = [], []
    for t in range(1, 6):
        M = 1 << t
        div = [l for l in sp if l % M == 0]
        ex = [l for l in sp if v2(l) == t]
        qs.append(min(div) // M if div else None)
        ms.append(min(ex) >> t if ex else None)
    DATA.append((name, G, sp, qs, ms))
    f = lambda xs: ' '.join(('%3s' % ('.' if x is None else x)) for x in xs)
    P('  %-19s %3d %2d | %-17s | %-17s'
      % (name, len(G), girth(G), f(qs), f(ms)), flush=True)

P()
P('=' * 78)
P('PART 2.  WHERE THE DESCENT STATEMENT IS TRUE', flush=True)
P('=' * 78)
P('  (DESC_t): the shortest cycle divisible by 2^t is a power of two,')
P('            i.e. q_t is a power of two.')
P('  (VAL_t) : m_t is 1 or 3, i.e. the shortest cycle of valuation')
P('            exactly t has odd part at most 3.')
P()
P('  %3s %10s %10s   %s' % ('t', 'DESC ok/n', 'VAL ok/n', 'counterexamples'),
  flush=True)
for t in range(1, 6):
    dok = dn = vok = vn = 0; bad = []; badv = []
    for (name, G, sp, qs, ms) in DATA:
        q = qs[t - 1]; m = ms[t - 1]
        if q is not None:
            dn += 1
            if ispow2(q): dok += 1
            else: bad.append('%s(q=%d)' % (name, q))
        if m is not None:
            vn += 1
            if m in (1, 3): vok += 1
            else: badv.append('%s(m=%d)' % (name, m))
    P('  %3d %10s %10s   DESC: %s' % (t, '%d/%d' % (dok, dn),
      '%d/%d' % (vok, vn), ', '.join(bad) if bad else 'none'), flush=True)
    if badv:
        P('  %3s %10s %10s   VAL : %s' % ('', '', '', ', '.join(badv)),
          flush=True)

P()
P('=' * 78)
P('PART 3.  THE BASE-CASE COLLAPSE AT m = 3', flush=True)
P('=' * 78)
P('  CLAIM.  For every t >= 0 the multiples of 2^t strictly below')
P('  3 . 2^t are exactly 2^t and 2 . 2^t, both powers of two.  Hence at')
P('  q = 3 the descent target set IS the target set: any strictly')
P('  shorter cycle divisible by 2^t is already a power of two.')
ok = True
for t in range(0, 20):
    below = [j * (1 << t) for j in range(1, 3)]
    if [x for x in below if not ispow2(x)]: ok = False
P('  checked t = 0..19: %s' % ('holds' if ok else 'FAILS'), flush=True)
P()
P('  CONTRAST at q = 5: targets are 2^t, 2.2^t, 3.2^t, 4.2^t; the')
P('  value 3 . 2^t is NOT a power of two, so shortening from q = 5 is')
P('  strictly weaker than winning.  The induction on m is therefore')
P('  easy above the base and equivalent to the problem AT the base.')
P()
P('  How often does the zoo actually sit at q_t = 3, for t >= 3?')
c3 = c5 = cp = 0
for (name, G, sp, qs, ms) in DATA:
    for t in (3, 4, 5):
        q = qs[t - 1]
        if q is None: continue
        if ispow2(q): cp += 1
        elif q == 3: c3 += 1
        else: c5 += 1
P('    q_t a power of two: %d      q_t = 3: %d      q_t >= 5: %d'
  % (cp, c3, c5), flush=True)

P()
P('=' * 78)
P('PART 4.  BASE-CASE WITNESS HUNT', flush=True)
P('=' * 78)
P('  Wanted: a cubic graph with NO 8-cycle, NO 16-cycle, and a')
P('  24-cycle.  Such a graph has q_3 = 3: it inhabits the base case')
P('  and refutes (DESC_3).  Search: 2-edge-swap descent on the number')
P('  of 8-cycles, then exact tests for 16 and 24.')
P()


def count_len(adj, L, cap=10 ** 9):
    """number of simple cycles of length exactly L (each counted once)."""
    n = len(adj); tot = 0
    for s in range(n):
        vis = [False] * n; vis[s] = True
        st = [(s, 1)]
        def dfs(v, ln):
            nonlocal tot
            if tot >= cap: return
            if ln == L:
                if s in adj[v]: tot += 1
                return
            for w in adj[v]:
                if w > s and not vis[w]:
                    vis[w] = True; dfs(w, ln + 1); vis[w] = False
        dfs(s, 1); vis[s] = False
    return tot // 2


def edges(adj):
    return [(u, w) for u in range(len(adj)) for w in adj[u] if u < w]


def swap(adj, rng):
    """random 2-edge swap preserving cubicity; returns new adj or None."""
    E = edges(adj)
    (a, b) = E[rng.randrange(len(E))]
    (c, d) = E[rng.randrange(len(E))]
    if len({a, b, c, d}) < 4: return None
    if rng.random() < 0.5: c, d = d, c
    if d in adj[a] or c in adj[b]: return None
    new = [set(x) for x in adj]
    new[a].discard(b); new[b].discard(a); new[c].discard(d); new[d].discard(c)
    new[a].add(d); new[d].add(a); new[b].add(c); new[c].add(b)
    A = [sorted(x) for x in new]
    if not all(len(x) == 3 for x in A): return None
    if not connected(A): return None
    return A


P('  %5s %8s %10s %10s   %s'
  % ('n', 'restarts', 'best #C8', 'C8 = 0 ?', 'verdict'), flush=True)
hits = []
for n in (24, 26, 28, 30):
    t0 = time.time(); best = 10 ** 9; bestG = None; restarts = 0
    while time.time() - t0 < 25.0:
        G = randcubic(n, rng, 3)
        if G is None: break
        restarts += 1
        cur = count_len(G, 8); stall = 0
        while stall < 60 and cur > 0 and time.time() - t0 < 25.0:
            H = swap(G, rng)
            if H is None: stall += 1; continue
            c = count_len(H, 8)
            if c <= cur:
                if c < cur: stall = 0
                else: stall += 1
                G, cur = H, c
            else: stall += 1
        if cur < best: best, bestG = cur, G
        if best == 0: break
    verdict = '-'
    if best == 0 and bestG is not None:
        h16 = has_cycle(bestG, 16); h24 = has_cycle(bestG, 24)
        verdict = ('16:%s 24:%s' % ('yes' if h16 else 'NO',
                                    'yes' if h24 else 'NO'))
        if (not h16) and h24: hits.append((n, bestG))
    P('  %5d %8d %10d %10s   %s'
      % (n, restarts, best, 'yes' if best == 0 else 'no', verdict), flush=True)

P()
for nm in ('trunc-K4', 'trunc-cube', 'trunc-Petersen', 'McGee'):
    for (name, G, sp, qs, ms) in DATA:
        if name != nm: continue
        P('  %-16s spectrum multiples of 8: %s'
          % (nm, sorted(l for l in sp if l % 8 == 0)), flush=True)

P()
P('=' * 78)
P('VERDICT', flush=True)
P('=' * 78)
if hits:
    P('  BASE CASE INHABITED: %d witness(es) with no 8-cycle, no' % len(hits))
    P('  16-cycle, and a 24-cycle.  (DESC_3) is FALSE.')
    for (n, G) in hits:
        P('    n = %d  adj = %s' % (n, G), flush=True)
        P('    re-verified: C8 %s, C16 %s, C24 %s'
          % (has_cycle(G, 8), has_cycle(G, 16), has_cycle(G, 24)), flush=True)
else:
    P('  No base-case witness found.  Every 8-cycle-free cubic graph')
    P('  reached carries a 16-cycle, so its q_3 is 2, not 3.  The base')
    P('  case m = 3 was NOT realised at t = 3 within this search.')
