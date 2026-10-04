"""base477.py -- THE BASE-CASE BOARD, WITH THE 24-CYCLE PINNED.

base476's search collapsed PAST the target: driving #C8 and #C16 to
zero it also destroyed the 24-cycle, landing on cubic graphs with NO
cycle length divisible by 8 at all.  Those are t <= 2 graphs and are
irrelevant to J7 - but they are worth one check of their own, below.

Here the 24-cycle is a HARD CONSTRAINT: every accepted swap must
preserve it.  The objective is #C8 + #C16.  A graph reaching zero is
a cubic graph whose shortest cycle divisible by 8 has length 24, i.e.
q_3 = 3: the base case of the m-descent, inhabited.  That would make
(DESC_3) false and kill the descent route at the level J7 needs.

PART 2 checks the t <= 2 escapees: a cubic graph with no 4-cycle AND
no cycle divisible by 8 has no power-of-two cycle at all, so it is an
Erdos-Gyarfas counterexample.  445.3 proves none exists.  Every
8-divisible-free graph found must therefore carry a 4-cycle; that is
a live consistency test of the archive's t <= 2 theorems.
"""
import sys, time, random
from collections import deque

sys.setrecursionlimit(100000)
P = print


def connected(adj):
    n = len(adj); seen = {0}; s = [0]
    while s:
        v = s.pop()
        for w in adj[v]:
            if w not in seen: seen.add(w); s.append(w)
    return len(seen) == n


def has_cycle(adj, L):
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


def count_len(adj, L, cap=10 ** 9):
    n = len(adj); tot = [0]
    for s in range(n):
        vis = [False] * n; vis[s] = True
        def dfs(v, ln):
            if tot[0] >= cap: return
            if ln == L:
                if s in adj[v]: tot[0] += 1
                return
            for w in adj[v]:
                if w > s and not vis[w]:
                    vis[w] = True; dfs(w, ln + 1); vis[w] = False
        dfs(s, 1); vis[s] = False
        if tot[0] >= cap: break
    return tot[0] // 2


def randcubic(n, rng, tries=60000):
    for _ in range(tries):
        st = [v for v in range(n) for _ in range(3)]; rng.shuffle(st)
        adj = [set() for _ in range(n)]; ok = True
        for i in range(0, len(st), 2):
            a, b = st[i], st[i + 1]
            if a == b or b in adj[a]: ok = False; break
            adj[a].add(b); adj[b].add(a)
        if not ok: continue
        A = [sorted(x) for x in adj]
        if connected(A): return A
    return None


def swap(adj, rng):
    E = [(u, w) for u in range(len(adj)) for w in adj[u] if u < w]
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
    return A if connected(A) else None


CAP = 300


def obj(G):
    return count_len(G, 8, CAP) + count_len(G, 16, CAP)


P('=' * 76)
P('PART 1.  BASE-CASE BOARD: no C8, no C16, C24 PINNED', flush=True)
P('=' * 76)
P('  objective = #C8 + #C16 (each capped at %d); every accepted swap' % CAP)
P('  must keep a 24-cycle.  Zero objective = q_3 = 3 = base case.')
P()
P('  %4s %8s %10s %10s   %s'
  % ('n', 'starts', 'best obj', 'seconds', 'verdict'), flush=True)
rng = random.Random(477477)
WIN = []
for n in (24, 26, 28, 30):
    t0 = time.time(); best = 10 ** 9; bestG = None; starts = 0
    while time.time() - t0 < 90.0 and best > 0:
        G = randcubic(n, rng)
        if G is None or not has_cycle(G, 24): continue
        starts += 1
        cur = obj(G); stall = 0
        while cur > 0 and stall < 120 and time.time() - t0 < 90.0:
            H = swap(G, rng)
            if H is None: stall += 1; continue
            if not has_cycle(H, 24): stall += 1; continue
            c = obj(H)
            if c <= cur:
                stall = 0 if c < cur else stall + 1
                G, cur = H, c
            else: stall += 1
        if cur < best: best, bestG = cur, G
    v = 'no'
    if best == 0 and bestG is not None:
        a8 = has_cycle(bestG, 8); a16 = has_cycle(bestG, 16)
        a24 = has_cycle(bestG, 24)
        if (not a8) and (not a16) and a24:
            v = 'BASE CASE FOUND'; WIN.append((n, bestG))
        else: v = 'stale (C8=%s C16=%s C24=%s)' % (a8, a16, a24)
    P('  %4d %8d %10s %10.1f   %s'
      % (n, starts, best if best < 10 ** 9 else '-', time.time() - t0, v),
      flush=True)

P()
if WIN:
    P('  WITNESSES, re-verified from stored adjacency:', flush=True)
    for (n, G) in WIN:
        P('    n = %d' % n, flush=True)
        P('      adj = %s' % G, flush=True)
        P('      C4=%s C8=%s C16=%s C24=%s C32=%s'
          % (has_cycle(G, 4), has_cycle(G, 8), has_cycle(G, 16),
             has_cycle(G, 24), has_cycle(G, 32) if n >= 32 else 'n/a'),
          flush=True)
        mult8 = [L for L in range(8, n + 1, 8) if has_cycle(G, L)]
        P('      cycle lengths divisible by 8: %s  =>  q_3 = %s'
          % (mult8, (min(mult8) // 8) if mult8 else 'stratum empty'),
          flush=True)

P()
P('=' * 76)
P('PART 2.  THE t <= 2 ESCAPEES MUST CARRY A 4-CYCLE', flush=True)
P('=' * 76)
P('  base476 produced cubic graphs with NO cycle divisible by 8.  If')
P('  such a graph also had no 4-cycle it would have no power-of-two')
P('  cycle at all - an Erdos-Gyarfas counterexample.  445.3 forbids')
P('  it.  Test: hunt 8-divisible-free cubic graphs and read off C4.')
P()
P('  %4s %8s %14s %10s   %s'
  % ('n', 'found', 'all have C4?', 'girth>=5?', 'note'), flush=True)
for n in (20, 22, 24, 26, 28):
    t0 = time.time(); got = 0; allc4 = True; g5 = 0
    while time.time() - t0 < 40.0 and got < 40:
        G = randcubic(n, rng)
        if G is None: break
        cur = sum(count_len(G, L, CAP) for L in range(8, n + 1, 8))
        stall = 0
        while cur > 0 and stall < 100 and time.time() - t0 < 40.0:
            H = swap(G, rng)
            if H is None: stall += 1; continue
            c = sum(count_len(H, L, CAP) for L in range(8, n + 1, 8))
            if c <= cur:
                stall = 0 if c < cur else stall + 1
                G, cur = H, c
            else: stall += 1
        if cur != 0: continue
        got += 1
        if not has_cycle(G, 4):
            allc4 = False
            P('    *** n=%d 8-divisible-free AND no 4-cycle: %s' % (n, G),
              flush=True)
        else:
            if not has_cycle(G, 3): g5 += 1
    P('  %4d %8d %14s %10s   %s'
      % (n, got, 'yes' if allc4 else 'NO  *** SEE ABOVE',
         '%d/%d' % (g5, got) if got else '-',
         'consistent with 445.3' if allc4 else 'CONTRADICTS 445.3'),
      flush=True)

P()
P('=' * 76)
P('VERDICT', flush=True)
P('=' * 76)
if WIN:
    P('  (DESC_3) IS FALSE: the base case q_3 = 3 is inhabited by an')
    P('  explicit cubic graph.  The descent route dies at t = 3, and')
    P('  the m-frame retains only its structural content.')
else:
    P('  Base case q_3 = 3 not realised under a pinned 24-cycle either.')
    P('  Evidence, not proof.')
