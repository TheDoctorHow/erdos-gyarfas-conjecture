"""couple480.py -- WHICH POWER OF TWO SURVIVES?

couple479 PART 3 measured the avoidance lattice and found: every
PROPER subset of {4, 8, 16} is avoidable by a cubic graph, the full
set never is.  Its final report line was WRONG - it called
count_len(G, L, 1), and count_len returns tot//2, so a cap of 1
always returns 0.  The lattice rows themselves used cap 30 and are
unaffected; only the "still carries" line was junk.  374.1.

This file redoes the decisive row properly and longer: minimise the
total number of power-of-two cycles over cubic graphs, and report
the exact surviving profile (#C4, #C8, #C16, #C32).  WHICH power of
two refuses to die, and how few copies of it remain, is the level
coupling in measured form.
"""
import sys, time, random
P = print
sys.setrecursionlimit(100000)


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
    """exact count of simple L-cycles, stopped once 2*cap raw hits seen."""
    n = len(adj); tot = [0]; lim = 2 * cap
    for s in range(n):
        vis = [False] * n; vis[s] = True
        def dfs(v, ln):
            if tot[0] >= lim: return
            if ln == L:
                if s in adj[v]: tot[0] += 1
                return
            for w in adj[v]:
                if w > s and not vis[w]:
                    vis[w] = True; dfs(w, ln + 1); vis[w] = False
        dfs(s, 1); vis[s] = False
        if tot[0] >= lim: break
    return tot[0] // 2


def girth(adj):
    for L in range(3, len(adj) + 1):
        if has_cycle(adj, L): return L
    return None


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


CAP = 30
rng = random.Random(480480)

P('=' * 78)
P('MINIMISING THE POWER-OF-TWO CYCLE COUNT', flush=True)
P('=' * 78)
P('  objective = #C4 + #C8 + #C16 (+ #C32 when n >= 32), each capped')
P('  at %d.  Zero would be an Erdos-Gyarfas counterexample.' % CAP)
P()
P('  %4s %8s %9s %8s %8s %8s %8s %7s'
  % ('n', 'restarts', 'best obj', '#C4', '#C8', '#C16', '#C32', 'girth'),
  flush=True)
BEST = {}
for n in (16, 20, 24, 28, 30, 32):
    live = [L for L in (4, 8, 16, 32) if L <= n]
    t0 = time.time(); best = 10 ** 9; bestG = None; starts = 0
    budget = 60.0 if n < 32 else 120.0
    while time.time() - t0 < budget and best > 0:
        G = randcubic(n, rng)
        if G is None: break
        starts += 1
        cur = sum(min(count_len(G, L, CAP), CAP) for L in live)
        stall = 0
        while cur > 0 and stall < 120 and time.time() - t0 < budget:
            H = swap(G, rng)
            if H is None: stall += 1; continue
            c = sum(min(count_len(H, L, CAP), CAP) for L in live)
            if c <= cur:
                stall = 0 if c < cur else stall + 1
                G, cur = H, c
            else: stall += 1
        if cur < best: best, bestG = cur, G
    BEST[n] = (best, bestG)
    prof = []
    for L in (4, 8, 16, 32):
        if L > n: prof.append('-')
        else: prof.append(str(count_len(bestG, L)))
    P('  %4d %8d %9d %8s %8s %8s %8s %7s'
      % (n, starts, best, prof[0], prof[1], prof[2], prof[3], girth(bestG)),
      flush=True)

P()
P('=' * 78)
P('VERIFICATION OF EACH BEST GRAPH (existence tests, not counts)',
  flush=True)
P('=' * 78)
P('  %4s %8s %8s %8s %8s   %s'
  % ('n', 'C4?', 'C8?', 'C16?', 'C32?', 'verdict'), flush=True)
ce = []
for n in (16, 20, 24, 28, 30, 32):
    best, G = BEST[n]
    if G is None: continue
    hits = {}
    for L in (4, 8, 16, 32):
        hits[L] = has_cycle(G, L) if L <= n else False
    surv = [L for L in (4, 8, 16, 32) if hits[L]]
    if not surv: ce.append((n, G))
    P('  %4d %8s %8s %8s %8s   %s'
      % (n, hits[4], hits[8], hits[16], hits[32] if n >= 32 else '-',
         ('*** NO POWER OF TWO ***' if not surv
          else 'survives: %s' % surv)), flush=True)

P()
P('=' * 78)
P('VERDICT', flush=True)
P('=' * 78)
if ce:
    P('  *** CANDIDATE COUNTEREXAMPLE - re-verify from stored adjacency ***')
    for (n, G) in ce: P('    n = %d  adj = %s' % (n, G), flush=True)
else:
    P('  Every optimised graph retains a power-of-two cycle.  The')
    P('  profile column says which one and how few copies survive: that')
    P('  minimum is the measured strength of the coupling.  Combined')
    P('  with couple479 PART 3 - every PROPER subset of {4,8,16} is')
    P('  avoidable, the full set never - the levels are pairwise free')
    P('  and jointly constrained.')
    for n in (16, 20, 24, 28, 30, 32):
        best, G = BEST[n]
        if G is None: continue
        P('    n = %2d  best adjacency = %s' % (n, G), flush=True)
