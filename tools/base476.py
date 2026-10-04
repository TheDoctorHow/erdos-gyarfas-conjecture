"""base476.py -- TWO DECISIVE CHECKS ON THE m-DESCENT.

CHECK A.  Is the descent step false at t = 1?  odd475 reports
trunc-Petersen with q_1 = 5: min even cycle 10, no 4-cycle, no
6-cycle, no 8-cycle.  If exact, then a cycle of length 2 . 5 has NO
shorter cycle divisible by 2, so (STEP) - "a cycle of length 2^t m,
m odd >= 3, has a strictly shorter companion divisible by 2^t" - is
FALSE at t = 1.  Verified here by targeted exact tests, not read off
a spectrum enumeration.

CHECK B.  Is the base case q_3 = 3 inhabited?  Wanted: a cubic graph
with no 8-cycle, no 16-cycle, and a 24-cycle.  Two-stage search:
anneal to 8-cycle-free, then anneal the (capped) 16-cycle count.
This is 451.2(C)'s m = 3 board in testable form.  A witness kills
the descent route at t = 3; persistent failure is evidence for it.
"""
import sys, time, random, itertools
from collections import deque

sys.setrecursionlimit(100000)
P = print


def gp(n, k):
    adj = {i: set() for i in range(2 * n)}
    def add(a, b):
        if a != b: adj[a].add(b); adj[b].add(a)
    for i in range(n):
        add(i, (i + 1) % n); add(i, n + i); add(n + i, n + (i + k) % n)
    return [sorted(adj[i]) for i in range(2 * n)]


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


P('=' * 76)
P('CHECK A.  IS THE DESCENT STEP FALSE AT t = 1?', flush=True)
P('=' * 76)
TP = truncate(gp(5, 2))
P('  truncated Petersen, n = %d' % len(TP))
for L in (3, 4, 5, 6, 7, 8, 9, 10, 12, 16, 24):
    P('    cycle of length %2d : %s' % (L, has_cycle(TP, L)), flush=True)
ev = [L for L in (4, 6, 8, 10) if has_cycle(TP, L)]
P()
P('  even cycle lengths among {4,6,8,10}: %s' % ev, flush=True)
if ev and min(ev) == 10:
    P('  => min even cycle is 10 = 2 . 5, q_1 = 5, and NO shorter even')
    P('     cycle exists.  (STEP) IS FALSE AT t = 1, and so is (DESC_1).')
else:
    P('  => claim NOT reproduced; odd475 reading was wrong.')

P()
P('=' * 76)
P('CHECK B.  IS THE BASE CASE q_3 = 3 INHABITED?', flush=True)
P('=' * 76)
P('  target: no 8-cycle, no 16-cycle, some 24-cycle.')
P()
P('  %4s %9s %9s %11s %9s   %s'
  % ('n', 'graphs', 'C8-free', 'best #C16', 'C24', 'witness'), flush=True)
rng = random.Random(476476)
CAP = 400
found = []
for n in (24, 26, 28, 30, 32):
    t0 = time.time(); tried = 0; c8free = 0; best16 = 10 ** 9; bestG = None
    while time.time() - t0 < 60.0:
        G = randcubic(n, rng)
        if G is None: break
        tried += 1
        # stage 1: drive #C8 to zero
        cur = count_len(G, 8); stall = 0
        while cur > 0 and stall < 80 and time.time() - t0 < 60.0:
            H = swap(G, rng)
            if H is None: stall += 1; continue
            c = count_len(H, 8)
            if c <= cur:
                stall = 0 if c < cur else stall + 1
                G, cur = H, c
            else: stall += 1
        if cur > 0: continue
        c8free += 1
        # stage 2: among C8-free graphs, drive the capped #C16 down
        c16 = count_len(G, 16, CAP); stall = 0
        while c16 > 0 and stall < 60 and time.time() - t0 < 60.0:
            H = swap(G, rng)
            if H is None: stall += 1; continue
            if count_len(H, 8) != 0: stall += 1; continue
            c = count_len(H, 16, CAP)
            if c <= c16:
                stall = 0 if c < c16 else stall + 1
                G, c16 = H, c
            else: stall += 1
        if c16 < best16: best16, bestG = c16, G
        if best16 == 0: break
    w = '-'
    c24 = '-'
    if bestG is not None and best16 == 0:
        a8 = has_cycle(bestG, 8); a16 = has_cycle(bestG, 16)
        a24 = has_cycle(bestG, 24)
        c24 = 'yes' if a24 else 'no'
        if (not a8) and (not a16) and a24:
            w = 'FOUND'; found.append((n, bestG))
        else: w = 'no (C8=%s C16=%s)' % (a8, a16)
    P('  %4d %9d %9d %11s %9s   %s'
      % (n, tried, c8free, best16 if best16 < 10 ** 9 else '-', c24, w),
      flush=True)

P()
P('=' * 76)
P('VERDICT', flush=True)
P('=' * 76)
if found:
    P('  BASE CASE INHABITED.  (DESC_3) is FALSE.')
    for (n, G) in found:
        P('    n = %d, adj = %s' % (n, G), flush=True)
        P('    re-verified from stored adjacency: C8=%s C16=%s C24=%s'
          % (has_cycle(G, 8), has_cycle(G, 16), has_cycle(G, 24)), flush=True)
else:
    P('  No cubic graph found with no 8-cycle and no 16-cycle.  Every')
    P('  8-cycle-free graph reached still carries a 16-cycle, so q_3 = 2')
    P('  and the base case is not realised.  This is EVIDENCE, not a')
    P('  proof: the search is local and the counts above say how hard')
    P('  the 16-cycle count was pushed.')
