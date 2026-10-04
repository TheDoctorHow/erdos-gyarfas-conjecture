"""cn481.py -- IS A CONGRUENCE THEOREM EVEN AVAILABLE AT MIN DEGREE 3?

The Combinatorial Nullstellensatz, Chevalley-Warning and every other
polynomial method over F_p conclude in RESIDUES: they can force a
subgraph whose degrees, or whose size, satisfy a congruence.  So
before asking whether such a method can help with Erdos-Gyarfas, ask
what congruence statements are even TRUE for cubic graphs:

    Q(k):  does every cubic graph have a cycle of length == 0 mod k ?

Q(3) is a classical theorem (Chen-Saito, to be verified verbatim
before use, 388.1).  Q(4) would say level 2 is never empty.  Q(8)
would be far stronger.  A single cubic graph avoiding all multiples
of k refutes Q(k) and removes that congruence target for good.

Method: 2-edge-swap descent on the total number of cycles of length
divisible by k, capped, with restarts.  A zero is a refutation and
is re-verified by exhaustive cycle enumeration.
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


def count_len(adj, L, cap=10 ** 9):
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


def spectrum_exact(adj):
    """complete cycle-length spectrum by exhaustive enumeration."""
    n = len(adj); found = set()
    for s in range(n):
        vis = [False] * n; vis[s] = True
        def dfs(v, ln):
            for w in adj[v]:
                if w < s: continue
                if w == s:
                    if ln >= 3: found.add(ln)
                    continue
                if vis[w]: continue
                vis[w] = True; dfs(w, ln + 1); vis[w] = False
        dfs(s, 1); vis[s] = False
    return found


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


CAP = 25
rng = random.Random(481481)
P('=' * 76)
P('Q(k): MUST A CUBIC GRAPH HAVE A CYCLE OF LENGTH DIVISIBLE BY k?',
  flush=True)
P('=' * 76)
P('  best = min over the search of the capped count of cycles whose')
P('  length is divisible by k.  0 means Q(k) is REFUTED at that n.')
P()
P('  %3s %9s %9s %9s %9s   %s'
  % ('k', 'n=12', 'n=16', 'n=20', 'n=24', 'refuted?'), flush=True)
WIT = {}
for k in (3, 4, 5, 6, 8):
    row = []; refuted = False
    for n in (12, 16, 20, 24):
        lens = [L for L in range(k if k >= 3 else 3, n + 1) if L % k == 0]
        lens = [L for L in lens if L >= 3]
        if not lens: row.append('n/a'); continue
        t0 = time.time(); best = 10 ** 9; bestG = None
        while time.time() - t0 < 20.0 and best > 0:
            G = randcubic(n, rng)
            if G is None: break
            cur = sum(min(count_len(G, L, CAP), CAP) for L in lens)
            stall = 0
            while cur > 0 and stall < 100 and time.time() - t0 < 20.0:
                H = swap(G, rng)
                if H is None: stall += 1; continue
                c = sum(min(count_len(H, L, CAP), CAP) for L in lens)
                if c <= cur:
                    stall = 0 if c < cur else stall + 1
                    G, cur = H, c
                else: stall += 1
            if cur < best: best, bestG = cur, G
        row.append(str(best) if best < 10 ** 9 else '-')
        if best == 0:
            refuted = True
            WIT.setdefault(k, (n, bestG))
    P('  %3d %9s %9s %9s %9s   %s'
      % (k, row[0], row[1], row[2], row[3], 'YES' if refuted else 'no'),
      flush=True)

P()
P('=' * 76)
P('RE-VERIFICATION OF EACH REFUTING WITNESS (exhaustive enumeration)',
  flush=True)
P('=' * 76)
for k in sorted(WIT):
    n, G = WIT[k]
    sp = sorted(spectrum_exact(G))
    div = [L for L in sp if L % k == 0]
    P('  k = %d, n = %d' % (k, n), flush=True)
    P('    spectrum   : %s' % sp, flush=True)
    P('    multiples of %d present : %s' % (k, div if div else 'NONE'),
      flush=True)
    P('    Q(%d) %s' % (k, 'REFUTED' if not div else 'NOT refuted - recheck'),
      flush=True)
    P('    adj = %s' % G, flush=True)
P()
P('=' * 76)
P('READING', flush=True)
P('=' * 76)
P('  A refuted Q(k) means no polynomial method over F_p can ever')
P('  deliver "some cycle has length divisible by k" for cubic graphs,')
P('  because the statement is false.  A surviving Q(k) is a candidate')
P('  target for such a method - and the only kind of conclusion those')
P('  methods produce.')
