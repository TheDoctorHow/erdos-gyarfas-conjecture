"""The deficiency-2 telescope, instrumented.

KEY IDENTITY (checked below): the re-presented object H after omitting
u, v is exactly the INDUCED SUBGRAPH G[V - {u,v}]: the new host is the
(M-2)-cycle, and every vertex's third edge survives iff its far end
survives.  So the telescope is a chain of induced subgraphs
    G = G_0 > G_1 > ... > G_k,   |G_{i+1}| = |G_i| - 2,
each Hamiltonian (the next host), ending at |G_k| = L where the host
itself is the L-cycle.  Girth is inherited because subgraphs only lose
cycles.

Fast deficiency-2 enumeration: for S with sorted endpoints e_1..e_{2k},
    length(S,0) = |S| + A,   length(S,1) = |S| + M - A,
    A = sum over odd j of (e_{j+1} - e_j),
so masks are filtered in O(k log k) before any edge set is built.
"""
import sys, time, random
from pushL import bigpow
from hitset import random_ladder

def build_check(M, ch, S, eps, target):
    tog = [0] * M
    for a, b in S:
        tog[a] ^= 1; tog[b] ^= 1
    t = eps
    edges = list(S)
    for v in range(M):
        t ^= tog[v]
        if t:
            edges.append((v, (v + 1) % M))
    if len(edges) != target:
        return None
    deg = {}
    for a, b in edges:
        deg[a] = deg.get(a, 0) + 1
        deg[b] = deg.get(b, 0) + 1
    if any(d != 2 for d in deg.values()):
        return None
    adj = {}
    for a, b in edges:
        adj.setdefault(a, []).append(b)
        adj.setdefault(b, []).append(a)
    st = [next(iter(adj))]; seen = {st[0]}
    while st:
        u = st.pop()
        for w in adj[u]:
            if w not in seen:
                seen.add(w); st.append(w)
    return edges if len(seen) == len(deg) else None

def cycles_of_len(M, ch, target, cap=100000):
    """all simple cycles of length exactly target, via the length filter"""
    k = len(ch)
    out = []
    for mask in range(1, 1 << k):
        S = [ch[i] for i in range(k) if mask >> i & 1]
        es = sorted(x for c in S for x in c)
        A = sum(es[j + 1] - es[j] for j in range(0, len(es), 2))
        for eps, ln in ((0, len(S) + A), (1, len(S) + M - A)):
            if ln != target:
                continue
            e = build_check(M, ch, S, eps, target)
            if e is not None:
                out.append(e)
                if len(out) >= cap:
                    return out
    if target == M:                      # host circle, S empty, eps=1
        out.append([(v, (v + 1) % M) for v in range(M)])
    return out

def has_len(M, ch, target):
    return bool(cycles_of_len(M, ch, target, cap=1))

def step(M, ch, cyc_edges):
    es = set(tuple(sorted(e)) for e in cyc_edges)
    verts = set()
    for a, b in es:
        verts.add(a); verts.add(b)
    adjfull = {v: set() for v in range(M)}
    for v in range(M):
        adjfull[v].add((v + 1) % M); adjfull[(v + 1) % M].add(v)
    for a, b in ch:
        adjfull[a].add(b); adjfull[b].add(a)
    cyc_adj = {}
    for a, b in es:
        cyc_adj.setdefault(a, []).append(b)
        cyc_adj.setdefault(b, []).append(a)
    start = min(verts)
    order = [start]
    prev, cur = None, start
    while True:
        nxt = next(w for w in cyc_adj[cur] if w != prev)
        if nxt == start:
            break
        order.append(nxt); prev, cur = cur, nxt
    assert len(order) == len(verts)
    pos = {v: i for i, v in enumerate(order)}
    newch = set()
    for v in order:
        for w in adjfull[v]:
            e = tuple(sorted((v, w)))
            if e in es:
                continue
            if w in pos:
                newch.add(tuple(sorted((pos[v], pos[w]))))
    return len(order), sorted(newch)

def parallel_pair(M, ch):
    s = set(ch)
    for a, b in ch:
        for da, db in (((a+2) % M, (b+2) % M), ((a-2) % M, (b-2) % M)):
            if tuple(sorted((da, db))) in s:
                return True
    return False

def telescope(M, ch, L, trail, budget, order_shuffle=None):
    if M == L:
        return True
    if budget[0] <= 0:
        return False
    budget[0] -= 1
    cands = cycles_of_len(M, ch, M - 2)
    if order_shuffle:
        order_shuffle.shuffle(cands)
    trail.append((M, len(ch), len(cands), parallel_pair(M, ch)))
    for e in cands:
        Mn, nch = step(M, ch, e)
        if telescope(Mn, nch, L, trail, budget, order_shuffle):
            return True
    return False

def uniform(n, s):
    """U(2N,s): chords (2j, 2j+s), s odd -- girth 6"""
    assert s % 2 == 1
    return sorted(tuple(sorted((2*j % n, (2*j+s) % n))) for j in range(n//2))

if __name__ == '__main__':
    mode = sys.argv[1]
    if mode == 'U':
        for n in [int(x) for x in sys.argv[2:]]:
            L = bigpow(n)
            if L == n:
                print("  U(%d,5) trivial" % n); continue
            ch = uniform(n, 5)
            trail = []; budget = [3000]
            t = time.time()
            r = telescope(n, ch, L, trail, budget)
            print("  U(%2d,5) -> L=%2d : %s  nodes=%d  [%.0fs]"
                  % (n, L, "REACHED" if r else "FAILED", 3000-budget[0],
                     time.time()-t))
            print("      trail", trail[:8])
            sys.stdout.flush()
    elif mode == 'R':
        rng = random.Random(int(sys.argv[2]))
        trials = int(sys.argv[3]); n = int(sys.argv[4])
        L = bigpow(n)
        ok = tot = 0
        for t in range(trials):
            ch = random_ladder(n, rng)
            if ch is None: continue
            trail = []; budget = [3000]
            t0 = time.time()
            r = telescope(n, ch, L, trail, budget)
            tot += 1; ok += bool(r)
            print("  ladder %2d: %s nodes=%d [%.0fs] trail %s"
                  % (t, "REACHED" if r else "FAILED", 3000-budget[0],
                     time.time()-t0, trail[:7]))
            sys.stdout.flush()
        print("  == %d/%d reached L=%d at 2N=%d" % (ok, tot, L, n))
    elif mode == 'HUNT':
        # exhaustive-DFS failure hunt: a ladder where the telescope fails
        # (with unlimited backtracking) despite having an L-cycle would
        # break the strategy.  Randomised order to vary the paths.
        rng = random.Random(int(sys.argv[2]))
        trials = int(sys.argv[3]); n = int(sys.argv[4])
        L = bigpow(n)
        fails = 0; tot = 0
        for t in range(trials):
            ch = random_ladder(n, rng)
            if ch is None: continue
            trail = []; budget = [100000]
            r = telescope(n, ch, L, trail, budget, order_shuffle=rng)
            tot += 1
            if not r:
                fails += 1
                print("  FAILURE at ladder %d: %s" % (t, ch))
            if tot % 50 == 0:
                print("  ... %d done, %d failures" % (tot, fails))
                sys.stdout.flush()
        print("  == %d/%d failures at 2N=%d (exhaustive DFS)"
              % (fails, tot, n))
