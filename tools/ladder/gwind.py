"""Winding-zero cycles of a GENERAL closed ladder.

Lift: each chord {u,v} gets a LEFT residue and a span (the short side).
In the cover Z the rung positions R and the rung ends mu(R) partition Z,
which is the only fact (WIND0)'s proof used.
"""
import itertools, sys
from collections import deque

def lift(mate, n):
    """isleft[r], span[r] for r in Z_n.  span[r] defined for left residues."""
    isleft = [False]*n; span = [0]*n
    for u in range(n):
        v = mate[u]
        if u >= v: continue
        d = v-u
        if d < n-d:   isleft[u] = True; span[u] = d
        elif d > n-d: isleft[v] = True; span[v] = n-d
        else:         isleft[u] = True; span[u] = d     # diameter: pick u
    assert sum(isleft) == n//2
    return isleft, span

def config(S, spans):
    """S sorted rung positions in Z, spans[i] the span of S[i].
    Returns (length, ok, vertices) with ok = pairing graph is one cycle."""
    k = len(S)
    pts = []
    for i, x in enumerate(S):
        pts.append((x, i)); pts.append((x+spans[i], i))
    pts.sort()
    if len(set(p for p, _ in pts)) != 2*k: return None, False, None
    adj = [[] for _ in range(k)]
    host = 0; verts = set()
    for j in range(0, 2*k, 2):
        (a, ia), (b, ib) = pts[j], pts[j+1]
        adj[ia].append(ib); adj[ib].append(ia)
        host += b-a
        for t in range(a, b+1): verts.add(t)
    if any(len(x) != 2 for x in adj): return None, False, None
    seen = {0}; prev = -1; cur = 0
    for _ in range(k):
        nxt = adj[cur][0] if adj[cur][0] != prev else adj[cur][1]
        prev, cur = cur, nxt; seen.add(cur)
    return k+host, len(seen) == k, verts

def wind0_lengths(mate, n, kmax=5, reps=2):
    """lengths of winding-zero cycles using at most kmax rungs, with rung
    positions taken from [0, reps*n)."""
    isleft, span = lift(mate, n)
    pos = [x for x in range(reps*n) if isleft[x % n]]
    out = {}
    for k in range(1, kmax+1):
        for S in itertools.combinations(pos, k):
            sp = [span[x % n] for x in S]
            L, ok, verts = config(list(S), sp)
            if not ok: continue
            if len(set(v % n for v in verts)) != len(verts): continue
            if L not in out: out[L] = list(S)
    return out

# ---------- ground truth ----------

def all_cycles(mate, n):
    A = [[(v+1) % n, (v-1) % n, mate[v]] for v in range(n)]
    out = []
    def rec(start, u, vis, path):
        for w in A[u]:
            if w == start and len(path) >= 3: out.append(list(path))
            elif w > start and not (vis >> w) & 1:
                path.append(w); rec(start, w, vis | (1 << w), path); path.pop()
    for s0 in range(n):
        rec(s0, s0, 1 << s0, [s0])
    return out

def wind(cyc, mate, n, isleft, span):
    d = 0
    for i in range(len(cyc)):
        u = cyc[i]; v = cyc[(i+1) % len(cyc)]
        if mate[u] == v and (v-u) % n not in (1, n-1):
            d += span[u] if isleft[u] else -span[v]
        else:
            d += 1 if (v-u) % n == 1 else -1
    assert d % n == 0, (cyc, d)
    return d // n

def girth(mate, n):
    A = [[(v+1) % n, (v-1) % n, mate[v]] for v in range(n)]
    g = 10**9
    for src in range(n):
        dist = {src: 0}; par = {src: -1}; q = deque([src])
        while q:
            u = q.popleft()
            for w in A[u]:
                if w not in dist:
                    dist[w] = dist[u]+1; par[w] = u; q.append(w)
                elif w != par[u]: g = min(g, dist[u]+dist[w]+1)
    return g

if __name__ == '__main__':
    from skip3 import all_ladders
    print("== (WIND0) on general ladders: criterion vs brute force ==")
    for n in (12, 14):
        Ls = [m for m in all_ladders(n) if girth(m, n) >= 5]
        bad = 0
        for mate in Ls:
            isleft, span = lift(mate, n)
            truth = set()
            for c in all_cycles(mate, n):
                if wind(c, mate, n, isleft, span) == 0: truth.add(len(c))
            got = set(wind0_lengths(mate, n, kmax=6, reps=2))
            if not got <= truth:
                bad += 1
                if bad == 1: print("  UNSOUND", mate, sorted(got-truth))
            if got != truth and bad == 0:
                pass
        print("  2N=%d : %d ladders, criterion never claims a false length: %s"
              % (n, len(Ls), bad == 0))
