"""Winding-zero cycles of U(2N,s).

A winding-zero cycle lifts to a cycle of the INFINITE ladder L(s):
vertices Z, host edges i~i+1, rungs {m, m+s} for m even.  In L(s) the
cycle space is spanned by the fundamental cycles of the rungs, so a
cycle is determined by its RUNG SET S alone: host edge (t,t+1) is used
iff an odd number of the intervals [m, m+s), m in S, contain t.
"""
import itertools

def cycle_from(S, s):
    """edges of the element of the cycle space with rung set S.
    Returns (rungs, host_edges) with host edges as left endpoints t."""
    S = sorted(S)
    ends = sorted([m for m in S] + [m+s for m in S])
    host = []
    for i in range(0, len(ends), 2):
        for t in range(ends[i], ends[i+1]):
            host.append(t)
    return S, host, ends

def pairing_is_one_cycle(S, s):
    """the k pairs of consecutive sorted endpoints, viewed as edges on the
    k rungs, must form a single cycle."""
    S = sorted(S); k = len(S)
    if k == 0: return False
    owner = {}
    for i, m in enumerate(S):
        owner[m] = i; owner[m+s] = i
    ends = sorted(owner)
    if len(ends) != 2*k: return False          # a start coinciding with an end
    adj = [[] for _ in range(k)]
    for i in range(0, 2*k, 2):
        a = owner[ends[i]]; b = owner[ends[i+1]]
        adj[a].append(b); adj[b].append(a)
    if any(len(a) != 2 for a in adj): return False
    seen = {0}; cur = 0; prev = -1
    for _ in range(k):
        nxt = adj[cur][0] if adj[cur][0] != prev else adj[cur][1]
        prev, cur = cur, nxt
        seen.add(cur)
    return len(seen) == k

def length_of(S, s):
    S = sorted(S)
    ends = sorted([m for m in S] + [m+s for m in S])
    return len(S) + sum(ends[i+1]-ends[i] for i in range(0, len(ends), 2))

def vertices_of(S, s):
    """the vertex set of the lifted cycle."""
    _, host, ends = cycle_from(S, s)
    vs = set()
    for i in range(0, len(ends), 2):
        for v in range(ends[i], ends[i+1]+1): vs.add(v)
    return vs

def embeds(S, s, n):
    """distinct residues mod n?"""
    vs = vertices_of(S, s)
    return len(set(v % n for v in vs)) == len(vs)

def realise(S, s, n):
    """build the cycle in U(n,s) and verify it is a simple cycle of the
    predicted length with winding zero.  Returns (ok, length)."""
    S = sorted(S)
    mate = [None]*n
    for i in range(0, n, 2):
        mate[i] = (i+s) % n; mate[(i+s) % n] = i
    _, host, ends = cycle_from(S, s)
    edges = []
    for m in S: edges.append(('c', m % n, (m+s) % n))
    for t in host: edges.append(('h', t % n, (t+1) % n))
    deg = {}
    for _, a, b in edges:
        deg[a] = deg.get(a, 0)+1; deg[b] = deg.get(b, 0)+1
    if any(d != 2 for d in deg.values()): return False, None
    if len(edges) != len(deg): return False, None
    nb = {}
    for _, a, b in edges:
        nb.setdefault(a, []).append(b); nb.setdefault(b, []).append(a)
    start = edges[0][1]; prev = None; cur = start; cnt = 0; disp = 0
    while True:
        w = nb[cur][0] if nb[cur][0] != prev else nb[cur][1]
        if mate[cur] == w and abs(((w-cur) % n)) != 1:
            disp += s if cur % 2 == 0 else -s
        else:
            disp += 1 if (w-cur) % n == 1 else -1
        prev, cur = cur, w; cnt += 1
        if cur == start: break
    return (cnt == len(edges) and disp == 0), cnt

# ---------- ground truth ----------

def all_cycles(n, s):
    mate = [None]*n
    for i in range(0, n, 2):
        mate[i] = (i+s) % n; mate[(i+s) % n] = i
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

def winding(cyc, n, s):
    mate = [None]*n
    for i in range(0, n, 2):
        mate[i] = (i+s) % n; mate[(i+s) % n] = i
    d = 0
    for i in range(len(cyc)):
        u = cyc[i]; v = cyc[(i+1) % len(cyc)]
        if mate[u] == v and (v-u) % n not in (1, n-1):
            d += s if u % 2 == 0 else -s
        else:
            d += 1 if (v-u) % n == 1 else -1
    assert d % n == 0
    return d // n

if __name__ == '__main__':
    print("== criterion vs brute force, winding-zero cycles ==")
    for n in (12, 14, 16, 18, 20):
        for s in range(5, n//2-1, 2):
            truth = set()
            for c in all_cycles(n, s):
                if winding(c, n, s) == 0: truth.add(len(c))
            got = set()
            W = n-1
            evens = [m for m in range(0, W+1, 2)]
            for k in range(1, min(len(evens), 8)+1):
                for S in itertools.combinations(evens, k):
                    if not pairing_is_one_cycle(S, s): continue
                    if not embeds(S, s, n): continue
                    ok, ln = realise(S, s, n)
                    if not ok: print("   REALISE FAILED", n, s, S); continue
                    assert ln == length_of(S, s), (S, s, ln, length_of(S, s))
                    got.add(ln)
            print("  2N=%2d s=%d  brute %s   criterion(window<2N) %s   extra %s"
                  % (n, s, sorted(truth), sorted(got), sorted(got-truth)))
