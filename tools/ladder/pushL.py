"""Enumerate, by pruning, the girth-legal closed ladders with NO cycle of
length exactly L.  Having such a cycle is monotone under adding chords.

Speedups over 588's version:
  - the search for an L-cycle through the new chord is pruned by distance
    to the target endpoint;
  - the first chord is fixed up to reflection.
"""
import sys, time
from collections import deque

def bfs(adj, n, src):
    d = [-1]*n; d[src] = 0; q = deque([src])
    while q:
        u = q.popleft()
        for w in adj[u]:
            if d[w] < 0: d[w] = d[u]+1; q.append(w)
    return d

def has_len_through(adj, n, a, b, L):
    """simple cycle of length exactly L using the edge (a,b)?
    i.e. a simple path of length L-1 from a to b avoiding that edge."""
    target = L-1
    db = bfs(adj, n, b)
    if db[a] < 0: return False
    seen = bytearray(n)
    found = [False]
    def rec(u, d):
        if found[0]: return
        rem = target-d
        if db[u] > rem: return
        if rem == 0:
            if u == b: found[0] = True
            return
        for w in adj[u]:
            if w == b:
                if rem == 1: found[0] = True; return
                continue
            if seen[w]: continue
            seen[w] = 1
            rec(w, d+1)
            seen[w] = 0
            if found[0]: return
    seen[a] = 1
    for w in adj[a]:
        if w == b: continue
        if seen[w]: continue
        seen[w] = 1
        rec(w, 1)
        seen[w] = 0
        if found[0]: return True
    return False

def short(adj, n, a, b, maxlen=4):
    d = {a: 0}; q = deque([a])
    while q:
        u = q.popleft()
        if d[u] >= maxlen-1: continue
        for w in adj[u]:
            if w not in d: d[w] = d[u]+1; q.append(w)
    return b in d

def enumerate_missing(n, L, report=None):
    adj = [[(v+1) % n, (v-1) % n] for v in range(n)]
    mate = [-1]*n
    out = []
    nodes = [0]
    def rec():
        nodes[0] += 1
        a = next((v for v in range(n) if mate[v] == -1), None)
        if a is None:
            out.append(list(mate)); return
        lo = a+1
        for b in range(lo, n):
            if mate[b] != -1: continue
            if short(adj, n, a, b): continue
            if a == 0 and b > n//2: break          # reflection
            adj[a].append(b); adj[b].append(a)
            mate[a] = b; mate[b] = a
            if not has_len_through(adj, n, a, b, L):
                rec()
            adj[a].pop(); adj[b].pop()
            mate[a] = -1; mate[b] = -1
    rec()
    return out, nodes[0]

def bigpow(n):
    L = 1
    while L*2 <= n: L *= 2
    return L

if __name__ == '__main__':
    sys.setrecursionlimit(100000)
    for n in [int(x) for x in sys.argv[1:]]:
        L = bigpow(n)
        if L == n:
            print("  2N=%2d : L=%2d = 2N, host cycle -- trivial" % (n, L)); continue
        t = time.time()
        out, nodes = enumerate_missing(n, L)
        print("  2N=%2d  L=%2d :  ladders with no L-cycle = %d   "
              "[search nodes %d, %.0fs]" % (n, L, len(out), nodes, time.time()-t))
        for m in out[:3]: print("        ", m)
        sys.stdout.flush()
