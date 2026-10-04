"""Closed ladders with NO 8-cycle: enumerate them by pruning (the presence
of an 8-cycle is monotone under adding chords), then ask whether they all
have a 16-cycle."""
import sys
from collections import deque

def has_cycle_through(adj, n, a, b, L):
    """is there a simple cycle of length exactly L using the edge (a,b)?
    i.e. a simple path of length L-1 from a to b avoiding that edge."""
    target = L-1
    found = [False]
    def dfs(u, vis, d):
        if found[0]: return
        if d == target:
            if u == b: found[0] = True
            return
        if d + 1 > target: return
        for w in adj[u]:
            if w == b and d+1 < target: continue
            if (vis >> w) & 1: continue
            if w == b and d+1 == target: found[0] = True; return
            dfs(w, vis | (1 << w), d+1)
            if found[0]: return
    # walk from a, forbidding the direct edge to b at step 1
    for w in adj[a]:
        if w == b: continue
        dfs(w, (1 << a) | (1 << w), 1)
        if found[0]: return True
    return False

def girth(adj, n):
    g = 10**9
    for src in range(n):
        d = {src: 0}; par = {src: -1}; q = deque([src])
        while q:
            u = q.popleft()
            if 2*d[u] >= g: break
            for w in adj[u]:
                if w not in d:
                    d[w] = d[u]+1; par[w] = u; q.append(w)
                elif w != par[u]: g = min(g, d[u]+d[w]+1)
    return g

def enumerate_no8(n, limit=None):
    """all girth>=5 closed ladders on Z_n with no 8-cycle."""
    adj = [[(v+1) % n, (v-1) % n] for v in range(n)]
    mate = [-1]*n
    out = []
    def short(a, b):
        """would (a,b) close a cycle of length <= 4?"""
        d = {a: 0}; q = deque([a])
        while q:
            u = q.popleft()
            if d[u] >= 3: continue
            for w in adj[u]:
                if w not in d: d[w] = d[u]+1; q.append(w)
        return b in d and d[b] <= 3
    def rec():
        if limit and len(out) >= limit: return
        a = next((v for v in range(n) if mate[v] == -1), None)
        if a is None:
            out.append(list(mate)); return
        for b in range(a+1, n):
            if mate[b] != -1: continue
            if short(a, b): continue
            adj[a].append(b); adj[b].append(a)
            mate[a] = b; mate[b] = a
            if not has_cycle_through(adj, n, a, b, 8):
                rec()
            adj[a].pop(); adj[b].pop()
            mate[a] = -1; mate[b] = -1
            if limit and len(out) >= limit: return
    rec()
    return out

def adj_of(mate, n):
    return [[(v+1) % n, (v-1) % n, mate[v]] for v in range(n)]

def has_len(mate, n, L):
    adj = adj_of(mate, n)
    for a in range(n):
        for b in adj[a]:
            if b > a and has_cycle_through(adj, n, a, b, L): return True
    return False

if __name__ == '__main__':
    sys.setrecursionlimit(100000)
    for n in (14, 16, 18, 20, 22):
        Ls = enumerate_no8(n)
        no16 = []
        gd = {}
        for m in Ls:
            g = girth(adj_of(m, n), n)
            gd[g] = gd.get(g, 0)+1
            if not has_len(m, n, 16): no16.append((m, g))
        print("2N=%2d : girth-legal ladders with NO 8-cycle: %5d   girths %s"
              "   of these, none with a 16-cycle: %d"
              % (n, len(Ls), dict(sorted(gd.items())), len(no16)))
        for m, g in no16[:3]:
            print("     girth %d  mate %s" % (g, m))
        sys.stdout.flush()
