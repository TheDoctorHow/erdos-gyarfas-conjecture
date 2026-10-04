import sys, itertools
from collections import deque

# ---------- girth-legal closed ladders ----------

def adj_of(mate, n):
    A = [[] for _ in range(n)]
    for v in range(n):
        A[v].append((v+1) % n)
        A[v].append((v-1) % n)
        A[v].append(mate[v])
    return A

def creates_short(mate, n, a, b, maxlen=4):
    """would adding chord (a,b) create a cycle of length <= maxlen?"""
    # BFS from a in current graph (host + chords placed) up to maxlen-1 steps
    A = [[] for _ in range(n)]
    for v in range(n):
        A[v].append((v+1) % n); A[v].append((v-1) % n)
        if mate[v] is not None: A[v].append(mate[v])
    dist = {a: 0}
    q = deque([a])
    while q:
        u = q.popleft()
        if dist[u] >= maxlen-1: continue
        for w in A[u]:
            if w not in dist:
                dist[w] = dist[u]+1
                q.append(w)
    return b in dist and dist[b] <= maxlen-1

def ladders(n):
    mate = [None]*n
    out = []
    def rec():
        a = None
        for v in range(n):
            if mate[v] is None: a = v; break
        if a is None:
            out.append(list(mate)); return
        for b in range(a+1, n):
            if mate[b] is not None: continue
            if not creates_short(mate, n, a, b):
                mate[a] = b; mate[b] = a
                rec()
                mate[a] = None; mate[b] = None
    rec()
    return out

def girth(mate, n):
    A = adj_of(mate, n)
    g = 10**9
    for src in range(n):
        dist = {src: 0}; par = {src: -1}
        q = deque([src])
        while q:
            u = q.popleft()
            for w in A[u]:
                if w not in dist:
                    dist[w] = dist[u]+1; par[w] = u
                    q.append(w)
                elif w != par[u]:
                    g = min(g, dist[u]+dist[w]+1)
    return g

# ---------- the one-lap family ----------

def onelap_lengths(mate, n, kmax=None):
    """all lengths of winding-1 all-forward simple cycles, by the criterion."""
    d = [ (mate[v]-v) % n for v in range(n) ]
    lens = set()
    for m in range(1, 1 << n):
        S = [v for v in range(n) if (m >> v) & 1]
        k = len(S)
        if kmax and k > kmax: continue
        ok = True
        b = [0]*k
        for j in range(k):
            x = S[j]; y = S[(j+1) % k]
            gap = n if k == 1 else (y-x) % n
            if gap == d[x]: ok = False; break
            b[j] = max(0, d[x]-gap)
        if not ok: continue
        for j in range(k):
            if b[j] + b[(j+1) % k] >= d[S[(j+1) % k]]:
                ok = False; break
        if not ok: continue
        L = n - sum(d[x]-1 for x in S) + 2*sum(b)
        lens.add(L)
    return lens

def all_cycle_lengths(mate, n):
    A = adj_of(mate, n)
    lens = set()
    # DFS from smallest vertex of the cycle
    def rec(start, u, visited, ln):
        for w in A[u]:
            if w == start and ln >= 3:
                lens.add(ln)
            elif w > start and w not in visited:
                visited.add(w)
                rec(start, w, visited, ln+1)
                visited.discard(w)
    for s0 in range(n):
        rec(s0, s0, {s0}, 1)
    return lens

def pow2(S):
    return sorted(x for x in S if x >= 4 and (x & (x-1)) == 0)

if __name__ == '__main__':
    for n in [12, 14, 16]:
        Ls = ladders(n)
        Ls = [m for m in Ls if girth(m, n) >= 5]
        print("2N = %d : %d girth-legal closed ladders" % (n, len(Ls)))
        nofull = 0; noone = 0; ex = []
        for mate in Ls:
            full = all_cycle_lengths(mate, n)
            one = onelap_lengths(mate, n)
            assert one <= full, (mate, sorted(one-full))
            if not pow2(full): nofull += 1
            if not pow2(one):
                noone += 1
                if len(ex) < 4: ex.append((list(mate), sorted(full), sorted(one)))
        print("   ladders with NO power-of-two cycle at all      : %d" % nofull)
        print("   ladders where the ONE-LAP family misses it     : %d" % noone)
        for e in ex:
            print("     mate=%s" % (e[0],))
            print("       full spectrum  %s" % (e[1],))
            print("       one-lap only   %s" % (e[2],))
        sys.stdout.flush()
