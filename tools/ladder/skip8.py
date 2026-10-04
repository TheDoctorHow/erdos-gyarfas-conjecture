"""The uniform ladders the corrected one-lap construction misses:
do they have a power-of-two cycle at all?  Mixed-direction one-lap search,
then a direct cycle search."""
import sys, itertools, random
from skip4 import uniform

def pows(n):
    out = []; L = 8
    while L <= n: out.append(L); L *= 2
    return out

def onelap_mixed(n, s, kmax):
    """lengths of one-lap simple cycles with at most kmax chords,
    chords in EITHER direction."""
    mate = uniform(n, s)
    d = [(mate[v]-v) % n for v in range(n)]
    lens = {}
    for k in range(1, kmax+1):
        for S in itertools.combinations(range(n), k):
            ok = True
            b = [0]*k
            for j in range(k):
                x = S[j]; y = S[(j+1) % k]
                g = n if k == 1 else (y-x) % n
                if g == d[x]: ok = False; break
                b[j] = d[x]-g if d[x] > g else 0
            if not ok: continue
            for j in range(k):
                if b[j] + b[(j+1) % k] >= d[S[(j+1) % k]]: ok = False; break
            if not ok: continue
            L = n - sum(d[x]-1 for x in S) + 2*sum(b)
            if L not in lens: lens[L] = (list(S), list(b))
    return lens

def find_cycle(n, s, target, seconds=20.0, seed=1):
    """randomised DFS for a simple cycle of exactly `target` edges."""
    import time
    mate = uniform(n, s)
    A = [[(v+1) % n, (v-1) % n, mate[v]] for v in range(n)]
    rng = random.Random(seed)
    t0 = time.time()
    # BFS distances from every vertex (cheap: n * n)
    from collections import deque
    dist = []
    for src in range(n):
        dd = [-1]*n; dd[src] = 0; q = deque([src])
        while q:
            u = q.popleft()
            for w in A[u]:
                if dd[w] < 0: dd[w] = dd[u]+1; q.append(w)
        dist.append(dd)
    start = 0
    def dfs(u, used, ln, path):
        if time.time() - t0 > seconds: raise TimeoutError
        rem = target - ln
        if rem == 0:
            return path if start in A[u] and ln == target else None
        nb = A[u][:]
        rng.shuffle(nb)
        for w in nb:
            if w == start:
                if ln+1 == target: return path+[w]
                continue
            if (used >> w) & 1: continue
            if dist[w][start] > rem-1: continue
            r = dfs(w, used | (1 << w), ln+1, path+[w])
            if r: return r
        return None
    try:
        return dfs(start, 1 << start, 0, [start])
    except TimeoutError:
        return None

if __name__ == '__main__':
    cases = [(41, 92), (43, 96), (45, 100), (73, 156)]
    for s, n in cases:
        print("== U(%d, %d)   powers of two available: %s ==" % (n, s, pows(n)))
        got = onelap_mixed(n, s, 3)
        hits = [L for L in pows(n) if L in got]
        print("   one-lap, mixed direction, k <= 3 :", hits,
              {L: got[L] for L in hits})
        sys.stdout.flush()
