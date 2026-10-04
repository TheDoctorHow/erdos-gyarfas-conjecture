"""Push the girth of a closed ladder as high as the machine allows, and
read off the smallest power-of-two cycle."""
import sys, random, time
from collections import deque

def build(n, g, rng, tries=200):
    """closed ladder (host cycle 0..n-1 plus perfect matching) of girth >= g."""
    for _ in range(tries):
        mate = [-1]*n
        nbr = [[(v+1) % n, (v-1) % n] for v in range(n)]
        order = list(range(n)); ok = True
        for a in order:
            if mate[a] != -1: continue
            # BFS from a to depth g-2 in the current graph
            d = {a: 0}; q = deque([a])
            while q:
                u = q.popleft()
                if d[u] >= g-2: continue
                for w in nbr[u]:
                    if w not in d: d[w] = d[u]+1; q.append(w)
            cand = [b for b in range(n)
                    if mate[b] == -1 and b != a and d.get(b, 10**9) >= g-1]
            if not cand: ok = False; break
            b = rng.choice(cand)
            mate[a] = b; mate[b] = a
            nbr[a].append(b); nbr[b].append(a)
        if ok and all(m != -1 for m in mate): return mate
    return None

def girth(mate, n):
    A = [[(v+1) % n, (v-1) % n, mate[v]] for v in range(n)]
    best = 10**9
    for src in range(n):
        d = {src: 0}; par = {src: -1}; q = deque([src])
        while q:
            u = q.popleft()
            if 2*d[u] >= best: break
            for w in A[u]:
                if w not in d:
                    d[w] = d[u]+1; par[w] = u; q.append(w)
                elif w != par[u]:
                    best = min(best, d[u]+d[w]+1)
    return best

def has_cycle(mate, n, T, seconds, seed=1):
    A = [[(v+1) % n, (v-1) % n, mate[v]] for v in range(n)]
    dist = []
    for src in range(n):
        dd = [-1]*n; dd[src] = 0; q = deque([src])
        while q:
            u = q.popleft()
            for w in A[u]:
                if dd[w] < 0: dd[w] = dd[u]+1; q.append(w)
        dist.append(dd)
    rng = random.Random(seed); t0 = time.time()
    def dfs(s0, u, used, ln):
        if time.time()-t0 > seconds: raise TimeoutError
        rem = T-ln
        if rem == 0: return s0 in A[u]
        nb = A[u][:]; rng.shuffle(nb)
        for w in nb:
            if w == s0:
                if ln+1 == T: return True
                continue
            if (used >> w) & 1: continue
            if dist[w][s0] > rem-1: continue
            if dfs(s0, w, used | (1 << w), ln+1): return True
        return False
    try:
        for s0 in range(0, n, max(1, n//8)):
            if dfs(s0, s0, 1 << s0, 0): return True
        return False
    except TimeoutError:
        return None

if __name__ == '__main__':
    sys.setrecursionlimit(200000)
    rng = random.Random(20260811)
    print("2N     target  girth   smallest power-of-two cycle")
    for n, g in [(70, 10), (100, 10), (140, 10), (200, 10), (200, 11),
                 (300, 11), (400, 11), (400, 12), (600, 12)]:
        t0 = time.time()
        m = build(n, g, rng, tries=30)
        if m is None:
            print("  %4d  %5d   -- not built (%.0fs)" % (n, g, time.time()-t0))
            sys.stdout.flush(); continue
        gg = girth(m, n)
        res = []
        for T in (8, 16, 32, 64):
            if T > n: break
            r = has_cycle(m, n, T, 25.0)
            res.append("%d:%s" % (T, {True: 'yes', False: 'NO', None: '?'}[r]))
            if r: break
        print("  %4d  %5d   %5d   %s   (%.0fs)"
              % (n, g, gg, "  ".join(res), time.time()-t0))
        sys.stdout.flush()
