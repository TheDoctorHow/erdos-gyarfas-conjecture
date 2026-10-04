"""How large can the girth of a girth-legal closed ladder be?
Any power-of-two cycle has length at least the girth, so a bounded-window
construction is possible only if the girth is bounded over the family.
"""
import random, sys, time
from collections import deque
from gwind import girth

def dist_from(a, mate, n, cap):
    """BFS distances from a in host + already placed chords, capped."""
    A = [[(v+1) % n, (v-1) % n] for v in range(n)]
    for v in range(n):
        if mate[v] is not None: A[v].append(mate[v])
    d = {a: 0}; q = deque([a])
    while q:
        u = q.popleft()
        if d[u] >= cap: continue
        for w in A[u]:
            if w not in d: d[w] = d[u]+1; q.append(w)
    return d

def build(n, rng, g, tries=400):
    """random closed ladder of girth >= g, or None."""
    for _ in range(tries):
        mate = [None]*n; ok = True
        while ok:
            a = next((v for v in range(n) if mate[v] is None), None)
            if a is None: break
            d = dist_from(a, mate, n, g-2)
            cand = [b for b in range(a+1, n)
                    if mate[b] is None and d.get(b, 10**9) >= g-1]
            if not cand: ok = False; break
            b = rng.choice(cand); mate[a] = b; mate[b] = a
        if ok and all(m is not None for m in mate):
            return mate
    return None

def has_cycle(mate, n, T, seconds=30.0, seed=5):
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
    def dfs(u, used, ln):
        if time.time()-t0 > seconds: raise TimeoutError
        rem = T-ln
        if rem == 0: return 0 in A[u]
        nb = A[u][:]; rng.shuffle(nb)
        for w in nb:
            if w == 0:
                if ln+1 == T: return True
                continue
            if (used >> w) & 1: continue
            if dist[w][0] > rem-1: continue
            if dfs(w, used | (1 << w), ln+1): return True
        return False
    try: return dfs(0, 1, 0)
    except TimeoutError: return None

if __name__ == '__main__':
    sys.setrecursionlimit(100000)
    rng = random.Random(2026)
    print("largest girth found, and which small powers of two survive:")
    for n in (20, 30, 40, 60, 80, 100, 140):
        best = None; bg = 0
        for g in range(6, 15):
            m = build(n, rng, g, tries=60)
            if m is None: break
            gg = girth(m, n)
            if gg > bg: bg = gg; best = m
        ps = []
        for T in (8, 16, 32, 64, 128):
            if T <= n:
                r = has_cycle(best, n, T)
                ps.append("%d:%s" % (T, {True: 'yes', False: 'NO', None: '?'}[r]))
        print("  2N=%3d   girth %2d   %s" % (n, bg, "   ".join(ps)))
        sys.stdout.flush()
