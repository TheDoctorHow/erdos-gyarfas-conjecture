"""Where does the length-15 exclusion set stop covering the whole ladder?"""
import sys, random, time

def endpoints_fast(adj, n, a, depth):
    """endpoints of simple paths of length exactly `depth` from a."""
    vis = bytearray(n)
    out = set()
    path = [a]
    vis[a] = 1
    def rec(u, d):
        if d == depth:
            out.add(u); return
        for w in adj[u]:
            if vis[w]: continue
            vis[w] = 1
            rec(w, d+1)
            vis[w] = 0
    sys.setrecursionlimit(10000)
    rec(a, 0)
    return out

def random_ladder(n, rng):
    """girth >= 5 closed ladder, greedy with restarts."""
    for _ in range(50):
        adj = [[(v+1) % n, (v-1) % n] for v in range(n)]
        mate = [-1]*n; ok = True
        for a in range(n):
            if mate[a] != -1: continue
            # forbid b within distance 3
            bad = {a}
            frontier = {a}
            for _ in range(3):
                nxt = set()
                for u in frontier:
                    for w in adj[u]:
                        if w not in bad: nxt.add(w); bad.add(w)
                frontier = nxt
            cand = [b for b in range(n) if mate[b] == -1 and b not in bad]
            if not cand: ok = False; break
            b = rng.choice(cand)
            adj[a].append(b); adj[b].append(a)
            mate[a] = b; mate[b] = a
        if ok and all(m != -1 for m in mate): return mate, adj
    return None, None

if __name__ == '__main__':
    rng = random.Random(20260811)
    print("   2N      |excl L=8|   |excl L=16|   fraction of 2N   time")
    for n in (1600, 6400, 12800, 25600, 51200, 102400):
        m, adj = random_ladder(n, rng)
        if m is None: print("  %6d : build failed" % n); continue
        t = time.time()
        e8 = len(endpoints_fast(adj, n, 0, 7))
        e16 = len(endpoints_fast(adj, n, 0, 15))
        print("  %6d      %6d       %6d        %.3f        %.0fs"
              % (n, e8, e16, e16/n, time.time()-t))
        sys.stdout.flush()
