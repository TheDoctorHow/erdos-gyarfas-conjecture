"""Where does the octagon first fail?  Randomised completion of a closed
ladder with no 8-cycle, then a check for a 16-cycle."""
import sys, random, time
from collections import deque
from no8 import has_cycle_through, girth, adj_of, has_len

def try_build(n, rng, forbid=(8,), tries=4000):
    for _ in range(tries):
        adj = [[(v+1) % n, (v-1) % n] for v in range(n)]
        mate = [-1]*n; ok = True
        while ok:
            a = next((v for v in range(n) if mate[v] == -1), None)
            if a is None: break
            cand = []
            for b in range(n):
                if b == a or mate[b] != -1: continue
                d = {a: 0}; q = deque([a])
                while q:
                    u = q.popleft()
                    if d[u] >= 3: continue
                    for w in adj[u]:
                        if w not in d: d[w] = d[u]+1; q.append(w)
                if b in d and d[b] <= 3: continue
                adj[a].append(b); adj[b].append(a)
                good = all(not has_cycle_through(adj, n, a, b, L) for L in forbid)
                adj[a].pop(); adj[b].pop()
                if good: cand.append(b)
            if not cand: ok = False; break
            b = rng.choice(cand)
            adj[a].append(b); adj[b].append(a)
            mate[a] = b; mate[b] = a
        if ok and all(m != -1 for m in mate):
            return mate
    return None

if __name__ == '__main__':
    sys.setrecursionlimit(100000)
    rng = random.Random(1234567)
    print("smallest closed ladder with NO 8-cycle, and whether it has a 16-cycle")
    first = None
    for n in range(24, 61, 2):
        t = time.time()
        m = try_build(n, rng, forbid=(8,), tries=300)
        if m is None:
            print("  2N=%2d : none found (%.0fs)" % (n, time.time()-t))
        else:
            g = girth(adj_of(m, n), n)
            h16 = has_len(m, n, 16)
            print("  2N=%2d : FOUND, girth %d, 16-cycle: %s   (%.0fs)"
                  % (n, g, h16, time.time()-t))
            print("        mate %s" % (m,))
            if first is None: first = (n, m, g, h16)
        sys.stdout.flush()
    print()
    print("first size with an 8-cycle-free ladder:", first[0] if first else None)
