"""Brute force: over ALL girth-legal partial matchings on the 2N-circle,
the largest number of chords carrying no cycle of length exactly L.

No ordering assumption at all -- every partial matching is generated,
including those that leave a low-numbered vertex unmatched.
"""
import sys
from pushL import short, bigpow

def has_cycle_len(adj, n, L):
    """any simple cycle of length exactly L?  start at its least vertex."""
    for s in range(n - L + 1):
        seen = bytearray(n)
        seen[s] = 1
        stack = [(s, 1)]
        # plain DFS over simple paths from s using only vertices > s
        def rec(u, d):
            if d == L:
                return s in adj[u]
            for w in adj[u]:
                if w < s or seen[w]:
                    continue
                seen[w] = 1
                if rec(w, d + 1):
                    seen[w] = 0
                    return True
                seen[w] = 0
            return False
        if rec(s, 1):
            return True
    return False

def brute(n, L):
    base = [[(v + 1) % n, (v - 1) % n] for v in range(n)]
    best = [0]
    wit = [None]
    mate = [-1] * n

    def rec(v, k):
        if v == n:
            if k > best[0]:
                adj = [list(x) for x in base]
                for a in range(n):
                    if mate[a] > a:
                        adj[a].append(mate[a])
                        adj[mate[a]].append(a)
                if not has_cycle_len(adj, n, L):
                    best[0] = k
                    wit[0] = list(mate)
            return
        if mate[v] != -1:
            rec(v + 1, k)
            return
        # skip v
        rec(v + 1, k)
        # match v
        adj = [list(base[x]) for x in range(n)]
        for a in range(n):
            if mate[a] > a:
                adj[a].append(mate[a])
                adj[mate[a]].append(a)
        for w in range(v + 1, n):
            if mate[w] != -1:
                continue
            if short(adj, n, v, w):
                continue
            mate[v] = w
            mate[w] = v
            rec(v + 1, k + 1)
            mate[v] = -1
            mate[w] = -1

    rec(0, 0)
    return best[0], wit[0]

if __name__ == '__main__':
    for n in [int(x) for x in sys.argv[1:]]:
        L = bigpow(n)
        if L == n:
            print("  2N=%2d  trivial" % n)
            continue
        b, w = brute(n, L)
        print("  2N=%2d  L=%2d  brute max chords = %d   margin = %d"
              % (n, L, b, n // 2 - b))
        print("        witness", [(a, w[a]) for a in range(n) if w[a] > a])
        sys.stdout.flush()
