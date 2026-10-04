"""How many chords can a girth-legal partial matching carry before a cycle
of length L is unavoidable?

The pruned search of 591 abandons a branch when an L-cycle APPEARS.  So
the maximum depth it reaches is the largest number of chords in a
girth-legal partial matching with no L-cycle.  If that maximum is below
N, no complete matching can avoid an L-cycle - which is exactly the zero
result, now read as a forcing statement.
"""
import sys, time
from collections import deque
from pushL import has_len_through, short, bigpow

def max_depth(n, L):
    adj = [[(v+1) % n, (v-1) % n] for v in range(n)]
    mate = [-1]*n
    best = [0]; witness = [None]
    def rec(depth):
        if depth > best[0]:
            best[0] = depth; witness[0] = list(mate)
        a = next((v for v in range(n) if mate[v] == -1), None)
        if a is None: return
        for b in range(a+1, n):
            if mate[b] != -1: continue
            if short(adj, n, a, b): continue
            if a == 0 and b > n//2: break
            adj[a].append(b); adj[b].append(a)
            mate[a] = b; mate[b] = a
            if not has_len_through(adj, n, a, b, L):
                rec(depth+1)
            adj[a].pop(); adj[b].pop()
            mate[a] = -1; mate[b] = -1
    rec(0)
    return best[0], witness[0]

if __name__ == '__main__':
    sys.setrecursionlimit(100000)
    print("  2N   N    L   max chords with no L-cycle   N - max")
    for n in [int(x) for x in sys.argv[1:]]:
        L = bigpow(n)
        if L == n:
            print("  %2d   -    %2d   trivial (L = 2N)" % (n, L)); continue
        t = time.time()
        d, w = max_depth(n, L)
        print("  %2d  %2d   %2d          %2d                  %2d     [%.0fs]"
              % (n, n//2, L, d, n//2-d, time.time()-t))
        sys.stdout.flush()
