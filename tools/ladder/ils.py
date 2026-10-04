"""Iterated local search for the deepest girth-legal partial matching of
the 2N-circle with NO cycle of length exactly L.

Every configuration it reports is a witness, so it certifies
    margin <= N - (chords found).
Ruin-and-recreate: drop a few chords at random, refill greedily, keep the
result if it is no worse.  This reproduces the exact optima known at
2N = 18, 20, 22 where plain random growth already reached them, and goes
much deeper than random growth at 2N = 34 and beyond.
"""
import sys, time, random
from pushL import short, has_len_through, bigpow

def build_adj(n, mate):
    adj = [[(v + 1) % n, (v - 1) % n] for v in range(n)]
    for v in range(n):
        if mate[v] > v:
            adj[v].append(mate[v]); adj[mate[v]].append(v)
    return adj

def extend(n, L, mate, adj, rng):
    """add chords in random order while no L-cycle appears"""
    k = sum(1 for v in range(n) if mate[v] > v)
    free = [v for v in range(n) if mate[v] == -1]
    rng.shuffle(free)
    for v in free:
        if mate[v] != -1:
            continue
        cand = [w for w in range(n) if w != v and mate[w] == -1
                and not short(adj, n, v, w)]
        rng.shuffle(cand)
        for w in cand:
            adj[v].append(w); adj[w].append(v)
            if has_len_through(adj, n, v, w, L):
                adj[v].pop(); adj[w].pop(); continue
            mate[v] = w; mate[w] = v; k += 1
            break
    return k

def ils(n, L, rng, seconds, ruin=2):
    t0 = time.time()
    mate = [-1] * n
    adj = build_adj(n, mate)
    best = extend(n, L, mate, adj, rng)
    bestmate = list(mate)
    cur, curmate = best, list(mate)
    while time.time() - t0 < seconds:
        mate = list(curmate)
        adj = build_adj(n, mate)
        ch = [(v, mate[v]) for v in range(n) if mate[v] > v]
        rng.shuffle(ch)
        for a, b in ch[:ruin]:
            adj[a].remove(b); adj[b].remove(a)
            mate[a] = -1; mate[b] = -1
        k = extend(n, L, mate, adj, rng)
        if k >= cur:
            cur, curmate = k, list(mate)
        if k > best:
            best, bestmate = k, list(mate)
            print("      %2d chords  (margin <= %d)  [%.0fs]"
                  % (k, n // 2 - k, time.time() - t0))
            sys.stdout.flush()
    return best, bestmate

if __name__ == '__main__':
    sys.setrecursionlimit(100000)
    rng = random.Random(int(sys.argv[1]))
    secs = float(sys.argv[2])
    for n in [int(x) for x in sys.argv[3:]]:
        L = bigpow(n)
        if L == n:
            print("  2N=%2d trivial" % n); continue
        print("  2N=%2d  N=%2d  L=%2d" % (n, n // 2, L))
        b, w = ils(n, L, rng, secs)
        print("  ==> best %d chords,  MARGIN <= %d" % (b, n // 2 - b))
        print("      ", [(a, w[a]) for a in range(n) if w[a] > a])
        sys.stdout.flush()
