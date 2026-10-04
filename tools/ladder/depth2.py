"""EXACT maximum number of chords in a girth-legal partial matching of the
2N-circle carrying NO cycle of length exactly L.   margin = N - that.

Corrects depth.py, which always matched the least unmatched vertex and so
could only reach matchings in which every chord's lower endpoint precedes
every unmatched vertex.  Here each vertex is either matched or SKIPPED.

Normalisation: any nonempty partial matching that is not perfect can be
rotated so that vertex 0 is matched and vertex n-1 is skipped.  (Perfect
matchings are the enumeration of 591/593 and are excluded separately.)

Bound: with `free` vertices still undecided, at most depth + free//2
chords can ever be reached, so a branch dies as soon as that is <= the
incumbent.  Near the top of the range the incumbent leaves room for only
a handful of skips, which is what makes the search finish.
"""
import sys, time, random
from pushL import short, has_len_through, bigpow

def max_chords(n, L, seed_best=0, verbose=False):
    adj = [[(v + 1) % n, (v - 1) % n] for v in range(n)]
    mate = [-1] * n
    state = [0] * n          # 0 undecided, 1 matched, 2 skipped
    best = [seed_best]
    wit = [None]
    nodes = [0]
    state[n - 1] = 2         # normalisation: n-1 skipped

    def rec(v, k, free):
        nodes[0] += 1
        if k > best[0]:
            best[0] = k
            wit[0] = list(mate)
            if verbose:
                print("      depth %d" % k, [(a, mate[a]) for a in range(n)
                                             if mate[a] > a])
                sys.stdout.flush()
        if k + free // 2 <= best[0]:
            return
        while v < n and state[v] != 0:
            v += 1
        if v >= n:
            return
        # branch 1: match v
        for w in range(v + 1, n):
            if state[w] != 0:
                continue
            if short(adj, n, v, w):
                continue
            # no reflection cut here: the rotation normalisation already
            # pins vertex n-1, and the reflection about 0 moves it.
            adj[v].append(w); adj[w].append(v)
            mate[v] = w; mate[w] = v
            state[v] = 1; state[w] = 1
            if not has_len_through(adj, n, v, w, L):
                rec(v + 1, k + 1, free - 2)
            adj[v].pop(); adj[w].pop()
            mate[v] = -1; mate[w] = -1
            state[v] = 0; state[w] = 0
        # branch 2: skip v   (not for v == 0, by the rotation normalisation)
        if v != 0:
            state[v] = 2
            rec(v + 1, k, free - 1)
            state[v] = 0

    rec(0, 0, n - 1)
    return best[0], wit[0], nodes[0]

def greedy(n, L, tries, rng):
    """random maximal girth-legal partial matchings with no L-cycle;
    returns the best size found (a certified lower bound)."""
    best = 0; wit = None
    for _ in range(tries):
        adj = [[(v + 1) % n, (v - 1) % n] for v in range(n)]
        mate = [-1] * n
        order = list(range(n)); rng.shuffle(order)
        k = 0
        for v in order:
            if mate[v] != -1:
                continue
            cand = [w for w in range(n) if w != v and mate[w] == -1
                    and not short(adj, n, v, w)]
            rng.shuffle(cand)
            for w in cand:
                adj[v].append(w); adj[w].append(v)
                if has_len_through(adj, n, v, w, L):
                    adj[v].pop(); adj[w].pop()
                    continue
                mate[v] = w; mate[w] = v; k += 1
                break
        if k > best:
            best = k; wit = list(mate)
    return best, wit

if __name__ == '__main__':
    sys.setrecursionlimit(100000)
    rng = random.Random(20260811)
    for n in [int(x) for x in sys.argv[1:]]:
        L = bigpow(n)
        if L == n:
            print("  2N=%2d  trivial (L = 2N)" % n); continue
        t = time.time()
        lo, lw = greedy(n, L, 300, rng)
        b, w, nd = max_chords(n, L, seed_best=lo)
        if w is None:
            w = lw
        print("  2N=%2d  L=%2d   max chords = %2d   MARGIN = %d"
              "   [greedy %d, nodes %d, %.0fs]"
              % (n, L, b, n // 2 - b, lo, nd, time.time() - t))
        print("        witness", [(a, w[a]) for a in range(n) if w[a] > a])
        sys.stdout.flush()
