"""Same maximisation as depth2, but the "does an L-cycle exist" oracle is
the GF(2) cycle-space sweep of verify_wit -- no path search anywhere.
Slow, but it shares no cycle-detection code with depth2/pushL, so an
agreement is a real check on the answer.
"""
import sys, time
from pushL import short, bigpow
from verify_wit import spectrum

def max_chords_cs(n, L, seed_best=0):
    mate = [-1] * n
    state = [0] * n
    state[n - 1] = 2
    best = [seed_best]; wit = [None]; nodes = [0]

    def adjnow():
        adj = [[(v + 1) % n, (v - 1) % n] for v in range(n)]
        for v in range(n):
            if mate[v] > v:
                adj[v].append(mate[v]); adj[mate[v]].append(v)
        return adj

    def rec(v, k, free):
        nodes[0] += 1
        if k > best[0]:
            best[0] = k; wit[0] = list(mate)
        if k + free // 2 <= best[0]:
            return
        while v < n and state[v] != 0:
            v += 1
        if v >= n:
            return
        adj = adjnow()
        for w in range(v + 1, n):
            if state[w] != 0:
                continue
            if short(adj, n, v, w):
                continue
            mate[v] = w; mate[w] = v
            state[v] = 1; state[w] = 1
            ch = [(a, mate[a]) for a in range(n) if mate[a] > a]
            if L not in spectrum(n, ch):
                rec(v + 1, k + 1, free - 2)
            mate[v] = -1; mate[w] = -1
            state[v] = 0; state[w] = 0
        if v != 0:
            state[v] = 2
            rec(v + 1, k, free - 1)
            state[v] = 0

    rec(0, 0, n - 1)
    return best[0], wit[0], nodes[0]

if __name__ == '__main__':
    sys.setrecursionlimit(100000)
    for n in [int(x) for x in sys.argv[1:]]:
        L = bigpow(n)
        t = time.time()
        b, w, nd = max_chords_cs(n, L)
        print("  2N=%2d  L=%2d  max chords = %d  MARGIN = %d"
              "   [cycle-space oracle, nodes %d, %.0fs]"
              % (n, L, b, n // 2 - b, nd, time.time() - t))
        print("        witness", [(a, w[a]) for a in range(n) if w[a] > a])
        sys.stdout.flush()
