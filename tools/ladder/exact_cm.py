"""EXACT completable margin: enumerate every girth-legal closed ladder on
2N vertices, and for each compute the minimum number of chords whose
deletion destroys all L-cycles (minimum hitting set over the L-cycles,
found by the (S,eps) dual).  The minimum over ladders is the completable
margin.
"""
import sys, time
from pushL import short, bigpow
from hitset import all_L_cycles, min_hitting_set

def ladders(n):
    adj = [[(v + 1) % n, (v - 1) % n] for v in range(n)]
    mate = [-1] * n
    out = []
    def rec():
        a = next((v for v in range(n) if mate[v] == -1), None)
        if a is None:
            out.append([(v, mate[v]) for v in range(n) if mate[v] > v])
            return
        for b in range(a + 1, n):
            if mate[b] != -1 or short(adj, n, a, b):
                continue
            if a == 0 and b > n // 2:
                break
            adj[a].append(b); adj[b].append(a)
            mate[a] = b; mate[b] = a
            rec()
            adj[a].pop(); adj[b].pop()
            mate[a] = -1; mate[b] = -1
    rec()
    return out

if __name__ == '__main__':
    for n in [int(x) for x in sys.argv[1:]]:
        L = bigpow(n)
        if L == n:
            print("  2N=%2d trivial" % n); continue
        t = time.time()
        ls = ladders(n)
        best = 10 ** 9; bw = None
        for ch in ls:
            cyc = all_L_cycles(n, ch, L)
            if not cyc:
                best = 0; bw = ch; break
            h = min_hitting_set(cyc, best)
            if h < best:
                best = h; bw = ch
        print("  2N=%2d  N=%2d  L=%2d   %d girth-legal ladders (up to the "
              "reflection cut)   COMPLETABLE MARGIN = %d   [%.0fs]"
              % (n, n // 2, L, len(ls), best, time.time() - t))
        print("        extremal ladder:", bw)
        sys.stdout.flush()
