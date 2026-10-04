"""(T1) enumeration over ALL closed ladders (span >= 2, ANY girth)."""
import sys, time
from collections import deque
from dualdet2 import detect2
from t1fast import path_budget

def enumerate_missing_free(n, L):
    adj = [[(v + 1) % n, (v - 1) % n] for v in range(n)]
    mate = [-1] * n
    chords = []
    out = []
    nodes = [0]

    def rec():
        nodes[0] += 1
        a = next((v for v in range(n) if mate[v] == -1), None)
        if a is None:
            out.append(list(chords)); return
        for b in range(a + 1, n):
            if mate[b] != -1: continue
            s = min(b - a, n - (b - a))
            if s < 2: continue
            if a == 0 and b > n // 2: break
            adj[a].append(b); adj[b].append(a)
            mate[a] = b; mate[b] = a
            chords.append((a, b))
            r = path_budget(adj, n, a, b, L, 1000)
            if r is not True:
                got, _ = detect2(n, chords, L, eps_list=(1,))
                r = bool(got)
            if not r:
                rec()
            chords.pop()
            adj[a].pop(); adj[b].pop()
            mate[a] = -1; mate[b] = -1
    rec()
    return out, nodes[0]

if __name__ == '__main__':
    sys.setrecursionlimit(200000)
    for n in [int(x) for x in sys.argv[1:]]:
        t = time.time()
        out, nodes = enumerate_missing_free(n, n - 2)
        real = [m for m in out
                if not detect2(n, m, n - 2, eps_list=(0, 1))[0]]
        print("  2N=%2d  L=%2d  ANY girth: survivors %d, truly lacking = %d"
              "   [nodes %d, %.0fs]"
              % (n, n - 2, len(out), len(real), nodes, time.time() - t))
        for m in real[:2]:
            print("       ", m)
        sys.stdout.flush()
