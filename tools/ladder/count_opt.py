"""How many maximum avoiding configurations are there?  Counts the
girth-legal partial matchings with exactly `target` chords and no
L-cycle, in the normalised frame (vertex 0 matched, vertex n-1 skipped).
"""
import sys, time
from pushL import short, has_len_through, bigpow

def count_at(n, L, target):
    adj = [[(v + 1) % n, (v - 1) % n] for v in range(n)]
    mate = [-1] * n
    state = [0] * n
    state[n - 1] = 2
    cnt = [0]; sample = []

    def rec(v, k, free):
        if k == target:
            cnt[0] += 1
            if len(sample) < 6:
                sample.append([(a, mate[a]) for a in range(n) if mate[a] > a])
            return
        if k + free // 2 < target:
            return
        while v < n and state[v] != 0:
            v += 1
        if v >= n:
            return
        for w in range(v + 1, n):
            if state[w] != 0 or short(adj, n, v, w):
                continue
            adj[v].append(w); adj[w].append(v)
            mate[v] = w; mate[w] = v
            state[v] = 1; state[w] = 1
            if not has_len_through(adj, n, v, w, L):
                rec(v + 1, k + 1, free - 2)
            adj[v].pop(); adj[w].pop()
            mate[v] = -1; mate[w] = -1
            state[v] = 0; state[w] = 0
        if v != 0:
            state[v] = 2
            rec(v + 1, k, free - 1)
            state[v] = 0

    rec(0, 0, n - 1)
    return cnt[0], sample

if __name__ == '__main__':
    n = int(sys.argv[1]); target = int(sys.argv[2])
    L = bigpow(n)
    t = time.time()
    c, s = count_at(n, L, target)
    print("  2N=%d  L=%d  configurations with exactly %d chords and no "
          "L-cycle (normalised frame): %d   [%.0fs]"
          % (n, L, target, c, time.time() - t))
    for x in s:
        print("      ", x)
