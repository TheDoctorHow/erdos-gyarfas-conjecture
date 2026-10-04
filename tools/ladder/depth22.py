"""Chain depth at 2N = 22, streaming: enumerate same-sign avoiders with
the relation prune and classify each on the fly (popcount-limited mask
scan), keeping only the deep ones."""
import sys, time
from itertools import combinations
from collections import Counter
sys.setrecursionlimit(100000)
from pushL import short
from descend2 import build_check

def relscan(m, n):
    for x0 in range(n):
        y0 = m.get(x0)
        if y0 is None:
            continue
        for d, dp in ((2, 2), (3, 1), (1, 3)):
            y1 = m.get((x0 + d) % n)
            if y1 is not None and y1 == (y0 + dp) % n:
                return True
    return False

def has_k(n, ch, k, target):
    for comb in combinations(range(len(ch)), k):
        S = [ch[i] for i in comb]
        es = sorted(x for c in S for x in c)
        if len(set(es)) < 2 * k:
            continue
        A = sum(es[j + 1] - es[j] for j in range(0, len(es), 2))
        for eps, ln in ((0, k + A), (1, k + n - A)):
            if ln == target and build_check(n, ch, S, eps, target):
                return True
    return False

def run(n):
    adj = [[(v + 1) % n, (v - 1) % n] for v in range(n)]
    mate = {}
    dist = Counter()
    deep = []
    cnt = [0]
    t0 = time.time()

    def classify():
        ch = sorted((x, y) for x, y in mate.items() if x < y)
        for k in (3, 4, 5, 6):
            if has_k(n, ch, k, n - 2):
                dist[k] += 1
                if k >= 5:
                    deep.append(ch)
                    print('  DEPTH-%d LADDER: %s' % (k, ch))
                    sys.stdout.flush()
                return
        dist['>6'] += 1
        deep.append(ch)
        print('  DEPTH->6 LADDER: %s' % (ch,))
        sys.stdout.flush()

    def rec():
        a = next((v for v in range(n) if v not in mate), None)
        if a is None:
            cnt[0] += 1
            classify()
            if cnt[0] % 100000 == 0:
                print('  ... %d avoiders  %s  [%.0fs]'
                      % (cnt[0], dict(sorted(dist.items(), key=str)),
                         time.time() - t0))
                sys.stdout.flush()
            return
        for b in range(a + 1, n):
            if b in mate or short(adj, n, a, b):
                continue
            if a == 0 and b > n // 2:
                break
            mate[a] = b; mate[b] = a
            if not relscan(mate, n):
                adj[a].append(b); adj[b].append(a)
                rec()
                adj[a].pop(); adj[b].pop()
            del mate[a]; del mate[b]
    rec()
    print('avoiders at 2N=%d: %d   min-k distribution %s   [%.0fs]'
          % (n, cnt[0], dict(sorted(dist.items(), key=str)), time.time() - t0))
    print('CHAIN DEPTH at 2N=%d = %s'
          % (n, max((k for k in dist if isinstance(k, int)), default=None)
             if not deep else 'SEE DEEP LIST'))

if __name__ == '__main__':
    run(int(sys.argv[1]) if len(sys.argv) > 1 else 22)
