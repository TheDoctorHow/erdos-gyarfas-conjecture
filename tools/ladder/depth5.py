"""Depth-5 hunt at 2N = 20: enumerate the same-sign relation avoiders
(pruned the moment a relation appears), then classify each by the least
k among its deficiency-2 cycles, scanning masks by popcount so the
k = 3 check costs C(10,3) masks, k = 4 costs C(10,4), and so on."""
import sys, time
from itertools import combinations
from collections import Counter
sys.setrecursionlimit(100000)
from pushL import short
from descend2 import build_check

def relscan_inc(m, n):
    for x0 in range(n):
        y0 = m.get(x0)
        if y0 is None:
            continue
        for d, dp in ((2, 2), (3, 1), (1, 3)):
            y1 = m.get((x0 + d) % n)
            if y1 is not None and y1 == (y0 + dp) % n:
                return True
    return False

def avoiders(n):
    adj = [[(v + 1) % n, (v - 1) % n] for v in range(n)]
    mate = {}
    out = []
    def rec():
        a = next((v for v in range(n) if v not in mate), None)
        if a is None:
            out.append(sorted((x, y) for x, y in mate.items() if x < y))
            return
        for b in range(a + 1, n):
            if b in mate or short(adj, n, a, b):
                continue
            if a == 0 and b > n // 2:
                break
            mate[a] = b; mate[b] = a
            if not relscan_inc(mate, n):
                adj[a].append(b); adj[b].append(a)
                rec()
                adj[a].pop(); adj[b].pop()
            del mate[a]; del mate[b]
    rec()
    return out

def has_k(n, ch, k, target):
    """a simple cycle of length target using exactly k chords?"""
    idx = list(range(len(ch)))
    for comb in combinations(idx, k):
        S = [ch[i] for i in comb]
        es = sorted(x for c in S for x in c)
        if len(set(es)) < 2 * k:
            continue
        A = sum(es[j + 1] - es[j] for j in range(0, len(es), 2))
        for eps, ln in ((0, k + A), (1, k + n - A)):
            if ln != target:
                continue
            if build_check(n, ch, S, eps, target) is not None:
                return True
    return False

if __name__ == '__main__':
    n = 20
    t0 = time.time()
    av = avoiders(n)
    print('same-sign avoiders at 2N=20: %d  [%.0fs]'
          % (len(av), time.time() - t0))
    sys.stdout.flush()
    dist = Counter()
    deep = []
    t0 = time.time()
    for i, ch in enumerate(av):
        mk = None
        for k in (3, 4, 5, 6, 7):
            if has_k(n, ch, k, n - 2):
                mk = k
                break
        dist[mk] += 1
        if mk is not None and mk >= 5:
            deep.append((mk, ch))
            print('  DEPTH-%d LADDER: %s' % (mk, ch))
            sys.stdout.flush()
        if (i + 1) % 20000 == 0:
            print('  ... %d done  %s  [%.0fs]'
                  % (i + 1, dict(sorted(dist.items(), key=str)), time.time() - t0))
            sys.stdout.flush()
    print('min-k distribution over avoiders:', dict(sorted(dist.items(), key=str)))
    print('CHAIN DEPTH at 2N=20 =', max(k for k in dist if k))
