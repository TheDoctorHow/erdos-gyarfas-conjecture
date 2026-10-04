"""The long-cycle dual of (WIND0-GEN).

In the circle the cycle space needs N+1 generators: the fundamental cycle
of each chord (w.r.t. the host path 0..2N-1) and the host cycle itself.
So an element is a pair (S, eps).  With H(S) the number of host edges
whose chord-coverage is odd,

    eps = 0 :  length = |S| + H(S)          (short cycles)
    eps = 1 :  length = |S| + 2N - H(S)     (long cycles)

so the two elements from one S have lengths summing to 2N + 2|S|, and

    DEFICIENCY  2N - length(S,1)  =  H(S) - |S|  =  sum (l_i - 1)

over the |S| pairing intervals of S.
"""
import sys, random, itertools

def chords_of(mate, n):
    return [(u, mate[u]) for u in range(n) if u < mate[u]]

def element(S, n, eps):
    """edge set of the cycle-space element (S, eps); S a list of (u,v), u<v."""
    cov = [0]*n
    for (u, v) in S:
        for t in range(u, v): cov[t] ^= 1
    host = [t for t in range(n) if (cov[t] == 1) != (eps == 1)]
    edges = [(u, v) for (u, v) in S] + [(t, (t+1) % n) for t in host]
    return edges

def is_cycle(edges, n):
    if not edges: return None
    nb = {}
    for u, v in edges:
        nb.setdefault(u, []).append(v); nb.setdefault(v, []).append(u)
    if any(len(x) != 2 for x in nb.values()): return None
    if len(nb) != len(edges): return None
    start = edges[0][0]; prev = None; cur = start; cnt = 0
    while True:
        w = nb[cur][0] if nb[cur][0] != prev else nb[cur][1]
        prev, cur = cur, w; cnt += 1
        if cur == start: break
    return cnt if cnt == len(edges) else None

def H(S, n):
    cov = [0]*n
    for (u, v) in S:
        for t in range(u, v): cov[t] ^= 1
    return sum(cov)

if __name__ == '__main__':
    from skip3 import random_ladder, girth
    rng = random.Random(8128)
    print("== the duality length(S,0) + length(S,1) = 2N + 2|S| ==")
    bad = 0; tested = 0; both = 0
    for n in (16, 20, 24):
        for _ in range(3):
            m = random_ladder(n, rng)
            if m is None or girth(m, n) < 5: continue
            C = chords_of(m, n)
            for k in range(1, 4):
                for S in itertools.combinations(C, k):
                    S = list(S)
                    h = H(S, n)
                    L0 = is_cycle(element(S, n, 0), n)
                    L1 = is_cycle(element(S, n, 1), n)
                    tested += 1
                    if L0 is not None and L0 != k+h: bad += 1
                    if L1 is not None and L1 != k+n-h: bad += 1
                    if L0 is not None and L1 is not None:
                        both += 1
                        if L0+L1 != n+2*k: bad += 1
    print("  subsets tested %d   both-are-cycles %d   formula failures %d"
          % (tested, both, bad))

    print()
    print("== which mechanisms give LONG cycles (eps = 1) ==")
    for n in (24, 30):
        m = random_ladder(n, rng)
        if m is None or girth(m, n) < 5: continue
        C = chords_of(m, n)
        best = {}
        for k in range(1, 5):
            for S in itertools.combinations(C, k):
                L1 = is_cycle(element(list(S), n, 1), n)
                if L1 is None: continue
                if L1 not in best or k < best[L1]: best[L1] = k
        P = [L for L in sorted(best) if L >= 4 and (L & (L-1)) == 0]
        print("  2N=%2d : long lengths %s" % (n, sorted(best)))
        print("          powers of two reachable with eps=1: %s   (chords: %s)"
              % (P, {L: best[L] for L in P}))
        sys.stdout.flush()
