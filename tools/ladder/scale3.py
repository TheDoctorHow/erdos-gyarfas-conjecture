"""How many rungs, and how wide a window, does a winding-zero power-of-two
cycle actually need?"""
import sys, random, itertools
from gwind import lift, girth, config
from scale import pows

def min_config(mate, n, kmax=4):
    """smallest k, then smallest window, over winding-zero configurations
    whose length is a power of two.  Returns (k, window, length, S)."""
    isleft, span = lift(mate, n)
    R = [(x, span[x % n]) for x in range(2*n) if isleft[x % n]]
    P = set(pows(n))
    best = None
    for k in range(1, kmax+1):
        found = None
        for idx in itertools.combinations(range(len(R)), k):
            xs = [R[i][0] for i in idx]
            if xs[-1] - xs[0] > n: continue
            sp = [R[i][1] for i in idx]
            L, ok, verts = config(xs, sp)
            if not ok or L not in P: continue
            if len(set(v % n for v in verts)) != len(verts): continue
            w = max(verts) - min(verts)
            if found is None or w < found[1]: found = (k, w, L, xs)
        if found: return found
    return None

def summarise(m, n, tag):
    g = girth(m, n)
    r = min_config(m, n)
    if r is None:
        print("  %-24s 2N=%3d g=%d :  NO winding-zero power of two with k<=4"
              % (tag, n, g))
        return None
    k, w, L, xs = r
    print("  %-24s 2N=%3d g=%d :  k=%d  window=%3d (%.2f of 2N)  length=%d"
          % (tag, n, g, k, w, w/n, L))
    return k, w, L

if __name__ == '__main__':
    from skip3 import random_ladder
    from highgirth import build
    from scale2 import ladder_with_spans
    from skip4 import uniform
    rng = random.Random(24601)
    print("== minimum rungs and window for a winding-zero power of two ==")
    rows = []
    for n in (20, 30, 40, 50, 60, 70):
        m = random_ladder(n, rng)
        if m and girth(m, n) >= 5:
            r = summarise(m, n, "random"); rows.append(('random', n, r))
        sys.stdout.flush()
    print()
    for n, tg in ((40, 7), (60, 7), (80, 8), (100, 8), (140, 9)):
        m = build(n, rng, tg, tries=40)
        if m:
            r = summarise(m, n, "high girth"); rows.append(('highg', n, r))
        sys.stdout.flush()
    print()
    for n in (24, 32, 40, 60, 80):
        for k in (2, 3):
            allowed = set(rng.sample(range(4, n//2+1), k))
            m = ladder_with_spans(n, rng, allowed)
            if m:
                r = summarise(m, n, "%d spans %s" % (k, sorted(allowed)))
                rows.append(('fewspan', n, r))
        sys.stdout.flush()
    print()
    for n, s in ((30, 11), (40, 17), (60, 25), (92, 41)):
        m = uniform(n, s)
        r = summarise(m, n, "uniform s=%d" % s); rows.append(('uniform', n, r))
        sys.stdout.flush()
    print()
    ks = [r[2][0] for r in rows if r[2]]
    ws = [r[2][1]/r[1] for r in rows if r[2]]
    print("ladders settled: %d of %d   rung counts used: %s   "
          "window as a fraction of 2N: min %.2f  max %.2f"
          % (len(ks), len(rows), sorted(set(ks)), min(ws), max(ws)))
