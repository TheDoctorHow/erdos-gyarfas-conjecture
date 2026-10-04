"""The repaired scaling law.

573 proved  window >= sigma_min , so no bound in g alone can hold.
The right parameter is therefore the SPAN.  Measure window / sigma.
"""
import sys, random, itertools
from gwind import lift, girth, config
from scale import pows

def spans_of(mate, n):
    isleft, span = lift(mate, n)
    return [span[r] for r in range(n) if isleft[r]]

def best_window(mate, n, kmax=4):
    """over winding-zero configurations of power-of-two length, the one with
    the smallest window.  Returns (window, k, L)."""
    isleft, span = lift(mate, n)
    R = [(x, span[x % n]) for x in range(2*n) if isleft[x % n]]
    P = set(pows(n))
    best = None
    for k in range(1, kmax+1):
        for idx in itertools.combinations(range(len(R)), k):
            xs = [R[i][0] for i in idx]
            if xs[-1] - xs[0] > n: continue
            sp = [R[i][1] for i in idx]
            L, ok, verts = config(xs, sp)
            if not ok or L not in P: continue
            if len(set(v % n for v in verts)) != len(verts): continue
            w = max(verts) - min(verts)
            if best is None or w < best[0]: best = (w, k, L)
    return best

def row(mate, n, tag):
    sp = spans_of(mate, n)
    smin, smax = min(sp), max(sp)
    b = best_window(mate, n)
    if b is None:
        print("  %-24s 2N=%3d g=%d  sigma %2d..%2d : none with k<=4"
              % (tag, n, girth(mate, n), smin, smax))
        return None
    w, k, L = b
    print("  %-24s 2N=%3d g=%d  sigma %2d..%2d :  k=%d L=%2d window=%3d"
          "   w/smin=%.2f  w/smax=%.2f"
          % (tag, n, girth(mate, n), smin, smax, k, L, w, w/smin, w/smax))
    return w/smin, w/smax, k

if __name__ == '__main__':
    from skip3 import random_ladder
    from highgirth import build
    from scale2 import ladder_with_spans
    from skip4 import uniform
    rng = random.Random(90210)
    rows = []
    print("== window measured against the span, not the girth ==")
    for n in (20, 30, 40, 50, 60, 70):
        m = random_ladder(n, rng)
        if m and girth(m, n) >= 5: rows.append(row(m, n, "random"))
    for n, tg in ((40, 7), (80, 8), (140, 9)):
        m = build(n, rng, tg, tries=40)
        if m: rows.append(row(m, n, "high girth"))
    P2M1 = set(2**j-1 for j in range(2, 9))
    for n in (32, 48, 60, 80):
        allowed = set(x for x in range(4, n//2+1) if x not in P2M1)
        sub = set(rng.sample(sorted(allowed), 3))
        m = ladder_with_spans(n, rng, sub)
        if m: rows.append(row(m, n, "no span 2^j-1"))
    for n, s in ((30, 11), (40, 17), (60, 25), (92, 41), (120, 55)):
        rows.append(row(uniform(n, s), n, "uniform s=%d" % s))
    rows = [r for r in rows if r]
    print()
    print("  ladders %d   max window/sigma_min %.2f   max window/sigma_max %.2f"
          "   max rungs %d"
          % (len(rows), max(r[0] for r in rows), max(r[1] for r in rows),
             max(r[2] for r in rows)))
