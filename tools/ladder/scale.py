"""Does the crossing-triple formula reach a power of two when girth forces
the spans to be large?

For a pairwise-crossing triple x1<x2<x3<y1<y2<y3 the length is
    3 + a + c + e,   a = x2-x1,  c = y1-x3,  e = y3-y2,
so the question is which values a+c+e takes.
"""
import sys, random
from gwind import lift, girth, config

def rungs_in_cover(mate, n, reps=2):
    isleft, span = lift(mate, n)
    return [(x, span[x % n]) for x in range(reps*n) if isleft[x % n]]

def descends(x1, s1, x2, s2, x3, s3, n):
    verts = set()
    for a, b in ((x1, x2), (x3, x1+s1), (x2+s2, x3+s3)):
        for t in range(a, b+1): verts.add(t)
    return len(set(v % n for v in verts)) == len(verts)

def triple_lengths(mate, n, reps=2, cap=None):
    """set of lengths of pairwise-crossing triples that descend."""
    R = rungs_in_cover(mate, n, reps)
    out = {}
    m = len(R)
    for i in range(m):
        x1, s1 = R[i]; y1 = x1+s1
        for j in range(i+1, m):
            x2, s2 = R[j]
            if x2 >= y1: break
            y2 = x2+s2
            if y2 <= y1: continue
            for l in range(j+1, m):
                x3, s3 = R[l]
                if x3 >= y1: break
                y3 = x3+s3
                if y3 <= y2: continue
                L = 3 + (x2-x1) + (y1-x3) + (y3-y2)
                if cap and L > cap: continue
                if L in out: continue
                if descends(x1, s1, x2, s2, x3, s3, n):
                    out[L] = (x1, s1, x2, s2, x3, s3)
    return out

def pows(n):
    out = []; L = 4
    while L <= n: out.append(L); L *= 2
    return out

def report(mate, n, tag):
    g = girth(mate, n)
    isleft, span = lift(mate, n)
    sp = sorted(span[r] for r in range(n) if isleft[r])
    T = triple_lengths(mate, n)
    ks = sorted(T)
    hits = [L for L in pows(n) if L in T]
    gaps = max((ks[i+1]-ks[i] for i in range(len(ks)-1)), default=0)
    print("  %-22s 2N=%3d girth %2d  spans %d..%d   triple lengths: %d values "
          "in [%d,%d], max gap %d   powers of two hit: %s"
          % (tag, n, g, sp[0], sp[-1], len(ks), ks[0] if ks else -1,
             ks[-1] if ks else -1, gaps, hits))
    return g, hits

if __name__ == '__main__':
    from highgirth import build
    from skip3 import random_ladder
    rng = random.Random(5150)
    print("== crossing-triple reach ==")
    for n, tgt in [(40, 7), (60, 8), (80, 8), (100, 8), (140, 9)]:
        m = build(n, rng, tgt, tries=60)
        if m is None: print("  2N=%d girth %d: not built" % (n, tgt)); continue
        report(m, n, "high-girth g>=%d" % tgt)
        sys.stdout.flush()
    print()
    for n in (40, 60, 80):
        m = random_ladder(n, rng)
        if m and girth(m, n) >= 5: report(m, n, "random")
        sys.stdout.flush()
