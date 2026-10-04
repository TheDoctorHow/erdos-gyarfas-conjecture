"""For each uniform ladder, find how L is actually realised in the one-lap
family: minimal k, and the overlap profile."""
import sys
from skip4 import uniform, biggest_pow2

def realising(mate, n, target, oddonly=True):
    """return (k, sorted b-profile) for a one-lap cycle of length target,
    minimal k, or None."""
    d = [(mate[v]-v) % n for v in range(n)]
    tails = [v for v in range(n) if v % 2 == 0] if oddonly else list(range(n))
    idx = {v: i for i, v in enumerate(tails)}
    best = None
    for ai in range(len(tails)):
        x1 = tails[ai]
        if n - d[x1] + 1 == target:
            cand = (1, [0])
            if best is None or cand[0] < best[0]: best = cand
        for bi in range(ai+1, len(tails)):
            x2 = tails[bi]
            g = x2 - x1
            if g == d[x1]: continue
            b1 = max(0, d[x1]-g)
            stack = [(x1, x2, b1, d[x1]-1-2*b1, [b1])]
            while stack:
                p, c, bp, tot, prof = stack.pop()
                dc = d[c]
                gcl = n - c + x1
                if gcl != dc:
                    bcl = max(0, dc-gcl)
                    if bp + bcl < dc and bcl + b1 < d[x1]:
                        if n - (tot + dc-1-2*bcl) == target:
                            cand = (len(prof)+1, sorted(prof+[bcl]))
                            if best is None or cand[0] < best[0]: best = cand
                if best is not None and len(prof)+1 >= best[0]: continue
                for ci in range(idx[c]+1, len(tails)):
                    nx = tails[ci]
                    g2 = nx - c
                    if g2 == dc: continue
                    bc = max(0, dc-g2)
                    if bp + bc >= dc: continue
                    stack.append((c, nx, bc, tot + dc-1-2*bc, prof+[bc]))
    return best

if __name__ == '__main__':
    print(" s   2N    L    D   k0   realised (k, overlap profile)")
    miss = 0
    for n in range(12, 45, 2):
        N = n//2
        L = biggest_pow2(n)
        D = n - L
        for s in range(5, N-1, 2):
            mate = uniform(n, s)
            r = realising(mate, n, L)
            k0 = -((-D)//(s-1)) if D > 0 else 0
            flag = ''
            if r is None:
                r2 = realising(mate, n, L, oddonly=False)
                flag = '   <-- odd-tail MISS, both-dir: %s' % (r2,)
                if r2 is None:
                    flag = '   <-- MISS ENTIRELY'
                    miss += 1
            print("%2d %4d %4d %4d %4d   %s%s" % (s, n, L, D, k0, r, flag))
        sys.stdout.flush()
    print("total complete misses:", miss)
