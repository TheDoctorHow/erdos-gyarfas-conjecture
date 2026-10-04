"""Validate the closed form for the SHORT-DIRECTION one-lap family of a
uniform ladder against a direct DFS, then use it at scale."""
import sys
from skip4 import uniform, biggest_pow2, Mmax, realisable

def onelap_odd(mate, n):
    """one-lap simple cycles whose chord tails are all ODD (0-based: even
    index), i.e. every chord taken in its short direction."""
    d = [(mate[v]-v) % n for v in range(n)]
    tails = [v for v in range(n) if v % 2 == 0]
    lens = set()
    for ai in range(len(tails)):
        x1 = tails[ai]
        lens.add(n - d[x1] + 1)
        for bi in range(ai+1, len(tails)):
            x2 = tails[bi]
            g = x2 - x1
            if g == d[x1]: continue
            b1 = max(0, d[x1]-g)
            t1 = d[x1]-1-2*b1
            stack = [(x1, x2, b1, t1)]
            while stack:
                p, c, bp, tot = stack.pop()
                dc = d[c]
                for ci in range(tails.index(c)+1, len(tails)):
                    nx = tails[ci]
                    g2 = nx - c
                    if g2 == dc: continue
                    bc = max(0, dc-g2)
                    if bp + bc >= dc: continue
                    stack.append((c, nx, bc, tot + dc-1-2*bc))
                gcl = n - c + x1
                if gcl != dc:
                    bcl = max(0, dc-gcl)
                    if bp + bcl < dc and bcl + b1 < d[x1]:
                        lens.add(n - (tot + dc-1-2*bcl))
    return lens

def realisable_lengths(s, n):
    """every length the closed form claims, for cross-checking."""
    out = set()
    for L in range(3, n+1):
        if realisable(s, n, L): out.add(L)
    return out

if __name__ == '__main__':
    print("== closed form vs DFS, short-direction one-lap family ==")
    bad = 0; rows = 0
    for n in range(12, 41, 2):
        N = n//2
        for s in range(5, N-1, 2):
            mate = uniform(n, s)
            dfs = onelap_odd(mate, n)
            cf = realisable_lengths(s, n)
            rows += 1
            if dfs != cf:
                bad += 1
                if bad <= 8:
                    print("  s=%2d 2N=%3d  dfs-only %s   cf-only %s"
                          % (s, n, sorted(dfs-cf), sorted(cf-dfs)))
    print("  rows %d  disagreements %d" % (rows, bad))
