"""How much of the non-backtracking trace is simple cycles?

For g <= m < 2g,  tr(B^m) = 2 m c_m  exactly (452, 555.2).  Above 2g the
trace counts closed non-backtracking walks that are not simple.  The
question is whether the simple part still dominates at m of order 2N,
where the target power of two lives.
"""
import sys

def hashimoto(mate, n):
    """directed edges and the non-backtracking matrix B."""
    E = []
    for v in range(n):
        for w in ((v+1) % n, (v-1) % n, mate[v]):
            E.append((v, w))
    idx = {e: i for i, e in enumerate(E)}
    m = len(E)
    B = [[0]*m for _ in range(m)]
    for i, (u, v) in enumerate(E):
        for w in ((v+1) % n, (v-1) % n, mate[v]):
            if w == u: continue          # no backtracking
            B[i][idx[(v, w)]] = 1
    return B, m

def trace_powers(B, m, upto):
    """tr(B^k) for k = 1..upto, exact integers."""
    cur = [row[:] for row in B]
    out = [None, sum(cur[i][i] for i in range(m))]
    for k in range(2, upto+1):
        nxt = [[0]*m for _ in range(m)]
        for i in range(m):
            Bi = B[i]
            ci = nxt[i]
            for t in range(m):
                if Bi[t]:
                    ct = cur[t]
                    for j in range(m):
                        if ct[j]: ci[j] += ct[j]
        cur = nxt
        out.append(sum(cur[i][i] for i in range(m)))
    return out

def cycle_counts(mate, n, upto):
    A = [[(v+1) % n, (v-1) % n, mate[v]] for v in range(n)]
    cnt = [0]*(upto+1)
    def rec(start, u, vis, d):
        if d > upto: return
        for w in A[u]:
            if w == start and d+1 >= 3:
                if d+1 <= upto: cnt[d+1] += 1
            elif w > start and not (vis >> w) & 1:
                rec(start, w, vis | (1 << w), d+1)
    for s0 in range(n):
        rec(s0, s0, 1 << s0, 0)
    return [c//2 for c in cnt]

if __name__ == '__main__':
    from skip3 import random_ladder, girth
    import random
    rng = random.Random(451)
    for n in (14, 16, 18):
        mate = None
        while mate is None or girth(mate, n) < 5:
            mate = random_ladder(n, rng)
        g = girth(mate, n)
        B, sz = hashimoto(mate, n)
        tr = trace_powers(B, sz, n)
        cc = cycle_counts(mate, n, n)
        print("2N=%d  girth %d   (2g = %d)" % (n, g, 2*g))
        print("   m   tr(B^m)        2m*c_m      simple share")
        for m in range(g, n+1):
            simple = 2*m*cc[m]
            share = simple/tr[m] if tr[m] else 0
            flag = ' <- m < 2g' if m < 2*g else ''
            print("  %2d  %12d  %12d      %.4f%s" % (m, tr[m], simple, share, flag))
        print()
        sys.stdout.flush()
