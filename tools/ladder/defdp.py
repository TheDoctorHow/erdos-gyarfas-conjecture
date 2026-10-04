"""Achievable (chord count, deficiency) pairs for a uniform ladder, by a
transfer matrix over the circle.

For U(2N,s) the chords are {2j, 2j+s}.  With a_j the indicator of chord j,
    cov(t) = sum of a_j over j with 2j <= t < 2j+s,
    H(S)   = #{t : cov(t) odd},   deficiency = H(S) - |S|.
Only the last ceil(s/2) indicators matter at any t, so the state is small.
"""
import sys

def reachable(n, s, kmax):
    """set of (k, deficiency) with k <= kmax, ignoring simplicity."""
    N = n//2
    w = (s+1)//2                       # how many chords can cover one t
    best = {}
    # brute force over the first w indicators (the wrap), then DP
    for head in range(1 << w):
        # state: tuple of the last w indicators, most recent last
        init = tuple((head >> i) & 1 for i in range(w))
        states = {init: {(sum(init), 0)}}
        for j in range(w, N):
            nxt = {}
            for st, pairs in states.items():
                for a in (0, 1):
                    ns = st[1:] + (a,)
                    # host edges t in [2j, 2j+2) get their coverage now
                    tot = 0
                    for t in (2*j, 2*j+1):
                        c = 0
                        for i in range(w):
                            jj = j - (w-1-i)
                            if 2*jj <= t < 2*jj+s: c += ns[i]
                        tot += c & 1
                    key = ns
                    acc = nxt.setdefault(key, set())
                    for (k, h) in pairs:
                        if k+a > kmax: continue
                        acc.add((k+a, h+tot))
            states = nxt
        for st, pairs in states.items():
            if st[-w:] != init[-w:] and False: continue
            for (k, h) in pairs:
                d = h-k
                if d >= 0 and (k not in best or d not in best.get(k, set())):
                    best.setdefault(k, set()).add(d)
    return best

def bigpow(n):
    L = 1
    while L*2 <= n: L *= 2
    return L if L < n else n//2

if __name__ == '__main__':
    print("U(2N,5): smallest k whose deficiency set contains D0 = 2N - L")
    print("   2N     L    D0    smallest k (arithmetic only)")
    for n in (30, 40, 50, 60, 70, 80, 100, 120):
        s = 5
        if not (4 <= s <= n//2 - 2): continue
        L = bigpow(n); D0 = n-L
        R = reachable(n, s, 8)
        k = None
        for kk in sorted(R):
            if D0 in R[kk]: k = kk; break
        print("  %4d  %4d  %4d      %s" % (n, L, D0, k))
        sys.stdout.flush()
