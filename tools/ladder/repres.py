"""Is the deficiency representable by PAIRWISE DISJOINT chords?

Arc replacement: shortcut a family of pairwise disjoint chord-arcs.  The
resulting cycle has length 2N - sum(sigma_i - 1), so the question is a
subset-sum over the savings with an interval-packing side condition.
"""
import sys, random

def arcs_of(mate, n):
    """each chord as its SHORT arc (a, b) with b - a = sigma <= n/2,
    given in the cover so that b may exceed n."""
    out = []
    for u in range(n):
        v = mate[u]
        if u >= v: continue
        d = v-u
        if d <= n-d: out.append((u, u+d))
        else: out.append((v, v+(n-d)))
    return out

def disjoint_sums(mate, n, cap=None):
    """set of sums of (sigma-1) over pairwise disjoint arc families."""
    A = arcs_of(mate, n)
    if cap is None: cap = n
    reach = set([0])
    for c in range(n):
        # linearise at cut c : keep arcs whose interior avoids c
        items = []
        for (a, b) in A:
            aa, bb = a % n, (a % n) + (b-a)
            # arc occupies [aa, bb] in the cover; shift into [c, c+n)
            start = aa if aa >= c else aa + n
            if start + (b-a) > c + n: continue
            items.append((start, start + (b-a), b-a-1))
        items.sort(key=lambda t: t[1])
        # interval DP: dp[j] = set of sums using arcs among the first j
        dp = [set([0])]
        ends = []
        for (s0, e0, w) in items:
            best = set()
            # arcs ending at or before s0
            j = 0
            for idx, (s1, e1, w1) in enumerate(items):
                if e1 <= s0: j = idx+1
                else: break
            best |= set(x+w for x in dp[j] if x+w <= cap)
            dp.append(dp[-1] | best)
            ends.append(e0)
        reach |= dp[-1]
    return reach

def bigpow(n):
    L = 1
    while L*2 <= n: L *= 2
    return L if L < n else n//2

if __name__ == '__main__':
    from skip3 import random_ladder, girth
    from skip4 import uniform
    rng = random.Random(2024)
    print("== is D0 = 2N - L a disjoint subset-sum of the savings? ==")
    ok = 0; tot = 0; fails = []
    for n in range(16, 49, 2):
        for _ in range(4):
            m = random_ladder(n, rng)
            if m is None or girth(m, n) < 5: continue
            L = bigpow(n); D0 = n-L
            R = disjoint_sums(m, n)
            tot += 1
            if D0 in R: ok += 1
            else: fails.append((n, D0, sorted(x for x in R if x <= D0+6)[-6:]))
    print("  random ladders: %d of %d representable by disjoint chords" % (ok, tot))
    for f in fails[:6]:
        print("     2N=%d  D0=%d  nearby reachable sums %s" % f)
    print()
    print("== uniform ladders: the prediction is (s-1) | D0 ==")
    for n, s in ((30, 11), (40, 17), (48, 21), (60, 25), (72, 31), (92, 41)):
        m = uniform(n, s); L = bigpow(n); D0 = n-L
        R = disjoint_sums(m, n)
        print("  U(%3d,%2d)  D0=%3d  (s-1)=%2d  divides: %-5s   D0 reachable: %s"
              % (n, s, D0, s-1, D0 % (s-1) == 0, D0 in R))
        sys.stdout.flush()
