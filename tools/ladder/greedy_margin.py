"""Certified LOWER bounds on the maximum chord count (hence UPPER bounds
on the margin), by randomised maximal partial matchings that avoid L.

Every value printed is witnessed, so margin <= N - value is a fact; the
exact margin can only be smaller or equal.
"""
import sys, time, random
from pushL import short, has_len_through, bigpow
from depth2 import greedy

if __name__ == '__main__':
    rng = random.Random(int(sys.argv[1]))
    tries = int(sys.argv[2])
    print("  2N    N    L    best found   margin <=   [s]")
    for n in [int(x) for x in sys.argv[3:]]:
        L = bigpow(n)
        if L == n:
            print("  %2d   trivial" % n); continue
        t = time.time()
        b, w = greedy(n, L, tries, rng)
        print("  %2d  %3d  %3d      %3d          %2d      [%.0f]"
              % (n, n // 2, L, b, n // 2 - b, time.time() - t))
        print("       ", [(a, w[a]) for a in range(n) if w[a] > a])
        sys.stdout.flush()
