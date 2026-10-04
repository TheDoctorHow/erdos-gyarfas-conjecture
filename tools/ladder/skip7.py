"""The corrected closed form for the short-direction one-lap family."""
import sys, itertools
from skip4 import uniform, biggest_pow2
from skip5 import onelap_odd

def capacity(k, m, s):
    """max total overlap with m positive odd overlaps (each <= s-2) among k
    cyclic slots, every adjacent pair summing to at most s-1."""
    if m == 0: return 0
    f = k - m
    if f == 0:
        base = m*(s-1)//2
        if m % 2 == 1 and ((s-1)//2) % 2 == 0: base -= 1
        return base
    odd = min(f, m)
    if (odd - m) % 2: odd -= 1
    return (m*(s-1) + (s-3)*odd)//2

def brute_capacity(k, m, s):
    """exhaustive: max and achievable set of sums."""
    best = -1; ach = set()
    vals = [v for v in range(1, s-1, 2)]
    for pos in itertools.combinations(range(k), m):
        posset = set(pos)
        for combo in itertools.product(vals, repeat=m):
            b = [0]*k
            for p, v in zip(pos, combo): b[p] = v
            ok = all(b[j] + b[(j+1) % k] < s for j in range(k))
            if ok:
                t = sum(combo); ach.add(t); best = max(best, t)
    return best, ach

def realisable(s, n, L):
    """is L reachable by a short-direction one-lap cycle of U(n,s)?"""
    D = n - L
    for k in range(1, L):
        num = k*(s-1) - D
        if num < 0 or num % 2: continue
        M = num//2
        P = L - k - M
        if P < 0: continue
        for m in range(0, k+1):
            if (M - m) % 2 or M < m: continue
            f = k - m
            if f == 0:
                if P != 0: continue
            elif P < f: continue
            if M > capacity(k, m, s): continue
            return (k, m, M, P)
    return None

if __name__ == '__main__':
    print("== capacity formula vs brute force ==")
    bad = 0; gaps = 0
    for s in range(5, 14, 2):
        for k in range(1, 8):
            for m in range(1, k+1):
                bf, ach = brute_capacity(k, m, s)
                cf = capacity(k, m, s)
                if bf != cf:
                    bad += 1
                    if bad <= 10:
                        print("  CAP MISMATCH s=%d k=%d m=%d brute=%d formula=%d"
                              % (s, k, m, bf, cf))
                want = set(range(m, bf+1, 2))
                if ach != want:
                    gaps += 1
                    if gaps <= 6:
                        print("  GAP s=%d k=%d m=%d achievable=%s expected=%s"
                              % (s, k, m, sorted(ach), sorted(want)))
    print("  capacity mismatches %d   interval gaps %d" % (bad, gaps))

    print()
    print("== closed form vs DFS on uniform ladders ==")
    bad = 0; rows = 0
    for n in range(12, 45, 2):
        N = n//2
        for s in range(5, N-1, 2):
            mate = uniform(n, s)
            dfs = onelap_odd(mate, n)
            cf = set(L for L in range(3, n+1) if realisable(s, n, L))
            rows += 1
            if dfs != cf:
                bad += 1
                if bad <= 10:
                    print("  s=%2d 2N=%3d  dfs-only %s   cf-only %s"
                          % (s, n, sorted(dfs-cf), sorted(cf-dfs)))
    print("  rows %d  disagreements %d" % (rows, bad))
