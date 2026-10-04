"""Does the one-lap construction really cover the uniform ladders?
Three independent computations:
  A. exact one-lap spectrum by DFS on the graph          (ground truth)
  B. the CORRECTED gap conditions (m = k forces P = 0, plus adjacency)
  C. 567's prescription, as recorded
"""
import sys
from skip3 import onelap_lengths

def uniform(n, s):
    mate = [None]*n
    for i in range(1, n+1, 2):
        a = i-1; b = (i+s-1) % n
        mate[a] = b; mate[b] = a
    return mate

def biggest_pow2(n):
    L = 1
    while L*2 <= n: L *= 2
    return L

# ---------- B : corrected conditions ----------

def Mmax(m, f, s):
    """largest total of m odd values in [1,s-2], split by f forward runs into
    groups on a path/cycle, with every adjacent pair summing to at most s-1."""
    if m == 0: return 0
    if f == 0:                      # one cyclic group of size m
        big = m//2                  # cannot exceed half, cyclically
    else:
        big = min(m, (m+f)//2)      # all groups odd-sized
    return big*(s-2) + (m-big)

def realisable(s, n, L):
    """is L reachable by the one-lap construction in U(n,s)?  corrected."""
    D = n - L
    for k in range(1, L):
        num = k*(s-1) - D
        if num < 0: continue
        if num % 2: continue
        M = num//2
        P = L - k - M
        if P < 0: continue
        for m in range(0, k+1):
            if (M - m) % 2: continue
            if M < m: continue
            f = k - m
            if f == 0:
                if P != 0: continue
            else:
                if P < f: continue
            if M > Mmax(m, f, s): continue
            return (k, m, M, P)
    return None

# ---------- C : 567's prescription, as recorded ----------

def prescription(s, n, L):
    D = n - L
    k0 = -((-D)//(s-1))
    if k0 < 1: k0 = 1
    for k in (k0, k0+1):
        M = (k*(s-1) - D)//2
        if (k*(s-1) - D) % 2 or M < 0: continue
        P = L - k - M
        if P < 0: continue
        if M <= k: m = M
        elif (k - M) % 2 == 0: m = k
        else: m = k-1
        if m < 0: continue
        return (k, m, M, P)
    return None

if __name__ == '__main__':
    print("== A vs B vs C on small uniform ladders ==")
    print("  s  2N    L   exact-one-lap?  corrected?  567-prescription")
    badAB = 0; presc_bad = 0; rows = 0
    for n in range(12, 49, 2):
        N = n//2
        L = biggest_pow2(n)
        for s in range(5, N-1, 2):
            mate = uniform(n, s)
            exact = L in onelap_lengths(mate, n)
            corr = realisable(s, n, L)
            pres = prescription(s, n, L)
            rows += 1
            if exact != (corr is not None):
                badAB += 1
                print("  MISMATCH A/B  s=%d 2N=%d L=%d exact=%s corr=%s"
                      % (s, n, L, exact, corr))
            if pres is not None:
                k, m, M, P = pres
                f = k - m
                ok = (M >= m and (M-m) % 2 == 0 and M <= Mmax(m, f, s)
                      and (P == 0 if f == 0 else P >= f))
                if not ok:
                    presc_bad += 1
                    if presc_bad <= 6:
                        print("  567 prescription UNREALISABLE  s=%2d 2N=%3d "
                              "k=%d m=%d M=%d P=%d  (exact says %s)"
                              % (s, n, k, m, M, P, exact))
    print("  rows %d   A/B mismatches %d   prescriptions unrealisable %d"
          % (rows, badAB, presc_bad))
