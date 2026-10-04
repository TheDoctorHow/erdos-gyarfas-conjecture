"""gdav445b.py -- the G-link of the quadratic Davenport recursion.

The descent E(t) <= 2 F(t-1), F(t) <= 3 G(2^t) reduces everything
to the G-PROBLEM at modulus 4: values V_a == 2 (mod 4) (a V == 0
would die as a singleton), couplings C_ab in Z/4 FREE, subset value
Q(T) = sum V + sum C; how large can an avoider (all nonempty T have
Q(T) =/= 0 mod 4) be?  B4 := the maximum; G4 = B4 + 1 units always
suffice.  Then E(3) <= 2 F(2) <= 6 G4 is a fully proven bound.

Machine: exhaustive DFS over C with incremental subset checks,
growing k until UNSAT; prints the maximum avoider and a witness.
"""
import sys, time
sys.setrecursionlimit(10000)

def max_avoider(MOD, V):
    """Grow k; return largest k admitting an avoider, with witness."""
    best = (0, None)
    k = 1
    while True:
        # DFS over C entries c[j][i] (j<i), values in 0..MOD-1
        C = [[0] * k for _ in range(k)]
        Q = [0] * (1 << k)
        found = [False]
        def place(i):
            """elements 0..i-1 placed with Q filled for masks < 2^i."""
            if i == k:
                found[0] = True; return True
            def place_c(j):
                if j == i:
                    # fill Q rows for element i, check all new subsets
                    for mm in range(1 << i):
                        q = V
                        for l in range(i):
                            if (mm >> l) & 1: q += C[l][i]
                        q = (Q[mm] + q) % MOD
                        if q == 0: return False
                        Q[mm | (1 << i)] = q
                    return place(i + 1)
                for cv in range(MOD):
                    C[j][i] = cv
                    # early check: subsets with max old element j
                    ok = True
                    for mm in range(1 << j):
                        base = mm | (1 << j)
                        q = V
                        for l in range(j + 1):
                            if (base >> l) & 1: q += C[l][i]
                        if (Q[base] + q) % MOD == 0: ok = False; break
                    if ok and place_c(j + 1): return True
                return False
            return place_c(0)
        Q[0] = 0
        if V % MOD == 0: return best      # singleton dies instantly
        Q[1] = V % MOD
        t0 = time.time()
        if k == 1:
            res = True                    # single element avoids
        else:
            res = place(1)
        if res:
            best = (k, [row[:] for row in C])
            print('  k = %d: avoider EXISTS  [%.1fs]'
                  % (k, time.time() - t0), flush=True)
            k += 1
            if k > 12:
                print('  stopping at k = 12 (no bound found)')
                return best
        else:
            print('  k = %d: NO avoider (exhaustive)  [%.1fs]'
                  % (k, time.time() - t0), flush=True)
            return best

print('G-problem mod 4, V == 2, free couplings:')
b, wit = max_avoider(4, 2)
print('B4 = %d  (=> G4 = %d, and E(3) <= 2*F(2) <= 6*G4 = %d)'
      % (b, b + 1, 6 * (b + 1)))
if wit:
    print('witness couplings (upper triangle):')
    for i in range(b):
        print('  ', [wit[j][i] for j in range(i)])

print()
print('G-problem mod 8, V == 2, free couplings (the next link):')
b8, wit8 = max_avoider(8, 2)
print('B8(V=2) = %d' % b8)
print()
print('G-problem mod 8, V == 4, free couplings:')
b84, _ = max_avoider(8, 4)
print('B8(V=4) = %d' % b84)
print()
print('G-problem mod 8, V == 6, free couplings:')
b86, _ = max_avoider(8, 6)
print('B8(V=6) = %d' % b86)
