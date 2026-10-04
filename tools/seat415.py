"""seat415.py -- the Seating Theorem: closed-form emptiness through A = 33.

Coverable intervals (chords span >= 4, pockets >= 3 atts with colour-diff
== 0 mod 3, colours alternating): exactly {0, 4, 6} union [7, oo).
Same-side window gaps even, intervals odd, odd lengths 1, 3, 5 dead: gaps
>= 8.  Window Parity Law fixes window parities per beta; window census
|WIN| == 2 (mod 3) at beta = 1, == 1 (mod 3) at beta = 0; s_1 even.
DP feasibility per A: infeasible for all A <= 33; first seat at A = 35
(beta = 0, s_1 = 8, one low + three high windows).
"""
from functools import lru_cache
COV = set([0, 4, 6]) | set(range(7, 400))

def side_ok(n, par, off, k):
    cand = [p for p in range(off, off + n) if p % 2 == par]
    @lru_cache(None)
    def go(i, placed, prev):
        if placed == k:
            return (off + n - 1 - prev) in COV
        if i >= len(cand): return False
        if go(i + 1, placed, prev): return True
        w = cand[i]
        return (w - prev - 1) in COV and go(i + 1, placed + 1, w)
    return go(0, 0, off - 1)

if __name__ == '__main__':
    feas = []
    for A in range(5, 120, 2):
        found = None
        for beta in (0, 1):
            need = 2 if beta == 1 else 1
            for s1 in range(2, A - 3, 2):
                for kL in range(1, 10):
                    if not side_ok(s1, (beta + 1) % 2, 0, kL): continue
                    for kR in range(1, 10):
                        if kL + kR < 3 or (kL + kR - need) % 3: continue
                        if side_ok(A - 3 - s1, beta % 2, s1 + 2, kR):
                            found = (beta, s1, kL, kR); break
                    if found: break
                if found: break
            if found: break
        if found: feas.append((A, found))
    dead = [a for a in range(5, 120, 2) if a not in [x[0] for x in feas]]
    print('infeasible A:', dead)
    print('first feasible:', feas[0] if feas else None)
