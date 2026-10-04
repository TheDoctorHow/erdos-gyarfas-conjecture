"""Delsarte LP feasibility for the cycle code of a cubic 3-connected girth->=5
mod-8-free graph.  Variables: the weight distribution A_0..A_m of the cycle code.
Forced constraints:
  A_0 = 1;  A_l = 0 for 1<=l<girth  (min weight = girth)
  A_l = 0 for l > n                 (max weight = n: codewords are unions of
                                     vertex-disjoint cycles, so |supp| <= n)
  A_8 = 0                           (WINDOW LEMMA: any codeword of weight < 2g is
                                     a single cycle, and 8 < 2*5, so a weight-8
                                     codeword would be an 8-cycle)
  sum A_l = 2^k, k = n/2+1
  dual distribution B_j = 2^{-k} sum_l A_l K_j(l) satisfies
     B_0 = 1, B_1 = B_2 = 0, B_3 = n (vertex stars; cyclically 4-edge-connected),
     B_j >= 0 for all j,  B_j = 0 for j odd where forced?  (not imposed)
If this LP is INFEASIBLE, no such graph exists and T(3) holds for that n.
"""
import sys
from math import comb
try:
    import numpy as np
    from scipy.optimize import linprog
except Exception as e:
    print("scipy/numpy unavailable:", e); sys.exit(0)

def krawtchouk(j, l, m):
    return sum((-1)**i * comb(l, i) * comb(m - l, j - i) for i in range(0, j + 1))

def test(n, girth=5, ban8=True, verbose=False):
    m = 3 * n // 2
    k = n // 2 + 1
    allowed = [0] + [l for l in range(girth, n + 1)]
    if ban8 and 8 in allowed and 8 < 2 * girth:
        allowed.remove(8)
    idx = {l: i for i, l in enumerate(allowed)}
    nv = len(allowed)
    Aeq, beq = [], []
    # A_0 = 1
    row = [0.0] * nv; row[idx[0]] = 1.0; Aeq.append(row); beq.append(1.0)
    # sum A_l = 2^k
    Aeq.append([1.0] * nv); beq.append(float(2 ** k))
    # B_1 = B_2 = 0 ; B_3 = n
    for j, val in [(1, 0.0), (2, 0.0), (3, float(n))]:
        row = [krawtchouk(j, l, m) / 2 ** k for l in allowed]
        Aeq.append(row); beq.append(val)
    # B_j >= 0  ->  -B_j <= 0
    Aub, bub = [], []
    for j in range(4, m + 1):
        row = [-krawtchouk(j, l, m) / 2 ** k for l in allowed]
        Aub.append(row); bub.append(0.0)
    res = linprog(c=[0.0] * nv, A_ub=np.array(Aub), b_ub=np.array(bub),
                  A_eq=np.array(Aeq), b_eq=np.array(beq),
                  bounds=[(0, None)] * nv, method="highs")
    return res

print("Delsarte LP feasibility for cycle codes of cubic girth->=5 mod-8-free graphs")
print("(INFEASIBLE would prove T(3) for that order)\n")
for n in [30, 32, 36, 40, 44, 50, 60]:
    r1 = test(n, ban8=True)
    r0 = test(n, ban8=False)
    print("n=%-3d m=%-3d k=%-3d   with A_8=0 : %-12s   control (A_8 free) : %s"
          % (n, 3*n//2, n//2+1, r1.message.split('.')[0], r0.message.split('.')[0]))

print()
print("=== higher-girth variants (window lemma bans every multiple of 8 below 2g) ===")
for n, g in [(30, 8), (60, 8), (70, 9), (90, 10), (126, 12)]:
    m = 3*n//2
    r = test(n, girth=g, ban8=True)
    print("n=%-4d girth=%-3d : %s" % (n, g, r.message.split('.')[0]))
