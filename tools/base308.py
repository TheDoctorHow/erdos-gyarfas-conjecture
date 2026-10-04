"""base308.py -- is the base case of the 307 descent EMPTY?

An ADMISSIBLE block is: 2-connected, bipartite, subcubic, girth >= 6, with at most FOUR
vertices of degree 2 (all others degree 3), and a rigid pair of two degree-2 vertices.

Necessary arithmetic.  Bipartition (L,R), |L|=p, |R|=q, degrees 2 or 3, both sides have
degree sum m.  Girth >= 6 means no two vertices of R share two neighbours in L, so the
number of cherries centred in R is at most the number of L-pairs:

        3 q_3 + q_2  <=  C(p,2)      and symmetrically   3 p_3 + p_2 <= C(q,2).

This enumerates all (n, p, q, p_2, q_2) with n <= 20 and reports for which n an admissible
degree sequence survives.
"""
import itertools
print('  n   feasible degree splits (p,q,p2,q2,m)')
first=None
for n in range(3, 21):
    sols=[]
    for p in range(1, n):
        q = n - p
        for p2 in range(0, min(p,4)+1):
            for q2 in range(0, min(q,4)+1):
                if p2 + q2 > 4: continue
                if p2 + q2 < 2: continue          # need two degree-2 vertices (the pair)
                p3 = p - p2; q3 = q - q2
                if p3 < 0 or q3 < 0: continue
                m1 = 3*p3 + 2*p2; m2 = 3*q3 + 2*q2
                if m1 != m2: continue
                m = m1
                if 3*q3 + q2 > p*(p-1)//2: continue      # girth >= 6, cherries in R
                if 3*p3 + p2 > q*(q-1)//2: continue      # girth >= 6, cherries in L
                if m < n: continue                        # 2-connected needs m >= n
                sols.append((p,q,p2,q2,m))
    if sols:
        if first is None: first = n
        print('  %2d  %s' % (n, sols[:6]))
    else:
        print('  %2d  -- NONE --' % n)
print()
print('smallest n admitting an admissible degree sequence:', first)
