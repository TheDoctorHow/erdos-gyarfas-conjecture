"""diam439.py -- the diameter board at t = 3: full-stall death and the
blocky-escape CSP.

Setting (438.2 + carrier forcing): at m = 3, t = 3 the supply is 12
diameter ears (u, u+12), u = 0..11, lengths T_u. Crossing-pair
composite length: |C'| = T_i + T_j + 24 - 2(j-i); progress <=> 16
divides |C'|. With even T (the cascade branch), T_u = 2 sigma_u mod
16, sigma in Z/8, and the survival system is
    sigma_i + sigma_j =/= (j - i) + 4  (mod 8)   for all i < j,
with the per-ear shallowness domain sigma_u =/= 2 (mod 4).

 (1) FULL STALL: all sigma equal -> pair sums fixed = 2 sigma; the
     gaps (j-i)+4 cover ALL residues mod 8 -> collision forced:
     DEAD. Verified.
 (2) BLOCKY ESCAPE: exhaustive CSP over sigma_0..sigma_11 in
     {0,1,3,4,5,7}: SAT (witness printed) or UNSAT (by exhaustion
     with pruning).
 (3) ODD BRANCH (cascade evaded): tau = 1 + 2u dodges all pairs -
     verified - the odd-diameter world survives THIS board and is
     handed to the odd-side weapons (437.3 extras), recorded.
"""
import sys

# ---------------- (1) full stall ----------------
gaps = sorted({(j - i + 4) % 8 for i in range(12) for j in range(i + 1, 12)})
print('(1) full stall: residues (j-i)+4 mod 8 realised by the gaps: %s'
      % gaps)
assert gaps == list(range(8))
print('    -> a fixed pair-sum 2 sigma collides for every sigma:'
      ' FULL-STALL EVEN DIAMETER WORLD DEAD.')

# ---------------- (2) blocky escape CSP ----------------
DOM = [0, 1, 3, 4, 5, 7]        # sigma =/= 2 (mod 4)
N = 12
sol = None
nodes = 0
def dfs(assign):
    global sol, nodes
    nodes += 1
    i = len(assign)
    if i == N:
        sol = list(assign); return True
    for v in DOM:
        ok = True
        for j, w in enumerate(assign):
            if (w + v) % 8 == ((i - j) + 4) % 8:
                ok = False; break
        if ok:
            assign.append(v)
            if dfs(assign): return True
            assign.pop()
    return False

sat = dfs([])
print('(2) blocky-escape CSP (12 vars, domain %s, 66 constraints):'
      ' %s  [%d nodes]' % (DOM, 'SAT' if sat else 'UNSAT', nodes))
if sat:
    print('    witness sigma = %s' % sol)
    # independent verification
    ok = all((sol[i] + sol[j]) % 8 != ((j - i) + 4) % 8
             for i in range(N) for j in range(i + 1, N))
    print('    verified: %s' % ok)
    assert ok

# ---------------- (3) odd branch ----------------
tau = [(1 + 2 * u) % 16 for u in range(12)]
ok = all((tau[i] + tau[j] + 24 - 2 * (j - i)) % 16 != 0
         for i in range(12) for j in range(i + 1, 12))
print('(3) odd branch tau = 1 + 2u: all 66 composites dodge 16: %s'
      % ok)
assert ok

# ---------------- (2b) if SAT: how blocky is every witness? ----------------
if sat:
    # count maximal constant runs in the witness and search for the
    # minimum number of distinct values over all solutions (bounded)
    best = [99]
    def dfs2(assign, distinct):
        if len(distinct) >= best[0]: return
        i = len(assign)
        if i == N:
            best[0] = len(distinct); return
        for v in DOM:
            if all((w + v) % 8 != ((i - j) + 4) % 8
                   for j, w in enumerate(assign)):
                assign.append(v)
                dfs2(assign, distinct | {v})
                assign.pop()
    dfs2([], set())
    print('(2b) minimum number of DISTINCT sigma-values over all'
          ' solutions: %d' % best[0])
    print('     -> every survivor carries at least that much'
          ' sigma-diversity: deviation fuel for the zigzag APs.')
