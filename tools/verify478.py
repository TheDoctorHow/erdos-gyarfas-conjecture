"""verify478.py -- INDEPENDENT RE-VERIFICATION OF THE q_3 = 3 WITNESS.

base477 reports a 24-vertex cubic graph with a 24-cycle, no 8-cycle
and no 16-cycle, refuting (DESC_3) and (STEP) at t = 3.  Per the
standing rule, a witness is never accepted from printed output: this
file re-reads the adjacency, re-checks that it IS a graph of the
claimed kind, and computes the COMPLETE cycle-length spectrum by
exhaustive enumeration - a different algorithm from the per-length
existence test that produced the claim.
"""
import sys
sys.setrecursionlimit(100000)
P = print

ADJ = [[2, 3, 22], [5, 17, 21], [0, 12, 18], [0, 20, 22], [7, 10, 21],
       [1, 11, 14], [9, 12, 18], [4, 10, 21], [10, 15, 23], [6, 14, 17],
       [4, 7, 8], [5, 16, 19], [2, 6, 18], [14, 17, 20], [5, 9, 13],
       [8, 16, 23], [11, 15, 19], [1, 9, 13], [2, 6, 12], [11, 16, 23],
       [3, 13, 22], [1, 4, 7], [0, 3, 20], [8, 15, 19]]
n = len(ADJ)

P('=' * 72)
P('STRUCTURAL CHECKS', flush=True)
P('=' * 72)
ok = True
P('  vertices                       : %d' % n)
deg = all(len(a) == 3 for a in ADJ)
P('  every degree exactly 3         : %s' % deg); ok &= deg
sym = all(v in ADJ[w] for v in range(n) for w in ADJ[v])
P('  adjacency symmetric            : %s' % sym); ok &= sym
simple = all(len(set(ADJ[v])) == 3 and v not in ADJ[v] for v in range(n))
P('  simple (no loops/multi-edges)  : %s' % simple); ok &= simple
seen = {0}; st = [0]
while st:
    v = st.pop()
    for w in ADJ[v]:
        if w not in seen: seen.add(w); st.append(w)
P('  connected                      : %s' % (len(seen) == n))
ok &= (len(seen) == n)

P()
P('=' * 72)
P('COMPLETE CYCLE SPECTRUM BY EXHAUSTIVE ENUMERATION', flush=True)
P('=' * 72)
counts = {}
for s in range(n):
    vis = [False] * n; vis[s] = True
    def dfs(v, ln):
        for w in ADJ[v]:
            if w < s: continue
            if w == s:
                if ln >= 3: counts[ln] = counts.get(ln, 0) + 1
                continue
            if vis[w]: continue
            vis[w] = True; dfs(w, ln + 1); vis[w] = False
    dfs(s, 1); vis[s] = False
counts = {L: c // 2 for L, c in counts.items()}
spec = sorted(counts)
P('  total simple cycles            : %d' % sum(counts.values()))
P('  spectrum                       : %s' % spec, flush=True)
P()
P('  %6s %10s' % ('length', 'count'), flush=True)
for L in spec:
    P('  %6d %10d' % (L, counts[L]), flush=True)

P()
P('=' * 72)
P('THE CLAIM', flush=True)
P('=' * 72)
mult8 = [L for L in spec if L % 8 == 0]
pows = [L for L in spec if L > 0 and (L & (L - 1)) == 0]
P('  cycle lengths divisible by 8   : %s' % mult8, flush=True)
P('  powers of two in the spectrum  : %s' % pows, flush=True)
P('  8-cycle present                : %s' % (8 in spec))
P('  16-cycle present               : %s' % (16 in spec))
P('  24-cycle present               : %s' % (24 in spec))
q3 = (min(mult8) // 8) if mult8 else None
P('  q_3 = min{ l/8 : 8 | l }       : %s' % q3)
P()
c1 = (mult8 == [24])
c2 = (q3 == 3)
c3 = bool(pows)
P('  (DESC_3) predicted min multiple of 8 to be a power of two.')
P('  Observed min multiple of 8 = %s.' % (min(mult8) if mult8 else '-'))
P()
if ok and c1 and c2:
    P('  CONFIRMED: (DESC_3) IS FALSE, and (STEP) is false at t = 3,')
    P('  m = 3 - the 24-cycle has no strictly shorter companion of')
    P('  length divisible by 8.')
else:
    P('  NOT CONFIRMED - the reported witness does not survive.')
P()
if c3:
    P('  Erdos-Gyarfas is NOT refuted: the spectrum contains the')
    P('  power(s) of two %s.  The witness escapes the descent by' % pows)
    P('  paying in a DIFFERENT stratum from the one being descended.')
else:
    P('  *** NO POWER OF TWO IN THE SPECTRUM - CHECK IMMEDIATELY ***')
