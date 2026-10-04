"""cyc337.py -- (CYC): a strictly weaker target than (BAL).

If Z is a cycle through terminals b_i and b_j, its two arcs are CONNECTING PATHS.  If both
arcs are EVEN then (SC) pins both to c_ij mod 4, so their lengths agree mod 4 and
     |Z| = arc1 + arc2 == 2*arc1 == 0 (mod 4),  arc1 being even.
Hence:
     (CYC)  a cycle through two terminals, BOTH ARCS EVEN, of length == 2 (mod 4)
            REFUTES (SC) outright.
Needs ONE cycle, TWO terminals, TWO prescribed parities - against (BAL)'s two vertices,
three internally disjoint paths and three prescribed parities.  Tested exhaustively.
"""
import sys, random, itertools, collections
from collections import deque
exec(open('tools/ecore2.py').read().split('def main()')[0].split(chr(34)*3,2)[2])

def all_simple(a, s, t, banned):
    out=[]; p=[s]; on=set(banned)|{s}
    def rec(u):
        for w in a[u]:
            if w==t: out.append(list(p)+[t])
            elif w not in on:
                on.add(w); p.append(w); rec(w); p.pop(); on.discard(w)
    rec(s); return out

rng=random.Random(int(sys.argv[1]) if len(sys.argv)>1 else 331)
st=collections.Counter()
for _ in range(int(sys.argv[2]) if len(sys.argv)>2 else 12):
    a,terms=make_bridge(10,rng)
    if a is None: continue
    st['bridges']+=1
    found=None
    for i,j in [(0,1),(0,2),(1,2)]:
        P=[p for p in all_simple(a,terms[i],terms[j],set()) if (len(p)-1)%2==0]
        for p1,p2 in itertools.combinations(P,2):
            if set(p1[1:-1]) & set(p2[1:-1]): continue      # internally disjoint -> a cycle
            L=(len(p1)-1)+(len(p2)-1)
            if L%4==2: found=(i,j,len(p1)-1,len(p2)-1,L); break
        if found: break
    st['(CYC) witness found' if found else '*** no (CYC) witness ***']+=1
    if found: st['  pair (%d,%d), arcs %d and %d, |Z|=%d'%found]+=1
for k in sorted(st): print('  %-48s %s'%(k,st[k]))
