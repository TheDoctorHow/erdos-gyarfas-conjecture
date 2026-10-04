"""bound351.py -- the Edge Odd Cycle bound against the truth.

LEMMA (proved): in a 3-connected non-bipartite graph, an edge e lies on an odd cycle of
length at most 1 + |C| + 2 d(e,C) for any odd cycle C.  Via the translation of 347.2, the
shortest EVEN cycle of B through terminal b is at most 2 + |C| + 2 d(b,C); with girth 5 and
C a 5-cycle that is 7 + 2 d(b,C).
Measured: d(b, nearest 5-cycle), the bound, and the true shortest even terminal cycle.
"""
import sys, random, itertools, collections
from collections import deque
exec(open('tools/canon345.py').read().split('rng=random.Random')[0].split('"""',2)[2])

def five_cycles(a):
    out=set()
    for s in a:
        for p1 in a[s]:
            for p2 in a[p1]:
                if p2==s: continue
                for p3 in a[p2]:
                    if p3 in (s,p1): continue
                    for p4 in a[p3]:
                        if p4 in (s,p1,p2): continue
                        if s in a[p4]: out.add(frozenset((s,p1,p2,p3,p4)))
    return [set(c) for c in out]
def shortest_even_terminal_cycle(a, terms):
    best=None
    for b in terms:
        nb=sorted(a[b],key=str)
        if len(nb)!=2: continue
        p,q=nb
        d={(p,0):0}; Q=deque([(p,0)])
        while Q:
            u,par=Q.popleft()
            for w in a[u]:
                if w==b: continue
                s=(w,1-par)
                if s not in d: d[s]=d[(u,par)]+1; Q.append(s)
        if (q,0) in d:
            L=d[(q,0)]+2
            if best is None or L<best[0]: best=(L,b)
    return best

rng=random.Random(int(sys.argv[1]) if len(sys.argv)>1 else 4001)
st=collections.Counter()
for NV in (30,40,50):
    for _ in range(int(sys.argv[2]) if len(sys.argv)>2 else 3):
        a,terms=make_real_bridge(NV,rng)
        if a is None: continue
        FC=five_cycles(a)
        if not FC: st['n=%-3d NO 5-cycle in B'%len(a)]+=1; continue
        r=shortest_even_terminal_cycle(a,terms)
        if r is None: continue
        L,b=r
        D=dist(a,b)
        dmin=min(min(D[v] for v in c if v in D) for c in FC)
        st['n=%-3d  d(b,5-cycle)=%d  bound %d  actual %d  %s'
           %(len(a),dmin,7+2*dmin,L,'OK' if L<=7+2*dmin else 'BOUND VIOLATED')]+=1
        st['  number of 5-cycles in B: %d'%len(FC)]+=1
for k in sorted(st): print('  %s  %s'%(k.ljust(56),st[k]))
