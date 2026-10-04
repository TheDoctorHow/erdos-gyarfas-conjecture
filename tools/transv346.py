"""transv346.py -- is some even cycle through a terminal NOT an odd-cycle transversal?

Choosing D through a terminal b_j makes the tail Q_2 empty, so the arc-shift hypothesis
reduces to: B - D connected (containing b_i and y') and NON-BIPARTITE.
B - D non-bipartite <=> some ODD CYCLE IS DISJOINT FROM D <=> D IS NOT AN ODD-CYCLE
TRANSVERSAL.  So the question is whether every bridge has an even cycle through a terminal
that fails to be a transversal.  Measured on genuine bridges.
"""
import sys, random, itertools, collections
from collections import deque
exec(open('tools/canon345.py').read().split('rng=random.Random')[0].split('"""',2)[2])

rng=random.Random(int(sys.argv[1]) if len(sys.argv)>1 else 2711)
st=collections.Counter()
for _ in range(int(sys.argv[2]) if len(sys.argv)>2 else 6):
    a,terms=make_real_bridge(30,rng)
    if a is None: continue
    st['genuine bridges (n=%d)'%len(a)]+=1
    seen=set(); thru=[]
    for _ in range(6000):
        c=rand_cycle(a,rng)
        if not c or len(c)%2: continue
        k=frozenset(c)
        if k in seen: continue
        seen.add(k)
        if any(t in c for t in terms): thru.append(c)
    good=[]; 
    for D in thru:
        R=comp_of(a,set(D),None) if False else None
        rest=set(a)-set(D)
        if not rest: continue
        s=next(iter(rest)); C=comp_of(a,set(D),s)
        if C is None or len(C)!=len(rest): continue        # B-D connected
        if nonbip(a,C): good.append(D)
    st['  even cycles through a terminal: %d'%len(thru)]+=1
    st['  of which B-D connected and non-bipartite: %d'%len(good)]+=1
    if good:
        st['  smallest such |D| = %d'%min(len(D) for D in good)]+=1
    st['SOME non-transversal even terminal cycle' if good
       else '*** every such cycle is a transversal ***']+=1
for k in sorted(st): print('  %-52s %s'%(k,st[k]))
