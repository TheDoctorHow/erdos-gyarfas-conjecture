"""cyctrans360.py -- is ANY cycle of a bridge an odd-cycle transversal?

359 showed tau_v is often 4 or 6, so SIZE does not forbid a hexagon from being a transversal.
Yet 351 found 12 of 12 terminal hexagons are not.  So the obstruction is structural.  Test it
in the sharpest form: over all short cycles of genuine bridges - terminal or not - is any one
of them an odd-cycle transversal?
256.1's transversal is a PATH INTERIOR, which induces a path; a cycle induces a cycle.  If no
cycle is ever a transversal, that shape difference is the mechanism.
"""
import sys, random, itertools, collections
from collections import deque
exec(open('tools/canon345.py').read().split('rng=random.Random')[0].split('"""',2)[2])

def bip_after(a,S):
    col={}
    for s in a:
        if s in S or s in col: continue
        col[s]=0; q=deque([s])
        while q:
            u=q.popleft()
            for w in a[u]:
                if w in S: continue
                if w not in col: col[w]=1-col[u]; q.append(w)
                elif col[w]==col[u]: return False
    return True

rng=random.Random(int(sys.argv[1]) if len(sys.argv)>1 else 6301)
st=collections.Counter()
for NV in (30,40):
    for _ in range(int(sys.argv[2]) if len(sys.argv)>2 else 3):
        a,terms=make_real_bridge(NV,rng)
        if a is None: continue
        seen=set(); cyc=[]
        for _ in range(9000):
            c=rand_cycle(a,rng)
            if not c or len(c)>12: continue
            k=frozenset(c)
            if k in seen: continue
            seen.add(k); cyc.append(c)
        hits=[c for c in cyc if bip_after(a,set(c))]
        st['n=%-3d  short cycles %3d, of which TRANSVERSALS %d'%(len(a),len(cyc),len(hits))]+=1
        if hits:
            st['  *** a cycle IS a transversal: |C|=%d ***'%min(len(c) for c in hits)]+=1
        # is any MINIMUM transversal cycle-shaped?  report smallest transversal size found
        st['  smallest transversal among these cycles: %s'
           %(min(len(c) for c in hits) if hits else 'none')]+=1
for k in sorted(st): print('  %s  %s'%(k.ljust(58),st[k]))
