"""diam365.py -- do minimum odd-cycle transversals have large diameter?

364.3 suggested: a minimum transversal of size k has diameter >= 2k - 2.  At k = 3 that gives
4, which would exclude any 3-subset of a hexagon (diameter at most 3) and settle
(NO-SPLIT-HEXAGON).  TESTED, not assumed: enumerate every minimum transversal and measure its
diameter in B.
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

rng=random.Random(int(sys.argv[1]) if len(sys.argv)>1 else 8101)
st=collections.Counter(); done=0
for _ in range(int(sys.argv[2]) if len(sys.argv)>2 else 60):
    a,terms=make_real_bridge(30,rng)
    if a is None: continue
    V=sorted(a,key=str); tv=None; mins=[]
    for m in range(0,5):
        cur=[S for S in itertools.combinations(V,m) if bip_after(a,set(S))]
        if cur: tv=m; mins=cur; break
    if tv is None or tv>4 or tv<2: continue
    done+=1
    D={v:dist(a,v) for v in a}
    diams=[max(D[x][y] for x,y in itertools.combinations(S,2)) for S in mins if len(S)>=2]
    if not diams: continue
    lo=min(diams); bound=2*tv-2
    st['tau_v=%d  min diameter over %2d minimum transversals = %d   bound 2k-2 = %d  %s'
       %(tv,len(mins),lo,bound,'OK' if lo>=bound else '*** VIOLATED ***')]+=1
    st['  diameter distribution: %s'%str(sorted(collections.Counter(diams).items()))]+=1
    if done>=6: break
print('  bridges examined: %d'%done)
for k in sorted(st): print('   %s  %s'%(k.ljust(66),st[k]))
