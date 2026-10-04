"""disj364.py -- does a bridge ever have TWO DISJOINT MINIMUM odd-cycle transversals?

360.1 forces a transversal hexagon at tau_v = 3 to split into two DISJOINT minimum
transversals.  So:
     no two disjoint minimum transversals  ==>  (NO-SPLIT-HEXAGON)  ==>  gap (B)  ==> ...
253.3 records that the path-interior transversals of a failing pair are PAIRWISE
INTERSECTING, which is suggestive but concerns only those transversals.  Test the general
statement directly: enumerate ALL minimum transversals and look for a disjoint pair.
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

rng=random.Random(int(sys.argv[1]) if len(sys.argv)>1 else 7901)
st=collections.Counter()
done=0
for _ in range(int(sys.argv[2]) if len(sys.argv)>2 else 60):
    a,terms=make_real_bridge(30,rng)
    if a is None: continue
    V=sorted(a,key=str); tv=None; mins=[]
    for m in range(0,5):
        cur=[S for S in itertools.combinations(V,m) if bip_after(a,set(S))]
        if cur: tv=m; mins=cur; break
    if tv is None or tv>4: continue
    done+=1
    disj=None
    for S,T in itertools.combinations(mins,2):
        if not (set(S)&set(T)): disj=(S,T); break
    st['tau_v=%d  minimum transversals %3d  DISJOINT PAIR: %s'
       %(tv,len(mins),'YES' if disj else 'no')]+=1
    if done>=8: break
print('  bridges examined: %d'%done)
for k in sorted(st): print('   %s  %s'%(k.ljust(56),st[k]))
