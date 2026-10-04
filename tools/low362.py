"""low362.py -- does gap (B) survive in bridges with tau_v <= 3?

361 left tau_v in {2,3} uncovered by the Transversal Cycle Bound.  360.1 says a transversal
hexagon there would need psi-classes exactly 3-3, both of them MINIMUM transversals - very
rigid.  Question: in bridges with tau_v <= 3, does some short even terminal cycle still have
a non-bipartite complement?  That is gap (B) holding without the bound certifying it.
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
def tau_v_upto(a,cap=3):
    V=sorted(a,key=str)
    for m in range(cap+1):
        for S in itertools.combinations(V,m):
            if bip_after(a,set(S)): return m
    return cap+1

rng=random.Random(int(sys.argv[1]) if len(sys.argv)>1 else 7101)
st=collections.Counter()
found=0
for _ in range(int(sys.argv[2]) if len(sys.argv)>2 else 60):
    a,terms=make_real_bridge(rng.choice([30,40]),rng)
    if a is None: continue
    tv=tau_v_upto(a,3)
    if tv>3: st['tau_v >= 4 (covered by 360.1)']+=1; continue
    found+=1
    st['*** tau_v = %d  (UNCOVERED case) ***'%tv]+=1
    seen=set(); best=None; anygood=False; ntrans=0; nt=0
    for _ in range(6000):
        c=rand_cycle(a,rng)
        if not c or len(c)%2 or len(c)>12: continue
        if not any(t in c for t in terms): continue
        k=frozenset(c)
        if k in seen: continue
        seen.add(k); nt+=1
        if bip_after(a,set(c)): ntrans+=1
        else:
            anygood=True
            if best is None or len(c)<best: best=len(c)
    st['   terminal even cycles <=12: %d, of which TRANSVERSALS: %d'%(nt,ntrans)]+=1
    st['   gap (B) HOLDS at |D|=%s'%best if anygood else '   *** gap (B) FAILS for all ***']+=1
    if found>=4: break
for k in sorted(st): print('  %s  %s'%(k.ljust(58),st[k]))
