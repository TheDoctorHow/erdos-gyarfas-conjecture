"""diam366.py -- (DIAM-4) on a wider and larger sample.

365 measured the minimum diameter of a minimum odd-cycle transversal as exactly 4, over six
bridges on 33 vertices.  Widen: more bridges, and n = 43 and 53, to see whether 4 is really
the floor or an artefact of size.
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

rng=random.Random(int(sys.argv[1]) if len(sys.argv)>1 else 8501)
st=collections.Counter(); allmin=[]
for NV in (30,40):
    done=0
    for _ in range(50):
        a,terms=make_real_bridge(NV,rng)
        if a is None: continue
        V=sorted(a,key=str); tv=None; mins=[]
        for m in range(0,5):
            cur=[S for S in itertools.combinations(V,m) if bip_after(a,set(S))]
            if cur: tv=m; mins=cur; break
        if tv is None or tv>4 or tv<2: continue
        done+=1
        D={v:dist(a,v) for v in a}
        diams=[max(D[x][y] for x,y in itertools.combinations(S,2)) for S in mins]
        lo=min(diams); allmin.append(lo)
        st['n=%-3d tau_v=%d  %3d minimum transversals, MIN DIAMETER %d'%(len(a),tv,len(mins),lo)]+=1
        if done>=3: break
for k in sorted(st): print('   %s  %s'%(k.ljust(58),st[k]))
print('  overall minimum diameter observed: %d  over %d bridges'%(min(allmin),len(allmin)))
print('  distribution of per-bridge minima: %s'%sorted(collections.Counter(allmin).items()))
