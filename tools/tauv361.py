"""tauv361.py -- is tau_v ever 2 or 3 in a genuine bridge with tau >= 2?

(TAU-4) needs tau_v >= 4.  It does NOT follow from tau_e >= 2: deleting all three edges at a
vertex gives only tau_e <= 3 tau_v, so tau_e >= 2 permits tau_v = 2.  So it must come from
girth and cubicity.  Direct test over many bridges: how small does tau_v get?
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

rng=random.Random(int(sys.argv[1]) if len(sys.argv)>1 else 6701)
st=collections.Counter()
for _ in range(int(sys.argv[2]) if len(sys.argv)>2 else 14):
    a,terms=make_real_bridge(rng.choice([30,40]),rng)
    if a is None: continue
    V=sorted(a,key=str); k=None
    for m in range(0,5):
        found=False
        for S in itertools.combinations(V,m):
            if bip_after(a,set(S)): found=True; break
        if found: k=m; break
    st['tau_v = %s'%(k if k is not None else '>= 5')]+=1
print('  distribution of tau_v over genuine bridges:')
for k in sorted(st): print('   %-14s %s'%(k,st[k]))
