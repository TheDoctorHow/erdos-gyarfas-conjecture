"""delta321.py -- the predicted pattern: firing is governed by delta = tau_3 - tau.

tau   = least |N| over ALL colourings (frustration index)
tau_3 = least |N| over colourings with the three terminals MONOCHROMATIC
delta = tau_3 - tau  = 0 exactly when some MINIMUM bipartition already has the terminals
        agreeing.
PREDICTION.  delta = 0  =>  fans avoiding N exist at vertices of both colours, so
sigma = c(v) takes both values and the Balanced Fan Theorem fires.
        delta > 0  =>  legs must pick up N-parity, balanced fans are rare, sigma is
constant, and the mechanism misses.
"""
import sys, random, itertools, collections
from collections import deque
exec(open('tools/yparity315.py').read().split('rng=random.Random')[0].split('"""',2)[2])

def tau_free(E, cap=4):
    El=sorted(E,key=lambda s:sorted(map(str,s)))
    for k in range(cap+1):
        for S in itertools.combinations(El,k):
            if bipn(adj_of(set(El)-set(S))) is not None: return k
    return cap+1

rng=random.Random(int(sys.argv[1]) if len(sys.argv)>1 else 77)
st=collections.Counter(); rows=[]
for _ in range(int(sys.argv[2]) if len(sys.argv)>2 else 40):
    a,terms=make_bridge(rng.choice([10,12,14]),rng)
    if a is None: continue
    E={frozenset((u,w)) for u in a for w in a[u]}
    t=tau_free(E)
    if t<2: continue
    best=None
    for trial in range(250):
        col={v:rng.randint(0,1) for v in a}
        for x in terms: col[x]=0
        imp=True
        while imp:
            imp=False
            for v in a:
                if v in terms: continue
                nc=sum(1 for w in a[v] if col[w]==col[v])
                if nc>len(a[v])-nc: col[v]^=1; imp=True
        N={e for e in E if col[tuple(e)[0]]==col[tuple(e)[1]]}
        if best is None or len(N)<len(best[1]): best=(dict(col),set(N))
    col,N=best; t3=len(N); d=t3-t
    P={v:y_parities(a,v,terms,N,rng,trials=160) for v in a if v not in terms}
    SIG={(col[v]+p)%2 for v in P for p in (0,1) if (p,p,p) in P[v]}
    fires=(len(SIG)==2)
    nbal=sum(1 for v in P if (0,0,0) in P[v] or (1,1,1) in P[v])
    st['bridges with tau>=2']+=1
    st['delta=%d : %s'%(d,'FIRES' if fires else 'MISS ')]+=1
    rows.append((d,fires,nbal,t,t3))
for k in sorted(st): print('  %-30s %s'%(k,st[k]))
print('  delta   fires/total   mean # balanced fans')
by=collections.defaultdict(list)
for d,f,nb,_,_ in rows: by[d].append((f,nb))
for d in sorted(by):
    v=by[d]
    print('    %2d      %2d/%-2d          %.1f'%(d,sum(1 for f,_ in v if f),len(v),
                                                 sum(nb for _,nb in v)/len(v)))
