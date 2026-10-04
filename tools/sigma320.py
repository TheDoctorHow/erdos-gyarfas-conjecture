"""sigma320.py -- THE UNIFORM REASON: a single invariant sigma.

With uniform leg parities p at w and p' at w', 295.1 collapses the shift parity to
     e_m == [c(w) =/= c(w')] + p + p'   (mod 2),   independent of m.
Step (a) of 199.2 forces every e_m EVEN.  So (SC) forces
     sigma(w,Y) := c(w) + p   (mod 2)
to be CONSTANT over every uniform-parity fan in B - including two fans at the SAME vertex,
since nothing in step (a) requires w and w' distinct.
Fires iff sigma takes both values.  Tested against the previous two forms.
"""
import sys, random, itertools, collections
from collections import deque
exec(open('tools/yparity315.py').read().split('rng=random.Random')[0].split('"""',2)[2])

rng=random.Random(int(sys.argv[1]) if len(sys.argv)>1 else 57)
st=collections.Counter(); miss=[]
for _ in range(int(sys.argv[2]) if len(sys.argv)>2 else 40):
    a,terms=make_bridge(rng.choice([10,12,14]),rng)
    if a is None: continue
    E={frozenset((u,w)) for u in a for w in a[u]}
    t=tau_edge(E)
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
    col,N=best; t3=len(N)
    P={v:y_parities(a,v,terms,N,rng,trials=160) for v in a if v not in terms}
    adjf=any((e not in N) and all(x in P for x in tuple(e))
             and any(val in P[tuple(e)[0]] and val in P[tuple(e)[1]]
                     for val in ((0,0,0),(1,1,1))) for e in E)
    oppf=any(len({col[v] for v in P if (p,p,p) in P[v]})==2 for p in (0,1))
    SIG={(col[v]+p)%2 for v in P for p in (0,1) if (p,p,p) in P[v]}
    sigf=(len(SIG)==2)
    st['bridges with tau>=2']+=1
    st['TOTAL adjacency']       +=1 if adjf else 0
    st['TOTAL opposite-colour'] +=1 if oppf else 0
    st['TOTAL sigma']           +=1 if sigf else 0
    st['tau_3=%d : sigma %s'%(t3,'FIRES' if sigf else 'MISS')]+=1
    if not sigf: miss.append((len(a),t,t3,sorted(SIG),
                              sum(1 for v in P if (0,0,0) in P[v] or (1,1,1) in P[v])))
for k in sorted(st): print('  %-38s %s'%(k,st[k]))
if miss:
    print('  sigma misses (n, tau, tau_3, sigma-values, #uniform-parity vertices):')
    for m in miss[:10]: print('    ',m)
