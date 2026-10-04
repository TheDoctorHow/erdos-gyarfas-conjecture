"""opp319.py -- ADJACENCY WAS NEVER NEEDED.

314.1 subtracts the three (SC) relations at w from those at w'.  The constants c_ij belong
to the terminal PAIR, not to w, so the subtraction is legitimate for ANY two vertices, and
gives every shift e_m even.  The other half needs only c(w') =/= c(w) to make every e_m
odd.  Adjacency merely guaranteed that.  So:

  THEOREM.  If two vertices of OPPOSITE COLOUR each carry a Y whose three legs have the
  same N-parity, and the two parities agree, then (SC) fails.

Tested against the adjacency version on the same bridges.
"""
import sys, random, itertools, collections
from collections import deque
exec(open('tools/yparity315.py').read().split('rng=random.Random')[0].split('"""',2)[2])

rng=random.Random(int(sys.argv[1]) if len(sys.argv)>1 else 57)
st=collections.Counter()
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
    # adjacency version
    adjfire=any((e not in N) and (tuple(e)[0] in P) and (tuple(e)[1] in P)
                and any(val in P[tuple(e)[0]] and val in P[tuple(e)[1]]
                        for val in ((0,0,0),(1,1,1)))
                for e in E)
    # opposite-colour version
    oppfire=False
    for p in (0,1):
        U=[v for v in P if (p,p,p) in P[v]]
        if len({col[v] for v in U})==2: oppfire=True; break
    st['bridges with tau>=2']+=1
    st['tau_3=%d  adjacency:%s  opposite-colour:%s'
       %(t3,'fires' if adjfire else 'MISS ','FIRES' if oppfire else 'MISS ')]+=1
    st['TOTAL adjacency fires']+=1 if adjfire else 0
    st['TOTAL opposite-colour fires']+=1 if oppfire else 0
for k in sorted(st): print('  %-46s %s'%(k,st[k]))
