"""table316.py -- the 2x2 that matters: over bridges with tau >= 2, cross
(does the parity Y-shift configuration exist?) against (does (SC) actually hold?).
The mechanism only has to cover bridges for which (SC) HOLDS.
"""
import sys, random, itertools, collections
from collections import deque
exec(open('tools/yparity315.py').read().split('rng=random.Random')[0].split('"""',2)[2])

def all_paths(a,s,t,cap=400000):
    out=[];p=[s];on={s}
    def rec(u):
        if len(out)>cap: return
        for w in a[u]:
            if w==t: out.append(len(p))
            elif w not in on:
                on.add(w);p.append(w);rec(w);p.pop();on.discard(w)
    rec(s);return out
def sc_holds(a,terms):
    for i,j in itertools.combinations(range(3),2):
        L=[x for x in all_paths(a,terms[i],terms[j]) if x%2==0]
        if L and len({x%4 for x in L})>1: return False
    return True

rng=random.Random(int(sys.argv[1]) if len(sys.argv)>1 else 23)
st=collections.Counter()
for _ in range(int(sys.argv[2]) if len(sys.argv)>2 else 60):
    a,terms=make_bridge(rng.choice([10,12,14]),rng)
    if a is None: continue
    E={frozenset((u,w)) for u in a for w in a[u]}
    t=tau_edge(E)
    if t<2: continue
    st['bridges with tau>=2']+=1
    best=None
    for trial in range(200):
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
    col,N=best
    P={v:y_parities(a,v,terms,N,rng) for v in a if v not in terms}
    fired=False
    for e in E:
        if e in N: continue
        u,v=tuple(e)
        if u in terms or v in terms: continue
        if any(val in P[u] and val in P[v] for val in ((0,0,0),(1,1,1))): fired=True;break
    sc=sc_holds(a,terms)
    st['Y-shift %s  /  (SC) %s'%('fires    ' if fired else 'NOT found',
                                 'HOLDS' if sc else 'fails')]+=1
for k in sorted(st): print('  %-40s %s'%(k,st[k]))
