"""multicol317.py -- the freedom not yet used: the CHOICE OF COLOURING.

314 and 315 tested a single colouring per bridge - the one minimising |N| subject to the
three terminals being monochromatic.  But (Y-PARITY) may hold for a DIFFERENT colouring:
every colouring with the terminals monochromatic gives a valid instance of 314.1.
This samples many such colourings per bridge and asks whether the mechanism fires for SOME
colouring.  It also diagnoses the misses: is the obstruction that no vertex carries a Y at
all, that Y's exist but never with uniform leg parity, or that uniform-parity Y's exist but
never at two ADJACENT vertices?
"""
import sys, random, itertools, collections
from collections import deque
exec(open('tools/yparity315.py').read().split('rng=random.Random')[0].split('"""',2)[2])

def fires(a, terms, E, N, rng):
    P={v:y_parities(a,v,terms,N,rng,trials=60) for v in a if v not in terms}
    anyY   = any(len(P[v])>0 for v in P)
    uni    = {v for v in P if (0,0,0) in P[v] or (1,1,1) in P[v]}
    hit=False
    for e in E:
        if e in N: continue
        u,v=tuple(e)
        if u in terms or v in terms: continue
        if any(val in P[u] and val in P[v] for val in ((0,0,0),(1,1,1))): hit=True;break
    return hit, anyY, len(uni)

def colourings(a, terms, rng, k):
    out=[]; seen=set()
    for _ in range(k*8):
        col={v:rng.randint(0,1) for v in a}
        for x in terms: col[x]=0
        if rng.random()<0.6:
            imp=True
            while imp:
                imp=False
                for v in a:
                    if v in terms: continue
                    nc=sum(1 for w in a[v] if col[w]==col[v])
                    if nc>len(a[v])-nc: col[v]^=1; imp=True
        key=tuple(sorted((str(v),col[v]) for v in col))
        if key in seen: continue
        seen.add(key); out.append(dict(col))
        if len(out)>=k: break
    return out

rng=random.Random(int(sys.argv[1]) if len(sys.argv)>1 else 31)
st=collections.Counter(); diag=[]
for _ in range(int(sys.argv[2]) if len(sys.argv)>2 else 45):
    a,terms=make_bridge(rng.choice([10,12,14]),rng)
    if a is None: continue
    E={frozenset((u,w)) for u in a for w in a[u]}
    t=tau_edge(E)
    if t<2: continue
    st['bridges with tau>=2']+=1
    ever=False; bestdiag=None; nc=0
    for col in colourings(a,terms,rng,24):
        N={e for e in E if col[tuple(e)[0]]==col[tuple(e)[1]]}
        if not N: continue
        nc+=1
        hit,anyY,nuni=fires(a,terms,E,N,rng)
        if hit: ever=True; break
        if bestdiag is None or nuni>bestdiag[1]: bestdiag=(anyY,nuni,len(N))
    st['FIRES for some colouring' if ever else '*** no colouring works ***']+=1
    st['colourings tried (mean bucket %d)'%min(nc,25)]+=0
    if not ever and bestdiag is not None:
        diag.append((len(a),t,bestdiag))
for k in sorted(st): print('  %-42s %s'%(k,st[k]))
if diag:
    print('  misses (n, tau, (someY, #uniform-parity vertices, |N|)):')
    for d in diag[:10]: print('    ',d)
