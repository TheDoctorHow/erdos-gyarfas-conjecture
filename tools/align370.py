"""align370.py -- can D be chosen DISJOINT from a transversal path alpha?

369.4's tension, made concrete.  Under (SC), 256.1 gives an even connecting path alpha
between the failing pair whose INTERIOR is a connected odd-cycle transversal of <= 7
vertices.  Label so the failing pair is (b1,b2) and put D through b3, taking b_i = b1.
If D avoids alpha then alpha lies in B - D and joins b1 to the transversal T, so the odd
cycles of B - D - which all meet T - sit in a block that the b1-to-T chain enters.
THE PLACEMENT QUESTION: can a shortest even cycle through b3 be chosen disjoint from a
b1-b2 path of length <= 8?  Measured on genuine bridges: |D| <= 10 and |alpha| <= 8 in a
graph of 33+ vertices, so it should be routine - but it is a hypothesis and is tested here.
"""
import sys, random, itertools, collections
from collections import deque
exec(open('tools/canon345.py').read().split('rng=random.Random')[0].split('"""',2)[2])

rng=random.Random(int(sys.argv[1]) if len(sys.argv)>1 else 8801)
st=collections.Counter()
for _ in range(int(sys.argv[2]) if len(sys.argv)>2 else 8):
    a,terms=make_real_bridge(rng.choice([30,40]),rng)
    if a is None: continue
    b1,b2,b3=terms
    # even b1-b2 connecting paths of length <= 8 (candidates for alpha)
    A=[p for p in rand_paths(a,b1,b2,rng,4000) if (len(p)-1)%2==0 and len(p)-1<=8]
    if not A: st['no short even b1-b2 path found']+=1; continue
    # even cycles through b3, shortest first
    seen=set(); Ds=[]
    for _ in range(6000):
        c=rand_cycle(a,rng)
        if c and len(c)%2==0 and b3 in c and frozenset(c) not in seen:
            seen.add(frozenset(c)); Ds.append(c)
    if not Ds: st['no even cycle through b3']+=1; continue
    Ds.sort(key=len)
    Lmin=len(Ds[0])
    ok=None
    for D in Ds:
        if len(D)>Lmin+2: break
        for al in A[:400]:
            if not (set(D) & set(al[1:-1])) and b1 not in D:
                ok=(len(D),len(al)-1); break
        if ok: break
    st['shortest |D| through b3 = %d, shortest even b1-b2 path = %d'
       %(Lmin,min(len(p)-1 for p in A))]+=1
    st['DISJOINT PAIR FOUND: |D|=%d |alpha|=%d'%ok if ok else '*** no disjoint pair ***']+=1
for k in sorted(st): print('   %s  %s'%(k.ljust(58),st[k]))
