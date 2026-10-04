"""hex2c353.py -- (HEX-2C): a terminal 6-cycle whose complement is 2-CONNECTED and
non-bipartite.

Derived theorem: if B satisfies (SC) then NO terminal 6-cycle has a 2-connected non-bipartite
complement - because Lemma B would then give both parities from b_i to EVERY outside
neighbour, any adjacent pair would work, and the Terminal Cycle Theorem would refute (SC).
So the unconditional statement
     (HEX-2C)  every bridge with tau >= 2 has a 6-cycle through a terminal whose complement
               is 2-connected and non-bipartite
implies (N2) at a stroke.  Here real bridges ARE the right population: the statement is
unconditional.  Measured across sizes.
"""
import sys, random, itertools, collections
from collections import deque
exec(open('tools/canon345.py').read().split('rng=random.Random')[0].split('"""',2)[2])

def six_cycles_through(a, terms):
    out=[]; seen=set()
    for b in terms:
        nb=sorted(a[b],key=str)
        if len(nb)!=2: continue
        p,q=nb
        # paths p -> q of length 4 avoiding b  => cycle of length 6 through b
        stack=[(p,[p])]
        while stack:
            u,path=stack.pop()
            if len(path)==5:
                if u==q:
                    c=[b]+path
                    k=frozenset(c)
                    if k not in seen and len(set(c))==6: seen.add(k); out.append(c)
                continue
            for w in a[u]:
                if w==b or w in path: continue
                stack.append((w,path+[w]))
    return out

rng=random.Random(int(sys.argv[1]) if len(sys.argv)>1 else 4507)
st=collections.Counter()
for NV in (30,40,50):
    for _ in range(int(sys.argv[2]) if len(sys.argv)>2 else 4):
        a,terms=make_real_bridge(NV,rng)
        if a is None: continue
        SC=six_cycles_through(a,terms)
        good=0
        for c in SC:
            Ds=set(c); rest=set(a)-Ds
            if len(rest)<3: continue
            s=next(iter(rest)); C=comp_of(a,Ds,s)
            if C is None or len(C)!=len(rest): continue
            if not nonbip(a,C): continue
            if two_conn(a,C): good+=1
        st['n=%-3d  terminal 6-cycles %2d, of which 2-conn+non-bip %2d'%(len(a),len(SC),good)]+=1
        st['(HEX-2C) HOLDS' if good else '*** (HEX-2C) FAILS for this bridge ***']+=1
for k in sorted(st): print('  %s  %s'%(k.ljust(58),st[k]))
